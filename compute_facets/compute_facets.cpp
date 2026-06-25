/*
  compute_facets.cpp
  -----------
  Compute PALP normal forms of the maximal cones (3D facets + origin) of a 4D
  reflexive lattice polytope, together with the original facet vertices.

  Self-contained: only depends on GlobalP.h (the PALP header bundled with the
  compute_facets module).  No Eigen, no RapidJSON, no MongoDB.

  Build:
    g++ -std=c++11 -O2 compute_facets.cpp -o compute_facets

  Usage (one JSON object per stdin line):
    echo '{"verts":[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-1,-1]]}' \
         | ./compute_facets

  Accepted input keys:
    "index"   -- (optional) integer, passed through to output unchanged
    "verts"   -- JSON array of 4D integer coordinate arrays

  Output (one JSON object per line):
    {
      "index":            <int>,               // passthrough from input 'index' (omitted if input has no 'index')
      "verts":            [[...], ...],        // original polytope vertices
      "dual_verts":       [[...], ...],        // vertices of the dual polytope
      "facets":           [[[...], ...], ...], // 3D facets in dual-vertex order (1:1)
      "facet_nfs":        [[[...], ...], ...], // GL(3,Z) NF of each 3D facet
      "maximal_cone_nfs": [[[...], ...], ...]  // GL(4,Z) NF of conv(facet ∪ {0})
    }
*/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cassert>

// Pull in the PALP implementation directly (header-only, C code).
namespace PALP {
#include "GlobalP.h"
}

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------
typedef std::vector<std::vector<int>> VMatrix;
typedef std::vector<int>              IVec;

// ---------------------------------------------------------------------------
// PALP helpers
// ---------------------------------------------------------------------------

/*
 * Initialise all PALP structures for a polytope given its vertex list.
 *
 *  verts  – rows are lattice points (each row is one vertex)
 *  ppl    – output: filled with ALL lattice points of the polytope
 *  vnl    – output: indices of true vertices within ppl
 *  eql    – output: facet equations
 *  vepm   – output: vertex-equation pairing matrix (caller allocates/zeroes)
 *  finf   – output: face incidence data; pass NULL to skip Make_Incidence
 *
 * Returns 1 if reflexive, 0 otherwise.
 */
int setup_palp(const VMatrix&        verts,
               PALP::PolyPointList*  ppl,
               PALP::VertexNumList*  vnl,
               PALP::EqList*         eql,
               PALP::PairMat         vepm,
               PALP::FaceInfo*       finf)
{
    int dim = (int)verts[0].size();
    int nin = (int)verts.size();

    // Scratch list: vertices only (before Complete_Poly)
    PALP::PolyPointList vlist;
    vlist.n  = dim;
    vlist.np = nin;
    for (int i = 0; i < nin; ++i)
        for (int j = 0; j < dim; ++j)
            vlist.x[i][j] = (long)verts[i][j];

    // Determine true vertices and facet equations
    int isref = PALP::Ref_Check(&vlist, vnl, eql);
    if (!isref)
        PALP::Find_Equations(&vlist, vnl, eql);

    // Deduplicate vertex indices (safety measure)
    std::vector<int> seen;
    for (int i = 0; i < vnl->nv; ++i) {
        int v = vnl->v[i];
        if (std::find(seen.begin(), seen.end(), v) == seen.end())
            seen.push_back(v);
    }
    vnl->nv = (int)seen.size();
    for (int i = 0; i < vnl->nv; ++i) vnl->v[i] = seen[i];

    // Build vertex-equation pairing matrix
    PALP::Make_VEPM(&vlist, vnl, eql, vepm);

    ppl->n  = dim;
    ppl->np = vnl->nv;
    for (int i = 0; i < vnl->nv; ++i)
        for (int j = 0; j < dim; ++j)
            ppl->x[i][j] = vlist.x[vnl->v[i]][j];

    if (finf)
        PALP::Make_Incidence(ppl, vnl, eql, finf);

    return isref;
}

/*
 * Extract all (dim-1)-dimensional facets directly from the VEPM.
 * vepm[e][v] == 0 iff vertex v lies on the facet defined by equation e.
 * This avoids Make_Incidence and its VERT_Nmax constraint on equation count.
 */
std::vector<VMatrix> get_facets(PALP::PolyPointList* ppl,
                                 PALP::VertexNumList* vnl,
                                 PALP::EqList*        eql,
                                 PALP::PairMat        vepm)
{
    int dim = ppl->n;
    std::vector<VMatrix> facets;
    for (int e = 0; e < eql->ne; ++e) {
        VMatrix fm;
        for (int v = 0; v < vnl->nv; ++v) {
            if (vepm[e][v] == 0) {
                IVec pt(dim);
                for (int j = 0; j < dim; ++j) pt[j] = (int)ppl->x[v][j];
                fm.push_back(pt);
            }
        }
        if (!fm.empty()) facets.push_back(fm);
    }
    return facets;
}

/*
 * Column-reduce a list of 4D integer vectors to a Z-basis for their span.
 * Uses the Euclidean (GCD) column reduction algorithm.
 * After reduction the basis is upper-triangular in pivot-row order:
 *   basis[k][pivot_rows[j]] = 0  for all j > k.
 * (Entries above the diagonal, j < k, may be non-zero when the GCD pivot
 * does not divide the corresponding basis entry.)
 * coords_in_basis recovers coordinates by forward substitution.
 */
static std::vector<IVec> col_basis(std::vector<IVec> cols,
                                    std::vector<int>&  pivot_rows)
{
    const int N = 4;
    std::vector<IVec> basis;
    pivot_rows.clear();

    for (int r = 0; r < N && !cols.empty(); ++r) {
        // Euclidean GCD reduction across all columns in row r
        bool changed = true;
        while (changed) {
            changed = false;
            // Move the column with smallest non-zero |entry| in row r to front
            int best = -1;
            for (int j = 0; j < (int)cols.size(); ++j) {
                if (cols[j][r] == 0) continue;
                if (best == -1 || std::abs(cols[j][r]) < std::abs(cols[best][r]))
                    best = j;
            }
            if (best == -1) break;
            if (best != 0) std::swap(cols[0], cols[best]);
            // Reduce all other columns modulo cols[0]
            for (int j = 1; j < (int)cols.size(); ++j) {
                if (cols[j][r] == 0) continue;
                long q = cols[j][r] / cols[0][r];
                if (q == 0) { changed = true; continue; } // will swap next iter
                for (int i = 0; i < N; ++i) cols[j][i] -= q * cols[0][i];
                changed = true;
            }
        }

        // Find the surviving non-zero column in row r
        int pivot = -1;
        for (int j = 0; j < (int)cols.size(); ++j)
            if (cols[j][r] != 0) { pivot = j; break; }
        if (pivot == -1) continue;

        // Canonical sign
        if (cols[pivot][r] < 0)
            for (int i = 0; i < N; ++i) cols[pivot][i] = -cols[pivot][i];

        // Eliminate row r from all existing basis vectors (RCEF)
        for (auto& b : basis) {
            if (b[r] == 0) continue;
            long q = b[r] / cols[pivot][r];
            for (int i = 0; i < N; ++i) b[i] -= q * cols[pivot][i];
        }
        // Eliminate row r from remaining candidate columns
        for (int j = 0; j < (int)cols.size(); ++j) {
            if (j == pivot || cols[j][r] == 0) continue;
            long q = cols[j][r] / cols[pivot][r];
            for (int i = 0; i < N; ++i) cols[j][i] -= q * cols[pivot][i];
        }

        pivot_rows.push_back(r);
        basis.push_back(cols[pivot]);
        cols.erase(cols.begin() + pivot);
    }
    return basis;
}

/*
 * Express a vector d as integer coordinates in the basis produced by col_basis.
 * The basis is upper-triangular in pivot-row order, so coordinates are
 * recovered by forward substitution:
 *   c[k] = (d[pivot_rows[k]] - sum_{j<k} c[j]*basis[j][pivot_rows[k]])
 *           / basis[k][pivot_rows[k]]
 */
static IVec coords_in_basis(const std::vector<IVec>& basis,
                              const std::vector<int>&  pivot_rows,
                              const IVec&              d)
{
    IVec c((int)basis.size(), 0);
    for (int k = 0; k < (int)basis.size(); ++k) {
        long num = d[pivot_rows[k]];
        for (int j = 0; j < k; ++j)
            num -= (long)c[j] * basis[j][pivot_rows[k]];
        c[k] = (int)(num / basis[k][pivot_rows[k]]);
    }
    return c;
}

/*
 * Compute the GL(3,Z) normal form of a 3D facet embedded in Z^4.
 *
 * Algorithm:
 *  1. Translate so facet_verts[0] is at the origin.
 *  2. Find a Z-basis for the 3D sublattice spanned by the difference vectors
 *     (using Euclidean column reduction → RCEF basis).
 *  3. Express all translated vertices as 3D integer coordinates in that basis.
 *  4. Compute the PALP GL(3,Z) NF of the 3D polytope.
 *
 * The result is canonical as an invariant of the facet as an abstract 3D
 * lattice polytope (in the lattice it generates).  It is coarser than the
 * 4D maximal-cone NF: same facet_nf ⟹ abstractly isomorphic 3D polytopes,
 * but the converse need not hold.
 */
VMatrix facet_normal_form(const VMatrix& facet_verts) {
    const int DIM4 = 4, DIM3 = 3;
    int n = (int)facet_verts.size();
    if (n < 2) return {};

    const IVec& v0 = facet_verts[0];

    // Build difference vectors (includes origin = v0 - v0)
    std::vector<IVec> diffs(n, IVec(DIM4, 0));
    std::vector<IVec> nonzero_diffs;
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < DIM4; ++j) diffs[i][j] = facet_verts[i][j] - v0[j];
        nonzero_diffs.push_back(diffs[i]);
    }

    // Find Z-basis for the 3D sublattice
    std::vector<int> pivot_rows;
    std::vector<IVec> basis = col_basis(nonzero_diffs, pivot_rows);
    if ((int)basis.size() != DIM3) return {}; // degenerate — shouldn't happen

    // Convert all translated vertices to 3D coordinates
    VMatrix pts3d;
    for (const IVec& d : diffs)
        pts3d.push_back(coords_in_basis(basis, pivot_rows, d));

    // Set up PALP for a 3D polytope and compute GL(3,Z) NF
    PALP::PolyPointList ppl;
    PALP::VertexNumList vnl;
    PALP::EqList        eql;
    PALP::PairMat       vepm;
    memset(&vepm, 0, sizeof(vepm));

    setup_palp(pts3d, &ppl, &vnl, &eql, vepm, nullptr);

    long pNF[POLY_Dmax][VERT_Nmax];
    PALP::Make_Poly_NF(&ppl, &vnl, &eql, pNF);

    VMatrix nff;
    for (int i = 0; i < vnl.nv; ++i) {
        IVec col(DIM3);
        for (int j = 0; j < DIM3; ++j) col[j] = (int)pNF[j][i];
        nff.push_back(col);
    }
    return nff;
}

/*
 * Compute the PALP normal form of conv(facet_verts ∪ {origin}).
 * The origin is invariant under GL(n,Z), so it remains at (0,...,0) in the NF
 * and is stripped from the returned result.
 */
VMatrix cone_normal_form(const VMatrix& facet_verts, int dim) {
    // Build cone = facet verts + origin
    VMatrix cone_verts = facet_verts;
    cone_verts.push_back(IVec(dim, 0));

    // Temporary PALP workspace
    PALP::PolyPointList ppl;
    PALP::VertexNumList vnl;
    PALP::EqList        eql;
    PALP::PairMat       vepm;
    memset(&vepm, 0, sizeof(vepm));

    setup_palp(cone_verts, &ppl, &vnl, &eql, vepm, nullptr);

    // Compute normal form
    long pNF[POLY_Dmax][VERT_Nmax];
    PALP::Make_Poly_NF(&ppl, &vnl, &eql, pNF);

    // Collect NF columns, skip origin
    IVec zero(dim, 0);
    VMatrix nff;
    for (int i = 0; i < vnl.nv; ++i) {
        IVec col(dim);
        for (int j = 0; j < dim; ++j) col[j] = (int)pNF[j][i];
        if (col != zero) nff.push_back(col);
    }
    return nff;
}

// ---------------------------------------------------------------------------
// Minimal JSON helpers (no external library)
// ---------------------------------------------------------------------------

std::string ivec_to_json(const IVec& v) {
    std::string s = "[";
    for (int i = 0; i < (int)v.size(); ++i) {
        if (i) s += ',';
        s += std::to_string(v[i]);
    }
    return s + ']';
}

std::string vmat_to_json(const VMatrix& m) {
    std::string s = "[";
    for (int i = 0; i < (int)m.size(); ++i) {
        if (i) s += ',';
        s += ivec_to_json(m[i]);
    }
    return s + ']';
}

/*
 * Parse either JSON  [[a,b,c,d],[e,f,g,h],...]
 * or Mathematica     {{a,b,c,d},{e,f,g,h},...}
 * into a VMatrix.  Handles both by normalising { → [ and } → ].
 */
VMatrix parse_verts(std::string s) {
    // Normalise to JSON array syntax
    for (char& c : s) {
        if (c == '{') c = '[';
        else if (c == '}') c = ']';
    }
    // Strip whitespace
    s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());

    VMatrix m;
    size_t pos = 0;
    while (true) {
        size_t open = s.find('[', pos + 1);    // skip the outermost [
        if (open == std::string::npos) break;
        size_t close = s.find(']', open + 1);
        if (close == std::string::npos) break;
        std::string row_s = s.substr(open + 1, close - open - 1);
        if (!row_s.empty() && row_s[0] != '[') {
            IVec row;
            std::stringstream ss(row_s);
            std::string tok;
            while (std::getline(ss, tok, ','))
                if (!tok.empty()) row.push_back(std::stoi(tok));
            if (!row.empty()) m.push_back(row);
        }
        pos = close;
    }
    return m;
}

/*
 * Extract the value associated with `key` in a flat JSON object string.
 * Works for string values (returned without quotes) and array/object values
 * (returned with their brackets, depth-matched).
 * Returns "" if the key is absent.
 */
std::string json_get(const std::string& line, const std::string& key) {
    std::string k = '"' + key + '"';
    size_t p = line.find(k);
    if (p == std::string::npos) return "";
    p = line.find(':', p + k.size());
    if (p == std::string::npos) return "";
    while (++p < line.size() && (line[p] == ' ' || line[p] == '\t')) {}
    if (p >= line.size()) return "";

    char delim = line[p];
    if (delim == '"') {
        // String value
        size_t start = p + 1;
        size_t end   = line.find('"', start);
        return (end == std::string::npos) ? "" : line.substr(start, end - start);
    }
    if (delim == '[' || delim == '{') {
        // Array / object: depth-match
        char open = delim, close = (delim == '[') ? ']' : '}';
        int depth = 0;
        for (size_t i = p; i < line.size(); ++i) {
            if (line[i] == open)  ++depth;
            else if (line[i] == close && --depth == 0)
                return line.substr(p, i - p + 1);
        }
    }
    // Bare integer or other scalar
    size_t end = line.find_first_of(",}", p);
    return line.substr(p, end - p);
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;

        // --- Parse vertices ---
        VMatrix verts;
        std::string vs = json_get(line, "verts");
        if (vs.empty()) vs = json_get(line, "NVERTS");
        if (!vs.empty()) verts = parse_verts(vs);
        if (verts.empty()) { std::cout << line << '\n'; continue; }

        std::string poly_idx = json_get(line, "index");
        int dim = (int)verts[0].size();

        // --- Set up main polytope ---
        PALP::PolyPointList ppl;
        PALP::VertexNumList vnl;
        PALP::EqList        eql;
        PALP::PairMat       vepm;
        memset(&vepm, 0, sizeof(vepm));

        setup_palp(verts, &ppl, &vnl, &eql, vepm, nullptr);

        // --- Dual vertices (= facet equation normals, 1:1 with facets) ---
        VMatrix dual_verts;
        for (int i = 0; i < eql.ne; ++i) {
            IVec dv(dim);
            for (int j = 0; j < dim; ++j) dv[j] = (int)eql.e[i].a[j];
            dual_verts.push_back(dv);
        }

        // --- Facets, facet NFs, and maximal cone NFs (in dual-vertex order) ---
        std::vector<VMatrix> facets = get_facets(&ppl, &vnl, &eql, vepm);
        std::vector<VMatrix> facet_nf_list, maxcone_nf_list;
        for (const VMatrix& fv : facets) {
            facet_nf_list.push_back(facet_normal_form(fv));
            maxcone_nf_list.push_back(cone_normal_form(fv, dim));
        }

        // --- Emit output ---
        std::cout << '{';
        if (!poly_idx.empty()) std::cout << "\"index\":" << poly_idx << ',';
        std::cout << "\"verts\":"      << vmat_to_json(verts)      << ',';
        std::cout << "\"dual_verts\":" << vmat_to_json(dual_verts) << ',';
        std::cout << "\"facets\":[";
        for (int i = 0; i < (int)facets.size(); ++i) {
            if (i) std::cout << ',';
            std::cout << vmat_to_json(facets[i]);
        }
        std::cout << "],\"facet_nfs\":[";
        for (int i = 0; i < (int)facet_nf_list.size(); ++i) {
            if (i) std::cout << ',';
            std::cout << vmat_to_json(facet_nf_list[i]);
        }
        std::cout << "],\"maximal_cone_nfs\":[";
        for (int i = 0; i < (int)maxcone_nf_list.size(); ++i) {
            if (i) std::cout << ',';
            std::cout << vmat_to_json(maxcone_nf_list[i]);
        }
        std::cout << "]}\n";
        std::cout.flush();
    }
    return 0;
}
