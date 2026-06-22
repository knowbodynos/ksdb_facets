/*
  cone_nf.cpp
  -----------
  Compute PALP normal forms of the maximal cones (3D facets + origin) of a 4D
  reflexive lattice polytope, together with the original facet vertices.

  Self-contained: only depends on GlobalP.h (the PALP header bundled with the
  facets module).  No Eigen, no RapidJSON, no MongoDB.

  Build:
    g++ -std=c++11 -O2 -I../modules/facets/code cone_nf.cpp -o cone_nf

  Usage (one JSON object per stdin line):
    echo '{"verts":[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-1,-1]]}' \
         | ./cone_nf
    echo '{"POLYID":1,"NVERTS":"{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},{-1,-1,-1,-1}}"}' \
         | ./cone_nf

  Accepted input keys:
    "verts"   -- JSON array of 4D integer coordinate arrays
    "NVERTS"  -- Mathematica-style string  {{a,b,c,d},{...},...}
    "POLYID"  -- (optional) integer, passed through to output

  Output (one JSON object per line):
    {
      "POLYID": <int>,              // if present in input
      "cones": [
        { "facet_verts": [[...], ...], "cone_nf": [[...], ...] },
        ...                          // one entry per facet
      ],
      "unique_cones": [
        { "cone_nf": [[...], ...], "multiplicity": N,
          "example_facet": [[...], ...] },
        ...                          // deduplicated by normal form
      ]
    }
*/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <bitset>
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

    // Fill ppl starting from the true vertices, then Complete_Poly adds all
    // interior/boundary lattice points.
    ppl->n  = dim;
    ppl->np = vnl->nv;
    for (int i = 0; i < vnl->nv; ++i)
        for (int j = 0; j < dim; ++j)
            ppl->x[i][j] = vlist.x[vnl->v[i]][j];
    PALP::Complete_Poly(vepm, eql, vnl->nv, ppl);

    // Face incidence (needed for extracting facets of the main polytope)
    if (finf)
        PALP::Make_Incidence(ppl, vnl, eql, finf);

    return isref;
}

/*
 * Decode an INCI bitmask into a sorted list of 0-based vertex indices.
 * Matches the encoding used by PALP's Make_Incidence.
 */
IVec inci_to_indices(PALP::INCI bitp, int nv) {
    const int BMAX = 8 * sizeof(unsigned long long);
    std::bitset<64> bits((unsigned long long)bitp);
    std::string s = bits.to_string();          // s[0] = MSB
    IVec inds;
    for (int i = 0; i < BMAX; ++i) {
        if (s[i] == '1') {
            int idx = nv - BMAX + i;
            if (idx >= 0) inds.push_back(idx);
        }
    }
    std::sort(inds.begin(), inds.end());
    return inds;
}

/*
 * Extract all (dim-1)-dimensional faces (facets) from precomputed PALP data.
 * Each returned VMatrix row is a vertex of that facet (original coordinates).
 */
std::vector<VMatrix> get_facets(PALP::PolyPointList* ppl,
                                 PALP::VertexNumList* vnl,
                                 PALP::FaceInfo*      finf,
                                 int                  dim)
{
    std::vector<VMatrix> facets;
    int nf = finf->nf[dim - 1];
    for (int i = 0; i < nf; ++i) {
        IVec inds = inci_to_indices(finf->v[dim - 1][i], vnl->nv);
        VMatrix fm;
        for (int idx : inds) {
            int vi = vnl->v[idx];
            IVec pt(dim);
            for (int j = 0; j < dim; ++j) pt[j] = (int)ppl->x[vi][j];
            fm.push_back(pt);
        }
        facets.push_back(fm);
    }
    return facets;
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

        std::string polyid = json_get(line, "POLYID");
        int dim = (int)verts[0].size();

        // --- Set up main polytope ---
        PALP::PolyPointList ppl;
        PALP::VertexNumList vnl;
        PALP::EqList        eql;
        PALP::PairMat       vepm;
        PALP::FaceInfo      finf;
        memset(&vepm, 0, sizeof(vepm));

        setup_palp(verts, &ppl, &vnl, &eql, vepm, &finf);

        // --- Facets and cone NFs ---
        std::vector<VMatrix> facets = get_facets(&ppl, &vnl, &finf, dim);

        // Parallel arrays: facet vertices and their cone NFs
        std::vector<VMatrix> facet_list, nf_list;
        // Deduplication: NF JSON string → (nf, multiplicity, example facet)
        std::map<std::string, std::pair<VMatrix, int>> unique;
        std::map<std::string, VMatrix>                 unique_ex;

        for (const VMatrix& fv : facets) {
            VMatrix nff = cone_normal_form(fv, dim);
            facet_list.push_back(fv);
            nf_list.push_back(nff);
            std::string key = vmat_to_json(nff);
            if (!unique.count(key)) {
                unique[key]    = {nff, 1};
                unique_ex[key] = fv;
            } else {
                ++unique[key].second;
            }
        }

        // --- Emit output ---
        std::cout << '{';
        if (!polyid.empty()) std::cout << "\"POLYID\":" << polyid << ',';

        std::cout << "\"cones\":[";
        for (int i = 0; i < (int)facet_list.size(); ++i) {
            if (i) std::cout << ',';
            std::cout << "{\"facet_verts\":" << vmat_to_json(facet_list[i])
                      << ",\"cone_nf\":"     << vmat_to_json(nf_list[i]) << '}';
        }

        std::cout << "],\"unique_cones\":[";
        bool first = true;
        for (const auto& kv : unique) {
            if (!first) std::cout << ',';
            first = false;
            std::cout << "{\"cone_nf\":"      << kv.first
                      << ",\"multiplicity\":" << kv.second.second
                      << ",\"example_facet\":" << vmat_to_json(unique_ex[kv.first])
                      << '}';
        }
        std::cout << "]}\n";
        std::cout.flush();
    }
    return 0;
}
