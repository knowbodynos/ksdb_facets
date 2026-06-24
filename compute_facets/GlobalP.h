#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

/*
These are include files that should exist in your C library.
*/




/*  ============	basic choice of PARAMETERS	      ============  */

#define	                Long            long
#define                 LLong           long long
/*
For reflexive polytopes in 4 or less dimensions, everything should work with
Long set to 32-bit-integer and LLong set to 64 bits.
Many applications will even work with LLong at 32 bits.
For higher dimensional or complicated non-reflexive polytopes it may be
necessary to set even Long to 64 bits.
*/

#define   		POLY_Dmax	4	/* max dim of polytope	    */
/*
POLY_Dmax should be set to the dimension of the polytopes that are analysed.
While the programs still work if POLY_Dmax is set to a higher value, they may
be considerably slowed down.
*/

#if	(POLY_Dmax <= 3)
#define   		POINT_Nmax	40	/* max number of points	    */
#define   		VERT_Nmax	16	/* max number of vertices   */
#define   		FACE_Nmax	30      /* max number of faces      */
#define	                SYM_Nmax	88	/* cube: 2^D*D! plus extra  */

#elif	(POLY_Dmax == 4)
#define   		POINT_Nmax	700	/* max number of points	    */
#define 		VERT_Nmax	34      /* max number of vertices   */
#define   		FACE_Nmax	824	/* max number of faces      */
#define	                SYM_Nmax	1200

#else
#define   		POINT_Nmax	2000000 
#define   		VERT_Nmax	64	/* !! use optimal value !!  */
#define   		FACE_Nmax	10000	/* max number of faces      */
#define                 SYM_Nmax        46080   /* symmetry (P_1)^6: 2^6*6! */
#define			EQUA_Nmax	1280    /* up to 20000 without alloc */
#endif

#ifndef			EQUA_Nmax			/* default setting */
#define                 EQUA_Nmax       VERT_Nmax
#endif
#ifndef	CEQ_Nmax	
#define CEQ_Nmax        EQUA_Nmax
#endif


/*
POINT_Nmax, VERT_Nmax and FACE_Nmax denote the maximal numbers of points,
vertices and faces, respectively.
SYM_Nmax is the maximal number of symmetries of a polytope, i.e. the order of
the finite subgroup S of the group GL(n,Z) of lattice automorphisms that leaves
a polytope invariant.
EQUA_Nmax denotes the maximal number of facets (given by equations) of a
polytope. By duality this is just the number of vertices of the dual polytope,
so it makes sense to have the default setting EQUA_Nmax = VERT_Nmax.
In applications not related to reflexive polytopes or in large dimensions a
larger value may be useful. While CPU-time is almost independent of EQUA_Nmax,
it strongly depends on VERT_Nmax/32 (i.e. use 32, 64, 96, ...).
Our settings for dimensions less than or equal to 4 are such that they work
for any reflexive polytope.
*/

#define                 AMBI_Dmax       (5 * POLY_Dmax)	/* default setting */
/*
If a polytope is determined by a combined weight system it is first realised
by an embeddeding in an ambient space of dimension (Poly-dim + number of
weight systems). AMBI_Dmax is the maximal dimension of this ambient space.
*/


#define                 FIB_Nmax	3000 /*NOW: 27/5/11 default setting*/
/*
Given a polytope P* it is possible to analyze the IP simplices among its 
points. These simplices are given in terms of weight relations among points 
of P*. FIB_Nmax is the maximal number of allowed relations.
*/


#define  CD2F_Nmax               FACE_Nmax
/*
Max number of codimension 2 faces.
*/


#define GL_Long		Long
/*
Uses W_to_GLZ like in Rat.c
*/


extern FILE *inFILE, *outFILE;
FILE * inFILE = stdin;
FILE * outFILE = stdout;
/*
Ascii-files for input and output. If not given in the parameter list they
default to stdin and stdout, respectively.
*/


/*  ==========         Global typedefs           		==========  */

typedef struct {int n, np; Long x[POINT_Nmax][POLY_Dmax];}   PolyPointList;
/*
A list (not necessarily complete) of lattice points of a polytope.
P.x[i][j] is the j'th coordinate of the i'th lattice point.
P.n is the dimension of the polytope and P.np the number of points in the list.
*/

typedef struct {int v[VERT_Nmax]; int nv;}                   VertexNumList;
/*
The list of vertices of a polytope, referring to some PolyPointList P.
The j'th coordinate of the i'th vertex is then given by P.x[V.v[i]][j].
V.nv is the number of vertices of P.
*/

typedef struct {Long a[POLY_Dmax], c;}                       Equation;
/*
This structure determines an equation of the type ax+c=0, explicitly:
sum_{i=1}^n E.a[i] x_i + E.c = 0.
*/

typedef struct {int ne; Equation e[CEQ_Nmax];}              CEqList;

/* A list of CEq; CEq.ne is the number of Equations in the list */

typedef struct {int ne; Equation e[EQUA_Nmax];}		     EqList;
/*
A list of equations; EL.ne is the number of equations in the list.
*/

typedef struct {EqList B; Long W[AMBI_Dmax][AMBI_Dmax], d[AMBI_Dmax]; 
  int nw, N, z[POLY_Dmax][AMBI_Dmax], m[POLY_Dmax], nz, index;}       CWS;
/*
Combined weight system: W[i][j] and d[i] are the j'th weight and the "degree"
of the i'th weight system, respectively; nw is the number of weight systems,
N is the dimension of the ambient space.
z[i][j]/m[i] are the phases of nz symmetries for polytopes on sublattices.
B describes the ambient space coordinate hyperplanes in terms of the new
(non-redundant) coordinates.
*/

typedef Long PairMat[EQUA_Nmax][VERT_Nmax];
/*
The matrix whose entries are the pairings av+c between the vertices v and
the equations (a,c).
*/

typedef struct {int mp, mv, np, nv, n, pic, cor, h22, h1[POLY_Dmax-1];}
                                                             BaHo;
/*
This structure is related to Batyrev's formulas for Hodge numbers.
n     ... dimension of the polytope
pic   ... Picard number
cor   ... sum of correction terms
h1[i] ... Hodge number h_{1i}
h22   ... Hodge number h_{22} (if n = 5)
mp, mv, np, nv denote the numbers of points/vertices in the M and N lattices,
repectively.
*/

typedef struct {
	Long W[FIB_Nmax][VERT_Nmax]; 
	int nw, PS, ZS, nv, f[VERT_Nmax],r[VERT_Nmax],nf,nz[FIB_Nmax], n0[FIB_Nmax],
         Z[FIB_Nmax][VERT_Nmax], M[FIB_Nmax];
	GL_Long G[VERT_Nmax][POLY_Dmax][POLY_Dmax];
	PolyPointList *P;
} FibW;
/*
This list is an extension of the PolyPointList with the combined weight system.
W[i][j] is the j'th weight; nw is the number of weight systems.
*/



/*  ==========         I/O functions (from Coord.c)		==========  */

int  Read_CWS_PP(CWS *C, PolyPointList *P);
/*
Reads either a CWS or a PolyPointList.
If *C is read, the PolyPointList *P determined by *C is calculated, otherwise
C->nw is set to 0 to indicate that no weight has been read.
CWS-input consists of a single line of the form
d1 w11 w12 ... d2 w21 w22 ...,
whereas PolyPointList-input begins with a line
#columns #lines
followed by #lines further lines. It reads P->x as a matrix such that
either P->n = #columns and P->np = #lines or vice versa (the result is
unique because of P->np > P->n).
*/

int  Read_CWS(CWS *_CW, PolyPointList *_P);
/*
 Reads CWS input *C, the PolyPointList *P determined by *C is calculated.
*/

int  Read_PP(PolyPointList *_P);
/*
Reads the PolyPointList input *P
*/

void Print_PPL(PolyPointList *P, const char *comment);
void Print_VL(PolyPointList *P, VertexNumList *V, const char *comment);
void Print_EL(EqList *EL, int *n, int suppress_c, const char *comment);
void Print_Matrix(Long Matrix[][VERT_Nmax], int n_lines, int n_columns,
		  const char *comment);
/*
Each of these routines prints a matrix in the format
#columns #lines  *comment
line_0
...
line_{#lines-1}.
With Print_PPL and Print_VL, points/vertices are displayed as column vectors
if there's enough space and as row vectors otherwise.
Print_EL always displays equations in line format.
If *suppress_c is zero line_i reads
EL->e[i].a[0] ... EL->e[i].a[*n-1]  EL->e[i].c,
otherwise the last entry EL->e[i].c is suppressed so that the
resulting output can be used as input for Read_CWS_PP.
*/

void Print_CWH(CWS *C, BaHo *BH);
/*
Writes a single line that reproduces *C (if C->nw isn't 0, i.e. if the
input was of CWS type), information on the numbers of points and
vertices of the polytope and its dual, and the Hodge numbers of the
corresponding Calabi-Yau variety.
*C is reproduced in the format
d1 w11 w12 ... d2 w21 w22 ...
Information on the polytope is given in the form
M:mp mv N:np nv
for reflexive polytopes.
Here mp and mv are the numbers of lattice points and vertices of the
polytope, respectively, and np and nv denote the corresponding numbers
for the dual polytope.
If a polytope is not reflexive, "N:np nv" is replaced by "F:ne" (the
number of facets/equations).
Hodge number information is given in the format
H: h11 h12 ... h1(n-2) [chi],
where the h1i are the corresponding Hodge numbers and chi is the Euler
number. This output is suppressed for polytopes that are not reflexive.
As an example, the complete output for the quintic threefold reads
5 1 1 1 1 1 M:126 5 N:6 5 H:1,101 [-200].
*/


/*  ==========              From Polynf.c                        ========== */

int  Make_Poly_Sym_NF(PolyPointList *P, VertexNumList *VNL, EqList *EL,
		      int *SymNum, int V_perm[][VERT_Nmax],
		      Long NF[POLY_Dmax][VERT_Nmax], int t, int S, int N);
/*
Given *P, *VNL and *EL, the following objects are determined:
the number *SymNum of GL(n,Z)-symmetries of the polytope,
the *SymNum vertex permutations V_perm realising these symmetries,
the normal form coordinates NF of the vertices,
the number of symmetries of the vertex pairing matrix
    (this number is the return value of Make_Poly_Sym_NF).
If t/S/N are non-zero, the output of the corresponding options of poly.x
is displayed.
*/

void IP_Simplex_Decomp(Long CM[][POLY_Dmax], int p, int d,
        int *nw, Long W[][VERT_Nmax], int Wmax, int codim);
/*
Given the matrix CM of coordinates of p points in Z^d, the list W[i] of *nw
weight systems corresponding to IP-simplices spanned by the points in CM is
created.
If codim!=0 only the IP-simplices with dimension > 1 and codimension
between 1 and codim are computed.
It is assumed that p<=VERT_Nmax and that W can hold at least Wmax sets of
coefficients.
*/

void IP_Simplices(PolyPointList *P, int nv, int PS, int VS, int CD);
/*
Realizes the -P,-V,-Z, and fibration options of poly (the results of this
routine are displayed as output; *P is not modified).
*/

int  Sublattice_Basis(int d, int p, Long *P[],     /* return index=det(D) */
	Long Z[][VERT_Nmax], Long *M, int *r, Long G[][POLY_Dmax], Long *D);
/*
Given a vector P[] of pointers at p points in N=Z^d that generate a 
(sub)lattice N' of the same dimension d, the following data are determined:
D[i] with 0 <= i < d  such that the lattice quotient N/N' is the product of 
cyclic groups Z_{D[i]} with D[i] dividing D[i+1], and a GL(d,Z) matrix G 
corresponding to a base change P->GxP such that the i'th new coordinate of 
each of the lattice points is divisible by D[i].
If p<=VERT_Nmax the program also computes *r coefficient vectors Z[i] for 
linear combinations of the points on P that are M[i]-fold multiples of
primitive lattice vectors, where M[i]=D[d-i] for i<*r.
If p>VERT_Nmax it is asserted that the index of the lattice quotient is 1.
*/

void Make_Poly_UTriang(PolyPointList *P);
/*
A coordinate change is performed that makes the matrix P->x upper triangular,
with minimal entries above the diagonal.
*/

void Make_ANF(PolyPointList *P, VertexNumList *V, EqList*E, 
	      Long ANF[][VERT_Nmax]);
/*
Given *P, *V and *E, the affine normal form ANF (i.e., a normal form
that also works for non-reflexive polytopes), is computed.
*/

int SimpUnimod(PolyPointList *P, VertexNumList *V, EqList *E, int vol);
/*
If vol is 0, the return value is 1 if all facets are simplicial, 0 otherwise
If vol is not 0, the return value is 1 if all facets are unimoular
(i.e. of volume 1) and 0 otherwise.
*/

int ConifoldSing(PolyPointList *P, VertexNumList *V, EqList *E,
		 PolyPointList *dP, EqList *dE, int CYorFANO);
/*
Realizes the -C1 or -C2 options of poly for CYorFANO being 1 or 2, respectively.
*/

int  Fano5d(PolyPointList *, VertexNumList *, EqList *);
/*
Realizes the -U5 option of poly.
*/

void Einstein_Metric(CWS *CW,PolyPointList *P,VertexNumList *V,EqList *E);
/*
Realizes the -E option of poly.
*/

int  Divisibility_Index(PolyPointList *P, VertexNumList *V);
/*
Returns the largest integer g for which *P is a g-fold multiple of some
other polytope.
*/

Long LatVol_Barycent(PolyPointList *P, VertexNumList *V, Long *B, Long *N);
/*
Given *P and *V, the coordinates of the barycenter of *P are computed (with 
the i'th coordinate as B[i] / *N) and the lattice volume of *P is returned.
*/

void IPs_degD(PolyPointList *P, VertexNumList *V, EqList *E, int l);
/*
 *P is interpreted as the origin and the first level of a Gorenstein cone. 
The points of the cone up to level l are computed and displayed together with 
information on the type of face of the cone they represent (option -B# of poly).
*/

void Make_Facet(PolyPointList *P, VertexNumList *V, EqList *E, int e, 
		Long vertices_of_facet[POLY_Dmax][VERT_Nmax], int *nv_of_facet);
/*
The e'th facet of *P is determined as a (P->n-1)-dimensional polytope:
*nv_of_facet vertices represented by vertices_of_facet.
*/

/*  ==========     General purpose functions from Vertex.c   	==========  */

void swap(int *i,int *j);
/*
Swaps *i and *j.
*/

void Sort_VL(VertexNumList *V);
/*
Sorts the entries _V->v[i] in ascending order.
*/

// Long Eval_Eq_on_V(Equation *E, Long *V, int n);
/*
Evaluates E on V, i.e. calculates \sum_{i=0}^{n-1} E->a[i] * V[i] + E->c.
*/

int  Span_Check(EqList *EL, EqList *HL, int *n);
/*
Returns 1 if every equation of *HL is contained in *EL and 0 otherwise.
*n is the dimension.
*/

int  Vec_Greater_Than(Long *X, Long *Y, int n);
/*
Returns 1 if *X > *Y in the sense that X[i] > Y[i] for the first i where
X[i] and Y[i] differ, returns 0 if *X < *Y and gives an error message if
X[i] equals Y[i] for all i in {0,...n-1}.
*/

int Vec_is_zero(Long *X, int n);
/*
Returns 1 if X[i]==0 for 0<=i<n; returns 0 otherwise.
*/

void Swap_Vecs(Long *X, Long *Y, int n);
/*
Exchanges the n-dimensional vectors X and Y.
*/

Equation EEV_To_Equation(Equation *E1, Equation *E2, Long *V, int n);
/*
Returns the equation describing the span of the vector V and the intersection
of the hyperplanes corresponding to E1 and E2; n is the dimension.
*/

void Make_VEPM(PolyPointList *P, VertexNumList *VNL, EqList *EL, PairMat PM);
/*
Calculates the matrix of pairings between the vertices in VNL and the
equations in EL.
*/

int EL_to_PPL(EqList *EL, PolyPointList *DP, int *n);
/*
Converts *EL to the incomplete PolyPointList *DP corresponding to the dual
polytope; *n is the dimension. Returns 1 if all equations of *EL are at
distance 1 from the origin and 0 otherwise.
*/

int VNL_to_DEL(PolyPointList *P, VertexNumList *V, EqList *DE);
/*
Converts *V, which refers to *P, into the list *DE of equations of the
dual polytope (assuming reflexivity).
Returns 0 if _V->nv exceeds EQUA_Nmax and 1 otherwise.
*/

int Transpose_PM(PairMat PM, PairMat DPM, int nv, int ne);
/*
Transposes PM into DPM; returns 1 if the dimensions nv, ne are within their
limits and 0 otherwise.
*/


/*  ==========   Polytope analysis functions (from Vertex.c)    ==========  */

int  Find_Equations(PolyPointList *P, VertexNumList *VNL, EqList *EL);
/*
For the polytope determined by P, *VNL and *EL are calculated.
*VNL is the complete list of vertices of P.
*EL is the complete list of equations determining the facets of P.
Find_Equations returns 1 if P has IP property (i.e., it has the
origin in its interior) and 0 otherwise.
*/

int  IP_Check(PolyPointList *P, VertexNumList *VNL, EqList *EL);
/*
Same as Find_Equations, but returns immediately without
calculating *VNL and *EL if P does not have the IP property.
*/

int  Ref_Check(PolyPointList *P, VertexNumList *VNL, EqList *EL);
/*
Returns 1 if P is reflexive and 0 otherwise.
Only in the reflexive case *VNL and *EL are calculated.
*/

void Make_Dual_Poly(PolyPointList *P, VertexNumList *VNL, EqList *EL,
		    PolyPointList *DP);
/*
Given P, VNL and EL for a reflexive polytope, the complete list *DP
of lattice points of the dual polytope is determined.
*/

void Complete_Poly(Long VPM[][VERT_Nmax],EqList *E,int nv,PolyPointList *P);
/*
Given the vertex pairing matrix VPM, the EqList *E and the number nv of
vertices, the complete list of lattice points *P is determined.
*/

void RC_Calc_BaHo(PolyPointList *P, VertexNumList *VNL, EqList *EL,
		  PolyPointList *DP, BaHo *BH);
/*
Given *P, *VNL, *EL and *DP (points of dual polytope) as input, the elements
of *BH are calculated. *P must be reflexive; *P and *DP must be complete.
*/


/*  ======  typedefs and functions (from Vertex.c) related to INCIs  ====  */

#define                 INT_Nbits            32
#define                 LONG_LONG_Nbits      64
/*
These numbers should be set to the actual numbers of bits occupied by the
structures "unsigned int" and "unsigned long long" in your version of C.
If they are set to lower values, everything still works but may be
considerably slowed down.
*/

#if (VERT_Nmax <= INT_Nbits)
typedef		        unsigned int            INCI;
#elif (VERT_Nmax <= LONG_LONG_Nbits)
typedef		        unsigned long long	INCI;
#else
#define I_NUI     ((VERT_Nmax-1)/INT_Nbits+1)
typedef struct {unsigned int ui[I_NUI];}   INCI;
#endif
/*
An INCI encodes the incidence relations between a face and a list of
vertices as a bit pattern (1 if a vertex lies on the face, 0 otherwise).
Depending on the allowed number VERT_Nmax of vertices, a single "unsigned int"
or "unsigned long long" may be sufficient.
If VERT_Nmax is larger than the number of bits in a "long long integer", an
array of unsigned integers is used to simulate an integer type of the required
size.
*/

typedef struct {int nf[POLY_Dmax+1];			  /* #(faces)[dim]  */
 	INCI v[POLY_Dmax+1][FACE_Nmax]; 		  /*  vertex info   */
 	INCI f[POLY_Dmax+1][FACE_Nmax]; 		  /* V-on-dual info */
 	Long nip[POLY_Dmax+1][FACE_Nmax];		   /* #IPs on face  */
 	Long dip[POLY_Dmax+1][FACE_Nmax];} 	FaceInfo;  /* #IPs on dual  */
/*
nf[i] denotes the number of faces of dimension i
   (the number of faces of dimension n-i-1 of the dual polytope).
v[i][j] encodes the incidence relation of the j'th dim-i face with the vertices
nip[i][j] is the number of interior points of the j'th dim-i face.
f[i][j] and dip[i][j] give the same informations for the dual (n-i-1
   dimensional) faces, with f[i][j] referring to the dual vertices.
*/

#if (VERT_Nmax <= LONG_LONG_Nbits)
#define INCI_M2(x)     ((x) % 2)              /* value of first bit      */
#define	INCI_AND(x,y)  ((x) & (y))            /* bitwise logical and     */
#define	INCI_OR(x,y)   ((x) | (y))            /* bitwise logical or      */
#define	INCI_XOR(x,y)  ((x) ^ (y))            /* bitwise exclusive or    */
#define	INCI_EQ(x,y)   ((x) == (y))           /* check on equality       */
#define INCI_LE(x,y)   INCI_EQ(INCI_OR(x,y),y)/* bitwise less or equal */
#define INCI_EQ_0(x)   INCI_EQ(x,INCI_0())    /* check if all bits = 0   */
#define INCI_0()       (0)                    /* set all bits to 0       */
#define INCI_1()       (1)                    /* set only first bit to 1 */
#define INCI_D2(x)     ((x) / 2)              /* shift by one bit        */
#define INCI_PN(x,y)   (2 * (x) + !(y))       /* shift and set first bit */
/*
For an INCI defined as a single unsigned (long long) integer whose bits are
regarded as representing incidences, these are useful definitions.
INCI_PN is particularly useful when a new vertex is added: if x represents
an equation E w.r.t. some vertex list and y is the result of evaluating E
on some new vertex V, then INCI_PN(x,y) represents x w.r.t. the vertex list
enhanced by V.
*/

#else
#define INCI_M2(x)      ((x).ui[0] % 2)
INCI INCI_AND(INCI x, INCI y);
INCI INCI_OR(INCI x, INCI y);
INCI INCI_XOR(INCI x, INCI y);
int  INCI_EQ(INCI x, INCI y);
int  INCI_LE(INCI x, INCI y);
int  INCI_EQ_0(INCI x);
INCI INCI_0();
INCI INCI_1();
INCI INCI_D2(INCI x);
INCI INCI_PN(INCI x, Long y);
#endif
/*
If we need more bits than can be represented by a single unsigned long long,
these routines are designed to simulate the above definitions.
*/

int  INCI_abs(INCI X);
/*
Returns the number of bits of X whose value is 1.
*/

int  Print_INCI(INCI X);
/*
Prints X as a pattern of 0's and 1's, omitting the 0's after the last 1.
*/

INCI Eq_To_INCI(Equation *E, PolyPointList *P, VertexNumList *VNL);
/*
Converts *E to an INCI.
*/

void Make_Incidence(PolyPointList *P, VertexNumList *VNL, EqList *EL,
                    FaceInfo *FI);
/*
Creates the structure FaceInfo *FI from *P, *VNL and *EL.
*/

void Print_FaceInfo(int n, FaceInfo *FI);
/*
Displays the information contained in the FaceInfo *FI.
*/

/* =========== Added functions from Vertex.c and Rat.h ============= */

void swap(int* i,int* j) {register int k; k=*i; *i=*j; *j=k;}

Long RoundQ(Long N,Long D)
{    Long F; if(D<0) {D=-D; N=-N;} F=N/D; return F+(2*(N-F*D))/D; 
}

LLong LFgcd(register LLong a, register LLong b)     /* Fast greatest common div */
{
     while( a %= b ) if ( !(b %= a) ) return a; return  b;
}

Long Egcd(register Long A0, register Long A1, Long *Vout0, Long *Vout1)  
{    register Long V0=A0, V1=A1, A2, X0=1, X1=0, X2=0;
     while((A2 = A0 % A1)) { X2=X0-X1*(A0/A1); A0=A1; A1=A2; X0=X1; X1=X2; }
     *Vout0=X1, *Vout1=(A1-(V0) * X1)/ (V1); return A1;
}

Long W_to_GLZ(Long *W, int *d, Long **GLZ)		
{    int i, j; Long G, *E=*GLZ, *B=GLZ[1]; for(i=0;i<*d;i++) assert(W[i]!=0);
     for(i=1;i<*d;i++)for(j=0;j<*d;j++)GLZ[i][j]=0;
     G=Egcd(W[0],W[1],&E[0],&E[1]); B[0]=-W[1]/G; B[1]=W[0]/G;
     for(i=2;i<*d;i++)                   
     {  Long a, b, g=Egcd(G,W[i],&a,&b); B=GLZ[i];
        B[i]= G/g; G=W[i]/g; for(j=0;j<i;j++) B[j]=-E[j]*G;  /* B=Base-Line */
        for(j=0;j<i;j++) E[j]*=a; E[j]=b;                     /* next Egcd */
        for(j=i-1;0<j;j--)                         /* I M P R O V E M E N T */
	{   int n; Long *Y=GLZ[j], rB=RoundQ(B[j],Y[j]), rE=RoundQ(E[j],Y[j]);
	/*  rB=CeilQ(B[j],Y[j]), rE=CeilQ(E[j],Y[j]);			*/
	/*  printf(" [%d/%d -> %d] ",(int)B[j],(int)Y[j],(int)rB);
	    printf(" [%d/%d -> %d] ",(int)E[j],(int)Y[j],(int)rE); 	*/
            for(n=0;n<=j;n++) { B[n] -= rB*Y[n]; E[n] -= rE*Y[n]; } 
	}   G=g;
     } 
     return G;
}

Long NNgcd(register Long a, register Long b)  /* NonNegative gcd handling 0 */
{
     a = (a<0) ? -a:a; b = (b<0) ? -b:b; if (!b) return a;
     while( a %= b ) if ( !(b %= a) ) return a; return  b;
}

LLong LNNgcd(register LLong a, register LLong b)  /* NonNegative gcd handling 0 */
{
     a = (a<0) ? -a:a; b = (b<0) ? -b:b; if (!b) return a;
     while( a %= b ) if ( !(b %= a) ) return a; return  b;
}

Long Eval_Eq_on_V(Equation *E, Long *V, int i){    
  Long p=E->c; while(i--) p+=V[i]*E->a[i]; return p;
}

int Vec_Greater_Than(Long *X, Long *Y, int i){	    /* return 1 iff `X > Y' */
  while(i--) {if(X[i]>Y[i]) return 1; if(X[i]<Y[i]) return 0;} 
  puts("Identical points in Vec_Greater_Than !!"); exit(0); return 0;
}

void Make_VEPM(PolyPointList *_P, VertexNumList *_V, EqList *_E, 
	       PairMat PM){
  int i, j;
  for (i=0;i<_E->ne;i++) for (j=0;j<_V->nv;j++) 
    PM[i][j]=Eval_Eq_on_V(&_E->e[i],_P->x[_V->v[j]],_P->n);
}

int Transpose_PM(PairMat PM, PairMat DPM, int nv, int ne){
  int i, j;
  if((nv>EQUA_Nmax)||(ne>VERT_Nmax)) return 0;
  for (i=0;i<ne;i++) for (j=0;j<nv;j++) DPM[j][i]=PM[i][j];
  return 1;
}

int VNL_to_DEL(PolyPointList *_P, VertexNumList *_V, EqList *_DE){
  int i, j;
  if (_V->nv>EQUA_Nmax) return 0;
  _DE->ne=_V->nv;
  for (i=0;i<_V->nv;i++){ 
    for (j=0;j<_P->n;j++) _DE->e[i].a[j]=_P->x[_V->v[i]][j];
    _DE->e[i].c=1;  }
  return 1;
}

#define	 LLong_EEV		(1)    /* 1 @ [4662 4 20 333 422 1554 2329] */
#define  TEST_EEV	      	(0)	       /* compare Long to LLong EEV */



Equation EEV_To_Equation(Equation *_E1, Equation *_E2, Long *_V, int n){
  /* Calculate the equation spanned by _V and the intersection of _E1, _E2  */
  int i; Long l, m, g; Equation Eq;
  l=Eval_Eq_on_V(_E2,_V,n);
  m=Eval_Eq_on_V(_E1,_V,n);
  g=NNgcd(l,m); assert(g); l/=g; m/=g;
#if ((!(LLong_EEV))||(TEST_EEV))			    /* Long version */
  for(i=0;i<n;i++) Eq.a[i]=l*_E1->a[i]-m*_E2->a[i];
  { int gcd=Eq.c=l*_E1->c-m*_E2->c;
    for(i=0;i<n;i++) gcd=NNgcd(gcd,Eq.a[i]); assert(gcd);
    if (gcd!=1) { for(i=0;i<n;i++) Eq.a[i]/=gcd; Eq.c/=gcd;}}
#endif
#if ((LLong_EEV)||(TEST_EEV))				   /* LLong version */
  { LLong A[POLY_Dmax], C, G; for(i=0;i<n;i++) 
	A[i]=((LLong) l)*((LLong)_E1->a[i])-((LLong)m)*((LLong)_E2->a[i]);
    G=C=((LLong) l)*((LLong)_E1->c)-((LLong)m)*((LLong)_E2->c);
    for(i=0;i<n;i++) G=LNNgcd(G,A[i]); assert(G);
    if(G!=1) {C/=G; for(i=0;i<n;i++) A[i]/=G;}
#if	(TEST_EEV)						 /* Compare */
    {	int e=(Eq.c!=C); for(i=0;i<n;i++) if(Eq.a[i]!=A[i]) e=1;     
	if(e) { printf("Error in EEV: l=%d m=%d g=%d\n",l,m,g);
	for(i=0;i<n;i++)printf("%d ",_E1->a[i]);printf("  %d = E1\n",_E1->c);
	for(i=0;i<n;i++)printf("%d ",_E2->a[i]);printf("  %d = E2\n",_E2->c);
	for(i=0;i<n;i++)printf("%d ",Eq.a[i]);printf("  %d = Eq\n",Eq.c);
	for(i=0;i<n;i++)printf("%d ",A[i]);printf("  %d = LL_Eq\n",C);
	exit(0); }
    }
#else
    Eq.c=C; for(i=0;i<n;i++) Eq.a[i]=A[i];
#endif
  }
#endif
  return Eq;
}


/*  ======================================================================  */
/*  ==========		     			  		==========  */
/*  ==========	   I N C I D E N C E S (as bit patterns)	==========  */
/*  ==========							==========  */ 
/*  ======================================================================  */

#if (VERT_Nmax > LONG_LONG_Nbits)
INCI INCI_AND(INCI x, INCI y){
  INCI z; int i; for (i=0;i<I_NUI;i++) z.ui[i]=(x.ui[i])&(y.ui[i]); return z;}
INCI INCI_OR(INCI x, INCI y){
  INCI z; int i; for (i=0;i<I_NUI;i++) z.ui[i]=(x.ui[i])|(y.ui[i]); return z;}
INCI INCI_XOR(INCI x, INCI y){
  INCI z; int i; for (i=0;i<I_NUI;i++) z.ui[i]=(x.ui[i])^(y.ui[i]); return z;}
int  INCI_EQ(INCI x, INCI y){
  int i; for (i=0;i<I_NUI;i++) if ((x.ui[i])!=(y.ui[i])) return 0; return 1;}
int  INCI_EQ_0(INCI x){
  int i; for (i=0;i<I_NUI;i++) if (x.ui[i]) return 0; return 1;}
int  INCI_LE(INCI x, INCI y){
  int i;
  /* for (i=0;i<I_NUI;i++) if ((x.ui[i]|y.ui[i])!=(y.ui[i])) return 0; */
  unsigned int *X=x.ui, *Y=y.ui;
  for (i=0;i<I_NUI;i++) if ((X[i]|Y[i])!=(Y[i])) return 0;
  return 1;}
INCI INCI_0(){
  INCI z; int i; for (i=0;i<I_NUI;i++) z.ui[i]=0; return z;}
INCI INCI_1(){
  INCI z; int i; z.ui[0]=1; for (i=1;i<I_NUI;i++) z.ui[i]=0; return z;}
INCI INCI_PN(INCI x, Long y){
  INCI z; int i; 
  z.ui[0]=(x.ui[0]<<1)|(!y); 
  for (i=1;i<I_NUI;i++) z.ui[i]=(x.ui[i]<<1)|(x.ui[i-1]>>(INT_Nbits-1)); 
  return z;}
INCI INCI_D2(INCI x){
  INCI z; int i; 
  for (i=0;i<I_NUI-1;i++) z.ui[i]=(x.ui[i]>>1)|(x.ui[i+1]<<(INT_Nbits-1));
  z.ui[I_NUI-1]=x.ui[I_NUI-1]>>1;
  return z;}
int  INCI_lex_GT(INCI *x, INCI *y){
  int i=I_NUI; while(i--) if(x->ui[i]>y->ui[i]) return 1; 
  else if(x->ui[i]<y->ui[i]) return 0; return 0; }
int  INCI_LmR(INCI *x, INCI *y){ puts("Implement INCI_LmR"); exit(0); }
#else
int  INCI_lex_GT(INCI *x, INCI *y){ return (*x > *y) ? 1 : 0 ; }
int  INCI_LmR(INCI *x, INCI *y){ return (*x>*y) ? 1 : (*x<*y) ? -1 : 0; }
/* int  Lead_Vert(INCI x){int i=0; while(!(x%2)) {i++; x/=2;} return i;} */
#endif

int INCI_abs(INCI X){
  int abs=0; while(!INCI_EQ_0(X)) {abs+=INCI_M2(X); X=INCI_D2(X);} return abs;
}

INCI Eq_To_INCI(Equation *_Eq, PolyPointList *_P, VertexNumList *_V){
  int j;  INCI X=INCI_0();
  for (j=0;j<_V->nv;j++) X=INCI_PN(X,Eval_Eq_on_V(_Eq,_P->x[_V->v[j]],_P->n));
  return X;
}

int Print_INCI(INCI X) { 
  int i=0;
  while(!INCI_EQ_0(X)) {printf("%d", (int) INCI_M2(X)); X=INCI_D2(X); i++;} 
  return i; /* number of printed bits */
}

void Print_FaceInfo(int M,FaceInfo *_I){
  int i,j, k, l;
  M--;
  printf("Incidences as binary numbers [F-vector=(%d",_I->nf[0]);
  for(i=1;i<=M;i++)printf(" %d",_I->nf[i]); puts(")]:");
  puts("v[d][i]: sum_j Incidence(i'th dim-d-face, j-th vertex) x 2^j");
  for(i=0;i<=M;i++)     {
    printf("v[%d]: ",i);
    for(j=0;j<_I->nf[i];j++){
      k=Print_INCI(_I->v[i][j]);
      for (l=k;l<_I->nf[0];l++) printf("0");
      printf(" ");}
    puts("");     }
  puts("f[d][i]: sum_j Incidence(i'th dim-d-face, j-th facet) x 2^j");
  for(i=0;i<=M;i++)     {
    printf("f[%d]: ",i);
    for(j=0;j<_I->nf[i];j++){
      k=Print_INCI(_I->f[i][j]);
      for (l=k;l<_I->nf[M];l++) printf("0");
      printf(" ");}
    puts("");     }
}

void Make_CD2Faces(PolyPointList *_P, VertexNumList *_V, EqList *_E,
		    FaceInfo *_I);



int  IsGoodCEq(Equation *_E, PolyPointList *_P, VertexNumList *_V){
  int i=_V->nv; 
  Long s;
  while(!(s=Eval_Eq_on_V(_E, _P->x[_V->v[--i]], _P->n))); 
  if(s < 0) { int j=_P->n; while(j--) _E->a[j]=-_E->a[j]; _E->c=-_E->c; }
  while(i) if(Eval_Eq_on_V(_E, _P->x[_V->v[--i]], _P->n) < 0) return 0;
  return 1;
}

void Make_New_CEqs(PolyPointList *_P, VertexNumList *_V, CEqList *_C, 
			EqList *_F, INCI *CEq_I, INCI *F_I){
  int i,j, Old_C_ne=_C->ne;
  static CEqList Bad_C;
  static INCI Bad_C_I[CEQ_Nmax];

#if     (SHOW_NEW_CEq)
  static int init; static clock_t CLOCK1; static time_t DATE1;
  if(!init){init=1; CLOCK1=clock(); DATE1=time(NULL);}
  printf("V=%d F=%d: Ceq=%d",_V->nv,_F->ne,_C->ne);fflush(stdout);
#endif

  Bad_C.ne=_C->ne=0;
  for (i=0;i<Old_C_ne;i++){
    Long dist = Eval_Eq_on_V(&_C->e[i],_P->x[_V->v[_V->nv-1]],_P->n);
    CEq_I[i]=INCI_PN(CEq_I[i],dist);
    if (dist<0) {Bad_C.e[Bad_C.ne]=_C->e[i]; Bad_C_I[Bad_C.ne++]=CEq_I[i];}
    else {_C->e[_C->ne]=_C->e[i]; CEq_I[_C->ne++]=CEq_I[i];}}
    
#if     (SHOW_NEW_CEq)
  printf("=%dg+%db",_C->ne,Bad_C.ne);fflush(stdout);
#endif

  Old_C_ne=_C->ne;
  for (i=0;i<_F->ne;i++) F_I[i]=
	INCI_PN(F_I[i],Eval_Eq_on_V(&_F->e[i],_P->x[_V->v[_V->nv-1]],_P->n));
  for (j=0;j<_F->ne;j++) if (!INCI_M2(F_I[j]))
    for (i=0;i<Bad_C.ne;i++){
      INCI New_Face=INCI_AND(Bad_C_I[i],F_I[j]);
      int k;
      if (INCI_abs(New_Face)<_P->n-1) continue;
      for (k=0;k<Bad_C.ne;k++) if (INCI_LE(New_Face,Bad_C_I[k])) 
	if (k!=i) break;
      if (k!=Bad_C.ne) continue;
      for (k=0;k<Old_C_ne;k++) if (INCI_LE(New_Face,CEq_I[k])) break;
      if (k!=Old_C_ne) continue;
      for (k=0;k<_F->ne;k++) if (INCI_LE(New_Face,F_I[k])) if (k!=j) break;
      if (k!=_F->ne) continue; 
      assert(_C->ne<CEQ_Nmax);
      CEq_I[_C->ne]=INCI_PN(INCI_D2(New_Face),0);
      _C->e[_C->ne]=EEV_To_Equation(&(Bad_C.e[i]),&(_F->e[j]),
				    _P->x[_V->v[_V->nv-1]],_P->n);
      assert(IsGoodCEq(&(_C->e[_C->ne++]),_P,_V));}
  for (j=0;j<Old_C_ne;j++) if (!INCI_M2(CEq_I[j])) 
    for (i=Bad_C.ne-1;i>=0;i--){
      INCI New_Face=INCI_AND(Bad_C_I[i],CEq_I[j]);
      int k;
      if (INCI_abs(New_Face)<_P->n-1) continue;
      for (k=0;k<Bad_C.ne;k++) if (INCI_LE(New_Face,Bad_C_I[k])) 
	if (k!=i) break;
      if (k!=Bad_C.ne) continue;
      for (k=0;k<Old_C_ne;k++) if (INCI_LE(New_Face,CEq_I[k]))
	if (k!=j) break;
      if (k!=Old_C_ne) continue;
      for (k=0;k<_F->ne;k++) if (INCI_LE(New_Face,F_I[k])) break;
      if (k!=_F->ne) continue;
      assert(_C->ne<CEQ_Nmax);
      CEq_I[_C->ne]=INCI_PN(INCI_D2(New_Face),0);
      _C->e[_C->ne]=EEV_To_Equation(&(Bad_C.e[i]),&(_C->e[j]),
				    _P->x[_V->v[_V->nv-1]],_P->n);
      assert(IsGoodCEq(&(_C->e[_C->ne++]),_P,_V));}

#if     (SHOW_NEW_CEq)
  {time_t DATE2=time(NULL); char sm[2]={'s',0}; 
  int Rs= (int)difftime(DATE2,DATE1);if(Rs>999){Rs/=60; *sm='m';}
  printf(" done: C.ne=%d  %d%s\n",_C->ne,Rs,sm);fflush(0);}
#endif
}

int  Search_New_Vertex(Equation *_E, PolyPointList *_P){
  int i, v=0; 
  Long *X=_P->x[0], x=Eval_Eq_on_V(_E,X,(_P->n));
  for(i=1;i<_P->np;i++)      {  
    Long *Y=_P->x[i], y=Eval_Eq_on_V(_E,Y,(_P->n));
    if(y>x) continue;
    if(y==x) if(Vec_Greater_Than(X,Y,_P->n)) continue;
    v=i; X=Y; x=y;     }
  return v;
}


Long VZ_to_Base(Long *V,int *d,Long M[POLY_Dmax][POLY_Dmax])  /* 0 iff V=0 */
{    int p[POLY_Dmax], i, j, J=0; Long g=0, W[POLY_Dmax], *G[POLY_Dmax]; 
     for(i=0;i<*d;i++) 	if(V[i]) {W[J]=V[i]; G[J]=M[i]; p[J++]=i;}
			else for(j=0;j<*d;j++) M[i][j]=(i==j);
     if(J) if(p[0]) { G[0]=M[0]; for(j=0;j<*d;j++) M[p[0]][j]=(j==0);}
     if(J>1) g=W_to_GLZ(W,&J,G); else if(J){g=*W; M[0][0]=0; M[0][p[0]]=1;}
     if(J>1)
     {  for(i=0;i<J;i++) { int I=J; 
	for(j=*d-1;j>=0;j--) G[i][j] = (V[j]) ? G[i][--I] : 0; assert(I==0);}
     }	return g;
}

int  OrthBase_red_by_V(Long *V, int *d, Long A[][POLY_Dmax], int *r,
	Long B[][POLY_Dmax])
{    int i, j, k; Long W[POLY_Dmax], G[POLY_Dmax][POLY_Dmax];
     for(i=0;i<*r;i++) {int j; W[i]=0; for(j=0;j<*d;j++) W[i]+=A[i][j]*V[j];}
     assert( VZ_to_Base(W,r,G) );
     for(i=0;i<*r-1;i++) for(k=0;k<*d;k++)
     {	B[i][k]=0; for(j=0;j<*r;j++) B[i][k]+=G[i+1][j]*A[j][k];
     }
#if	(TEST_GLZ_EQ)
	printf("A -> B ... V = "); for(k=0;k<*d;k++) printf(" %5d",V[k]);
	printf("  W=");for(k=0;k<*r;k++)printf(" %5d",W[k]);puts(""); {int 
	a,b; for(a=0;a<*r-1;a++){for(b=0;b<*d;b++)printf(" %5d",A[a][b]);
	printf("  =A  B=  ");for(b=0;b<*d;b++)printf(" %5d",B[a][b]);
	puts("");}for(b=0;b<*d;b++)printf(" %5d",A[a][b]);printf("  =A\n");}
#endif
     	return (*r)--;
}

int  New_Start_Vertex(Long *V0,Long *Ea, PolyPointList *P,int *v) /* P.x[v] */
{    Equation E; int i, n=0, p=0; Long d, dn=0, dp=0, *Xn=P->x[0], *Xp=Xn; 
     for(i=0;i<P->n;i++) E.a[i]=Ea[i]; E.c=0; E.c=-Eval_Eq_on_V(&E,V0,P->n);
     d=Eval_Eq_on_V(&E,P->x[0],P->n); if(d>0) dp=d; if(d<0) dn=d;
     for(i=1;i<P->np;i++)
     {	d=Eval_Eq_on_V(&E,P->x[i],P->n); if(d==0) continue; 
	if(d==dp) if(Vec_Greater_Than(P->x[i],Xp,P->n)) Xp=P->x[p=i];
        if(d>dp)  {dp=d; Xp=P->x[p=i];}
	if(d==dn) if(Vec_Greater_Than(P->x[i],Xn,P->n)) Xn=P->x[n=i];
        if(d<dn)  {dn=d; Xn=P->x[n=i];}
     }
     if(dp) if(dn) 				 /* points on both sides */
#if	(VERT_WITH_MAX_DISTANCE)
			{if(dp+dn>0) *v=p; else *v=n;}
#else
			{if(dp+dn>0) *v=n; else *v=p;}
#endif
	  else   *v=p;				/* d >=0 */
     else if(dn) *v=n;				/* d <=0 */
          else return 0;
/*	for(i=0;i<P->n;i++) printf(" %d ",Xp[i]); printf(" = Xp  Xn =");
	for(i=0;i<P->n;i++) printf(" %d ",Xn[i]); printf(" \n");
*/
     return 1;
}

int  GLZ_Start_Simplex(PolyPointList *_P, VertexNumList *_V, CEqList *_C)
{    int i, x=0, y=0, *VN=_V->v, *d=&_P->n, r=*d, b[POLY_Dmax]; 
     Long *X=_P->x[x], *Y=_P->x[y], XX=0, YY=0, 
	B[(POLY_Dmax*(POLY_Dmax+1))/2][POLY_Dmax], W[POLY_Dmax];  if(_P->np<2) 
     {	for(x=0;x<_P->n;x++) for(y=0;y<_P->n;y++) _C->e[x].a[y]=(x==y);
	assert(_P->np>0); for(x=0;x<_P->n;x++) _C->e[x].c=-_P->x[0][x];
	return _C->ne=_P->n;
     }
     for(i=1; i<_P->np; i++) 
     {	Long *Z=_P->x[i]; 			    
	if(Vec_Greater_Than(X,Z,_P->n)) X=_P->x[x=i];	/* (x_n)-max: VN[0] */
	if(Vec_Greater_Than(Z,Y,_P->n)) Y=_P->x[y=i];	/* (x_n)-min: VN[1] */
     }	assert(x!=y);	     /* at this point I need two different vertices */
     for(i=0;i<*d;i++) { Long Xi=(X[i]>0) ? X[i]: -X[i],
	Yi=(Y[i]>0) ? Y[i]: -Y[i]; if(Xi>XX)XX=Xi; if(Yi>YY)YY=Yi;}
     if(YY<XX) {VN[0]=y;VN[1]=x;} else {VN[0]=x;VN[1]=y;} _V->nv=2; y=VN[1];
     X=_P->x[VN[0]];Y=_P->x[VN[1]]; for(i=0;i<*d;i++) b[i]=(i*(2*(*d)-i+1))/2;
     for(x=0;x<*d;x++) for(i=0;i<*d;i++) B[x][i]=(x==i); /* b[i+1]-b[i]=d-i */
     for(x=1;x<*d;x++)
     {	for(i=0;i<*d;i++) W[i]=_P->x[y][i]-X[i];
	OrthBase_red_by_V(W,d,&B[b[x-1]],&r,&B[b[x]]); for(i=0;i<r;i++) 
#if	(LONG_EQ_FIRST)
	if(New_Start_Vertex(X,B[b[x]+r-i-1],_P,&y)) break;
#else
	if(New_Start_Vertex(X,B[b[x]+i],_P,&y)) break;
#endif
	if(i==r) break;	_V->v[_V->nv++]=y;	       /* x = dim(span) < d */
     }
     if(x<*d)
     {  for(y=0;y<r;y++)
	{   Equation *E=&_C->e[y]; Long *Z=B[b[x]+y]; 
	    E->c=0; for(i=0;i<*d;i++) E->a[i]=Z[i]; 
	    E->c=-Eval_Eq_on_V(E,X,_P->n); 
	}   return _C->ne=r;
     }
     else
     {	Equation *E=_C->e; Long *Z=B[b[*d-1]]; E->c=0; _C->ne=2;
	for(i=0;i<*d;i++)E->a[i]=Z[i]; E->c=-Eval_Eq_on_V(E,X,_P->n); 
	if(Eval_Eq_on_V(E,_P->x[_V->v[*d]],_P->n)<0) {
		for(i=0;i<*d;i++)E->a[i]=-Z[i]; E->c*=-1;}
	X=_P->x[_V->v[r=*d]]; 
        for(x=1;x<*d;x++)			/* now the 2nd equation */
	{   Y=_P->x[_V->v[x-1]]; for(i=0;i<*d;i++) W[i]=X[i]-Y[i];
	    OrthBase_red_by_V(W,d,&B[b[x-1]],&r,&B[b[x]]);
	}
	E=&_C->e[1]; E->c=0;
	for(i=0;i<*d;i++)E->a[i]=Z[i]; E->c=-Eval_Eq_on_V(E,X,_P->n); 
	assert(XX=Eval_Eq_on_V(E,_P->x[_V->v[*d-1]],_P->n)); 
	if(XX<0) {for(i=0;i<*d;i++)E->a[i]=-Z[i]; E->c*=-1;}
        for(x=*d-2;x>=0;x--)			/* omit vertex #x */
	{   r=*d-x; for(y=x+1;y<*d;y++)
	    {	Y=_P->x[_V->v[y]]; for(i=0;i<*d;i++) W[i]=X[i]-Y[i];
	    	OrthBase_red_by_V(W,d,&B[b[y-1]],&r,&B[b[y]]);
	    }
	E=&_C->e[(_C->ne)++]; E->c=0;
	for(i=0;i<*d;i++)E->a[i]=Z[i]; E->c=-Eval_Eq_on_V(E,X,_P->n); 
	assert(XX=Eval_Eq_on_V(E,_P->x[_V->v[x]],_P->n)); 
	if(XX<0) {for(i=0;i<*d;i++)E->a[i]=-Z[i]; E->c*=-1;}
	}
     } 
     assert(*d+1==_C->ne); for(x=0;x<_C->ne;x++) for(i=0;i<=*d;i++)
     assert((x==i)==(0!=Eval_Eq_on_V(&_C->e[x],_P->x[_V->v[*d-i]],_P->n)));
     return 0;
}

int  REF_Search_Bad_Eq(CEqList *_C, EqList *_F, INCI *CEq_I, INCI *F_I,
	       PolyPointList *_P, int *_REF){   /* return 0 :: no bad eq. */
  while(_C->ne--)  {	
    int j; 
    for(j=0;j<_P->np;j++)			
      if(Eval_Eq_on_V(&(_C->e[_C->ne]),_P->x[j],_P->n) < 0) return ++_C->ne;
    if(_C->e[_C->ne].c != 1) { *_REF=0; return 1;}
    assert(_F->ne<EQUA_Nmax);
    _F->e[_F->ne]=_C->e[_C->ne];
    F_I[_F->ne++]=CEq_I[_C->ne];}
  return 0;
}

int  Finish_REF_Check(PolyPointList *_P, VertexNumList *_V, EqList *_F,
		     CEqList *_CEq, INCI *F_I, INCI *CEq_I){
  int REF=1;
  while(0<=_CEq->ne) if(REF_Search_Bad_Eq(_CEq,_F,CEq_I,F_I,_P,&REF)){
    if(!REF) return 0;	/* found d!=1 */
    assert(_V->nv<VERT_Nmax);
    _V->v[_V->nv++]=Search_New_Vertex(&(_CEq->e[_CEq->ne-1]),_P);
    Make_New_CEqs(_P,_V,_CEq,_F,CEq_I,F_I); }
  return 1;					   
}

int  Ref_Check(PolyPointList *_P, VertexNumList *_V, EqList *_F){
  int i; 
  CEqList *CEq = (CEqList *) malloc(sizeof(CEqList)); 
  INCI *CEq_I = (INCI *) malloc(sizeof(INCI)*CEQ_Nmax);
  INCI *F_I = (INCI *) malloc(sizeof(INCI)*EQUA_Nmax);  
  if((CEq==NULL)||(CEq_I==NULL)||(F_I==NULL)) {
    printf("Allocation failure in Ref_Check\n"); exit(0);}
  if (GLZ_Start_Simplex(_P, _V, CEq)) {
    free(CEq); free(CEq_I); free(F_I); return 0;}
  for (i=0;i<CEq->ne;i++) CEq_I[i]=Eq_To_INCI(&(CEq->e[i]),_P,_V);
  _F->ne=0;
  i=Finish_REF_Check(_P, _V, _F, CEq, F_I, CEq_I);
  free(CEq); free(CEq_I); free(F_I);
  return i;
}

int  FE_Search_Bad_Eq(CEqList *_C, EqList *_F, INCI *CEq_I, INCI *F_I,
	       PolyPointList *_P, int *_IP){   /* return 0 :: no bad eq. */
  while(_C->ne--)  {	
    int j; 
    for(j=0;j<_P->np;j++)			
      if(Eval_Eq_on_V(&(_C->e[_C->ne]),_P->x[j],_P->n) < 0) return ++_C->ne;
    if(_C->e[_C->ne].c < 1) *_IP=0;
    assert(_F->ne<EQUA_Nmax);
    _F->e[_F->ne]=_C->e[_C->ne];
    F_I[_F->ne++]=CEq_I[_C->ne];}
  return 0;
}

int  Finish_Find_Equations(PolyPointList *_P, VertexNumList *_V, 
		     EqList *_F, CEqList *_CEq, INCI *F_I, INCI *CEq_I){
  int IP=1;
  while(0<=_CEq->ne) if (FE_Search_Bad_Eq(_CEq,_F,CEq_I,F_I,_P,&IP)){
    assert(_V->nv<VERT_Nmax);
    _V->v[_V->nv++]=Search_New_Vertex(&(_CEq->e[_CEq->ne-1]),_P);
    Make_New_CEqs(_P,_V,_CEq,_F,CEq_I,F_I); }
  return IP;					   
}

int  Find_Equations(PolyPointList *_P, VertexNumList *_V, EqList *_F){
  /* return: IP, finds Vertices and Equations for _P even if not IP */
  int i; 
  CEqList *CEq = (CEqList *) malloc(sizeof(CEqList)); 
  INCI *CEq_I = (INCI *) malloc(sizeof(INCI)*CEQ_Nmax);
  INCI *F_I = (INCI *) malloc(sizeof(INCI)*EQUA_Nmax); 
  CEq->ne=0;
  if((CEq==NULL)||(CEq_I==NULL)||(F_I==NULL)) {
     printf("Allocation failure in Find_Equations\n"); exit(0);}
  if (GLZ_Start_Simplex(_P, _V, CEq)) {
    _F->ne=CEq->ne; 
    for(i=0;i<_F->ne;i++) _F->e[i]=CEq->e[i]; 
    free(CEq); free(CEq_I); free(F_I); 
    return 0;}
  _F->ne=0;
  for (i=0;i<CEq->ne;i++) 
    if(INCI_abs(CEq_I[i]=Eq_To_INCI(&(CEq->e[i]),_P,_V))<_P->n)
      {fprintf(outFILE,"Bad CEq in Find_Equations"); exit(0);}
  i=Finish_Find_Equations(_P, _V, _F, CEq, F_I, CEq_I);
  free(CEq); free(CEq_I); free(F_I);
  return i;
}

/*  ======================================================================  */
/*  ==========		     			  		==========  */
/*  ==========	  D U A L   P O L Y   &   C O M P L E T I O N 	==========  */
/*  ==========							==========  */ 
/*  ======================================================================  */

typedef	struct {LLong N; LLong D;} 		             LRat;  /* = N/D */

LRat  LrI(LLong a) 			       		      /*  a -> a/1  */
{    LRat c; c.N=a; c.D=1; return c; 
}

LRat  LrP(LRat a, LRat b)                                            /*  a * b  */
{    register LLong g=LFgcd(a.N,b.D); register LLong h=LFgcd(b.N,a.D);
     LRat c; c.N=(a.N/g)*(b.N/h); 
     if((c.D=(a.D/h)*(b.D/g)) < 0) {c.D=-c.D; c.N=-c.N;} return c;
}
LRat  LrQ(LRat a, LRat b)                                            /*  a / b  */
{    register LLong g=LNNgcd(a.N,b.N); register LLong h=LFgcd(b.D,a.D);
     LRat c; c.N=(a.N/g)*(b.D/h);
     if((c.D=(a.D/h)*(b.N/g)) < 0) {c.N=-c.N; c.D=-c.D;} return c;
}

LRat  LrD(LRat a, LRat b)                                            /*  a - b  */
{    LRat c; register LLong g=LFgcd(a.D,b.D); 
     g = LFgcd(c.N=a.N*(b.D/g)-b.N*(a.D/g), c.D=a.D*(b.D/g));
     if(g < 0) g=-g; c.N/=g; c.D/=g; /* LLong in line above ... */
#ifdef	TEST
     if(c.D<=0) {LRpr(c); puts(" *** c.D<=0 in rD! ***"); exit(0);}
#endif
     return c;
}

int Vec_Equal(Long *X, Long *Y, int i){	    /* return 1 iff `X == Y' */
  while(i--) if(X[i]!=Y[i]) return 0; return 1;
}

void Make_Dual_Poly(PolyPointList *_P, VertexNumList *_V, EqList *_E,
                  PolyPointList *_DP){
  EqList DE;
  PairMat PM, DPM;
  Make_VEPM(_P,_V,_E,PM);
  Transpose_PM(PM, DPM, _V->nv, _E->ne);
  VNL_to_DEL(_P,_V,&DE);
  _DP->n=_P->n;
  _DP->np=0;
  Complete_Poly(DPM, &DE, _E->ne, _DP);
}

void add_for_completion(Long *yDen, Long Den,
    EqList *_E, PolyPointList *_CP, int *old_np){
  int i,n=_CP->n;
  Long yold[POLY_Dmax];

  if(Den>1) for(i=0;i<n;i++) {
    if(yDen[i]%Den) return;
    yold[i]=yDen[i]/Den;}
  else for(i=0;i<n;i++) yold[i]=yDen[i];
  for (i=0;i<_E->ne;i++) if (Eval_Eq_on_V(&(_E->e[i]), yold, n) < 0) return;
  for (i=0;i<*old_np;i++) if (Vec_Equal(_CP->x[i],yold,n)) return;
  assert(_CP->np<POINT_Nmax);
  for(i=0;i<n;i++) _CP->x[_CP->np][i]=yold[i];
  _CP->np++;
}

void Complete_Poly(PairMat VPM, EqList *_E, int nv, 
			     PolyPointList *_CP){
  int i,j,k,l,InsPoint,rank=0,n=_CP->n,old_np=_CP->np;
  Long MaxDist[EQUA_Nmax], InvMat[POLY_Dmax][POLY_Dmax], Den=1;
  Long yDen[POLY_Dmax];
  int OrdFac[VERT_Nmax], 
    BasFac[POLY_Dmax], one[POLY_Dmax], position[POLY_Dmax];
  LRat ind[POLY_Dmax][POLY_Dmax], x[POLY_Dmax], y[POLY_Dmax], f, 
    PInvMat[POLY_Dmax][POLY_Dmax];
  /*_CP->np=0;*/

  /* Calculate maximal distances from facets of Delta^* (Vertices of Delta) */

  for (i=0;i<_E->ne;i++) {  
    MaxDist[i]=0;
    for (j=0;j<nv;j++) 
    if (MaxDist[i]<VPM[i][j]) MaxDist[i]=VPM[i][j];}
	
  /* Order facets of Delta^* (Vertices of Delta) w.r.t. MaxDist   */

  OrdFac[0]=0;
  for (i=1;i<_E->ne;i++){
    InsPoint=i; 
    while (InsPoint&&(MaxDist[i]<MaxDist[OrdFac[InsPoint-1]])) InsPoint--;
    for (j=i;j>InsPoint;j--) OrdFac[j]=OrdFac[j-1];
    OrdFac[InsPoint]=i; }
	
  /* Find first POLY_Dmax linearly independent facets + Inverse Matrix */

  for (i=0;i<n;i++) for (j=0;j<n;j++) PInvMat[i][j]=LrI(0);
  for (i=0;i<n;i++) PInvMat[i][i]=LrI(1);
  i=0;
  while (rank<n){
    for (j=0;j<n;j++) x[j]=LrI(_E->e[OrdFac[i]].a[j]);
    for (j=0;j<n;j++) y[j]=LrI(0);
    y[rank]=LrI(1);
    for (j=0;j<rank;j++) {
      f=x[one[j]];
      for (k=0;k<n;k++) {
        x[k]=LrD(x[k],LrP(f,ind[j][k])); 
        y[k]=LrD(y[k],LrP(f,PInvMat[j][k]));  } }
    one[rank]=-1;
    for (l=0;(l<n)&&(one[rank]==-1);l++) if (x[l].N) one[rank]=l;
    if(one[rank]>-1){
      for (k=0;k<n;k++) {
        ind[rank][k]=LrQ(x[k],x[one[rank]]);
        PInvMat[rank][k]=LrQ(y[k],x[one[rank]]); }
      for (j=0;j<rank;j++) {
        f=ind[j][one[rank]];
        for (k=0;k<n;k++)         {
          ind[j][k]=LrD(ind[j][k],LrP(ind[rank][k],f));   
          PInvMat[j][k]=LrD(PInvMat[j][k],LrP(PInvMat[rank][k],f));  }     }
      BasFac[rank]=OrdFac[i];
      rank++; }  
    i++; }
  for (i=0;i<n;i++) for (j=0;j<n;j++) 
    Den=(Den/LFgcd(Den,PInvMat[i][j].D))*PInvMat[i][j].D;
  for (i=0;i<n;i++) for (j=0;j<n;j++) 
    InvMat[one[i]][j]=(Den/PInvMat[i][j].D)*PInvMat[i][j].N;

  for (i=0;i<n;i++){
    for (j=0;j<n;j++) {
      long long s=0;
      for(k=0;k<n;k++) s+=((long long) (InvMat[k][i]))*
			 ((long long) (_E->e[BasFac[j]].a[k]));
      if (s!=Den*(i==j)) {
	puts("something wrong in Make_Dual_Poly");
	exit(0);}}} 

  /* Examine all integer points of parallelogram:                         */
  /* The basic structure of the algorithm is:
  for (k=0;k<n-1;k++) position[k]=-1;      / * sets k=n-1; important!      *
  position[n-1]=-2;  / * starting point just outside the parallelogram     *
  while(k>=0){
    position[k]++;
    DO AT position;
    for(k=n-1;((position[k]==MaxDist[BasFac[k]]-1)&&(k>=0));k--) 
       position[k]=-1;  }
         / * sets k to the highest value where pos.[k] wasn't the max value; 
            resets the following max values to min values                 */
  /* Quantities linear in position can be changed with every change of
     position (here: yDen)                                                */

  for(i=0;i<n;i++) yDen[i]=0;
  for (k=0;k<n-1;k++) {   /* sets k=n-1; important!   */
    position[k]=-_E->e[BasFac[k]].c;   
    for(i=0;i<n;i++) yDen[i]-=_E->e[BasFac[k]].c*InvMat[i][k]; }
  position[n-1]=-_E->e[BasFac[n-1]].c-1;
  for(i=0;i<n;i++) yDen[i]-=(_E->e[BasFac[k]].c+1)*InvMat[i][n-1];
  while(k>=0){
    position[k]++;
    for(i=0;i<n;i++) yDen[i]+=InvMat[i][k];
    add_for_completion(yDen, Den, _E, _CP, &old_np);
    for(k=n-1;(k>=0);k--){
      if (position[k]!=MaxDist[BasFac[k]]-_E->e[BasFac[k]].c) break;
      position[k]=-_E->e[BasFac[k]].c;
      for (i=0;i<n;i++) yDen[i]-=MaxDist[BasFac[k]]*InvMat[i][k]; }}
}

void Print_PPL(PolyPointList *_P, const char *comment){
  int i,j;
  if(_P->np>20){
    fprintf(outFILE,"%d %d  %s\n",_P->np,_P->n,comment);
    for(i=0;i<_P->np;i++) {
      for(j=0;j<_P->n;j++) fprintf(outFILE,"%d ",(int) _P->x[i][j]); 
      fprintf(outFILE,"\n");}}
  else {
    fprintf(outFILE,"%d %d  %s\n",_P->n,_P->np,comment);
    for(i=0;i<_P->n;i++) {
      for(j=0;j<_P->np;j++) fprintf(outFILE," %4d",(int) _P->x[j][i]); 
      fprintf(outFILE,"\n");}}
}

void Make_Incidence(PolyPointList *_P, VertexNumList *_V, EqList *_E,
		    FaceInfo *_I)
/*    The incidence relations for faces are stored on the structure FaceInfo:
 *
 *    	int  nf[d]      ==  #faces(dim.=d)  ==  #dual faces[dim.=n-d-1]
 *    	INCI v[d][i]    ::  vertices on i-th dim=d face
 *    	INCI f[d][i]    ::  dual vertices on face dual to i-th dim=n-d-1 face
 *    	Long nip[d][i]  ::  #(IPs of i-th dim=d face)
 *    	Long dip[d][i]  ::  #(IPs of i-th dim=n-d-1 face on dual)
 *
 * .v: compute vertices on facets; make intersections `&' of faces with facets
 *     while keeping only maximal intersections;
 * .f: take analogous union; if same intersection again make union `|' of f's
 */
{    int i, j, M=_P->n-1, d=M, D;  
     assert(_E->ne<=VERT_Nmax);
     _I->nf[M]=_E->ne; 
     for(i=0;i<_E->ne;i++)_I->v[M][i]=Eq_To_INCI(&_E->e[i],_P,_V); /* init .v*/
     assert(i>0);
     _I->f[M][--i]=INCI_1();
     while(i--) _I->f[M][i]=INCI_PN(_I->f[M][i+1],1);    	  /* init .f */
     while((D=d--))
     {	int *n=&_I->nf[d];  *n=0;			       /* #(d-faces) */
	for(i=0; i<_I->nf[D]; i++) for(j=0; j<_I->nf[M]; j++)
	{   int k; INCI x=INCI_AND(_I->v[D][i],_I->v[M][j]); /* x=candidate */
	    INCI u=INCI_OR(_I->f[D][i],_I->f[M][j]);
	    if( (!INCI_EQ(x,_I->v[D][i]))&&(INCI_abs(x)>d) )/* x!=vD & |x|>d */
	    {	for(k=0;k<*n;k++)
	    	{   INCI *y=&_I->v[d][k],*v=&_I->f[d][k]; /* (*y)==v[d][k] */
		    if(INCI_LE(*y,x))
		    {	if(INCI_EQ(x,*y)) 
		    	{   *v=INCI_OR(*v,u); break;     /* x=y :: .f&=... */
		    	} 
		    	else {int l=k; *y=x; _I->f[d][k]=u;/* x>y :: y=x;f= */
	while((++l) < (*n))
if(INCI_LE(_I->v[d][l],x)) {(*n)--; 
	if(l<(*n)){_I->v[d][l]=_I->v[d][*n]; _I->f[d][l]=_I->f[d][*n];}}
else assert(!INCI_LE(x,_I->v[d][l]));

/* for(k++;k<*n;k++) if(!INCI_LE(x,_I->v[d][k])&&!INCI_LE(_I->v[d][k],x))
   {Print_PPL(_P,"FACE_Info trouble");assert(0);} */
			    break;}
		    }
		    else if(INCI_LE(x,*y))   break;	   /* x<y :: break */
	        }
	        if(*n==k) 	      	  /* non-comparable => new face */
	        {   assert(k<FACE_Nmax); _I->v[d][k]=x; _I->f[d][k]=u; (*n)++;
	        }
	    }
	}
     }	d=_P->n; M=0; for(i=0;i<d;i++) M += _I->nf[i] * (1-2*(i%2));
     if(M!=2*(d%2)){for(i=0;i<d;i++)printf("%3d",_I->nf[i]);puts("=F-vector");
     	Print_PPL(_P,"PPL for incidence error");Print_FaceInfo(_P->n,_I);
	printf("d=%d  euler=%d\n",d,M); assert(M==2*(d%2));}
}

/* -------- Functions added from Polynf.c ---------*/

#define SL_Long		LLong		    /* has same problems as REgcd */

#if	(POLY_Dmax < 5)

#define NFX_Limit       903		/* 138b->255  153e->279  165c->327 */
#define X_Limit         9999		   /* 178c->375  218->399,462,483 */
#define VPM_Limit	9999

#else

#define NFX_Limit  1631721   /* 1631721 1 903 37947 233103 543907 815860    */
#define X_Limit    3263441   /* 3263442 1 1806 75894 466206 1087814 1631721 */
#define VPM_Limit  3263442   /* 1631721 1 903 37947 233103 543907 815860    */

#endif

typedef struct {int C[VERT_Nmax], L[VERT_Nmax], s;}             PERM;
typedef struct {int nv, nf, ns;}  				vNF;

#define	Fputs(S)	{fputs(S,outFILE);fputs("\n",outFILE);}

/*   ------	forward declarations	------ */

void NF_Coordinates(PolyPointList *_P, VertexNumList *_V, EqList *_F);

int  GLZ_Make_Trian_NF(Long X[][VERT_Nmax], int *n, int *nv,
		       GL_Long G[POLY_Dmax][POLY_Dmax]);    /* current=best */

int  SL2Z_Make_Poly_NF(Long X[][VERT_Nmax], int *n, int *nv,
		       SL_Long S[POLY_Dmax][POLY_Dmax]);    /* previous=bad */

int  Init_rVM_VPM(PolyPointList *P, VertexNumList *_V,EqList *_F,/* in */
	    	int *d,int *v,int *f, Long VM[POLY_Dmax][VERT_Nmax], /* out */
	    	Long VPM[VERT_Nmax][VERT_Nmax]);	/* return reflexive */

void Eval_Poly_NF(int *d,int *v,int *f, Long VM[POLY_Dmax][VERT_Nmax],
		Long VPM[VERT_Nmax][VERT_Nmax],			      /* in */
		Long pNF[POLY_Dmax][VERT_Nmax], int t);		     /* out */

void Make_VPM_NF(int *v, int *f, Long VPM[VERT_Nmax][VERT_Nmax],      /* in */
		PERM *CL,int *ns,Long VPM_NF[VERT_Nmax][VERT_Nmax]); /* out */

void Aux_pNF_from_vNF(PERM *CL,int *ns,int *v,int *d,
		Long VM[POLY_Dmax][VERT_Nmax],			      /* in */
		Long pNF[POLY_Dmax][VERT_Nmax],int *t);		     /* out */

void New_pNF_Order(int *v,int *f,PERM *CL,int *ns,Long VPM_NF[][VERT_Nmax]);

int  Make_Poly_NF(PolyPointList *_P, VertexNumList *_V, EqList *_F,
		Long pNF[POLY_Dmax][VERT_Nmax]);	  /* 1 if reflexive */

/*   ------	definitions	------ */

void New_pNF_Order(int *v,int *f,PERM *CL,int *ns,Long VPM_NF[][VERT_Nmax])
{    int i, j, pi[VERT_Nmax], c[VERT_Nmax]; 
     Long maxP[VERT_Nmax], sumP[VERT_Nmax];
     for(i=0;i<*v;i++)
     {	pi[i]=i; maxP[i]=sumP[i]=0; for(j=0;j<*f;j++)
	{   sumP[i]+=VPM_NF[j][i];
	    if(VPM_NF[j][i]>maxP[i])maxP[i]=VPM_NF[j][i];
	}
     }
     for(i=0;i<*v-1;i++)
     {	int n=i; for(j=i+1;j<*v;j++) 
        {   if(maxP[j]<maxP[n]) n=j; else 
	    if(maxP[j]==maxP[n]) if(sumP[j]<sumP[n]) n=j;
	}
        if(n!=i)	
     	{   Long aP=maxP[i]; int a=pi[i]; 
	    maxP[i]=maxP[n]; maxP[n]=aP; pi[i]=pi[n]; pi[n]=a;
	    aP=sumP[i]; sumP[i]=sumP[n]; sumP[n]=aP; 
	}
     }
     for(i=0;i<*ns;i++) 
     {	int *C=CL[i].C; for(j=0;j<*v;j++) c[j]=C[pi[j]]; 
	for(j=0;j<*v;j++) C[j]=c[j];
     }
}

GL_Long GL_Egcd(GL_Long A0, GL_Long A1, GL_Long *Vout0, GL_Long *Vout1)  
{    GL_Long V0=A0, V1=A1, A2, X0=1, X1=0, X2=0;
     while((A2 = A0 % A1)) { X2=X0-X1*(A0/A1); A0=A1; A1=A2; X0=X1; X1=X2; }
     *Vout0=X1, *Vout1=(A1-(V0) * X1)/ (V1); return A1;
}
GL_Long GL_RoundQ(GL_Long N,GL_Long D)
{    GL_Long F; if(D<0) {D=-D; N=-N;} F=N/D; return F+(2*(N-F*D))/D; 
}

GL_Long GL_W_to_GLZ(GL_Long *W, int d, GL_Long **GLZ)		
{    int i, j; GL_Long G, *E=*GLZ, *B=GLZ[1];for(i=0;i<d;i++)assert(W[i]!=0);
     for(i=1;i<d;i++)for(j=0;j<d;j++)GLZ[i][j]=0;
     G=GL_Egcd(W[0],W[1],&E[0],&E[1]); B[0]=-W[1]/G; B[1]=W[0]/G;
     for(i=2;i<d;i++)                   
     {  GL_Long a, b, g=GL_Egcd(G,W[i],&a,&b); B=GLZ[i];
        B[i]= G/g; G=W[i]/g; for(j=0;j<i;j++) B[j]=-E[j]*G;  /* B=Base-Line */
        for(j=0;j<i;j++) E[j]*=a; E[j]=b;                     /* next Egcd */
        for(j=i-1;0<j;j--)                         /* I M P R O V E M E N T */
	{   GL_Long *Y=GLZ[j],rB=GL_RoundQ(B[j],Y[j]),rE=GL_RoundQ(E[j],Y[j]);
            int n; for(n=0;n<=j;n++) { B[n] -= rB*Y[n]; E[n] -= rE*Y[n]; } 
	}   G=g;
     } 
     return G;
}


int  PermChar(int n)
{    if(n<10) return '0'+n; else if(n<36) return 'a'+n-10; else
     if(n<62) return 'A'+n-36; else
     {puts("Printing permutations only for #Vert<=62 !!");exit(0);} return 0;
}

int  Perm_String(int *p,int v,char *s)
{    int i=0; if(v<62) for(i=0;i<v;i++) s[i]=PermChar(p[i]); s[i]=0;return i;
}

void Print_Perm(int *p,int v,const char *s)
{    int i; for(i=0;i<v;i++) fprintf(outFILE,"%c",PermChar(p[i]));
     /*puts("");for(i=48;i<128;i++)printf("%c",i);*/ fprintf(outFILE,"%s",s);
}

int  Aux_XltY_Poly_NF(Long X[][VERT_Nmax],Long Y[][VERT_Nmax], int *n,int *nv)
{    int i, j; Long d;					  /* return "X < Y" */
     for(i=0;i<*n;i++)for(j=0;j<*nv;j++) if((d=X[i][j]-Y[i][j])) 
     {	if(d<0) return 1; else return 0;
     }
     return 0;
}	

int  GLZ_Make_Trian_NF(Long X[][VERT_Nmax], int *n, int *nv,
		       			       GL_Long G[POLY_Dmax][POLY_Dmax])
{    int i, j, C=-1, L; GL_Long g, W[POLY_Dmax], NF[POLY_Dmax][VERT_Nmax], 
	*_G[POLY_Dmax], NG[POLY_Dmax][POLY_Dmax];for(i=0;i<*n;i++)_G[i]=NG[i];
     for(i=0;i<*n;i++)for(j=0;j<*n;j++)G[i][j]=(i==j);		  /* init G */
     for(i=0;i<*n;i++)for(j=0;j<*nv;j++)NF[i][j]=0; 
     for(L=0;L<*n;L++)
     {  int N=0, p[POLY_Dmax];
        while(N==0)               /* find column C with N>0 non-zero entries */
        {   ++C; for(i=0;i<*n;i++) 
	    {	for(j=0;j<*n;j++) NF[i][C]+=G[i][j]*X[j][C];
	    }
	    for(i=L;i<*n;i++) if(NF[i][C]) {W[N]=NF[i][C]; p[N++]=i;}
	}
        assert(N); if(N==1) {g=W[0]; _G[0][0]=1;} else g=GL_W_to_GLZ(W,N,_G);
	if(g<0) { g *= -1; for(i=0;i<N;i++)_G[0][i] *= -1; }
	NF[L][C]=g; for(i=L+1;i<*n;i++) NF[i][C]=0;
        for(i=0;i<*n;i++)
	{   GL_Long Cp[POLY_Dmax]; for(j=0;j<N;j++) Cp[j]=G[p[j]][i];
	    for(j=0;j<N;j++) 
	    {	int k; G[p[j]][i]=0; 
		for(k=0;k<N;k++) G[p[j]][i] += _G[j][k]*Cp[k];
	    }
	}
        if(L!=p[0])for(i=0;i<*n;i++)	     /* swap lines G[L] <-> G[p[0]] */
	{   GL_Long A=G[L][i]; G[L][i]=G[p[0]][i]; G[p[0]][i]=A; 
	}
	for(i=0;i<L;i++)	 /* make upper diag minimal nonneg. */
	{   GL_Long R=NF[i][C]/NF[L][C];
	    if((NF[i][C]-R*NF[L][C])<0) R-=1;
	    NF[i][C]-=R*NF[L][C]; for(j=0;j<*n;j++)G[i][j]-=R*G[L][j];
	}
     }
     while(++C<*nv)for(i=0;i<*n;i++)for(j=0;j<*n;j++)NF[i][C]+=G[i][j]*X[j][C];
     for(i=0;i<*n;i++)for(j=0;j<*nv;j++) { 
#ifdef	SHOW_NFX_LIMIT
	g=NF[i][j]; if(g<0) g=-g; if(g>NFX_Limit) { fprintf(stderr,
	    "NFX_Limit in GL -> %lld !!\n",(long long) g); return 0; } else 
#endif
	X[i][j]=NF[i][j]; }
     return 1;
}	

int  Aux_Make_Poly_NF(Long X[][VERT_Nmax], int *n, int *nv)
{    GL_Long G[POLY_Dmax][POLY_Dmax];
#if	(TEST_GLZ_VS_SL)
     int i,j,x; SL_Long S[POLY_Dmax][POLY_Dmax];Long XS[POLY_Dmax][VERT_Nmax];
     for(i=0;i<*n;i++)for(j=0;j<*nv;j++)XS[i][j]=X[i][j];
     x=GLZ_Make_Trian_NF(X,n,nv,G); SL2Z_Make_Poly_NF(XS,n,nv,S);
     for(i=0;i<*n;i++)for(j=0;j<*n;j++)  assert( S[i][j]==G[i][j]);}
     for(i=0;i<*n;i++)for(j=0;j<*nv;j++) assert(XS[i][j]==X[i][j]);
	return x;
#else
	return GLZ_Make_Trian_NF(X,n,nv,G);
#endif
}	

void Aux_vNF_Line(int l,vNF *_X,Long x[][VERT_Nmax], PERM *CL,int *S,int *_ns)
{    int n=(*_ns), cf=0;	/*  cf=CompareFlag (ref. exists & o.k.      */
     Long *y, r[VERT_Nmax];	/*  r=ReferenceLine; y->X[line]		    */
     while(n--)      /*  go over CL (n from  *_ns-1  to  0), ns_* changes!  */
     {	PERM nP[VERT_Nmax];	        
	int c=0, L=l-1, np=0, *C, ccf=cf;	/*  ccf=column compare flag */
	*nP=CL[n];
	while(++L<_X->nf)			/*  init nP (from 1st col.) */
	{   int j=0; C=nP[np].C; y=x[nP[np].L[L]];
	    while(++j<=*S) if(y[C[c]]<y[C[j]]) swap(&C[c],&C[j]);
	    if(ccf)
	    {   Long d=y[*C]-*r;
		if(d<0) ;					/* BAD line */
		else if(d)				     /* BETTER line */
		{   *r=y[*C]; cf=0; *nP=nP[np]; nP[np=1]=CL[n]; *_ns=n+1;
		    swap(&(nP->L[l]),&(nP->L[L]));
		}
		else					      /* EQUAL line */
		{   swap(&(nP[np].L[l]),&(nP[np].L[L])); nP[++np]=CL[n];
		}
	    }
	    else 						/* NEW line */
	    {	*r=y[*C]; swap(&(nP[np].L[l]),&(nP[np].L[L]));
		nP[++np]=CL[n]; ccf=1;
	    }
	}
	while(++c<_X->nv)			       /* check/complete nP */
	{   int s=S[c]; L=np; ccf=cf;
	    if(s<c) s=S[s];	
	    while(L--)
	    {	int j=c; C=nP[L].C; y=x[nP[L].L[l]];
		while(++j<=s) if(y[C[c]]<y[C[j]]) swap(&C[c],&C[j]);
		if(ccf)
		{   Long d=y[C[c]]-r[c];
		    if(d<0) {if(--np>L) nP[L]=nP[np];}		/* BAD line */
		    else if(d) 				     /* BETTER line */
		    {	r[c]=y[C[c]]; cf=0; np=L+1; *_ns=n+1;
		    }	/* else	; */			      /* EQUAL line */
		}
		else { r[c]=y[C[c]]; ccf=1; }
	    }
	}
	cf=1;
	if(--(*_ns) > n) CL[n]=CL[*_ns]; 		/*  write nP to CL  */
	if(SYM_Nmax < (cf=(*_ns+np)))
	{   printf("Need SYM_Nmax > %d !!\n",cf);exit(0);
	}
	for(L=0;L<np;L++) CL[(*_ns)++]=nP[L];
     }
     y=x[CL->L[l]];					       /* compute S */
     {	int c=0, *C=CL->C;
	while(c<_X->nv)	
	{   int s=S[c]+1; S[c]=c; while(++c<s)
	    {   if(y[C[c]]==y[C[c-1]]) ++S[ S[c]=S[c-1] ]; 
	        else S[c]=c;
	    }
	}
     }
}
void Aux_vNF_Init(vNF *_X, Long x[][VERT_Nmax], PERM *CL, int *S, int *_ns)
{    int i, j, nn; Long *b, *y;
     PERM P, *q, *p;		       /* b=x[nb] -> best;  y=x[nn] -> next */
     for(i=0;i<_X->nf;i++) P.L[i]=i;
     for(j=0;j<_X->nv;j++) P.C[j]=j; /* init P */
     q=CL; *q=P; b=*x;          /* P=CL[ns-1] StartPerm; maximize 0-th line */
     for(j=1;j<_X->nv;j++)if(b[q->C[0]]<b[q->C[j]])swap(&q->C[0],&q->C[j]);
     for(i=1;i<_X->nv;i++) 
     {  for(j=i+1;j<_X->nv;j++) if(b[q->C[i]]<b[q->C[j]]) 
        swap(&q->C[i],&q->C[j]);
     }
     for(nn=1;nn<_X->nf;nn++)			     /* maximize nn-th line */
     {  Long d; p=&CL[*_ns]; *p=P; y=x[nn]; /* nb=*q=*b=best, nn=*p=*y=next */
        {   int m=0; for(j=1;j<_X->nv;j++) if(y[p->C[m]]<y[p->C[j]]) m=j;
            if(m) swap(&p->C[0],&p->C[m]);
        }
        if((d=y[p->C[0]]-b[q->C[0]]) < 0) continue;   /* d<0 => forget this */
        for(i=1;i<_X->nv;i++) 
        {   int m=i; for(j=i+1;j<_X->nv;j++) if(y[p->C[m]]<y[p->C[j]]) m=j; 
            if(m>i) swap(&p->C[i],&p->C[m]);
            if(d==0) if((d=y[p->C[i]]-b[q->C[i]]) <0) break;
        }
        if(d<0) continue;
        swap(&p->L[0],&p->L[nn]);		 /* p->L[nn]=0; p->L[0]=nn; */ 
        if(d==0) (*_ns)++; 
        else {*q=*p; *_ns=1; b=y;}                    /* d>0 => forget prev */
     }
     y=x[CL->L[0]]; S[0]=0; for(i=1;i<_X->nv;i++)	       /* compute S */
     /*std::cout << "i: " << i << std::endl;
     std::cout << "CL->C[i]: " << CL->C[i] << std::endl;
     std::cout << "CL->C[i-1]: " << CL->C[i-1] << std::endl;
     std::cout << "y[CL->C[i]]: " << y[CL->C[i]] << std::endl;
     std::cout << "y[CL->C[i-1]]: " << y[CL->C[i-1]] << std::endl;
     std::cout << "S[i]: " << S[i] << std::endl << std::endl;*/
     if(y[CL->C[i]]==y[CL->C[i-1]]) ++S[ S[i]=S[i-1] ]; else S[i]=i;
}	

void Print_vNF(int *v, int *f, Long VPM[][VERT_Nmax], Long VPM_NF[][VERT_Nmax])
{    int i,j; fprintf(outFILE,"\nVPM NF (v=%d f=%d):\n",*v,*f); fflush(stdout);
     for(i=0;i<*f;i++)
     {  for(j=0;j<*v;j++)fprintf(outFILE,"%3d",(int)VPM[i][j]);
	fprintf(outFILE," =>");
	fflush(stdout);
        for(j=0;j<*v;j++)fprintf(outFILE,"%3d",(int)VPM_NF[i][j]);Fputs(""); 
	fflush(stdout);
     }	Fputs("");
}

void TEST_pNF(int* C,Long V[][VERT_Nmax],Long X[][VERT_Nmax],
	      int* n,int* nv,int* try_)
{    int i,j; fprintf(outFILE,"Poly NF try[%d]:   C=",*try_);
     Print_Perm(C,*nv,"\n");
     for(i=0;i<*n;i++)
     {	for(j=0;j<*nv;j++) fprintf(outFILE," %3d",(int) V[i][j]); 
	fprintf(outFILE," =>");
	for(j=0;j<*nv;j++) fprintf(outFILE," %3d",(int) X[i][j]); Fputs("");
     }
}

void Aux_Make_Triang(PERM *CL,int ns,Long V[][VERT_Nmax],int*n,int*nv,int *t)
{    int i, j, s, x=0, g=0, ps=1;		   /* x :: make X :: if X>Y */
     Long X[POLY_Dmax][VERT_Nmax], Y[POLY_Dmax][VERT_Nmax];
     for(i=0;i<*n;i++) for(j=0;j<*nv;j++) X[i][j]=V[i][CL->C[j]];
     if(!Aux_Make_Poly_NF(X,n,nv)) exit(0); 		  /* t>0: print NFs */
							  /*  -1: calc CL.s */
     if(*t) { if(*t>0) TEST_pNF(CL->C,V,X,n,nv,&g); else 
              { CL->s=1; if(*t+1){puts("t<-1 in Aux_Make_Triang");exit(0);} }
            } for(s=1;s<ns;s++)CL[s].s=0;

     for(s=1;s<ns;s++)
     if(x)
     {	for(i=0;i<*n;i++) for(j=0;j<*nv;j++) X[i][j]=V[i][CL[s].C[j]];
     	if(!Aux_Make_Poly_NF(X,n,nv)) exit(0); 
	if(Aux_XltY_Poly_NF(X,Y,n,nv)) x=0;

	if(*t)
	{   if(*t>0) TEST_pNF(CL[s].C,V,X,n,nv,&s); 
	    if(x==0)
	    {	if(*t<0) { int k; for(k=g;k<s;k++) CL[k].s=0; CL[s].s=1;*t=-1;}
		g=s; ps=1; 
	    } 
	    else if(!Aux_XltY_Poly_NF(Y,X,n,nv))
		{ if(*t<0) {CL[s].s=1;(*t)--;} ps++;}
        }
     }
     else
     {	for(i=0;i<*n;i++) for(j=0;j<*nv;j++) Y[i][j]=V[i][CL[s].C[j]];
     	if(!Aux_Make_Poly_NF(Y,n,nv)) exit(0); 
	if(Aux_XltY_Poly_NF(Y,X,n,nv)) x=1;

	if(*t)
	{   if(*t>0) TEST_pNF(CL[s].C,V,Y,n,nv,&s); 
	    if(x==1)
	    {	if(*t<0) { int k; for(k=g;k<s;k++) CL[k].s=0; CL[s].s=1;*t=-1;}
		g=s; ps=1;
	    } 
	    else if(!Aux_XltY_Poly_NF(X,Y,n,nv))
		{ if(*t<0) {CL[s].s=1;(*t)--;} ps++;}
	}
     }
     if(*t>0)
     fprintf(outFILE,
      "\nPoly NF:  NormalForm=try[%d]  #Sym(VPM)=%d  #Sym(Poly)=%d\n",g,ns,ps);
     if(x) for(i=0;i<*n;i++)for(j=0;j<*nv;j++) V[i][j]=Y[i][j];
     else  for(i=0;i<*n;i++)for(j=0;j<*nv;j++) V[i][j]=X[i][j];
}

void TEST_rVM_VPM(int *d,int *v,int *f, Long X[POLY_Dmax][VERT_Nmax],
	    	Long x[VERT_Nmax][VERT_Nmax])
{    int i,j,err=0; for(i=0;i<*v;i++)
     {	for(j=0;j<*d;j++) if(abs(X[j][i])>X_Limit) err=X[j][i];
	for(j=0;j<*f;j++) if(abs(x[j][i])>VPM_Limit) err=x[j][i]; 
     }	if(err)
     {	/*printf("TEST_VM_VPM: limits exceeded %d\n",err);
	printf("%d %d VM[%d][%d]:\n",*v,*d,*d,*v);
	for(j=0;j<*d;j++)
	{   for(i=0;i<*v;i++)printf("%3d ",(int) X[j][i]);puts("");
	}   puts("");
	printf("VPM[%d][%d]:\n",*f,*v);
	for(j=0;j<*f;j++)
	{   for(i=0;i<*v;i++)printf("%3d ",(int) x[j][i]);puts("");
	}   puts("");*/
	//exit(0);
     }
}

int  Init_rVM_VPM(PolyPointList *_P,VertexNumList *_V,EqList *_F,/* in */
	    	int *d,int *v,int *f, Long X[POLY_Dmax][VERT_Nmax],  /* out */
	    	Long x[VERT_Nmax][VERT_Nmax])		/* return reflexive */
{    int i,j, ref=1; 
     *v=_V->nv; *f=_F->ne; *d=_P->n;
     for(j=0;j<_F->ne;j++)                            	     /* compute VPM */
     {
	if(_F->e[j].c!=1) ref=0; 
	for(i=0;i<_V->nv;i++)
        x[j][i]=Eval_Eq_on_V(&_F->e[j],_P->x[_V->v[i]],_P->n);
     }
     for(i=0;i<_V->nv;i++)
     {  Long *pv=_P->x[_V->v[i]];
	for(j=0;j<_P->n;j++) X[j][i]=pv[j];
     }
     TEST_rVM_VPM(d,v,f,X,x);
     return ref;
}

void Eval_Poly_NF(int *d,int *v,int *f, Long VM[POLY_Dmax][VERT_Nmax],
		Long VPM[VERT_Nmax][VERT_Nmax],			      /* in */
		Long pNF[POLY_Dmax][VERT_Nmax],int t)		     /* out */
{    PERM *CL=(PERM *) malloc((SYM_Nmax+1) * sizeof(PERM)); 
     Long VPM_NF[VERT_Nmax][VERT_Nmax]; int ns; assert(CL!=NULL);
     Make_VPM_NF(v,f,VPM,CL,&ns,VPM_NF);	if(t)Print_vNF(v,f,VPM,VPM_NF);
     New_pNF_Order(v,f,CL,&ns,VPM_NF);
     Aux_pNF_from_vNF(CL,&ns,v,d,VM,pNF,&t);	free(CL);
#ifdef WARN_BIG_NS
#ifdef SHOW_BIG_NS
     if(SHOW_BIG_NS<=ns) {int i,j;    printf("ns=%d VM:\n",ns);fflush(stdout);
	for(i=0;i<*d;i++){for(j=0;j<*v;j++)printf("%2d ",(int)VM[i][j]);
	puts("");}	printf("ns=%d VPM:\n",ns);fflush(stdout);
	for(i=0;i<*f;i++){for(j=0;j<*v;j++)printf("%2d ",(int)VPM[i][j]);
	puts("");}	printf("ns=%d VPM_NF:\n",ns);fflush(stdout);
	for(i=0;i<*f;i++){for(j=0;j<*v;j++)printf("%2d ",(int)VPM_NF[i][j]);
	puts("");}	printf("ns=%d pNF:\n",ns);fflush(stdout);
	for(i=0;i<*d;i++){for(j=0;j<*v;j++)printf("%2d ",(int)pNF[i][j]);
	puts("");}	exit(0);}
#endif
#endif
}

void Make_VPM_NF(int *v, int *f, Long x[VERT_Nmax][VERT_Nmax],      /* in */
		PERM *CL,int *ns,Long VPM_NF[VERT_Nmax][VERT_Nmax])  /* out */
{    int i, j, S[VERT_Nmax]; int nsF=0, nsM=0; 		     /* make VPM NF */

     volatile vNF auX; vNF *_X= (vNF*) &auX;  _X->nv=*v;_X->nf=*f; /* x=VPM */
	 /*for (int& gh : CL->C) {
	 	std::cout << gh << endl;
	 }*/
     *ns=1; Aux_vNF_Init(_X, x, CL, S, ns);             /* init = 1st line */
     /*std::cout << "_X:" << endl;
     std::cout << _X->nf << endl;
     std::cout << _X->nv << endl;
     std::cout << _X->ns << endl;*/
     for(i=1;i<_X->nf-1;i++){Aux_vNF_Line(i,_X,x,CL,S,ns);  /* lines of NF */
#ifdef	WARN_BIG_NS
	if((WARN_BIG_NS<=(*ns))||nsF){ nsF=1;/*printf("ns[%d]=%d\n",i,*ns);*/}
#endif
	if(*ns>nsM) nsM=*ns; }
     _X->ns=*ns; for(i=0;i<_X->nv;i++)                /* write VPM-NF to _X */
     {  for(j=0;j<_X->nf;j++) /* _X->x */ VPM_NF[j][i]=x[CL->L[j]][CL->C[i]];
     }
     if(nsF)printf("WARNing: ns_max=%d -> ns=%d\n",nsM,*ns);
}

void Aux_pNF_from_vNF(PERM *CL,int *ns,int *v,int *d,
		Long VM[POLY_Dmax][VERT_Nmax],			      /* in */
		Long pNF[POLY_Dmax][VERT_Nmax],int *t)		     /* out */
{    int i,j;
     for(i=0;i<*d;i++)for(j=0;j<*v;j++) pNF[i][j]=VM[i][j];
     Aux_Make_Triang(CL,*ns,pNF,d,v,t);
}

int  Make_Poly_Sym_NF(PolyPointList *_P, VertexNumList *_V, EqList *_F, 
		      int *SymNum, int V_perm[][VERT_Nmax], 
		      Long NF[POLY_Dmax][VERT_Nmax], int traced, int S, int N)
{    int i, j, ns, t=-1, *d=&_P->n, *v=&_V->nv, *f=&_F->ne, *C; 
     PERM *CL = (PERM *) malloc ( sizeof(PERM) *(SYM_Nmax+1));
     Long VM[POLY_Dmax][VERT_Nmax], VPM[VERT_Nmax][VERT_Nmax];
     Long VPM_NF[VERT_Nmax][VERT_Nmax];

     Init_rVM_VPM(_P,_V,_F,d,v,f,VM,VPM);
     if (traced) Eval_Poly_NF(&_P->n,&_V->nv,&_F->ne,VM,VPM,NF,1);
     Make_VPM_NF(v,f,VPM,CL,&ns,VPM_NF);
     New_pNF_Order(v,f,CL,&ns,VPM_NF);
     Aux_pNF_from_vNF(CL,&ns,v,d,VM,NF,&t);
     *SymNum=-t;i=0; while(0==CL[i].s) i++; C=CL[i].C;
     for(t=0;i<ns;i++) if(CL[i].s)		/* inv Perm: C[c0[i]]=C[i] */
       {for(j=0;j<*v;j++) V_perm[t][C[j]]=CL[i].C[j]; t++;}
     if(*SymNum<SYM_Nmax){int *s=V_perm[*SymNum];for(i=0;i<*v;i++)s[i]=C[i];}
     if(t!=*SymNum) { puts("Error in Poly_Sym!!"); exit(0);} 
     if(traced)
     {	fprintf(outFILE,
	    "\nV_perm made by Poly_Sym (order refers to VertNumList):\n");
	for(i=0;i<*SymNum;i++) Print_Perm(V_perm[i],_V->nv,"\n");
      /*{for(j=0;j<_V->nv;j++)printf("%c",PermChar(V_perm[i][j]));puts("");}*/
     } 
      if (S) fprintf(outFILE,"#GL(Z,%d)-Symmetries=%d, #VPM-Symmetries=%d\n",
		     _P->n, *SymNum, ns);
      if (N) {
	char c[VERT_Nmax+38]="Normal form of vertices of P";
	if(*SymNum<SYM_Nmax) 
	  if(Perm_String(V_perm[*SymNum],_V->nv,&c[37]))
	    {strcpy(&c[28],"    perm");c[36]='=';}
	Print_Matrix(NF, _P->n, _V->nv,c);} 
     free(CL); return ns;
 }

int  Make_Poly_NF(PolyPointList *_P, VertexNumList *_V, EqList *_F,
		Long pNF[POLY_Dmax][VERT_Nmax])		  /* 1 if reflexive */
{    int d, v, f;
     Long VM[POLY_Dmax][VERT_Nmax], VPM[VERT_Nmax][VERT_Nmax];
     int ref=Init_rVM_VPM(_P,_V,_F,&d,&v,&f,VM,VPM);
     Eval_Poly_NF(&d,&v,&f,VM,VPM,pNF,0); return ref;
}

/* -------- Functions added from Coord.c ---------*/

void Print_Matrix(Long Matrix[][VERT_Nmax], int n_lines, int n_columns, 
		  const char *comment){
  int i,j;
  fprintf(outFILE,"%d %d  %s\n",n_lines, n_columns, comment);
  for(i=0;i<n_lines;i++) {
    for(j=0;j<n_columns;j++) fprintf(outFILE," %3d",(int) Matrix[i][j]); 
    fprintf(outFILE,"\n");}
}

/*  ======================================================================  */
/*  ==========		     			  		==========  */
/*  ==========	  B A T Y R E V ' S    F O R M U L A S       	==========  */
/*  ==========							==========  */ 
/*  ======================================================================  */

Long DualBraP1(Long *X, Long *Y, int n){    
  Long p=1; while(n--) p+=X[n] * Y[n]; return (Long) p;
}


void Make_FaceIPs(PolyPointList *_P, VertexNumList *_V, EqList *_E, 
		  PolyPointList *_DP, FaceInfo *_I){    
  /*   compute IP's of faces by computing Incidences for all points and
   *   comparing with Incidences of dual faces                             */
  int i, j, k;
  INCI x;
  for(i=0;i<_P->n;i++) for(j=0;j<_I->nf[i];j++) {
    _I->nip[i][j]=0; 
    _I->dip[i][j]=0; }
  for(k=0;k<_P->np;k++){	
    x=INCI_0(); 
    for(i=0;i<_E->ne;i++) 
      x=INCI_PN(x,Eval_Eq_on_V(&(_E->e[i]),_P->x[k],_P->n)); 
    for(i=0;i<_P->n;i++) for(j=0;j<_I->nf[i];j++) 
      if(INCI_EQ(x,_I->f[i][j])) _I->nip[i][j]++; }	
  for (k=0;k<_DP->np;k++){
    x=INCI_0(); 
    for(i=0;i<_V->nv;i++) 
      x=INCI_PN(x,DualBraP1(_P->x[_V->v[i]],_DP->x[k],_P->n));
    for(i=0;i<_P->n;i++) for(j=0;j<_I->nf[i];j++) 
      if(INCI_EQ(x,_I->v[i][j])) _I->dip[i][j]++;}
}

void PrintFaceIPs(PolyPointList *_P,FaceInfo *_I){
  int i,j, M=_P->n-1;
  for(i=0;i<=M;i++)     {	
    printf("ip[%d]:",i);
    for(j=0;j<_I->nf[i];j++) printf(" %ld",(long) _I->nip[i][j]);
    puts("");     }
  for(i=0;i<=M;i++)     {	
    printf("dip[%d]:",i);
    for(j=0;j<_I->nf[i];j++)printf(" %ld",(long) _I->dip[i][j]);
    puts("");     }
}

void Eval_BaHo(FaceInfo *_I, BaHo *_BH){
  /* Calculate Hodge/Picard numbers from FaceInfo */
  int i,j, n=_BH->n;
  int *h1;
  _BH->cor=0;
  h1=_BH->h1;
  for(i=0;i<n-1;i++) h1[i]=0;
  h1[1]+=_BH->np-n-1;					 
  for (i=0;i<_I->nf[0];i++) h1[1]-=_I->dip[0][i];
  for (j=1;j<n-1;j++) for (i=0;i<_I->nf[j];i++) {
    h1[j]+=_I->dip[j][i]*_I->nip[j][i];
    _BH->cor+=_I->dip[j][i]*_I->nip[j][i];}
  if (n==3) _BH->pic=h1[1];
  for (i=0;i<_I->nf[n-1];i++) h1[n-2]-=_I->nip[n-1][i];
  h1[n-2]+=_BH->mp-n-1;
  if (n==5) _BH->h22=44+4*h1[1]+4*h1[3]-2*h1[2];
}

void RC_Calc_BaHo(PolyPointList *_P, VertexNumList *_V, EqList *_E, 
                   PolyPointList *_DP, BaHo *_BH){
  /* Needs reflexive and complete _P and _DP */
  FaceInfo *_FI=(FaceInfo *) malloc(sizeof(FaceInfo)); 
  if(_FI==NULL) {printf("RC_Calc_BaHo: Unable to allocate _FI\n"); exit(0);}
  _BH->mp=_P->np; _BH->mv=_V->nv; _BH->nv=_E->ne; _BH->np=_DP->np; 
  _BH->n=_P->n;
  Make_Incidence(_P, _V, _E, _FI);
  Make_FaceIPs(_P, _V, _E, _DP, _FI);
  Eval_BaHo(_FI, _BH);
  free(_FI);
}