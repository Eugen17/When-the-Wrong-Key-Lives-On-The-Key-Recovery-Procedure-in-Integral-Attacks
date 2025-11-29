// alpha_admissible_present.c
// Fast C translation of alpha_admissible_present.py with micro-opts.
// Build: cc -O3 -march=native -DNDEBUG alpha_admissible_present.c -o alpha_admissible_present
// Optional: -funroll-loops

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <stdbool.h>
#include <time.h>
#include <errno.h>
#include <assert.h>

#ifndef likely
#  define likely(x)   __builtin_expect(!!(x),1)
#  define unlikely(x) __builtin_expect(!!(x),0)
#endif

enum { BLOCK_SIZE = 64, SBOX_SIZE = 4, NUM_SBOXES = 16 };
// GIFT-64 parameters / helpers
// Inverse GIFT S-box GS^{-1}
static const uint8_t GIFT_S_INV[16] = {
    0xD,0x0,0x8,0x6, 0x2,0xC,0x4,0xB,
    0xE,0x7,0x1,0xA, 0x3,0x9,0xF,0x5
};

// S^{-1} on 16-bit chunks: table of 65536 u16
static uint16_t SBOX_INV_LUT[1u<<16];

// p-layer permutation maps and LUTs (4 chunks × 65536)
static uint64_t P_LUT[4][1u<<16];
static uint64_t P_INV_LUT[4][1u<<16];

// bit permutation map
static uint8_t P_MAP[64];
static uint8_t P_INV_MAP[64];

// ---------- fast popcount parity ----------
static inline uint32_t parity64(uint64_t x) {
    return (uint32_t)(__builtin_popcountll(x) & 1);
}

// ---------- init LUTs ----------
static void init_sbox_inv_lut(void) {
    for (uint32_t val=0; val < (1u<<16); ++val) {
        uint32_t out = 0;
        uint32_t v = val;
        // process 4 nibbles
        out |= (uint32_t)GIFT_S_INV[(v >>  0) & 0xF] <<  0;
        out |= (uint32_t)GIFT_S_INV[(v >>  4) & 0xF] <<  4;
        out |= (uint32_t)GIFT_S_INV[(v >>  8) & 0xF] <<  8;
        out |= (uint32_t)GIFT_S_INV[(v >> 12) & 0xF] << 12;
        SBOX_INV_LUT[val] = (uint16_t)out;
    }
}

static void init_p_maps_and_luts(void) {
    // GIFT-64 bit permutation P64(i)
    static const uint8_t GIFT_P64[64] = {
        48, 1, 18, 35, 32, 49, 2, 19, 16, 33, 50, 3, 0, 17, 34, 51, 52,
             5, 22, 39, 36, 53, 6, 23, 20, 37, 54, 7, 4, 21, 38, 55, 56, 
             9, 26, 43, 40, 57, 10, 27, 24, 41, 58, 11, 8, 25, 42, 59, 60,
              13, 30, 47, 44, 61, 14, 31, 28, 45, 62, 15, 12, 29, 46, 63

    };

    for (int i=0;i<64;++i) {
        P_MAP[i] = GIFT_P64[i];
        P_INV_MAP[P_MAP[i]] = (uint8_t)i;
    }

    // Build P_LUT and P_INV_LUT chunkwise (unchanged)
    for (int chunk=0; chunk<4; ++chunk) {
        const int shift_in = 16*chunk;
        for (uint32_t val=0; val < (1u<<16); ++val) {
            uint64_t perm = 0, permi = 0;
            uint32_t v = val;
            for (int bit=0; bit<16; ++bit) {
                if (v & (1u<<bit)) {
                    int src = shift_in + bit;
                    int dst = P_MAP[src];
                    perm |= (1ULL<<dst);

                    int dsti = P_INV_MAP[src];
                    permi |= (1ULL<<dsti);
                }
            }
            P_LUT[chunk][val]     = perm;
            P_INV_LUT[chunk][val] = permi;
        }
    }
}

static inline uint64_t sbox_layer_inv(uint64_t x) {
    // apply 16-bit LUT to each 16-bit lane
    uint64_t a = SBOX_INV_LUT[ (uint16_t)( x        & 0xFFFFu) ];
    uint64_t b = SBOX_INV_LUT[ (uint16_t)((x>>16)   & 0xFFFFu) ];
    uint64_t c = SBOX_INV_LUT[ (uint16_t)((x>>32)   & 0xFFFFu) ];
    uint64_t d = SBOX_INV_LUT[ (uint16_t)((x>>48)   & 0xFFFFu) ];
    return (a) | (b<<16) | (c<<32) | (d<<48);
}

static inline uint64_t p_layer(uint64_t x) {
    return P_LUT[0][ (uint16_t)(x      & 0xFFFFu) ] |
           P_LUT[1][ (uint16_t)((x>>16)& 0xFFFFu) ] |
           P_LUT[2][ (uint16_t)((x>>32)& 0xFFFFu) ] |
           P_LUT[3][ (uint16_t)((x>>48)& 0xFFFFu) ];
}

static inline uint64_t p_layer_inv(uint64_t x) {
    return P_INV_LUT[0][ (uint16_t)(x      & 0xFFFFu) ] |
           P_INV_LUT[1][ (uint16_t)((x>>16)& 0xFFFFu) ] |
           P_INV_LUT[2][ (uint16_t)((x>>32)& 0xFFFFu) ] |
           P_INV_LUT[3][ (uint16_t)((x>>48)& 0xFFFFu) ];
}

static inline uint64_t F_inverse_round_core(uint64_t x) {
    // F = S^{-1} ∘ L^{-1}
    return sbox_layer_inv(p_layer_inv(x));
}

// ============================================================
// RNG (xoshiro256**): fast, good quality
// ============================================================
typedef struct { uint64_t s[4]; } rng_t;

static inline uint64_t rotl64(uint64_t x, int k){ return (x<<k) | (x>>(64-k)); }

static inline uint64_t rng_next(rng_t* r){
    const uint64_t result = rotl64(r->s[1] * 5u, 7) * 9u;
    const uint64_t t = r->s[1] << 17;
    r->s[2] ^= r->s[0];
    r->s[3] ^= r->s[1];
    r->s[1] ^= r->s[2];
    r->s[0] ^= r->s[3];
    r->s[2] ^= t;
    r->s[3] = rotl64(r->s[3], 45);
    return result;
}

static void rng_seed(rng_t* r, uint64_t seed){
    // SplitMix64 to expand seed
    uint64_t x = seed ? seed : (uint64_t)time(NULL);
    for (int i=0;i<4;++i){
        x += 0x9e3779b97f4a7c15ULL;
        uint64_t z = x;
        z = (z ^ (z>>30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z>>27)) * 0x94d049bb133111ebULL;
        r->s[i] = z ^ (z>>31);
    }
}

// ============================================================
// GF(2) helpers on 64-bit
// ============================================================
static inline int bitlen64(uint64_t x){
    return x ? (64 - __builtin_clzll(x)) : 0;
}

// standard row-reduction, basis sorted by decreasing bitlen
typedef struct {
    uint64_t v[64]; // safe upper bound; realistic dims are small
    int n;
} basis64_t;

static void basis64_init(basis64_t* B){ B->n = 0; }

static void basis64_add(basis64_t* B, uint64_t x){
    uint64_t t = x;
    for (int i=0;i<B->n;++i){
        uint64_t b = B->v[i];
        if ((t ^ b) < t) t ^= b;
    }
    if (!t) return;
    // insert sorted by decreasing bitlen (descending numeric works under this construction)
    int pos = 0;
    while (pos < B->n && B->v[pos] > t) ++pos;
    // shift right
    for (int i=B->n; i>pos; --i) B->v[i] = B->v[i-1];
    B->v[pos] = t;
    B->n++;
}

static inline int in_span64(uint64_t x, const basis64_t* B){
    uint64_t t = x;
    for (int i=0;i<B->n;++i){
        uint64_t b = B->v[i];
        if (bitlen64(t) == bitlen64(b)) t ^= b;
        // else skip; second pass as in python isn't really needed with our construction
    }
    for (int i=0;i<B->n;++i){
        uint64_t b = B->v[i];
        if (bitlen64(t) == bitlen64(b)) t ^= b;
    }
    return t==0;
}

// ============================================================
// Affine space
// ============================================================
typedef struct {
    uint64_t offset;
    uint64_t basis[64];
    int dim;
} AffineSpace;

static inline void affine_init(AffineSpace* A, uint64_t offset){
    A->offset = offset;
    A->dim = 0;
}

static inline void affine_copy(AffineSpace* dst, const AffineSpace* src){
    dst->offset = src->offset;
    dst->dim = src->dim;
    for (int i=0;i<src->dim;++i) dst->basis[i] = src->basis[i];
}

static inline void affine_add_vector(AffineSpace* out, const AffineSpace* A, uint64_t v){
    affine_copy(out, A);
    out->offset ^= v;
}

static void affine_span_from_points(AffineSpace* out, const uint64_t* pts, int npts){
    assert(npts>0);
    uint64_t o0 = pts[0];
    basis64_t B; basis64_init(&B);
    for (int i=1;i<npts;++i){
        basis64_add(&B, pts[i]^o0);
    }
    out->offset = o0;
    out->dim = B.n;
    for (int i=0;i<B.n;++i) out->basis[i] = B.v[i];
}

static void affine_sample(const AffineSpace* A, uint64_t* out, int num, rng_t* rng){
    if (A->dim == 0){
        for (int i=0;i<num;++i) out[i] = A->offset;
        return;
    }
    for (int i=0;i<num;++i){
        uint64_t comb = 0;
        // pick random combination
        for (int j=0;j<A->dim;++j){
            if (rng_next(rng) & 1u) comb ^= A->basis[j];
        }
        out[i] = A->offset ^ comb;
    }
}

// ============================================================
// Propagation for affine subspaces via sampling
// ============================================================
static void propagate_affine(
    const AffineSpace* U,
    const uint64_t* xs, int Z,
    AffineSpace* V_out
){
    // O = {F(x) ^ F(x ^ u)}
    // B_i = {F(x) ^ F(x ^ u_i ^ u)}
    // gather points: size = Z*(1 + dim)
    const int total = Z * (1 + U->dim);
    uint64_t* pts = (uint64_t*)malloc((size_t)total * sizeof(uint64_t));
    int idx = 0;

    for (int k=0;k<Z;++k){
        uint64_t x = xs[k];
        pts[idx++] = F_inverse_round_core(x) ^ F_inverse_round_core(x ^ U->offset);
    }
    for (int j=0;j<U->dim;++j){
        uint64_t uij = U->basis[j] ^ U->offset;
        for (int k=0;k<Z;++k){
            uint64_t x = xs[k];
            pts[idx++] = F_inverse_round_core(x) ^ F_inverse_round_core(x ^ uij);
        }
    }
    affine_span_from_points(V_out, pts, total);
    free(pts);
}

// ============================================================
// Candidates: <= 1 active nibble after L
// ============================================================
typedef struct {
    // We’ll iterate on the fly—no storage needed—but keep helpers here if desired.
    int include_zero;
} CandGen;

static inline void candgen_init(CandGen* G, int include_zero){ G->include_zero = include_zero; }

// Iterate via nested loops in-place where needed.

// ============================================================
// Inner product <alpha,x> mod 2 and constancy test
// ============================================================
static inline uint32_t inner_product_mod2(uint64_t a, uint64_t b){
    return parity64(a & b);
}

static int is_inner_product_constant_on_affine(uint64_t alpha, const AffineSpace* V, int Z, rng_t* rng){
    if (V->dim == 0) {
        uint32_t v = inner_product_mod2(alpha, V->offset);
        (void)Z;
        return 1; // trivially constant
    }
    // sample Z points
    uint64_t x0;
    // first sample
    {
        uint64_t comb = 0;
        for (int j=0;j<V->dim;++j) if (rng_next(rng)&1u) comb ^= V->basis[j];
        x0 = V->offset ^ comb;
    }
    uint32_t v0 = inner_product_mod2(alpha, x0);
    for (int i=1;i<Z;++i){
        uint64_t comb = 0;
        for (int j=0;j<V->dim;++j) if (rng_next(rng)&1u) comb ^= V->basis[j];
        uint64_t xi = V->offset ^ comb;
        if (inner_product_mod2(alpha, xi) != v0) return 0;
    }
    return 1;
}

// ============================================================
// Recursive search
// ============================================================
typedef struct {
    // dynamic array of solutions; each solution is r words (w1..wr)
    uint64_t* data;
    int r;
    size_t n, cap;
} SolVec;

static void solvec_init(SolVec* S, int r){
    S->data = NULL; S->r = r; S->n = 0; S->cap = 0;
}
static void solvec_push(SolVec* S, const uint64_t* tuple_r){
    if (S->n == S->cap){
        S->cap = S->cap ? (S->cap*2) : 64;
        S->data = (uint64_t*)realloc(S->data, (size_t)S->cap * (size_t)S->r * sizeof(uint64_t));
    }
    memcpy(&S->data[S->n * (size_t)S->r], tuple_r, (size_t)S->r * sizeof(uint64_t));
    S->n++;
}
static void solvec_free(SolVec* S){ free(S->data); }

// propagate recursion
typedef struct {
    int r;
    uint64_t alpha;
    int Z;
    const uint64_t* xs; // size Z
    CandGen G;
    rng_t* rng;
} SearchCtx;

static void recurse_level(int i, const AffineSpace* V_i, uint64_t* w_suffix_rev, int wlen,
                          const SearchCtx* C, SolVec* out)
{
    if (i == 0){
        if (is_inner_product_constant_on_affine(C->alpha, V_i, C->Z, C->rng)) {
            // reverse the w_suffix_rev (which holds w_{i+1}..w_r) to (w1..wr)
            uint64_t* wtuple = (uint64_t*)alloca((size_t)C->r * sizeof(uint64_t));
            for (int k=0;k<C->r;++k) wtuple[k] = w_suffix_rev[C->r-1-k];
            solvec_push(out, wtuple);
        }
        return;
    }

    const uint64_t v_i = V_i->offset;

    // include_zero?
    if (C->G.include_zero){
        uint64_t c_i = 0;
        uint64_t w_i = v_i ^ p_layer(c_i);
        AffineSpace U_i; affine_add_vector(&U_i, V_i, w_i);
        AffineSpace V_im1; propagate_affine(&U_i, C->xs, C->Z, &V_im1);
        w_suffix_rev[wlen] = w_i;
        recurse_level(i-1, &V_im1, w_suffix_rev, wlen+1, C, out);
    }

    // for each nibble position and nonzero value
    for (int nib=0; nib<NUM_SBOXES; ++nib){
        uint64_t base = (uint64_t)nib * SBOX_SIZE;
        for (uint64_t val=1; val<(1u<<SBOX_SIZE); ++val){
            uint64_t c_i = (val << base);
            uint64_t w_i = v_i ^ p_layer(c_i);
            AffineSpace U_i; affine_add_vector(&U_i, V_i, w_i);
            AffineSpace V_im1; propagate_affine(&U_i, C->xs, C->Z, &V_im1);
            w_suffix_rev[wlen] = w_i;
            recurse_level(i-1, &V_im1, w_suffix_rev, wlen+1, C, out);
        }
    }
}

// ============================================================
// Big-int over GF(2) for basis over concatenated tuples
// ============================================================
typedef struct {
    int nlimbs;          // = r
    uint64_t* limb;      // array of length nlimbs, big-endian (limb 0 most significant)
} Big;

static inline Big big_make(int nlimbs){
    Big b; b.nlimbs = nlimbs;
    b.limb = (uint64_t*)calloc((size_t)nlimbs, sizeof(uint64_t));
    return b;
}

static inline void big_free(Big* b){ free(b->limb); b->limb=NULL; b->nlimbs=0; }

static inline void big_copy(Big* dst, const Big* src){
    dst->nlimbs = src->nlimbs;
    memcpy(dst->limb, src->limb, (size_t)src->nlimbs * sizeof(uint64_t));
}

static inline int big_is_zero(const Big* b){
    for (int i=0;i<b->nlimbs;++i) if (b->limb[i]) return 0;
    return 1;
}

static inline int big_cmp(const Big* a, const Big* b){
    // lexicographic on limbs
    for (int i=0;i<a->nlimbs;++i){
        if (a->limb[i] != b->limb[i]) return (a->limb[i] > b->limb[i]) ? 1 : -1;
    }
    return 0;
}

static inline void big_xor(Big* dst, const Big* src){
    for (int i=0;i<dst->nlimbs;++i) dst->limb[i] ^= src->limb[i];
}

static int big_bitlen(const Big* b){
    for (int i=0;i<b->nlimbs;++i){
        uint64_t w = b->limb[i];
        if (w){
            int rem = 64 - __builtin_clzll(w);
            return (b->nlimbs-1-i)*64 + rem;
        }
    }
    return 0;
}

static Big concat_tuple_to_big(const uint64_t* wtuple, int r){
    // big-endian: limb 0 is the most significant w, i.e., w1 at limb 0, … wr at limb r-1
    Big B = big_make(r);
    for (int i=0;i<r;++i) B.limb[i] = wtuple[i];
    return B;
}

typedef struct {
    Big* vec;
    int n, cap;
    int nlimbs;
} BigVec;

static void bigvec_init(BigVec* V, int nlimbs){
    V->vec=NULL; V->n=0; V->cap=0; V->nlimbs=nlimbs;
}
static void bigvec_push(BigVec* V, const Big* b){
    if (V->n == V->cap){
        V->cap = V->cap ? V->cap*2 : 32;
        V->vec = (Big*)realloc(V->vec, (size_t)V->cap*sizeof(Big));
    }
    Big c = big_make(V->nlimbs);
    big_copy(&c, b);
    V->vec[V->n++] = c;
}
static void bigvec_free(BigVec* V){
    for (int i=0;i<V->n;++i) big_free(&V->vec[i]);
    free(V->vec);
}

// eliminate basis over GF(2), descending by bitlen
static void gf2_eliminate_basis_big(BigVec* V){
    // simple insertion with reduction
    BigVec out; bigvec_init(&out, V->nlimbs);

    for (int i=0;i<V->n;++i){
        Big x = big_make(V->nlimbs);
        big_copy(&x, &V->vec[i]);

        // reduce by current out basis
        for (int j=0;j<out.n;++j){
            Big* b = &out.vec[j];
            // if (x ^ b) < x  in lex order
            Big tmp = big_make(V->nlimbs);
            big_copy(&tmp, &x);
            big_xor(&tmp, b);
            if (big_cmp(&tmp, &x) < 0){
                big_free(&x); x = tmp; // x = tmp
            } else {
                big_free(&tmp);
            }
        }
        if (!big_is_zero(&x)){
            // insert x to keep descending "numeric" (lex) order
            int pos = 0;
            while (pos<out.n && big_cmp(&out.vec[pos], &x) > 0) ++pos;
            // shift
            out.vec = (Big*)realloc(out.vec, (size_t)(out.n+1)*sizeof(Big));
            for (int k=out.n; k>pos; --k) out.vec[k] = out.vec[k-1];
            out.vec[pos] = x;
            out.n++;
        } else {
            big_free(&x);
        }
    }
    // move out back to V
    for (int i=0;i<V->n;++i) big_free(&V->vec[i]);
    free(V->vec);
    *V = out;
}

// ============================================================
// CLI and main
// ============================================================
static void usage(const char* prog){
    fprintf(stderr,
        "Usage: %s [-r ROUNDS] [-Z SAMPLES] [-a ALPHA_HEX] [-s SEED] [--no-zero]\n"
        "  -r ROUNDS     number of rounds (default 3)\n"
        "  -Z SAMPLES    sampling size Z (default 100)\n"
        "  -a ALPHA      64-bit mask in hex (default 0x1)\n"
        "  -s SEED       RNG seed (default time)\n"
        "  --no-zero     exclude c=0 candidate (default include)\n"
        "  -h            this help\n", prog);
}

int main(int argc, char** argv){
    int r = 3;
    int Z = 50;
    uint64_t alpha = 0x1ULL;
    uint64_t seed = 0;
    int include_zero = 1;

    // parse CLI
    for (int i=1;i<argc;++i){
        if (!strcmp(argv[i], "-r") && i+1<argc){ r = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "-Z") && i+1<argc){ Z = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "-a") && i+1<argc){
            const char* s = argv[++i];
            if (!strncmp(s,"0x",2) || !strncmp(s,"0X",2)) alpha = strtoull(s, NULL, 16);
            else alpha = strtoull(s, NULL, 10);
        }
        else if (!strcmp(argv[i], "-s") && i+1<argc){
            seed = (uint64_t)strtoull(argv[++i], NULL, 10);
        }
        else if (!strcmp(argv[i], "--no-zero")){ include_zero = 0; }
        else if (!strcmp(argv[i], "-h")){ usage(argv[0]); return 0; }
        else { fprintf(stderr,"Unknown arg: %s\n", argv[i]); usage(argv[0]); return 1; }
    }

    if (r <= 0){ fprintf(stderr,"r must be >=1\n"); return 1; }
    if (Z <= 0){ fprintf(stderr,"Z must be >=1\n"); return 1; }

    // init LUTs
    init_sbox_inv_lut();
    init_p_maps_and_luts();

    // RNG & sample xs
    rng_t R; rng_seed(&R, seed);
    uint64_t* XS = (uint64_t*)malloc((size_t)Z*sizeof(uint64_t));
    for (int i=0;i<Z;++i) XS[i] = rng_next(&R);

    // Search setup
    SearchCtx Ctx;
    Ctx.r = r; Ctx.alpha = alpha; Ctx.Z = Z; Ctx.xs = XS;
    candgen_init(&Ctx.G, include_zero);
    Ctx.rng = &R;

    // Start with V_r = {0}
    AffineSpace V_r; affine_init(&V_r, 0);

    // run recursion
    SolVec sols; solvec_init(&sols, r);
    uint64_t* w_suffix_rev = (uint64_t*)alloca((size_t)r * sizeof(uint64_t));
    recurse_level(r, &V_r, w_suffix_rev, 0, &Ctx, &sols);

    // Build GF(2) basis over concatenated tuples
    BigVec bigs; bigvec_init(&bigs, r);
    // de-duplicate by a tiny open-addressing set (optional). For speed, skip strict dedup;
    // elimination + XOR identities make duplicates harmless. If needed, we could add a hash set.

    for (size_t i=0;i<sols.n; ++i){
        const uint64_t* tuple_i = &sols.data[i*(size_t)r];
        Big B = concat_tuple_to_big(tuple_i, r);
        bigvec_push(&bigs, &B);
        big_free(&B);
    }
    gf2_eliminate_basis_big(&bigs);

    printf("Found a basis of %d alpha-admissible key differences for alpha = 0x%016" PRIx64 ", r=%d:\n",
           bigs.n, alpha, r);
    for (int i=0;i<bigs.n;++i){
        // print tuple (hex)
        printf("(");
        for (int k=0;k<r;++k){
            printf("0x%016" PRIx64, bigs.vec[i].limb[k]);
            if (k+1<r) printf(", ");
        }
        printf(")\n");
    }

    // cleanup
    bigvec_free(&bigs);
    solvec_free(&sols);
    free(XS);
    return 0;
}
