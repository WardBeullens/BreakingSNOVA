#include "snova.h"

#include "gf16_matrix_inline.h"
#include "snova_kernel.h"
#include "snova_plasma/snova_plasms_option.h"
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <sys/time.h>


static inline int64_t cpucycles(void) {
    unsigned int hi, lo;

    asm volatile ("rdtsc" : "=a" (lo), "=d"(hi));
    return ((int64_t) lo) | (((int64_t) hi) << 32);
}

#define TIC printf("\n"); \
        int64_t tic_toc_cycles = cpucycles();

#define TOC(name) printf(" %-30s cycles: %lu \n", name, cpucycles() - tic_toc_cycles); \
        tic_toc_cycles = cpucycles();



/**
 * SNOVA init
 */
void snova_init() {
    static int first_time = 1;
    if (first_time) {
        first_time = 0;
        init_gf16_tables();
        gen_S_array();

#if OPTIMISATION == 2 && rank == 4
        init_4x4();
#elif OPTIMISATION == 2 && rank == 2
        init_2x2();
#elif OPTIMISATION == 2
        init_avx_table();
#endif
    }
}

/**
 * generate snova key elements.
 * @param key_elems - pointer to output snova key elements.
 * @param pk_seed - pointer to input public key seed.
 * @param sk_seed - pointer to input private key elements.
 */
void generate_keys_core(snova_key_elems* key_elems, const uint8_t* pk_seed, const uint8_t* sk_seed) {
    gen_seeds_and_T12(key_elems->T12, sk_seed);
    memcpy(key_elems->pk.pt_public_key_seed, pk_seed, seed_length_public);
    gen_A_B_Q_P(&(key_elems->map1), pk_seed);
    gen_F(&(key_elems->map2), &(key_elems->map1), key_elems->T12);
    gen_P22(key_elems->T12, key_elems->map1.P21, key_elems->map2.F12, key_elems->pk.P22);
}

/**
 * Generates public and private key. where private key is the seed of private
 * key.
 * @param pkseed - pointer to input public key seed.
 * @param skseed - pointer to input private key seed.
 * @param pk - pointer to output public key.
 * @param ssk - pointer to output private key.
 */
void generate_keys_ssk(const uint8_t* pkseed, const uint8_t* skseed, uint8_t* pk, uint8_t* ssk) {
    snova_init();
    snova_key_elems key_elems;
    generate_keys_core(&key_elems, pkseed, skseed);
    pk_pack(pk, &key_elems);
    memcpy(ssk, pkseed, seed_length_public);
    memcpy(ssk + seed_length_public, skseed, seed_length_private);
}

/**
 * Generates public and private key. where private key is the expanded version.
 * @param pkseed - pointer to input public key seed.
 * @param skseed - pointer to input private key seed.
 * @param pk - pointer to output public key.
 * @param esk - pointer to output private key. (expanded)
 */
void generate_keys_esk(const uint8_t* pkseed, const uint8_t* skseed, uint8_t* pk, uint8_t* esk) {
    snova_init();
    snova_key_elems key_elems;
    generate_keys_core(&key_elems, pkseed, skseed);
    pk_pack(pk, &key_elems);
    sk_pack(esk, &key_elems, skseed);
}

/**
 * Compute the signature using ssk (private key seed). some preparatory work
 * before using sign_digest_core()
 * @param pt_signature - pointer to output signature.
 * @param digest - pointer to input digest.
 * @param array_salt - pointer to input salt.
 * @param ssk - pointer to input private key (seed).
 */
void sign_digest_ssk(uint8_t* pt_signature, const uint8_t* digest, uint64_t bytes_digest, uint8_t* array_salt,
                     const uint8_t* ssk) {
    snova_init();
    snova_key_elems key_elems;
    const uint8_t* pk_seed = ssk;
    const uint8_t* sk_seed = ssk + seed_length_public;

    gen_seeds_and_T12(key_elems.T12, sk_seed);
    gen_A_B_Q_P(&(key_elems.map1), pk_seed);
    gen_F(&(key_elems.map2), &(key_elems.map1), key_elems.T12);
    sign_digest_core(pt_signature, digest, bytes_digest, array_salt, key_elems.map1.Aalpha, key_elems.map1.Balpha,
                     key_elems.map1.Qalpha1, key_elems.map1.Qalpha2, key_elems.T12, key_elems.map2.F11, key_elems.map2.F12,
                     key_elems.map2.F21, pk_seed, sk_seed);
}

/**
 * Compute the signature using esk (). some preparatory work before using
 * sign_digest_core()
 * @param pt_signature - pointer to output signature.
 * @param digest - pointer to input digest.
 * @param array_salt - pointer to input salt.
 * @param esk - pointer to input private key (expanded).
 */
void sign_digest_esk(uint8_t* pt_signature, const uint8_t* digest, uint64_t bytes_digest, uint8_t* array_salt,
                     const uint8_t* esk) {
    snova_init();
    sk_gf16 sk_upk;
    sk_unpack(&sk_upk, esk);
    sign_digest_core(pt_signature, digest, bytes_digest, array_salt, sk_upk.Aalpha, sk_upk.Balpha, sk_upk.Qalpha1,
                     sk_upk.Qalpha2, sk_upk.T12, sk_upk.F11, sk_upk.F12, sk_upk.F21, sk_upk.pt_public_key_seed,
                     sk_upk.pt_private_key_seed);
}

/**
 * Verifies signature.
 * @param pt_digest - pointer to input digest.
 * @param pt_signature - pointer to output signature.
 * @param pk - pointer to output public key.
 * @returns - 0 if signature could be verified correctly and -1 otherwise
 */
int verify_signture(const uint8_t* pt_digest, uint64_t bytes_digest, const uint8_t* pt_signature, const uint8_t* pk) {
    snova_init();
    return verify_core(pt_digest, bytes_digest, pt_signature, pk);
}

//// Attack stuff added by Ward Beullens

// given M = a_0 S^0 + a_1 S^1 + ...  outputs a_0, a_1, ... 
void get_coefss_from_S_matrix(gf16_t *coeffs, gf16m_t M){
    #if rank > 4
        only implemented for l=2, l=3, l=4 so far
    #endif

    switch (rank)
    {
    case 2:
        coeffs[1] = gf16_get_mul(M[1],inv(7));
        coeffs[0] = gf16_get_add(M[0],gf16_get_mul(coeffs[1],8));
        break;
    case 3:
        coeffs[0] = M[0] ^ gf16_get_mul(M[1],8) ^ gf16_get_mul(M[2],8);
        coeffs[1] = gf16_get_mul(M[1],7) ^ gf16_get_mul(M[2],6);
        coeffs[2] = gf16_get_mul(M[1],7) ^ M[2];
        break;
    case 4:
        coeffs[0] = M[0] ^ gf16_get_mul(M[1],3) ^ gf16_get_mul(M[2],12) ^ gf16_get_mul(M[3],3);
        coeffs[1] = gf16_get_mul(M[1],6) ^ gf16_get_mul(M[2],6) ^ gf16_get_mul(M[3],4);
        coeffs[2] = gf16_get_mul(M[1],1) ^ gf16_get_mul(M[2],8) ^ gf16_get_mul(M[3],9);
        coeffs[3] = gf16_get_mul(M[1],4) ^ gf16_get_mul(M[2],13) ^ gf16_get_mul(M[3],4);
        break;
    }
}

void compute_matrices(uint8_t *pk_seed, gf16_t *matrices){
    map_group1 randomness;
    gen_A_B_Q_P(&randomness, pk_seed);

    gf16_t Lc[rank],Rc[rank],U[sq_rank],V[sq_rank];

    // prepare powers of S
    gf16m_t powers_of_S[rank];
    be_aI(powers_of_S[0], 1);
    be_the_S(powers_of_S[1]);
    for (size_t i = 2; i < rank; i++)
    {
        gf16m_mul( powers_of_S[i-1], powers_of_S[1], powers_of_S[i] );
    }

    for (size_t alpha = 0; alpha < sq_rank; alpha++)
    {
        for (size_t i = 0; i < sq_rank; i++)
        {
            for (size_t j = 0; j < sq_rank; j++)
            {
                gf16m_t Q1S,Q2S;
                
                gf16m_mul(randomness.Qalpha1[alpha], powers_of_S[i%rank], Q1S);
                gf16m_mul(randomness.Qalpha2[alpha], powers_of_S[j%rank], Q2S);

                get_coefss_from_S_matrix(Lc, Q1S);
                get_coefss_from_S_matrix(Rc, Q2S);

                for (size_t a = 0; a < rank; a++)
                {
                    for (size_t b = 0; b < rank; b++)
                    {
                        V[a + rank*b] = gf16_get_mul(Lc[b],Rc[a]);
                        U[a + rank*b] = gf16_get_mul(get_gf16m(randomness.Aalpha[alpha], b, i/rank),get_gf16m(randomness.Balpha[alpha], j/rank, a));
                    }
                }

                int ii,jj;

                if(i <= j){
                    ii = i;
                    jj = j;
                }
                else{
                    ii = j;
                    jj = i;
                }

                for (size_t x = 0; x < sq_rank; x++)
                {
                    for (size_t y = 0; y < sq_rank; y++)
                    {
                        matrices[(ii*sq_rank + jj)*sq_rank*sq_rank + x*sq_rank + y] ^= gf16_get_mul(V[y],U[x]);
                    }
                }
            }
        }
    }
}

// check if a matrix with l^2 columns has rank 1, assumes that top row is nonzero
int _has_rank_1( const gf16_t *matrix, const size_t rows, const size_t nonzero_col ){
    
    gf16_t inverse = inv(matrix[nonzero_col]);

    for (size_t row = 1; row < rows; row++)
    {
        gf16_t multiplier = gf16_get_mul(inverse, matrix[row*sq_rank + nonzero_col]);

        for (size_t col = 0; col < sq_rank; col++)
        {
            if (gf16_get_mul(matrix[col], multiplier) != matrix[row*sq_rank + col]){
                return 0;
            }
        }
    }
    
    return 1;
}

// check if a matrix with l^2 columns has rank at most 1
int has_rank_1( const gf16_t *matrix, const size_t rows ){
   
    // find nonzero element in top row
    for (size_t i = 0; i < sq_rank; i++)
    {
        if (matrix[i] != 0){
            return _has_rank_1(matrix, rows, i);
        }
    }
    
    // top row is zero -> matrix has rank <= 1 iff remaining rows have rank <=1
    return has_rank_1(matrix + sq_rank, rows-1);
}

static inline uint64_t mul_fx8(unsigned char a, uint64_t b) {
    // carryless multiply
    uint64_t p;
    p  = (a & 1)*b;
    p ^= (a & 2)*b;
    p ^= (a & 4)*b;
    p ^= (a & 8)*b;

    // reduce mod x^4 + x + 1
    uint64_t top_p = p & 0xf0f0f0f0f0f0f0f0;
    uint64_t out = (p ^ (top_p >> 4) ^ (top_p >> 3)) & 0x0f0f0f0f0f0f0f0f;
    return out;
}

static inline uint64_t cl_mul_fx8(unsigned char a, uint64_t b) {
    // carryless multiply
    return (a & 1)*b ^ (a & 2)*b ^ (a & 4)*b ^ (a & 8)*b;
}

static inline void redx8(uint64_t *p) {
    // reduce mod x^4 + x + 1
    uint64_t top_p = *p & 0xf0f0f0f0f0f0f0f0;
    *p = (*p ^ (top_p >> 4) ^ (top_p >> 3)) & 0x0f0f0f0f0f0f0f0f;
}

# if rank == 3
void load_mat(gf16_t* in, uint64_t* out){
    out[0] = in[0*9];
    out[1] = in[1*9];
    out[2] = in[2*9];
    out[3] = in[3*9];
    out[4] = in[4*9];
    out[5] = in[5*9];
    out[6] = in[6*9];
    out[7] = in[7*9];
    out[8] = in[8*9];
    memcpy(out+9 , in + 1, 8);
    memcpy(out+10, in + 1+ 1*9, 8);
    memcpy(out+11, in + 1+ 2*9, 8);
    memcpy(out+12, in + 1+ 3*9, 8);
    memcpy(out+13, in + 1+ 4*9, 8);
    memcpy(out+14, in + 1+ 5*9, 8);
    memcpy(out+15, in + 1+ 6*9, 8);
    memcpy(out+16, in + 1+ 7*9, 8);
    memcpy(out+17, in + 1+ 8*9, 8);
}

void mul_mat(const uint64_t *A, const uint64_t*B, uint64_t *C){

    #if rank != 3
        not implemented;
    #endif

    #pragma GCC unroll 9
    for (size_t i = 0; i < 9; i++)
    {
        C[i] = cl_mul_fx8(A[i],B[0]);
        C[9+i] = cl_mul_fx8(A[i],B[9]);
    }

    #pragma GCC unroll 8
    for (size_t col = 0; col < 8; col++)
    {
        #pragma GCC unroll 9
        for (size_t i = 0; i < 9; i++)
        {
            uint8_t c = A[9+i] >> (col)*8;
            C[i] ^= cl_mul_fx8(c,B[col+1]);
            C[9+i] ^= cl_mul_fx8(c,B[10+col]);
        }
    }

    #pragma GCC unroll 18
    for (size_t i = 0; i < 18; i++)
    {
        redx8(C + i);
    }
}

int ef_mat(uint64_t *mat, int tres){

    #if rank != 3
        not implemented;
    #endif

    uint64_t *first_col = mat;
    mat += 9;

    int row = 0;

    // first row
    // force nonzero pivot
    int pivot_row = row;
    while ((pivot_row < (sq_rank-1)) && first_col[ pivot_row ] == 0)
    {
        pivot_row += 1;
    }

    if(first_col[pivot_row] != 0){
        if(pivot_row > row){
            first_col[row] ^= first_col[pivot_row];
            mat[row] ^= mat[pivot_row];
        }

        gf16_t inverse = inv(first_col[row]);

        for (size_t i = row+1; i < sq_rank; i++)
        {
            gf16_t scalar = gf16_get_mul(inverse, first_col[i]);
            mat[i] ^= mul_fx8(scalar, mat[row]);
        }

        row ++;        
        if (row > tres){
            return row;
        }
    }

    for (uint64_t pivot_col = 0; pivot_col < 8; pivot_col ++)
    {
        // force nonzero pivot
        pivot_row = row;
        while ((pivot_row < (sq_rank-1)) && ((mat[pivot_row] >> (pivot_col*8)) & 0xf )== 0)
        {
            pivot_row += 1;
        }

        if((mat[pivot_row] >> (pivot_col*8) & 0xf ) == 0){
            continue;
        }

        if(pivot_row > row){
            mat[row] ^= mat[pivot_row];
        }

        gf16_t inverse = inv((mat[row] >> (pivot_col*8)) & 0xf );

        for (size_t i = row+1; i < sq_rank; i++)
        {
            gf16_t scalar = gf16_get_mul(inverse, (mat[i] >> (pivot_col*8)) & 0xf);
            mat[i] ^= mul_fx8(scalar, mat[row]);
        }

        row ++;        
        if (row > tres){
            return row;
        }
    }
    
    return row;
}

#endif

int ef(gf16_t* mat, int tres){
    int row = 0;
    for (int pivot_col = 0; pivot_col < sq_rank; pivot_col ++)
    {
        // force nonzero pivot
        int pivot_row = row;
        while ((pivot_row < (sq_rank-1)) && mat[ pivot_col + pivot_row*sq_rank ] == 0)
        {
            pivot_row += 1;
        }

        if(mat[ pivot_col + pivot_row*sq_rank ] == 0){
            continue;
        }

        if(pivot_row > row){
            for (size_t i = pivot_col; i < sq_rank; i++)
            {
                mat[ i + row*sq_rank ] ^= mat[ i + pivot_row*sq_rank ];
            }
        }

        gf16_t inverse = inv(mat[ pivot_col + row*sq_rank ]);

        for (size_t i = row+1; i < sq_rank; i++)
        {
            gf16_t scalar = gf16_get_mul(inverse, mat[ pivot_col + i*sq_rank ] );
            for (size_t j = pivot_col; j < sq_rank; j++)
            {
                mat[ j + i*sq_rank ] = gf16_get_add(mat[ j + i*sq_rank ], gf16_get_mul(scalar, mat[j + row*sq_rank]));
            }
        }

        row ++;        
        if (row > tres){
            return row;
        }
    }

    /*
    for (int i = 0; i < sq_rank; i++)
    {
        for (int j = 0; j < sq_rank; j++)
        {
            printf("%2d ", mat[j +i*sq_rank]);
        }
        printf("\n");
    }
    printf("\n");*/
    
    return row;
}

int compute_rank(const gf16_t *mat){
    gf16_t m[sq_rank*sq_rank];
    memcpy(m,mat, sizeof(gf16_t[sq_rank*sq_rank]));
    return ef(m,sq_rank);
}

int check_rank(const gf16_t *mat, int target){
    gf16_t m[sq_rank*sq_rank];
    memcpy(m,mat, sizeof(gf16_t[sq_rank*sq_rank]));
    return ef(m,target) <= target;
}

void combine_matrices(const gf16_t* mul_matrices, const gf16_t* coefs, gf16_t *out){
    size_t offset;
    // constant term
    memcpy(out, mul_matrices + sq_rank*sq_rank*sq_rank*sq_rank, sq_rank*sq_rank);

    for (int i = 0; i < sq_rank-rank; i++)
    {
        // linear terms
        offset = coefs[i]*sq_rank*sq_rank*sq_rank*sq_rank + (rank+i)*sq_rank*sq_rank;
        for (size_t k = 0; k < sq_rank*sq_rank; k++){
            out[k] ^= mul_matrices[offset + k ];
        }

        //quadratic terms
        for (int j = i; j < sq_rank-rank; j++)
        {
            gf16_t scalar = gf16_get_mul(coefs[i],coefs[j]);
            offset = scalar*sq_rank*sq_rank*sq_rank*sq_rank + ((rank+i)*sq_rank+rank+j)*sq_rank*sq_rank;
            for (size_t k = 0; k < sq_rank*sq_rank; k++){
                out[k] ^= mul_matrices[ offset + k ];
            }
        }
    }
}

void combine_matrices_q(const gf16_t* mul_matrices, const gf16_t* coefs, gf16_t *out){
    size_t offset;
    // constant term
    memcpy(out, mul_matrices + sq_rank*sq_rank*sq_rank*sq_rank, sq_rank*sq_rank);

    for (int i = 0; i < sq_rank-rank-1; i++)
    {
        // linear terms
        offset = coefs[i]*sq_rank*sq_rank*sq_rank*sq_rank + (rank+i)*sq_rank*sq_rank;
        for (size_t k = 0; k < sq_rank*sq_rank; k++){
            out[k] ^= mul_matrices[offset + k ];
        }

        //quadratic terms
        for (int j = i; j < sq_rank-rank-1; j++)
        {
            gf16_t scalar = gf16_get_mul(coefs[i],coefs[j]);
            offset = scalar*sq_rank*sq_rank*sq_rank*sq_rank + ((rank+i)*sq_rank+rank+j)*sq_rank*sq_rank;
            for (size_t k = 0; k < sq_rank*sq_rank; k++){
                out[k] ^= mul_matrices[ offset + k ];
            }
        }
    }

    for (size_t last_coef = 1; last_coef < 16; last_coef++)
    {
        memcpy(out + last_coef*sq_rank*sq_rank, out, sq_rank*sq_rank);

        offset = last_coef*sq_rank*sq_rank*sq_rank*sq_rank + (rank+sq_rank-rank-1)*sq_rank*sq_rank;
        for (size_t k = 0; k < sq_rank*sq_rank; k++){
            out[k + last_coef*sq_rank*sq_rank] ^= mul_matrices[offset + k ];
        }

        for (int i = 0; i < sq_rank-rank-1; i++)
        {
            gf16_t scalar = gf16_get_mul(coefs[i],last_coef);
            offset = scalar*sq_rank*sq_rank*sq_rank*sq_rank + ((rank+i)*sq_rank+rank+sq_rank-rank-1)*sq_rank*sq_rank;
            for (size_t k = 0; k < sq_rank*sq_rank; k++){
                out[k + last_coef*sq_rank*sq_rank] ^= mul_matrices[ offset + k ];
            }
        }

        gf16_t scalar = gf16_get_mul(last_coef,last_coef);
        offset = scalar*sq_rank*sq_rank*sq_rank*sq_rank + ((rank+sq_rank-rank-1)*sq_rank+rank+sq_rank-rank-1)*sq_rank*sq_rank;
        for (size_t k = 0; k < sq_rank*sq_rank; k++){
            out[k + last_coef*sq_rank*sq_rank] ^= mul_matrices[ offset + k ];
        }
    }
    
}


void combine_matrices_q_mat(const uint64_t* mul_matrices, const gf16_t* coefs, uint64_t *out){
    size_t offset;
    // constant term
    memcpy(out, mul_matrices + sq_rank*sq_rank*2*9, sizeof(uint64_t[2*9]));

    for (int i = 0; i < sq_rank-rank-1; i++)
    {
        // linear terms
        offset = coefs[i]*sq_rank*sq_rank*2*9 + (rank+i)*2*9;
        for (size_t k = 0; k < 2*9; k++){
            out[k] ^= mul_matrices[offset + k ];
        }

        //quadratic terms
        for (int j = i; j < sq_rank-rank-1; j++)
        {
            gf16_t scalar = gf16_get_mul(coefs[i],coefs[j]);
            offset = scalar*sq_rank*sq_rank*2*9 + ((rank+i)*sq_rank+rank+j)*2*9;
            for (size_t k = 0; k < 2*9; k++){
                out[k] ^= mul_matrices[ offset + k ];
            }
        }
    }

    for (size_t last_coef = 1; last_coef < 16; last_coef++)
    {
        memcpy(out + last_coef*2*9, out, sizeof(uint64_t[2*9]));

        offset = last_coef*sq_rank*sq_rank*2*9 + (rank+sq_rank-rank-1)*2*9;
        for (size_t k = 0; k < 2*9; k++){
            out[k + last_coef*2*9] ^= mul_matrices[offset + k ];
        }

        for (int i = 0; i < sq_rank-rank-1; i++)
        {
            gf16_t scalar = gf16_get_mul(coefs[i],last_coef);
            offset = scalar*sq_rank*sq_rank*2*9 + ((rank+i)*sq_rank+rank+sq_rank-rank-1)*2*9;
            for (size_t k = 0; k < 2*9; k++){
                out[k + last_coef*2*9] ^= mul_matrices[ offset + k ];
            }
        }

        gf16_t scalar = gf16_get_mul(last_coef,last_coef);
        offset = scalar*sq_rank*sq_rank*2*9 + ((rank+sq_rank-rank-1)*sq_rank+rank+sq_rank-rank-1)*2*9;
        for (size_t k = 0; k < 2*9; k++){
            out[k + last_coef*2*9] ^= mul_matrices[ offset + k ];
        }
    }
}

#define COEFS (sq_rank-rank)

int lowest_rank(uint8_t* pk){
    gf16_t matrices[sq_rank*sq_rank*sq_rank*sq_rank] = {0};

    compute_matrices(pk, matrices);
#if rank == 3 // optimization for rank 3
    uint64_t combined[16*2*9];
    uint64_t mul_matrices[16*sq_rank*sq_rank*2*9] = {0};

    for (size_t i = 0; i < sq_rank*sq_rank; i++)
    {
        load_mat(matrices + i*sq_rank*sq_rank, mul_matrices + i*2*9 + sq_rank*sq_rank*2*9);
    }

    for (size_t c = 2; c < 16; c++)
    {
        for (size_t i = 0; i < sq_rank*sq_rank*2*9; i++)
        {
            mul_matrices[(c*sq_rank*sq_rank)*2*9 + i] = mul_fx8(c,mul_matrices[(sq_rank*sq_rank)*2*9 + i]);
        }
    }
#else
    gf16_t combined[16*sq_rank*sq_rank];
    gf16_t mul_matrices[sq_rank*sq_rank*sq_rank*sq_rank*16] = {0};
    for (size_t i = 0; i < 16; i++)
    {
        for (size_t j = 0; j < sq_rank*sq_rank*sq_rank*sq_rank; j++)
        {
            mul_matrices[sq_rank*sq_rank*sq_rank*sq_rank*i + j] = gf16_get_mul(matrices[j],i);
        }
    }
#endif
    
    gf16_t coefs[COEFS-1];

    int lr = sq_rank;

    for (int c = 0; c < (1 << (4*(COEFS-1))); c++)
    {
        int cc = c;
        for (size_t i = 0; i < COEFS-1; i++)
        {
            coefs[i] = cc % 16;
            cc >>= 4;
        }

#if rank == 3
        combine_matrices_q_mat(mul_matrices,coefs,combined);
#else
        combine_matrices_q(mul_matrices,coefs,combined);
#endif
        for (int lc = 0; lc < 16; lc++)
        {
#if rank == 3
            int r = ef_mat(combined + lc*2*9,9);
#else
            int r = ef(combined + lc*sq_rank*sq_rank,9);
#endif
            if (r < lr){
                lr = r;
            }
        }
    }
    return lr;
}

void print_matrices(uint8_t* pk){

    gf16_t matrices[sq_rank*sq_rank*sq_rank*sq_rank] = {0};
    compute_matrices(pk, matrices);

    for (size_t a = 0; a < 4; a++)
    {
        for (size_t b = 0; b < 4; b++)
        {
            printf("%d ", matrices[a*4 + b]);
        }
        printf("\n");
    }
    printf("\n");
}

// write public key in format expected by our SAGE script
void print_public_key(uint8_t* pk){

    map_group1 map1;
    public_key* pk_stru = (public_key*)pk;
    gf16m_t P22[m_SNOVA][o_SNOVA][o_SNOVA];

    // generate PRNG part of public key
    gen_A_B_Q_P(&map1, pk_stru->pt_public_key_seed);
    // read  P22
    input_P22((uint8_t*)P22, pk_stru->P22);

    // print A_alpha
    printf("A: \n");
    for (int alpha = 0; alpha < sq_rank; alpha++)
    {
        for (size_t i = 0; i < sq_rank; i++)
        {
            printf("%d ", map1.Aalpha[alpha][i]);
        }
    }
    printf("\n");

    // print B_alpha
    printf("B: \n");
    for (int alpha = 0; alpha < sq_rank; alpha++)
    {
        for (size_t i = 0; i < sq_rank; i++)
        {
            printf("%d ", map1.Balpha[alpha][i]);
        }
    }
    printf("\n");

    // print Q_1alpha
    printf("Q1: \n");

    for (int alpha = 0; alpha < sq_rank; alpha++)
    {
        for (size_t i = 0; i < sq_rank; i++)
        {
            printf("%d ", map1.Qalpha1[alpha][i]);
        }
    }
    printf("\n");

    // print Q_2alpha
    printf("Q2: \n");

    for (int alpha = 0; alpha < sq_rank; alpha++)
    {
        for (size_t i = 0; i < sq_rank; i++)
        {
            printf("%d ", map1.Qalpha2[alpha][i]);
        }
    }
    printf("\n");

    // print P_i each on a separate line
    printf("P: \n");
    for (size_t i = 0; i < m_SNOVA; i++)
    {
        // P11
        for (size_t row = 0; row < v_SNOVA; row++)
        {
            for ( size_t col = 0; col < v_SNOVA; col ++){
                for (size_t j = 0; j < sq_rank; j++){
                    printf("%d ", map1.P11[i][row][col][j]);
                }
            } 
        }

        // P12
        for (size_t row = 0; row < v_SNOVA; row++)
        {
            for ( size_t col = 0; col < m_SNOVA; col ++){
                for (size_t j = 0; j < sq_rank; j++){
                    printf("%d ", map1.P12[i][row][col][j]);
                }
            } 
        }

        // P21
        for (size_t row = 0; row < o_SNOVA; row++)
        {
            for ( size_t col = 0; col < v_SNOVA; col ++){
                for (size_t j = 0; j < sq_rank; j++){
                    printf("%d ", map1.P21[i][row][col][j]);
                }
            } 
        }

        // P22
        for (size_t row = 0; row < o_SNOVA; row++)
        {
            for ( size_t col = 0; col < o_SNOVA; col ++){
                for (size_t j = 0; j < sq_rank; j++){
                    printf("%d ", P22[i][row][col][j]);
                }
            } 
        }
        printf("\n");
    }

    // print input-output pair
    print_input_output_pair(pk);
    
}

void print_input_output_pair(uint8_t *pk){

    gf16m_t Left[lsq_SNOVA][n_SNOVA], Right[lsq_SNOVA][n_SNOVA];
    gf16m_t temp1, temp2;
    map_group1 map1;
    public_key* pk_stru = (public_key*)pk;
    gf16m_t P22[m_SNOVA][o_SNOVA][o_SNOVA];

    // generate PRNG part of public key
    gen_A_B_Q_P(&map1, pk_stru->pt_public_key_seed);
    // read  P22
    input_P22((uint8_t*)P22, pk_stru->P22);

    gf16m_t hash_in_GF16Matrix[m_SNOVA];
    gf16m_t signature_in_GF16Matrix[n_SNOVA];

    gf16_t *input = (gf16_t*) &signature_in_GF16Matrix;
    gf16_t *output =(gf16_t*) &hash_in_GF16Matrix;

    // set input
    for (size_t i = 0; i < sq_rank*n_SNOVA; i++)
    {
        input[i] = rand()%16;
        printf("%d ", input[i]);
    }
    printf("\n");
    
    // compute output
    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
        for (int index = 0; index < n_SNOVA; ++index) {
            // Left[alpha][index]
            gf16m_transpose(signature_in_GF16Matrix[index], temp1);
            gf16m_mul(temp1, map1.Qalpha1[alpha], temp2);
            gf16m_mul(map1.Aalpha[alpha], temp2, Left[alpha][index]);
            // Right[alpha][index]
            gf16m_mul(map1.Qalpha2[alpha], signature_in_GF16Matrix[index], temp2);
            gf16m_mul(temp2, map1.Balpha[alpha], Right[alpha][index]);

            /*
            Left[alpha][index] = Aalpha[alpha] *
            (signature_in_GF16Matrix[index].transpose()) * Qalpha1[alpha];
            Right[alpha][index] = Qalpha2[alpha] *
            signature_in_GF16Matrix[index] * Balpha[alpha];
            */
        }
    }

    for (int i = 0; i < m_SNOVA; ++i) {
        // gf16m_add(hash_in_GF16Matrix[i], hash_in_GF16Matrix[i],
        // hash_in_GF16Matrix[i]);
        gf16m_set_zero(hash_in_GF16Matrix[i]);
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
            for (int dj = 0; dj < v_SNOVA; ++dj) {
                for (int dk = 0; dk < v_SNOVA; ++dk) {
                    gf16m_mul(Left[alpha][dj], map1.P11[i][dj][dk], temp1);
                    gf16m_mul(temp1, Right[alpha][dk], temp2);
                    gf16m_add(hash_in_GF16Matrix[i], temp2, hash_in_GF16Matrix[i]);
                    // hash_in_GF16Matrix[i] = hash_in_GF16Matrix[i] +
                    // Left[alpha][dj] * P11[i][dj][dk] * Right[alpha][dk];
                }
            }
            for (int dj = 0; dj < v_SNOVA; ++dj) {
                for (int dk = 0; dk < o_SNOVA; ++dk) {
                    gf16m_mul(Left[alpha][dj], map1.P12[i][dj][dk], temp1);
                    gf16m_mul(temp1, Right[alpha][dk + v_SNOVA], temp2);
                    gf16m_add(hash_in_GF16Matrix[i], temp2, hash_in_GF16Matrix[i]);
                    // hash_in_GF16Matrix[i] = hash_in_GF16Matrix[i] +
                    // Left[alpha][dj] * P12[i][dj][dk] *
                    // Right[alpha][dk+v_SNOVA];
                }
            }
            for (int dj = 0; dj < o_SNOVA; ++dj) {
                for (int dk = 0; dk < v_SNOVA; ++dk) {
                    gf16m_mul(Left[alpha][dj + v_SNOVA], map1.P21[i][dj][dk], temp1);
                    gf16m_mul(temp1, Right[alpha][dk], temp2);
                    gf16m_add(hash_in_GF16Matrix[i], temp2, hash_in_GF16Matrix[i]);
                    // hash_in_GF16Matrix[i] = hash_in_GF16Matrix[i] +
                    // Left[alpha][dj + v_SNOVA] * P21[i][dj][dk] *
                    // Right[alpha][dk];
                }
            }
            for (int dj = 0; dj < o_SNOVA; ++dj) {
                for (int dk = 0; dk < o_SNOVA; ++dk) {
                    gf16m_mul(Left[alpha][dj + v_SNOVA], P22[i][dj][dk], temp1);
                    gf16m_mul(temp1, Right[alpha][dk + v_SNOVA], temp2);
                    gf16m_add(hash_in_GF16Matrix[i], temp2, hash_in_GF16Matrix[i]);
                    // hash_in_GF16Matrix[i] = hash_in_GF16Matrix[i] +
                    // Left[alpha][dj + v_SNOVA] * P22[i][dj][dk] *
                    // Right[alpha][dk + v_SNOVA];
                }
            }
        }
    }

    // print output
    for (size_t i = 0; i < sq_rank*m_SNOVA; i++)
    {
        printf("%d ", output[i]);
    }
    printf("\n");
}