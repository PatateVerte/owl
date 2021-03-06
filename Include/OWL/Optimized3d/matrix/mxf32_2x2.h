#ifndef OWL_MXF32_2X2_H_INCLUDED
#define OWL_MXF32_2X2_H_INCLUDED

#include <OWL/owl.h>

//SSE / 128-bits AVX
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>

//store as followed : (a00, a10, a01, a11)
typedef __m128 owl_mxf32_2x2;

//M = 0
//Return M
static inline owl_mxf32_2x2* owl_mxf32_2x2_zero(owl_mxf32_2x2* M)
{
	*M = _mm_setzero_ps();
	return M;
}

//
static inline owl_mxf32_2x2* owl_mxf32_2x2_diag(owl_mxf32_2x2* M, float diag_val)
{
    __m128 tmp = _mm_set_ss(diag_val);
    tmp = _mm_shuffle_ps(tmp, tmp, 0b00111100);
    *M = tmp;

    return M;
}

//M is given by column
static inline owl_mxf32_2x2* owl_mxf32_2x2_set(owl_mxf32_2x2* M, float a00, float a10, float a01, float a11)
{
    *M = _mm_set_ps(a11, a01, a10, a00);
    return M;
}

//
static inline float* owl_mxf32_2x2_store(float* dst, owl_mxf32_2x2 const* M)
{
    _mm_storeu_ps(dst, *M);
    return dst;
}

//
static inline owl_mxf32_2x2* owl_mxf32_2x2_load(owl_mxf32_2x2* M, float const* src)
{
    *M = _mm_loadu_ps(src);
    return M;
}

//
OWL_DLL_EXPORT float owl_mxf32_2x2_get_element(owl_mxf32_2x2* M, unsigned int i, unsigned int j);

//
OWL_DLL_EXPORT owl_mxf32_2x2* owl_mxf32_2x2_set_element(owl_mxf32_2x2* M, float value, unsigned int i, unsigned int j);

//M = A + B
static inline owl_mxf32_2x2* owl_mxf32_2x2_add(owl_mxf32_2x2* M, owl_mxf32_2x2 const* A, owl_mxf32_2x2 const* B)
{
    *M = _mm_add_ps(*A, *B);
    return M;
}

//M = A - B
static inline owl_mxf32_2x2* owl_mxf32_2x2_sub(owl_mxf32_2x2* M, owl_mxf32_2x2 const* A, owl_mxf32_2x2 const* B)
{
    *M = _mm_sub_ps(*A, *B);
    return M;
}

//M = a * A
static inline owl_mxf32_2x2* owl_mxf32_2x2_scalar_mul(owl_mxf32_2x2* M, owl_mxf32_2x2 const* A, float a)
{
    __m128 broadcast = _mm_set1_ps(a);
    *M = _mm_mul_ps(*A, broadcast);
    return M;
}

//M = tA
static inline owl_mxf32_2x2* owl_mxf32_2x2_transp(owl_mxf32_2x2* M, owl_mxf32_2x2 const* A)
{
    *M = _mm_shuffle_ps(*A, *A, 0b11011000);
    return M;
}


//M = A * B
static inline owl_mxf32_2x2* owl_mxf32_2x2_mul(owl_mxf32_2x2* M, owl_mxf32_2x2 const* A, owl_mxf32_2x2 const* B)
{
    __m128 C1 = _mm_shuffle_ps(*A, *A, 0b01000100);
    __m128 tmp = _mm_mul_ps(
                                C1,
                                _mm_shuffle_ps(*B, *B, 0b10100000)
                            );

    __m128 C2 = _mm_shuffle_ps(*A, *A, 0b11101110);
    tmp = _mm_add_ps(
                        tmp,
                        _mm_mul_ps(
                                    C2,
                                    _mm_shuffle_ps(*B, *B, 0b11110101)
                                   )
                     );
    *M = tmp;
    return M;
}

//M = A
static inline owl_mxf32_2x2* owl_mxf32_2x2_copy(owl_mxf32_2x2* M, owl_mxf32_2x2 const* A)
{
    *M = *A;
    return M;
}

//norm2(A)
static inline float owl_mxf32_2x2_norm2(owl_mxf32_2x2 const* A)
{
    return _mm_cvtss_f32( _mm_sqrt_ss( _mm_dp_ps(*A, *A, 0b11110001) ) );
}

//norm_inf(A)
static inline float owl_mxf32_2x2_norminf(owl_mxf32_2x2 const* A)
{
    __m128 abs_A = _mm_max_ps(*A, _mm_sub_ps(_mm_setzero_ps(), *A));

    __m128 abs_A_C2 = _mm_shuffle_ps(abs_A, abs_A, 0b00001110);
    __m128 tmp = _mm_max_ps(abs_A, abs_A_C2);

    return _mm_cvtss_f32(
                            _mm_max_ss(
                                        tmp,
                                        _mm_insert_ps(tmp, tmp, 0b01001110)
                                       )
                         );
}

//Tr(A)
static inline float owl_mxf32_2x2_trace(owl_mxf32_2x2 const* A)
{
    return _mm_cvtss_f32(
                            _mm_add_ps(
                                        *A,
                                        _mm_insert_ps(*A, *A, 0b11001110)
                                       )
                         );
}

//det(A)
static inline float owl_mxf32_2x2_det(owl_mxf32_2x2 const* A)
{
    float det = _mm_cvtss_f32(*A) * _mm_cvtss_f32( _mm_insert_ps(*A, *A, 0b11001110) );
    det -= _mm_cvtss_f32( _mm_insert_ps(*A, *A, 0b01001110) ) * _mm_cvtss_f32( _mm_insert_ps(*A, *A, 0b10001110) );
    return det;
}

//M = A ^(-1)
OWL_DLL_EXPORT owl_mxf32_2x2* owl_mxf32_2x2_Inv(owl_mxf32_2x2* M, owl_mxf32_2x2 const* A);

//A = P * D * tP with A symmetric and D=diag(eigenvalue_list)
//Parameter P is optional
//P is an orthogonal direct base of R^2
//Return eigenvalue_list
OWL_DLL_EXPORT float* owl_mxf32_2x2_diagonalize_sym(float* eigenvalue_list, owl_mxf32_2x2* P, owl_mxf32_2x2 const* A);

#endif // OWL_MXF32_2X2_H_INCLUDED
