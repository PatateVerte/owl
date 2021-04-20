#ifndef OWL_Q32_H_INCLUDED
#define OWL_Q32_H_INCLUDED

#include <OWL/owl.h>
#include <OWL/Optimized3d/vector/v3f32.h>
#include <OWL/Optimized3d/matrix/mxf32_3x3.h>

#include <xmmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>

typedef __m128 owl_q32;

static inline owl_q32 owl_q32_zero()
{
    return _mm_setzero_ps();
}

//Sets a quaternion (w = Real part, (x, y, z) = vect part)
static inline owl_q32 owl_q32_set(float w, float x, float y, float z)
{
    return _mm_set_ps(z, y, x, w);
}

//q = w
static inline owl_q32 owl_q32_from_real(float w)
{
    return _mm_set_ss(w);
}

//Load quaternion
static inline owl_q32 owl_q32_load4(float const* src)
{
    return _mm_loadu_ps(src);
}

//Store quaternion
static inline float* owl_q32_store4(float* dst, owl_q32 q)
{
    _mm_storeu_ps(dst, q);
    return dst;
}

//Re(q) as a float
static inline float owl_q32_Ref(owl_q32 q)
{
    return _mm_cvtss_f32(q);
}

//Im(q) as a vector
static inline owl_v3f32 owl_q32_Imv(owl_q32 q)
{
    owl_v3f32 v = _mm_shuffle_ps(q, q, 0b00111001);
    return _mm_insert_ps(v, v, 0b00001000);
}

//Re(q) as a q32
static inline owl_q32 owl_q32_Re(owl_q32 q)
{
    return _mm_insert_ps(q, q, 0b00001110);
}

//Im(q)
static inline owl_q32 owl_q32_Im(owl_q32 q)
{
    return _mm_insert_ps(q, q, 0b00000001);
}

//|q1|^2
static inline float owl_q32_square_mod(owl_q32 q)
{
    return _mm_cvtss_f32( _mm_dp_ps(q, q, 0b11110001) );
}

//|q1|
static inline float owl_q32_mod(owl_q32 q)
{
    return _mm_cvtss_f32( _mm_sqrt_ss( _mm_dp_ps(q, q, 0b11110001) ) );
}

//q1 + q2
static inline owl_q32 owl_q32_add(owl_q32 q1, owl_q32 q2)
{
    return _mm_add_ps(q1, q2);
}

//q1 - q2
static inline owl_q32 owl_q32_sub(owl_q32 q1, owl_q32 q2)
{
    return _mm_sub_ps(q1, q2);
}

//Conj(q)
static inline owl_q32 owl_q32_conj(owl_q32 q)
{
    return owl_q32_sub(
                            owl_q32_Re(q),
                            owl_q32_Im(q)
                         );
}

//a * q
static inline owl_q32 owl_q32_real_mul(owl_q32 q, float a)
{
    return _mm_mul_ps(q, _mm_set1_ps(a));
}

//q / a
static inline owl_q32 owl_q32_real_div(owl_q32 q, float a)
{
    return _mm_div_ps(q, _mm_set1_ps(a));
}

//Inv(q)
static inline owl_q32 owl_q32_inv(owl_q32 q)
{
    return owl_q32_real_div(
                                owl_q32_conj(q),
                                owl_q32_square_mod(q)
                              );
}

//q1 * q2
OWL_DLL_IMPORT owl_q32 OWL_VECTORCALL owl_q32_mul(owl_q32 q1, owl_q32 q2);
OWL_DLL_IMPORT owl_q32 owl_q32_extern_mul(owl_q32 q1, owl_q32 q2);

//||v|| = 1
//Rotation of alpha around the unitary vector v
//||v|| can be null if alpha == 0
OWL_DLL_IMPORT owl_q32 OWL_VECTORCALL owl_q32_from_rotation(owl_v3f32 v, float alpha);
OWL_DLL_IMPORT owl_q32 owl_q32_extern_from_rotation(owl_v3f32 v, float alpha);

//|q| = 1
//(q) * u * (q^-1)
OWL_DLL_IMPORT owl_v3f32 OWL_VECTORCALL owl_q32_transform_v3f32(owl_q32 q, owl_v3f32 u);
OWL_DLL_IMPORT owl_v3f32 owl_q32_extern_transform_v3f32(owl_q32 q, owl_v3f32 u);

//Transform a rotation matrix to a quaternion
OWL_DLL_IMPORT owl_q32 owl_q32_from_rotation_matrix(owl_mxf32_3x3 const* O);

#endif // OWL_Q32_H_INCLUDED
