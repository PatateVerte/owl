#ifndef OWL_V3F32_H_INCLUDED
#define OWL_V3F32_H_INCLUDED

#include <OWL/owl.h>

#include <math.h>

//SSE / 128-bits AVX
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>

typedef __m128 owl_v3f32;

//null vector
static inline owl_v3f32 owl_v3f32_zero()
{
    return _mm_setzero_ps();
}

//(x, y, z) vector
static inline owl_v3f32 owl_v3f32_set(float x, float y, float z)
{
    return _mm_set_ps(0.0, z, y, x);
}

//Vector from a float[4]
static inline owl_v3f32 owl_v3f32_load4(const float* src)
{
    return _mm_load_ps(src);
}

//Vector to float[4]
static inline float* owl_v3f32_store4(float* dst, owl_v3f32 v)
{
    _mm_store_ps(dst, v);
    return dst;
}

//Fill base[3] with the basic trihedral multiplied by f
static inline owl_v3f32* owl_v3f32_base_xyz(owl_v3f32* base, float f)
{
    __m128 tmp = _mm_set_ss(f);
    base[0] = _mm_insert_ps(tmp, tmp, 0b00001110);
    base[1] = _mm_insert_ps(tmp, tmp, 0b00011101);
    base[2] = _mm_insert_ps(tmp, tmp, 0b00101011);

    return base;
}

//Broadcast f within a vector
static inline owl_v3f32 owl_v3f32_broadcast(float f)
{
    __m128 tmp = _mm_set1_ps(f);
    return _mm_insert_ps(tmp, tmp, 0b00001000);
}

//a + b
static inline owl_v3f32 owl_v3f32_add(owl_v3f32 a, owl_v3f32 b)
{
    return _mm_add_ps(a, b);
}

//a - b
static inline owl_v3f32 owl_v3f32_sub(owl_v3f32 a, owl_v3f32 b)
{
    return _mm_sub_ps(a, b);
}

//
#define owl_v3f32_set_component(v, i, val) \
    _mm_insert_ps( (v) , _mm_set_ss((val)) , (0b1000 | (0b00 << 6) | ((i) << 4)) )


#define owl_v3f32_get_component(v, i) \
    _mm_cvtss_f32( _mm_insert_ps((v) , (v) , (0b1110 | ((i) << 6) | (0b00 << 4)) ) )

//a * v
static inline owl_v3f32 owl_v3f32_scalar_mul(owl_v3f32 v, float a)
{
    return _mm_mul_ps(v, _mm_set1_ps(a));
}

//v / a
static inline owl_v3f32 owl_v3f32_scalar_div(owl_v3f32 v, float a)
{
    return _mm_div_ps(v, _mm_set1_ps(a));
}

//v1 + a*v2
static inline owl_v3f32 owl_v3f32_add_scalar_mul(owl_v3f32 v1, owl_v3f32 v2, float a)
{
    return _mm_add_ps(
                        v1,
                        _mm_mul_ps( _mm_set1_ps(a), v2 )
                      );
}

//a.b
static inline float owl_v3f32_dot(owl_v3f32 a, owl_v3f32 b)
{
    return _mm_cvtss_f32( _mm_dp_ps(a, b, 0b01110001) );
}

//||a||
static inline float owl_v3f32_norm(owl_v3f32 a)
{
    return _mm_cvtss_f32( _mm_sqrt_ss( _mm_dp_ps(a, a, 0b01110001) ) );
}

//
static inline owl_v3f32 owl_v3f32_normalize(owl_v3f32 a)
{
    return owl_v3f32_scalar_div(a, owl_v3f32_norm(a) );
}

//a^b
static inline owl_v3f32 owl_v3f32_cross(owl_v3f32 a, owl_v3f32 b)
{
    __m128 tmp = _mm_mul_ps(
                                _mm_shuffle_ps(a, a, 0b11001001),
                                _mm_shuffle_ps(b, b, 0b11010010)
                            );

    tmp = _mm_sub_ps(
                     tmp,
                     _mm_mul_ps(
                                    _mm_shuffle_ps(a, a, 0b11010010),
                                    _mm_shuffle_ps(b, b, 0b11001001)
                                )
                    );

    return tmp;
}

//[a, b, c]
static inline float owl_v3f32_triple(owl_v3f32 a, owl_v3f32 b, owl_v3f32 c)
{
    return owl_v3f32_dot(owl_v3f32_cross(a, b), c);
}

//
static inline owl_v3f32 owl_v3f32_rotate_comp(owl_v3f32 v)
{
    return _mm_shuffle_ps(v, v, 0b11010010);
}

//Infinity norm
static inline float owl_v3f32_inf_norm(owl_v3f32 v)
{
    __m128 v_abs = _mm_max_ps(
                                v,
                                _mm_sub_ps(_mm_setzero_ps(), v)
                              );

    __m128 tmp = _mm_insert_ps(v_abs, v_abs, 0b00001110);
    tmp = _mm_max_ss(tmp, _mm_insert_ps(v_abs, v_abs, 0b01001110));
    tmp = _mm_max_ss(tmp, _mm_insert_ps(v_abs, v_abs, 0b10001110));

    return _mm_cvtss_f32(tmp);
}

//Sign mask in int[2:0]
//0 : < 0
//1 : >= 0
static inline int owl_v3f32_sign_mask(owl_v3f32 v)
{
    __m128i tmp = _mm_castps_si128( _mm_cmpge_ps(v, _mm_setzero_ps()) );
    tmp = _mm_packs_epi32(tmp, tmp);
    tmp = _mm_packs_epi16(tmp, tmp);

    int tmp_mask = _mm_cvtsi128_si32(tmp);
    tmp_mask &= (1<<0) | (1<<8) | (1<<16);
    int sign_mask = (tmp_mask>>0) | (tmp_mask >> 7) | (tmp_mask >> 14);
    sign_mask &= 0b111;

    return sign_mask;
}

#endif // OWL_V3F32_H_INCLUDED
