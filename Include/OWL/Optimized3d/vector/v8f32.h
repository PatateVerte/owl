#ifndef OWL_V8F32_H_INCLUDED
#define OWL_V8F32_H_INCLUDED

#include <OWL/owl.h>

#include <malloc.h>

//256-bits AVX
#include <immintrin.h>

typedef __m256 owl_v8f32;

//null vector
static inline owl_v8f32 owl_v8f32_zero()
{
    return _mm256_setzero_ps();
}

//Vector
static inline owl_v8f32 owl_v8f32_set(float x0, float x1, float x2, float x3, float x4, float x5, float x6, float x7)
{
    return _mm256_set_ps(x7, x6, x5, x4, x3, x2, x1, x0);
}

//Allocate an owl_v3f32-aligned array with 8*len float
static inline float* owl_v8f32_array_alloc(size_t len)
{
    return _aligned_malloc(len * sizeof(owl_v8f32), 32);
}

static inline void owl_v8f32_array_free(float* ptr)
{
    _aligned_free(ptr);
}

//Vector from a float[8]
static inline owl_v8f32 owl_v8f32_load8(const float* src)
{
    return _mm256_loadu_ps(src);
}

//Vector to float[8]
static inline float* owl_v8f32_store8(float* dst, owl_v8f32 v)
{
    _mm256_storeu_ps(dst, v);
    return dst;
}

//Broadcast f within a vector
static inline owl_v8f32 owl_v8f32_broadcast(float f)
{
    return _mm256_broadcastss_ps(_mm_set_ss(f));
}

//-a
static inline owl_v8f32 owl_v8f32_negate(owl_v8f32 a)
{
    return _mm256_sub_ps(_mm256_setzero_ps(), a);
}

//a + b
static inline owl_v8f32 owl_v8f32_add(owl_v8f32 a, owl_v8f32 b)
{
    return _mm256_add_ps(a, b);
}

//a - b
static inline owl_v8f32 owl_v8f32_sub(owl_v8f32 a, owl_v8f32 b)
{
    return _mm256_sub_ps(a, b);
}

//v1 + v2 * v3
static inline owl_v8f32 owl_v8f32_addmul(owl_v8f32 v1, owl_v8f32 v2, owl_v8f32 v3)
{
    return _mm256_fmadd_ps(v2, v3, v1);
}

//a * v
static inline owl_v8f32 owl_v8f32_scalar_mul(owl_v8f32 v, float a)
{
    return _mm256_mul_ps(v, owl_v8f32_broadcast(a));
}

//v1 + a*v2
static inline owl_v8f32 owl_v8f32_add_scalar_mul(owl_v8f32 v1, owl_v8f32 v2, float a)
{
    return _mm256_fmadd_ps(owl_v8f32_broadcast(a), v2, v1);
}

//Component by component mul
static inline owl_v8f32 owl_v8f32_mul(owl_v8f32 v1, owl_v8f32 v2)
{
    return _mm256_mul_ps(v1, v2);
}

//Component by component div
static inline owl_v8f32 owl_v8f32_comp_div(owl_v8f32 v1, owl_v8f32 v2)
{
    return _mm256_div_ps(v1, v2);
}

//v / a
static inline owl_v8f32 owl_v8f32_scalar_div(owl_v8f32 v, float a)
{
    return _mm256_div_ps(v, owl_v8f32_broadcast(a));
}

#endif // OWL_V8F32_H_INCLUDED
