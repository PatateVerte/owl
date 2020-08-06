#include <OWL/Optimized3d/vector/v3f32.h>

//Infinity norm
//
//
float OWL_VECTORCALL owl_v3f32_norminf(owl_v3f32 v)
{
    __m128 v_abs = _mm_max_ps(
                                v,
                                _mm_sub_ps(_mm_setzero_ps(), v)
                              );

    __m128 tmp = v_abs;
    tmp = _mm_max_ss(tmp, _mm_insert_ps(v_abs, v_abs, 0b01001110));
    tmp = _mm_max_ss(tmp, _mm_insert_ps(v_abs, v_abs, 0b10001110));

    return _mm_cvtss_f32(tmp);
}

//Sign mask in int[2:0]
//0 : < 0
//1 : >= 0
int OWL_VECTORCALL owl_v3f32_sign_mask(owl_v3f32 v)
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
