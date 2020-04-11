#include <OWL/mxf32_2x2.h>

//
//
//
float owl_mxf32_2x2_get_element(owl_mxf32_2x2* M, int i, int j)
{
    int i_ = i % 2;
    int j_ = j % 2;

    float const* flat_matrix = (float const*)M;
    return flat_matrix[2 * j_ + i_];
}

//
//
//
owl_mxf32_2x2* owl_mxf32_2x2_set_element(owl_mxf32_2x2* M, float value, int i, int j)
{
    int i_ = i % 2;
    int j_ = j % 2;

    float* flat_matrix = (float*)M;
    flat_matrix[2 * j_ + i_] = value;
    *M = _mm_load_ps(flat_matrix);

    return M;
}

//M = A ^(-1)
//
//
owl_mxf32_2x2* owl_mxf32_2x2_Inv(owl_mxf32_2x2* M, owl_mxf32_2x2 const* A)
{
    __m128 broadcast_inv_det = _mm_set1_ps( 1.0 / owl_mxf32_2x2_det(A) );

    __m128 tmp = _mm_shuffle_ps(*A, *A, 0b11100001);
    tmp = _mm_addsub_ps(
                            _mm_setzero_ps(),
                            _mm_mul_ps(broadcast_inv_det, tmp)
                        );
    *M = _mm_shuffle_ps(tmp, tmp, 0b01100011);

    return M;
}

//A = P * D * tP with A symmetric
//A is considered invertible
//Parameter P is optional
//Return D
owl_mxf32_2x2* owl_mxf32_2x2_diagonalize_sym(owl_mxf32_2x2* D, owl_mxf32_2x2* P, owl_mxf32_2x2 const* A)
{
    if(owl_mxf32_2x2_norminf(A) == 0.0)
    {
        owl_mxf32_2x2_zero(D);
        owl_mxf32_2x2_diag(P, 1.0);
    }
    else
    {
        float vp_list[4] OWL_ALIGN16;
        {
            float flat_A[4] OWL_ALIGN16;
            _mm_store_ps(flat_A, A);

            float delta = (flat_A[0] - flat_A[3]) * (flat_A[0] - flat_A[3]) + 4 * flat_A[1] * flat_A[1];
            vp_list[0] = 0.5f * (flat_A[0] + flat_A[3]) + sqrtf(delta);
            vp_list[1] = 0.5f * (flat_A[0] + flat_A[3]) - sqrtf(delta);
        }
    }

    return D;
}
