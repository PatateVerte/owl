#include <OWL/Optimized3d/matrix/mxf32_2x2.h>

#include <OWL/owl.h>

#include <math.h>

//
//
//
OWL_DLL_EXPORT float owl_mxf32_2x2_get_element(owl_mxf32_2x2* M, unsigned int i, unsigned int j)
{
    unsigned int i_ = i % 2;
    unsigned int j_ = j % 2;

    float const* flat_matrix = (float const*)M;
    return flat_matrix[2 * j_ + i_];
}

//
//
//
OWL_DLL_EXPORT owl_mxf32_2x2* owl_mxf32_2x2_set_element(owl_mxf32_2x2* M, float value, unsigned int i, unsigned int j)
{
    unsigned int i_ = i % 2;
    unsigned int j_ = j % 2;

    float* flat_matrix = (float*)M;
    flat_matrix[2 * j_ + i_] = value;
    *M = _mm_load_ps(flat_matrix);

    return M;
}

//M = A ^(-1)
//
//
OWL_DLL_EXPORT owl_mxf32_2x2* owl_mxf32_2x2_Inv(owl_mxf32_2x2* M, owl_mxf32_2x2 const* A)
{
    __m128 broadcast_inv_det = _mm_set1_ps( 1.0f / owl_mxf32_2x2_det(A) );

    __m128 tmp = _mm_shuffle_ps(*A, *A, 0b11100001);
    tmp = _mm_addsub_ps(
                            _mm_setzero_ps(),
                            _mm_mul_ps(broadcast_inv_det, tmp)
                        );
    *M = _mm_shuffle_ps(tmp, tmp, 0b01100011);

    return M;
}

//A = P * D * tP with A symmetric and D=diag(eigenvalue_list)
//Parameter P is optional
//Return eigenvalue_list
OWL_DLL_EXPORT float* owl_mxf32_2x2_diagonalize_sym(float* eigenvalue_list, owl_mxf32_2x2* P, owl_mxf32_2x2 const* A)
{
    //Detection of diagonal matrix
    float const non_diag_term = _mm_cvtss_f32( _mm_insert_ps(*A, *A, 0b01001110) );

    if(non_diag_term == 0.0f)
    {
        eigenvalue_list[0] = _mm_cvtss_f32(_mm_insert_ps(*A, *A, 0b00001110));
        eigenvalue_list[1] = _mm_cvtss_f32(_mm_insert_ps(*A, *A, 0b11001110));

        if(P != NULL)
        {
            owl_mxf32_2x2_diag(P, 1.0);
        }
    }
    else
    {
        float dominant_eigenvalue;
        {
            float OWL_ALIGN16 flat_A[4];
            owl_mxf32_2x2_store(flat_A, A);

            float const delta = (flat_A[0] - flat_A[3]) * (flat_A[0] - flat_A[3]) + 4 * flat_A[1] * flat_A[1];

            if(flat_A[0] + flat_A[3] >= 0.0f)
            {
                dominant_eigenvalue = 0.5f * ((flat_A[0] + flat_A[3]) + sqrtf(delta));
            }
            else
            {
                dominant_eigenvalue = 0.5f * ((flat_A[0] + flat_A[3]) - sqrtf(delta));
            }
        }

        owl_mxf32_2x2 H;
        owl_mxf32_2x2_diag(&H, dominant_eigenvalue);
        owl_mxf32_2x2_sub(&H, A, &H);

        float a, b;

        float n1 = _mm_cvtss_f32(_mm_dp_ps(H, H, 0b01010001));
        float n2 = _mm_cvtss_f32(_mm_dp_ps(H, H, 0b10100001));
        float OWL_ALIGN16 flat_H[4];
        owl_mxf32_2x2_store(flat_H, &H);
        if(n1 >= n2)
        {
            a = flat_H[0];
            b = flat_H[2];
        }
        else
        {
            a = flat_H[1];
            b = flat_H[3];
        }

        float x, y;
        if(fabsf(a) < fabsf(b))
        {
            x = 1.0f;
            y = - a / b;
        }
        else
        {
            y = 1.0f;
            x = - b / a;
        }

        float inv_norm = 1.0f / sqrtf(x*x + y*y);
        x *= inv_norm;
        y *= inv_norm;

        owl_mxf32_2x2 P_;
        owl_mxf32_2x2_set(&P_, x, y, -y, x);

        owl_mxf32_2x2 M;
        owl_mxf32_2x2_mul(&M, A, &P_);

        eigenvalue_list[0] = dominant_eigenvalue;
        eigenvalue_list[1] = _mm_cvtss_f32(_mm_dp_ps(M, P_, 0b11000001));

        if(P != NULL)
        {
            *P = P_;
        }
    }

    return eigenvalue_list;
}
