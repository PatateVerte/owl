#include <OWL/mxf32_3x3.h>
#include <OWL/mxf32_2x2.h>

#include <math.h>

//
//
//
owl_mxf32_3x3* owl_mxf32_3x3_diag(owl_mxf32_3x3* M, float diag_val)
{
    owl_v3f32_base_xyz(M->column, diag_val);

    return M;
}

//
//
//
float owl_mxf32_3x3_get_element(owl_mxf32_3x3* M, unsigned int i, unsigned int j)
{
    unsigned int i_ = i % 3;
    unsigned int j_ = j % 3;

    float flat_column[4] OWL_ALIGN16;
    owl_v3f32_store4(flat_column, M->column[j_]);
    return flat_column[i_];
}

//
//
//
owl_mxf32_3x3* owl_mxf32_3x3_set_element(owl_mxf32_3x3* M, float value, unsigned int i, unsigned int j)
{
    unsigned int i_ = i % 3;
    unsigned int j_ = j % 3;

    float flat_column[4] OWL_ALIGN16;
    owl_v3f32_store4(flat_column, M->column[j_]);
    flat_column[i_] = value;
    M->column[j_] = owl_v3f32_load4(flat_column);

    return M;
}

//M = A + B
//
//
owl_mxf32_3x3* owl_mxf32_3x3_add(owl_mxf32_3x3* M, owl_mxf32_3x3 const* A, owl_mxf32_3x3 const* B)
{
    for(int j = 0 ; j < 3 ; j++)
    {
        M->column[j] = owl_v3f32_add(A->column[j], B->column[j]);
    }

    return M;
}

//M = A - B
//
//
owl_mxf32_3x3* owl_mxf32_3x3_sub(owl_mxf32_3x3* M, owl_mxf32_3x3 const* A, owl_mxf32_3x3 const* B)
{
    for(int j = 0 ; j < 3 ; j++)
    {
        M->column[j] = owl_v3f32_sub(A->column[j], B->column[j]);
    }

    return M;
}

//
//
//
owl_mxf32_3x3* owl_mxf32_3x3_scalar_mul(owl_mxf32_3x3* M, owl_mxf32_3x3 const* A, float a)
{
    for(int j = 0 ; j < 3 ; j++)
    {
        M->column[j] = owl_v3f32_scalar_mul(A->column[j], a);
    }

    return M;
}

//M = tA
//
//
owl_mxf32_3x3* owl_mxf32_3x3_transp(owl_mxf32_3x3* M, owl_mxf32_3x3 const* A)
{
    float const* base_index = (float const*) &(A->column);
    __m128i ind = _mm_set_epi32(3, 8, 4, 0);
    __m128i one = _mm_set_epi32(0, 1, 1, 1);

    for(int j = 0 ; j < 3 ; j++)
    {
        M->column[j] = _mm_i32gather_ps(base_index, ind, sizeof(*base_index));
        ind = _mm_add_epi32(ind, one);
    }

    return M;
}

//A * v
//
//
owl_v3f32 owl_mxf32_3x3_transform_v3f32(owl_mxf32_3x3 const* A, owl_v3f32 v)
{
    owl_v3f32 vr = owl_v3f32_zero();

    vr = owl_v3f32_filter_apply(owl_v3f32_unsafe_broadcast_comp(v, 0), A->column[0]);
    vr = owl_v3f32_add(
                        vr,
                        owl_v3f32_filter_apply(owl_v3f32_unsafe_broadcast_comp(v, 1), A->column[1])
                       );
    vr = owl_v3f32_add(
                        vr,
                        owl_v3f32_filter_apply(owl_v3f32_unsafe_broadcast_comp(v, 2), A->column[2])
                       );

    return vr;
}

//M = A * B
//
//
owl_mxf32_3x3* owl_mxf32_3x3_mul(owl_mxf32_3x3* M, owl_mxf32_3x3 const* A, owl_mxf32_3x3 const* B)
{
    owl_mxf32_3x3 S;

    for(int j = 0 ; j < 3 ; j++)
    {
        S.column[j] = owl_mxf32_3x3_transform_v3f32(A, B->column[j]);
    }

    owl_mxf32_3x3_copy(M, &S);

    return M;
}

//norm2(A)
//
//
float owl_mxf32_3x3_norm2(owl_mxf32_3x3 const* A)
{
    owl_v3f32 buffer = owl_v3f32_zero();
    for(int j = 0 ; j < 3 ; j++)
    {
        buffer = owl_v3f32_add(
                                buffer,
                                owl_v3f32_comp_mul(A->column[j], A->column[j])
                               );
    }

    return sqrtf(owl_v3f32_dot(buffer, owl_v3f32_broadcast(1.0)));
}

//norminf(A)
//
//
float owl_mxf32_3x3_norminf(owl_mxf32_3x3 const* A)
{
    float norminf = 0.0;

    for(int j = 0 ; j < 3 ; j++)
    {
        norminf = fmaxf(norminf, owl_v3f32_norminf(A->column[j]));
    }

    return norminf;
}

//M = A ^(-1)
//
//
owl_mxf32_3x3* owl_mxf32_3x3_Inv(owl_mxf32_3x3* M, owl_mxf32_3x3 const* A)
{
    float inv_det = 1.0 / owl_mxf32_3x3_det(A);
    owl_mxf32_3x3 inv_det_comA;

    for(int j = 0 ; j < 3 ; j++)
    {
        owl_v3f32 v = owl_v3f32_cross(A->column[(j + 1) % 3], A->column[(j + 2) % 3]);
        inv_det_comA.column[j] = owl_v3f32_scalar_mul(v, inv_det);
    }

    return owl_mxf32_3x3_transp(M, &inv_det_comA);
}

//Return dominant eigenvalue
//
//
float owl_mxf32_3x3_sym_dominant_eigenvalue(owl_v3f32* eigenvector_ptr, owl_mxf32_3x3 const* A)
{
    float const delta = 1.0 / ((float)(1<<24));
    float const eps = delta / (2.0 * sqrtf(3.0));
    float const sqrt_eps = sqrtf(eps);
    float const f_sqrt_eps = 2.0 * sqrt_eps / ((1.0 + sqrt_eps) * (1.0 + sqrt_eps));

    float eigenvalue = 0.0;
    owl_v3f32 eigenvector = owl_v3f32_zero();

    unsigned int k = 1;
    owl_mxf32_3x3 A1;
    owl_mxf32_3x3_mul(&A1, A, A);
    owl_mxf32_3x3 Ak;
    owl_mxf32_3x3_copy(&Ak, &A1);
    float Tr = owl_mxf32_3x3_trace(&Ak);

    if(Tr > 0.0)
    {
        do
        {
            owl_mxf32_3x3_scalar_mul(&Ak, &Ak, 1.0 / Tr);
            owl_mxf32_3x3_mul(&Ak, &Ak, &Ak);
            Tr = owl_mxf32_3x3_trace(&Ak);

            k++;

        } while(k <= 27 && Tr > 0.0 && 1.0 - Tr > f_sqrt_eps);
    }

    float square_norm = 0.0;
    owl_v3f32 k_eigenvector;
    for(unsigned int i = 0 ; i < 3 ; i++)
    {
        owl_v3f32 v = Ak.column[i];
        float tmp = owl_v3f32_dot(v, v);
        if(tmp > square_norm)
        {
            k_eigenvector = v;
            square_norm = tmp;
        }
    }


    if(square_norm > 0.0)
    {
        owl_v3f32 A1_k_eigenvector = owl_mxf32_3x3_transform_v3f32(&A1, k_eigenvector);
        float abs_eigenvalue = sqrtf(owl_v3f32_dot(k_eigenvector, A1_k_eigenvector) / square_norm);

        owl_v3f32 A_k_eigenvector = owl_mxf32_3x3_transform_v3f32(A, k_eigenvector);

        owl_v3f32 v1 = owl_v3f32_add_scalar_mul(A_k_eigenvector, k_eigenvector, abs_eigenvalue);
        float n1 = owl_v3f32_dot(v1, v1);
        owl_v3f32 v2 = owl_v3f32_add_scalar_mul(A_k_eigenvector, k_eigenvector, -abs_eigenvalue);
        float n2 = owl_v3f32_dot(v2, v2);

        eigenvector = owl_v3f32_normalize((n1 >= n2) ? v1 : v2);
        eigenvalue = owl_v3f32_dot(
                                    eigenvector,
                                    owl_mxf32_3x3_transform_v3f32(A, eigenvector)
                                   );
    }
    else
    {
        eigenvalue = 0.0;
        eigenvector = owl_v3f32_set(1.0, 0.0, 0.0);
    }

    if(eigenvector_ptr != NULL)
    {
        *eigenvector_ptr = eigenvector;
    }

    return eigenvalue;
}

//A = P * D * tP
//
//
float* owl_mxf32_3x3_sym_diagonalize(float* eigenvalue_list, owl_mxf32_3x3* P, owl_mxf32_3x3 const* A)
{
    //Computes an eigenvector associated with the dominant eigenvalue
    owl_v3f32 V0;
    eigenvalue_list[0] = owl_mxf32_3x3_sym_dominant_eigenvalue(&V0, A);

    owl_mxf32_3x3 B;
    owl_mxf32_3x3_diag(&B, 1.0);

    owl_v3f32 B_sev[2];
    {
        if( owl_v3f32_unsafe_get_component(V0, 0) < owl_v3f32_unsafe_get_component(V0, 1) )
        {
            B_sev[0] = owl_v3f32_cross(B.column[0], V0);
        }
        else
        {
            B_sev[0] = owl_v3f32_cross(B.column[1], V0);
        }

        B_sev[0] = owl_v3f32_normalize(B_sev[0]);
        B_sev[1] = owl_v3f32_cross(V0, B_sev[0]);
    }

    owl_mxf32_2x2 As;
    {
        owl_v3f32 Im_Bs[2] =
        {
            owl_mxf32_3x3_transform_v3f32(A, B_sev[0]),
            owl_mxf32_3x3_transform_v3f32(A, B_sev[1])
        };

        owl_mxf32_2x2_set(
                            &As,
                            owl_v3f32_dot(Im_Bs[0], B_sev[0]),
                            owl_v3f32_dot(Im_Bs[0], B_sev[1]),
                            owl_v3f32_dot(Im_Bs[0], B_sev[1]),
                            owl_v3f32_dot(Im_Bs[1], B_sev[1])
                          );
    }

    owl_mxf32_2x2 Ps;
    owl_mxf32_2x2_diagonalize_sym(eigenvalue_list + 1, &Ps, &As);

    if(P != NULL)
    {
        float flat_Ps[4] OWL_ALIGN16;
        owl_mxf32_2x2_store(flat_Ps, &Ps);

        P->column[0] = V0;
        P->column[1] = owl_v3f32_add_scalar_mul(
                                                    owl_v3f32_scalar_mul(B_sev[0], flat_Ps[0]),
                                                    B_sev[1], flat_Ps[1]
                                                );
        P->column[2] = owl_v3f32_add_scalar_mul(
                                                    owl_v3f32_scalar_mul(B_sev[0], flat_Ps[2]),
                                                    B_sev[1], flat_Ps[3]
                                                );
    }

    return eigenvalue_list;
}
