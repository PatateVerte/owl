#include <OWL/mxf32_3x3.h>

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
float owl_mxf32_3x3_get_element(owl_mxf32_3x3* M, int i, int j)
{
    int i_ = i % 3;
    int j_ = j % 3;

    float flat_column[4] OWL_ALIGN16;
    owl_v3f32_store4(flat_column, M->column[j_]);
    return flat_column[i_];
}

//
//
//
owl_mxf32_3x3* owl_mxf32_3x3_set_element(owl_mxf32_3x3* M, float value, int i, int j)
{
    int i_ = i % 3;
    int j_ = j % 3;

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
owl_v3f32 owl_mxf32_3x3_transform(owl_mxf32_3x3 const* A, owl_v3f32 v)
{
    owl_v3f32 v_ = v;
    owl_v3f32 vr = owl_v3f32_zero();

    for(int j = 2 ; j >= 0 ; j--)
    {
        v_ = owl_v3f32_rotate_comp(v_);
        vr = owl_v3f32_add(
                                vr,
                                owl_v3f32_scalar_mul(A->column[j], owl_v3f32_get_component(v_, 0))
                            );

    }

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
        S.column[j] = owl_mxf32_3x3_transform(A, B->column[j]);
    }

    owl_mxf32_3x3_copy(M, &S);

    return M;
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

//A = P * D * tP
//
//
owl_mxf32_3x3* owl_mxf32_3x3_diagonalize_sym(owl_mxf32_3x3* D, owl_mxf32_3x3* P, owl_mxf32_3x3 const* A)
{
    #define diag_nb_iter 25

    float vp_list[3] = {0.0, 0.0, 0.0};

    owl_mxf32_3x3 M;
    owl_mxf32_3x3_copy(&M, A);

    owl_mxf32_3x3 B;
    owl_mxf32_3x3_diag(&B, 1.0);

    owl_mxf32_3x3 VP_Base;
    owl_mxf32_3x3_zero(&VP_Base);

    for(int base_j = 0 ; base_j < 2 ; base_j++)
    {
        for(int k = 0 ; k < diag_nb_iter ; k++)
        {
            owl_mxf32_3x3_mul(&M, &M, &M);
            float norm2_M =  owl_mxf32_3x3_norm2(&M);
            owl_mxf32_3x3_scalar_mul(&M, &M, 1.0 / norm2_M);
        }

        owl_mxf32_3x3 T;
        owl_mxf32_3x3_mul(&T, &M, &B);

        int r = 0;
        float vp_abs = 0.0;
        for(int j = 0 ; j < 3 ; j++)
        {
            float tmp = owl_v3f32_norm( T.column[j] );
            if(tmp >= vp_abs)
            {
                r = j;
                vp_abs = tmp;
            }
        }

        VP_Base.column[base_j] = owl_v3f32_normalize( T.column[r] );
        vp_list[base_j] = owl_v3f32_dot(
                                            VP_Base.column[base_j],
                                            owl_mxf32_3x3_transform(A, VP_Base.column[base_j])
                                          );

        if(base_j == 0)
        {
            owl_mxf32_3x3_copy(&M, A);
            for(int j = 0 ; j < 3 ; j++)
            {
                M.column[j] = owl_v3f32_sub(
                                                M.column[j],
                                                owl_v3f32_scalar_mul(
                                                                        VP_Base.column[0],
                                                                        vp_list[0] * owl_v3f32_dot(VP_Base.column[0], B.column[j])
                                                                       )
                                               );
            }
        }
    }

    vp_list[2] = owl_mxf32_3x3_det(A) / (vp_list[0] * vp_list[1]);
    VP_Base.column[2] = owl_v3f32_cross(VP_Base.column[0], VP_Base.column[1]);

    for(int j = 0 ; j < 3 ; j++)
    {
        D->column[j] = owl_v3f32_scalar_mul(B.column[j], vp_list[j]);
    }

    if(P != NULL)
    {
        owl_mxf32_3x3_copy(P, &VP_Base);
    }

    return D;

    #undef diag_nb_iter
}
