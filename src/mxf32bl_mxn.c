#include <OWL/mxf32bl_mxn.h>

#include <malloc.h>

owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Create(unsigned int m, unsigned int n, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;
    owl_Mxf32bl_mxn* M = NULL;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        M = malloc(sizeof(*M));
        if(M != NULL)
        {
            M->m = m;
            M->n = n;

            M->nb_m_blocks = (m % 4 == 0) ? m / 4 : (m / 4) + 1;
            M->nb_n_blocks = (n % 4 == 0) ? n / 4 : (n / 4) + 1;
            M->nb_blocks = (size_t)M->nb_m_blocks * (size_t)M->nb_n_blocks;

            M->blocks = _aligned_malloc(M->nb_blocks * sizeof(*M->blocks), 16);
            if(M->blocks != NULL)
            {
                owl_Mxf32bl_mxn_Zero(M, &error, &error);
            }
            else
            {
                error = OWL_MEMORY_ERROR;
            }
        }
        else
        {
            error = OWL_MEMORY_ERROR;
        }

    }

    if(error != OWL_SUCCESS)
    {
        owl_Mxf32bl_mxn_Destroy(M);
        M = NULL;
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return M;
}

//
//
//
void owl_Mxf32bl_mxn_Destroy(owl_Mxf32bl_mxn* M)
{
    if(M != NULL)
    {
        _aligned_free(M->blocks);
        free(M);
    }
}

//M = 0
//
//
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Zero(owl_Mxf32bl_mxn* M, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        for(size_t k = 0 ; k < M->nb_blocks ; k++)
        {
            owl_mxf32_4x4_zero(M->blocks + k);
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return M;
}

//M = "Diag_mn(a)"
//
//
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Diag(owl_Mxf32bl_mxn* M, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        owl_Mxf32bl_mxn_Zero(M, &error, &error);

        if(error == OWL_SUCCESS)
        {
            unsigned int r = M->m;
            if(r > M->n)
            {
                r = M->n;
            }

            for(unsigned int jb = 0 ; jb < r / 4 ; jb++)
            {
                owl_mxf32_4x4_diag(M->blocks + jb * M->m + jb, a);
            }

            for(unsigned int j = 0 ; j < (r % 4) ; j++)
            {
                owl_mxf32_4x4_set_element(M->blocks + (r / 4) * M->m + (r / 4), a, j, j);
            }
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return M;
}

//A[i, j]
//
//
float owl_Mxf32bl_mxn_GetElement(owl_Mxf32bl_mxn const* A, unsigned int i, unsigned int j, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    float a = 0.0;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(i < A->m && j < A->n)
        {
            a = owl_mxf32_4x4_get_element(A->blocks + (j / 4) * A->m + (i / 4), i % 4, j % 4);
        }
        else
        {
            error = OWL_OUT_OF_BOUND_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return a;
}

//M[i, j] = a
//
//
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_SetElement(owl_Mxf32bl_mxn* M, unsigned int i, unsigned int j, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(i < M->m && j < M->n)
        {
            owl_mxf32_4x4_set_element(M->blocks + (j / 4) * M->m + (i / 4), a, i % 4, j % 4);
        }
        else
        {
            error = OWL_OUT_OF_BOUND_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return M;
}

//M = A
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Copy(owl_Mxf32bl_mxn* M, owl_Mxf32bl_mxn const* A, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(M->m == A->m && M->n == A->n)
        {
            for(unsigned int k = 0 ; k < M->nb_blocks ; k++)
            {
                owl_mxf32_4x4_copy(M->blocks + k, A->blocks + k);
            }
        }
        else
        {
            error = OWL_DIMENSION_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return M;
}

//M = A + B
//
//
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Add(owl_Mxf32bl_mxn* M, owl_Mxf32bl_mxn const* A, owl_Mxf32bl_mxn const* B, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(M->m == A->m && M->m == B->m && M->n == A->n && M->n == B->n)
        {
            for(unsigned int k = 0 ; k < M->nb_blocks ; k++)
            {
                owl_mxf32_4x4_add(M->blocks + k, A->blocks + k, B->blocks + k);
            }
        }
        else
        {
            error = OWL_DIMENSION_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return M;
}

//M = A - B
//
//
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Sub(owl_Mxf32bl_mxn* M, owl_Mxf32bl_mxn const* A, owl_Mxf32bl_mxn const* B, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(M->m == A->m && M->m == B->m && M->n == A->n && M->n == B->n)
        {
            for(unsigned int k = 0 ; k < M->nb_blocks ; k++)
            {
                owl_mxf32_4x4_sub(M->blocks + k, A->blocks + k, B->blocks + k);
            }
        }
        else
        {
            error = OWL_DIMENSION_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return M;
}

//M = a * A
//
//
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_ScalarMul(owl_Mxf32bl_mxn* M, owl_Mxf32bl_mxn const* A, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(M->m == A->m && M->n == A->n)
        {
            for(unsigned int k = 0 ; k < M->nb_blocks ; k++)
            {
                owl_mxf32_4x4_scalar_mul(M->blocks + k, A->blocks + k, a);
            }
        }
        else
        {
            error = OWL_DIMENSION_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return M;
}

//M = A + a * B
//
//
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_AddScalarMul(owl_Mxf32bl_mxn* M, owl_Mxf32bl_mxn const* A, owl_Mxf32bl_mxn const* B, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(M->m == A->m && M->m == B->m && M->n == A->n && M->n == B->n)
        {
            for(unsigned int k = 0 ; k < M->nb_blocks ; k++)
            {
                owl_mxf32_4x4_add_scalar_mul(M->blocks + k, A->blocks + k, B->blocks + k, a);
            }
        }
        else
        {
            error = OWL_DIMENSION_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return M;
}

//M = A * B
//
//
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Mul(owl_Mxf32bl_mxn* M, owl_Mxf32bl_mxn const* A, owl_Mxf32bl_mxn const* B, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(A->n == B->m && M->m == A->m && M->n == B->n)
        {
            owl_Mxf32bl_mxn* S = NULL;
            owl_mxf32_4x4 Macc;

            owl_Mxf32bl_mxn* F;

            if(M != A)
            {
                F = M;
            }
            else
            {
                S = owl_Mxf32bl_mxn_Create(M->m, M->n, &error, &error);
                F = S;
            }

            if(error == OWL_SUCCESS)
            {
                for(unsigned int jbl = 0 ; jbl < M->n ; jbl++)
                {
                    for(unsigned int ibl = 0 ; ibl < M->m ; ibl++)
                    {
                        owl_mxf32_4x4_zero(&Macc);

                        for(unsigned int k = 0 ; k < A->nb_n_blocks ; k++)
                        {
                            owl_mxf32_4x4 Mprod;
                            owl_mxf32_4x4_mul(&Mprod, A->blocks + A->m * k + ibl, B->blocks + B->m * jbl + k);
                            owl_mxf32_4x4_add(&Macc, &Macc, &Mprod);
                        }

                        owl_mxf32_4x4_copy(F->blocks + F->m * jbl + ibl, &Macc);
                    }
                }

                if(F != M)
                {
                    owl_Mxf32bl_mxn_Copy(M, A, &error, &error);
                }
            }

            owl_Mxf32bl_mxn_Destroy(S);
        }
        else
        {
            error = OWL_DIMENSION_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return M;
}
