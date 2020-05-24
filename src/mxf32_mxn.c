#include <OWL/mxf32_mxn.h>

owl_Mxf32_mxn* owl_Mxf32_mxn_Create(unsigned int m, unsigned int n, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;
    owl_Mxf32_mxn* M = NULL;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        M = malloc(sizeof(*M));
        if(M != NULL)
        {
            M->m = m;
            M->n = n;
            M->column = malloc(n * sizeof(*M->column));
            if(M->column != NULL)
            {
                for(int j = 0 ; j < n ; j++)
                {
                    M->column[j] = NULL;
                }

                for(int j = 0 ; j < n && error == OWL_SUCCESS ; j++)
                {
                    M->column[j] = owl_Vnf32_Create(m, &error, &error);
                }
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
        owl_Mxf32_mxn_Destroy(M);
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
void owl_Mxf32_mxn_Destroy(owl_Mxf32_mxn* M)
{
    if(M != NULL)
    {
        if(M->column != NULL)
        {
            for(unsigned int j = 0 ; j < M->n ; j++)
            {
                owl_Vnf32_Destroy(M->column[j]);
            }

            free(M->column);
        }

        free(M);
    }
}

//M = 0
//
//
owl_Mxf32_mxn* owl_Mxf32_mxn_Zero(owl_Mxf32_mxn* M, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        for(unsigned int j = 0 ; j < M->n && error == OWL_SUCCESS ; j++)
        {
            owl_Vnf32_Zero(M->column[j], &error, &error);
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
owl_Mxf32_mxn* owl_Mxf32_mxn_Diag(owl_Mxf32_mxn* M, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        owl_Mxf32_mxn_Zero(M, &error, &error);

        if(error == OWL_SUCCESS)
        {
            unsigned int r = M->m;
            if(r > M->n)
            {
                r = M->n;
            }

            for(unsigned int j = 0 ; j < r ; j++)
            {
                owl_Mxf32_mxn_SetElement(M, j, j, a, &error, &error);
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
float owl_Mxf32_mxn_GetElement(owl_Mxf32_mxn const* A, unsigned int i, unsigned int j, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    float a = 0.0;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(i < A->m && j < A->n)
        {
            a = owl_Vnf32_GetComponent(A->column[j], i, &error, &error);
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
owl_Mxf32_mxn* owl_Mxf32_mxn_SetElement(owl_Mxf32_mxn* M, unsigned int i, unsigned int j, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(i < M->m && j < M->n)
        {
            owl_Vnf32_SetComponent(M->column[j], i, a, &error, &error);
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

//M = A + B
//
//
owl_Mxf32_mxn* owl_Mxf32_mxn_Add(owl_Mxf32_mxn* M, owl_Mxf32_mxn const* A, owl_Mxf32_mxn const* B, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(M->m == A->m && M->m == B->m && M->n == A->n && M->n == B->n)
        {
            for(unsigned int j = 0 ; j < M->n && error == OWL_SUCCESS ; j++)
            {
                owl_Vnf32_Add(M->column[j], A->column[j], B->column[j], &error, &error);
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
owl_Mxf32_mxn* owl_Mxf32_mxn_Sub(owl_Mxf32_mxn* M, owl_Mxf32_mxn const* A, owl_Mxf32_mxn const* B, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(M->m == A->m && M->m == B->m && M->n == A->n && M->n == B->n)
        {
            for(unsigned int j = 0 ; j < M->n && error == OWL_SUCCESS ; j++)
            {
                owl_Vnf32_Sub(M->column[j], A->column[j], B->column[j], &error, &error);
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
owl_Mxf32_mxn* owl_Mxf32_mxn_ScalarMul(owl_Mxf32_mxn* M, owl_Mxf32_mxn const* A, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(M->m == A->m && M->n == A->n)
        {
            for(unsigned int j = 0 ; j < M->n && error == OWL_SUCCESS ; j++)
            {
                owl_Vnf32_ScalarMul(M->column[j], A->column[j], a, &error, &error);
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
owl_Mxf32_mxn* owl_Mxf32_mxn_AddScalarMul(owl_Mxf32_mxn* M, owl_Mxf32_mxn const* A, owl_Mxf32_mxn const* B, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(M->m == A->m && M->m == B->m && M->n == A->n && M->n == B->n)
        {
            for(unsigned int j = 0 ; j < M->n && error == OWL_SUCCESS ; j++)
            {
                owl_Vnf32_AddScalarMul(M->column[j], A->column[j], B->column[j], a, &error, &error);
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
owl_Mxf32_mxn* owl_Mxf32_mxn_Mul(owl_Mxf32_mxn* M, owl_Mxf32_mxn const* A, owl_Mxf32_mxn const* B, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(A->n == B->m && M->m == A->m && M->n == B->n)
        {
            owl_Vnf32* Uacc = owl_Vnf32_Create(M->m, &error, &error);
            owl_Mxf32_mxn* S = NULL;
            owl_Mxf32_mxn* F;

            if(M != A)
            {
                F = M;
            }
            else
            {
                S = owl_Mxf32_mxn_Create(M->m, M->n, &error, &error);
                F = S;
            }

            if(error == OWL_SUCCESS)
            {
                for(unsigned int j = 0 ; j < M->n && error == OWL_SUCCESS ; j++)
                {
                    owl_Vnf32_Zero(Uacc, &error, &error);

                    for(unsigned int i = 0 ; i < A->n && error == OWL_SUCCESS ; i++)
                    {
                        float bij = owl_Vnf32_GetComponent(B->column[j], i, &error, &error);
                        owl_Vnf32_AddScalarMul(Uacc, Uacc, A->column[i], bij, &error, &error);
                    }

                    owl_Vnf32_ScalarMul(F->column[j], Uacc, 1.0, &error, &error);
                }

                if(error == OWL_SUCCESS)
                {
                    if(F != M)
                    {
                        owl_Mxf32_mxn_ScalarMul(M, F, 1.0, &error, &error);
                    }
                }
            }

            owl_Mxf32_mxn_Destroy(S);
            owl_Vnf32_Destroy(Uacc);
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
