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
