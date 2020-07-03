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
            if(n == 0 || M->column != NULL)
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

//M = A
//
//
owl_Mxf32_mxn* owl_Mxf32_mxn_Copy(owl_Mxf32_mxn* M, owl_Mxf32_mxn const* A, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(M->m == A->m && M->n == A->n)
        {
            for(unsigned int j = 0 ; j < M->n && error == OWL_SUCCESS ; j++)
            {
                owl_Vnf32_Copy(M->column[j], A->column[j], &error, &error);
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

//Tr(A)
//
//
float owl_Mxf32_mxn_Tr(owl_Mxf32_mxn const* A, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    float Tr = 0.0;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(A->m == A->n)
        {
            unsigned int n = A->n;
            for(unsigned int j = 0 ; j < n ; j++)
            {
                Tr += owl_Vnf32_GetComponent(A->column[j], j, &error, &error);
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

    return Tr;
}

//Return the dominant eigenvector of a symetric matrix
//
//
float owl_Mxf32_mxn_sym_DominantEigenvalue(owl_Vnf32* Eigenvector, owl_Mxf32_mxn const* A, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    float eigenvalue = 0.0;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(A->m == A->n && (Eigenvector == NULL || Eigenvector->n == A->n))
        {
            unsigned int n = A->n;

            if(n == 1)
            {
                eigenvalue = owl_Vnf32_GetComponent(A->column[0], 0, &error, &error);

                if(Eigenvector != NULL)
                {
                    owl_Vnf32_SetComponent(Eigenvector, 0, 1.0, &error, &error);
                }
            }
            else if(n > 1)
            {
                float const delta = 1.0 / ((float)(1<<24));
                float const eps = delta / (2.0 * sqrtf((float)n));
                float const sqrt_eps = sqrtf(eps);
                float const f_sqrt_eps = 2.0 * sqrt_eps / ((1.0 + sqrt_eps) * (1.0 + sqrt_eps));

                unsigned int k = 1;
                owl_Mxf32_mxn* A1 = owl_Mxf32_mxn_Create(n, n, &error, &error);
                owl_Mxf32_mxn_Mul(A1, A, A, &error, &error);
                owl_Mxf32_mxn* Ak = owl_Mxf32_mxn_Create(n, n, &error, &error);
                owl_Mxf32_mxn_Copy(Ak, A1, &error, &error);
                float Tr = owl_Mxf32_mxn_Tr(Ak, &error, &error);
                owl_Mxf32_mxn* Tk = owl_Mxf32_mxn_Create(n, n, &error, &error);

                owl_Mxf32_mxn* A1_k_MatEigenvector = owl_Mxf32_mxn_Create(n, 1, &error, &error);
                owl_Mxf32_mxn* k_MatEigenvector = owl_Mxf32_mxn_Create(n, 1, &error, &error);

                owl_Vnf32* V1 = owl_Vnf32_Create(n, &error, &error);
                owl_Vnf32* V2 = owl_Vnf32_Create(n, &error, &error);

                if(error == OWL_SUCCESS)
                {
                    if(Tr > 0.0)
                    {
                        do
                        {
                            owl_Mxf32_mxn_ScalarMul(Tk, Ak, 1.0 / Tr, &error, &error);
                            owl_Mxf32_mxn_Mul(Ak, Tk, Tk, &error, &error);
                            Tr = owl_Mxf32_mxn_Tr(Ak, &error, &error);

                            k++;

                        } while(k <= 27 && Tr > 0.0 && 1.0 - Tr > f_sqrt_eps);
                    }

                    unsigned int j0 = 0;
                    float square_norm = 0.0;
                    for(unsigned int j = 0 ; j < n ; j++)
                    {
                        float tmp_square_norm = owl_Vnf32_Dot(Ak->column[j], Ak->column[j], &error, &error);
                        if(tmp_square_norm > square_norm)
                        {
                            j0 = j;
                            square_norm = tmp_square_norm;
                        }
                    }

                    owl_Vnf32* k_Eigenvector = k_MatEigenvector->column[0];
                    owl_Vnf32_ScalarMul(k_Eigenvector, Ak->column[j0], 1.0 / sqrtf(square_norm), &error, &error);
                    owl_Vnf32* A1_k_Eigenvector = A1_k_MatEigenvector->column[0];

                    if(square_norm > 0.0)
                    {
                        owl_Mxf32_mxn_Mul(A1_k_MatEigenvector, A1, k_MatEigenvector, &error, &error);
                        float abs_eigenvalue = sqrtf(owl_Vnf32_Dot(k_Eigenvector, A1_k_Eigenvector, &error, &error));

                        owl_Vnf32_AddScalarMul(V1, A1_k_Eigenvector, k_Eigenvector, abs_eigenvalue, &error, &error);
                        float n1 = owl_Vnf32_Dot(V1, V1, &error, &error);
                        owl_Vnf32_AddScalarMul(V2, A1_k_Eigenvector, k_Eigenvector, -abs_eigenvalue, &error, &error);
                        float n2 = owl_Vnf32_Dot(V2, V2, &error, &error);

                        if(n1 >= n2)
                        {
                            owl_Vnf32_ScalarMul(k_Eigenvector, V1, 1.0 / sqrtf(n1), &error, &error);
                        }
                        else
                        {
                            owl_Vnf32_ScalarMul(k_Eigenvector, V2, 1.0 / sqrtf(n2), &error, &error);
                        }

                        owl_Mxf32_mxn_Mul(A1_k_MatEigenvector, A, k_MatEigenvector, &error, &error);
                        eigenvalue = owl_Vnf32_Dot(k_Eigenvector, A1_k_Eigenvector, &error, &error);

                        if(Eigenvector != NULL)
                        {
                            owl_Vnf32_Copy(Eigenvector, k_Eigenvector, &error, &error);
                        }
                    }
                    else
                    {
                        eigenvalue = 0.0;

                        if(Eigenvector != NULL)
                        {
                            owl_Vnf32_Zero(Eigenvector, &error, &error);
                            owl_Vnf32_SetComponent(Eigenvector, 0, 1.0, &error, &error);
                        }
                    }
                }

                owl_Mxf32_mxn_Destroy(A1);
                owl_Mxf32_mxn_Destroy(Ak);
                owl_Mxf32_mxn_Destroy(Tk);
                owl_Mxf32_mxn_Destroy(A1_k_MatEigenvector);
                owl_Mxf32_mxn_Destroy(k_MatEigenvector);

                owl_Vnf32_Destroy(V1);
                owl_Vnf32_Destroy(V2);
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

    return eigenvalue;
}

//A = P * D * tP with A symmetric and D = diag(eigenvalue_list)
//The length of eigenvalue_list must be at least n
//Parameter P is optional
//Return eigenvalue_list
float* owl_Mxf32_mxn_sym_Diagonalize(float* eigenvalue_list, owl_Mxf32_mxn* P, owl_Mxf32_mxn const* A, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(A->m == A->n)
        {
            unsigned int n = A->n;

            if(n > 0)
            {
                owl_Vnf32* DominantEigenvector = owl_Vnf32_Create(n, &error, &error);
                owl_Mxf32_mxn* Q = owl_Mxf32_mxn_Create(n, n, &error, &error);
                owl_Mxf32_mxn* Q_ = owl_Mxf32_mxn_Create(n, n - 1, &error, &error);
                owl_Mxf32_mxn* Asub = owl_Mxf32_mxn_Create(n - 1, n - 1, &error, &error);

                if(error == OWL_SUCCESS)
                {
                    eigenvalue_list[0] = owl_Mxf32_mxn_sym_DominantEigenvalue(DominantEigenvector, A, &error, &error);


                    //Fills Q=(Q1|...|Qn) with an orthonormal base and Q1=DominantEigenVector
                    unsigned int i0 = 0;
                    float max_abs_component = 0.0;
                    for(unsigned int i = 0 ; i < n ; i++)
                    {
                        float tmp_abs_component = fabsf(owl_Vnf32_GetComponent(DominantEigenvector, i, &error, &error));
                        if(tmp_abs_component > max_abs_component)
                        {
                            i0 = i;
                            max_abs_component = tmp_abs_component;
                        }
                    }
                    owl_Vnf32_Copy(Q->column[0], DominantEigenvector, &error, &error);
                    for(unsigned int j = 1 ; j < n ; j++)
                    {
                        owl_Vnf32_Zero(Q->column[j], &error, &error);
                        unsigned int i = (j > i0) ? j : j - 1;
                        owl_Vnf32_SetComponent(Q->column[j], i, 1.0, &error, &error);
                    }
                    owl_Vnf32_GramSchmidt(Q->column, (owl_Vnf32 const**)Q->column, n, &error, &error);
                    //Fills Q_=(Q2|...|Qn)
                    for(unsigned int j = 1 ; j < n ; j++)
                    {
                        owl_Vnf32_Copy(Q_->column[j - 1], Q->column[j], &error, &error);
                    }

                    //Fills Asub
                    owl_Mxf32_mxn_Mul(Q, A, Q, &error, &error);
                    for(unsigned int i = 0 ; i < n - 1 ; i++)
                    {
                        for(unsigned int j = 0 ; j < n - 1 ; j++)
                        {
                            float a = owl_Vnf32_Dot(Q_->column[i], Q->column[j + 1], &error, &error);
                            owl_Mxf32_mxn_SetElement(Asub, i, j, a, &error, &error);
                        }
                    }

                    if(P == NULL)
                    {
                        owl_Mxf32_mxn_sym_Diagonalize(eigenvalue_list + 1, NULL, Asub, &error, &error);
                    }
                    else
                    {
                        owl_Mxf32_mxn* Psub = owl_Mxf32_mxn_Create(n - 1, n - 1, &error, &error);

                        if(error == OWL_SUCCESS)
                        {
                            owl_Mxf32_mxn_sym_Diagonalize(eigenvalue_list + 1, Psub, Asub, &error, &error);

                            owl_Mxf32_mxn_Mul(Q_, Q_, Psub, &error, &error);

                            owl_Vnf32_Copy(P->column[0], DominantEigenvector, &error, &error);
                            for(unsigned int j = 1 ; j < n ; j++)
                            {
                                owl_Vnf32_Copy(P->column[j], Q_->column[j - 1], &error, &error);
                            }
                        }

                        owl_Mxf32_mxn_Destroy(Psub);
                    }

                }

                owl_Vnf32_Destroy(DominantEigenvector);
                owl_Mxf32_mxn_Destroy(Q);
                owl_Mxf32_mxn_Destroy(Q_);
                owl_Mxf32_mxn_Destroy(Asub);
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

    return eigenvalue_list;
}
