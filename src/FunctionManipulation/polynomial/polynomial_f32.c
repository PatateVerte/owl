#include <OWL/FunctionManipulation/polynomial/polynomial_f32.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

//Create a polynomial
//
//
owl_Polynomial_f32* owl_Polynomial_f32_Create(unsigned int nb_coeff, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;
    owl_Polynomial_f32* P = NULL;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        P = malloc(sizeof(*P));
        if(P != NULL)
        {
            P->nb_coeff = 0;
            P->coeff = NULL;

            owl_Polynomial_f32_AdjustNbCoeff(P, nb_coeff, &error, &error);
        }
        else
        {
            error = OWL_MEMORY_ERROR;
        }
    }

    if(error != OWL_SUCCESS)
    {
        owl_Polynomial_f32_Destroy(P);
        P = NULL;
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return P;
}

//Destroy a polynomial
//
//
void owl_Polynomial_f32_Destroy(owl_Polynomial_f32* P)
{
    if(P != NULL)
    {
        free(P->coeff);
        free(P);
    }
}

//Copy
//
//
owl_Polynomial_f32* owl_Polynomial_f32_Copy(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        owl_Polynomial_f32_AdjustNbCoeff(P, P1->nb_coeff, &error, &error);
        if(error == OWL_SUCCESS)
        {
            memcpy(P->coeff, P1->coeff, (size_t)P->nb_coeff * sizeof(*P->coeff));
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return P;
}

//P = 0
//
//
owl_Polynomial_f32* owl_Polynomial_f32_ZeroCoeff(owl_Polynomial_f32* P, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        for(unsigned int i = 0 ; i < P->nb_coeff ; i++)
        {
            P->coeff[i] = 0.0;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return P;
}

//Adjust number of coefficients
//
//
owl_Polynomial_f32* owl_Polynomial_f32_AdjustNbCoeff(owl_Polynomial_f32* P, unsigned int new_nb_coeff, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        void* tmp_ptr = realloc(P->coeff, (size_t)new_nb_coeff * sizeof(*P->coeff));
        if(new_nb_coeff == 0 || tmp_ptr != NULL)
        {
            P->coeff = tmp_ptr;

            unsigned int old_nb_coeff = P->nb_coeff;

            if(new_nb_coeff > old_nb_coeff)
            {
                for(unsigned int i = old_nb_coeff ; i < new_nb_coeff ; i++)
                {
                    P->coeff[i] = 0.0;
                }
            }

            P->nb_coeff = new_nb_coeff;
        }
        else
        {
            error = OWL_MEMORY_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return P;
}

//Cuts null high degree coefficients
//
//
owl_Polynomial_f32* owl_Polynomial_f32_Cut(owl_Polynomial_f32* P, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        int deg = owl_Polynomial_f32_Deg(P, &error, &error);
        owl_Polynomial_f32_AdjustNbCoeff(P, (unsigned int)(deg + 1), &error, &error);
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return P;
}

//deg(P1)
//deg(0) = -1
//
int owl_Polynomial_f32_Deg(owl_Polynomial_f32 const* P1, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;
    int deg = (int)P1->nb_coeff - 1;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        while(deg >= 0 && P1->coeff[deg] == 0.0)
        {
            deg--;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return deg;
}

//P(x)
float owl_Polynomial_f32_Evaluate(owl_Polynomial_f32 const* P1, float x, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;
    float P_x = 0.0;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        for(unsigned int i = 0 ; i < P1->nb_coeff ; i++)
        {
            P_x = x * P_x + P1->coeff[P1->nb_coeff - 1 - i];
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return P_x;
}

//P = a * P1
owl_Polynomial_f32* owl_Polynomial_f32_ScalarMul(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(a == 0.0)
        {
            owl_Polynomial_f32_AdjustNbCoeff(P, 0, &error, &error);
        }
        else
        {
            unsigned int r = P1->nb_coeff;
            owl_Polynomial_f32_AdjustNbCoeff(P, r, &error, &error);

            if(error == OWL_SUCCESS)
            {
                for(unsigned int i = 0 ; i < r ; i++)
                {
                    P->coeff[i] = a * P1->coeff[i];
                }
            }
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return P;
}

//P = P1 + a * P2
owl_Polynomial_f32* owl_Polynomial_f32_AddScalarMul(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, owl_Polynomial_f32 const* P2, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        unsigned int r = (P1->nb_coeff >= P2->nb_coeff) ? P1->nb_coeff : P2->nb_coeff;
        owl_Polynomial_f32_AdjustNbCoeff(P, r, &error, &error);

        if(error == OWL_SUCCESS)
        {
            for(unsigned int i = 0 ; i < r ; i++)
            {
                P->coeff[i] = P1->coeff[i] + a * P2->coeff[i];
            }

            owl_Polynomial_f32_Cut(P, &error, &error);
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return P;
}

//P = P1 * P2
//
//
owl_Polynomial_f32* owl_Polynomial_f32_Mul(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, owl_Polynomial_f32 const* P2, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        unsigned int r = (P1->nb_coeff + P2->nb_coeff > 0) ? P1->nb_coeff + P2->nb_coeff - 1 : 0;

        owl_Polynomial_f32* Q = owl_Polynomial_f32_Create(r, &error, &error);

        if(error == OWL_SUCCESS)
        {
            for(unsigned int i = 0 ; i < P1->nb_coeff ; i++)
            {
                for(unsigned int j = 0 ; j < P2->nb_coeff ; j++)
                {
                    Q->coeff[i + j] += P1->coeff[i] * P2->coeff[j];
                }
            }

            owl_Polynomial_f32_Cut(Q, &error, &error);
            owl_Polynomial_f32_Copy(P, Q, &error, &error);
        }

        owl_Polynomial_f32_Destroy(Q);
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return P;
}
