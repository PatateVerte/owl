#include <OWL/polynomial_f32.h>

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
void owl_Polynomial_f32_Destroy(owl_Polynomial_f32* P)
{
    if(P != NULL)
    {
        free(P->coeff);
        free(P);
    }
}

//Copy
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

//Adjust number of coefficients
owl_Polynomial_f32* owl_Polynomial_f32_AdjustNbCoeff(owl_Polynomial_f32* P, unsigned int new_nb_coeff, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        void* tmp_ptr = realloc(P->coeff, (size_t)new_nb_coeff * sizeof(*P->coeff));
        if(new_nb_coeff == 0 || tmp_ptr != NULL)
        {
            P->coeff = tmp_ptr;
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
owl_Polynomial_f32* owl_Polynomial_f32_Cut(owl_Polynomial_f32* P, owl_error* ret_error, owl_error const* pass_through_error);

//deg(P1)
//deg(0) = 0
unsigned int owl_Polynomial_f32_Deg(owl_Polynomial_f32 const* P1, owl_error* ret_error, owl_error const* pass_through_error);

//P(x)
float owl_Polynomial_f32_Evaluate(owl_Polynomial_f32 const* P1, float x, owl_error* ret_error, owl_error const* pass_through_error);

//P = P1 + a * P2
owl_Polynomial_f32* owl_Polynomial_f32_ScalarMulAdd(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, owl_Polynomial_f32 const* P2, owl_error* ret_error, owl_error const* pass_through_error);

//P = P1 * P2
owl_Polynomial_f32* owl_Polynomial_f32_Mul(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, owl_Polynomial_f32 const* P2, owl_error* ret_error, owl_error const* pass_through_error);
