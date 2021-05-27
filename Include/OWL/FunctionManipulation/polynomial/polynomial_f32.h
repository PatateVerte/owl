#ifndef OWL_POLYNOMIAL_F32_H_INCLUDED
#define OWL_POLYNOMIAL_F32_H_INCLUDED

#include <OWL/owl.h>

//
typedef struct
{
    unsigned int nb_coeff;
    float* coeff;

} owl_Polynomial_f32;

//Create a polynomial
OWL_DLL_EXPORT owl_Polynomial_f32* owl_Polynomial_f32_Create(unsigned int nb_coeff, owl_error* ret_error, owl_error const* pass_through_error);

//Destroy a polynomial
OWL_DLL_EXPORT void owl_Polynomial_f32_Destroy(owl_Polynomial_f32* P);

//Copy
OWL_DLL_EXPORT owl_Polynomial_f32* owl_Polynomial_f32_Copy(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, owl_error* ret_error, owl_error const* pass_through_error);

//P = 0
OWL_DLL_EXPORT owl_Polynomial_f32* owl_Polynomial_f32_ZeroCoeff(owl_Polynomial_f32* P, owl_error* ret_error, owl_error const* pass_through_error);

//Adjust number of coefficients
OWL_DLL_EXPORT owl_Polynomial_f32* owl_Polynomial_f32_AdjustNbCoeff(owl_Polynomial_f32* P, unsigned int new_nb_coeff, owl_error* ret_error, owl_error const* pass_through_error);

//Cuts null high degree coefficients
OWL_DLL_EXPORT owl_Polynomial_f32* owl_Polynomial_f32_Cut(owl_Polynomial_f32* P, owl_error* ret_error, owl_error const* pass_through_error);

//deg(P1)
//deg(0) = -1
OWL_DLL_EXPORT int owl_Polynomial_f32_Deg(owl_Polynomial_f32 const* P1, owl_error* ret_error, owl_error const* pass_through_error);

//P(x)
OWL_DLL_EXPORT float owl_Polynomial_f32_Evaluate(owl_Polynomial_f32 const* P1, float x, owl_error* ret_error, owl_error const* pass_through_error);

//P = a * P1
OWL_DLL_EXPORT owl_Polynomial_f32* owl_Polynomial_f32_ScalarMul(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, float a, owl_error* ret_error, owl_error const* pass_through_error);

//P = P1 + a * P2
OWL_DLL_EXPORT owl_Polynomial_f32* owl_Polynomial_f32_AddScalarMul(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, owl_Polynomial_f32 const* P2, float a, owl_error* ret_error, owl_error const* pass_through_error);

//P = P1 * P2
OWL_DLL_EXPORT owl_Polynomial_f32* owl_Polynomial_f32_Mul(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, owl_Polynomial_f32 const* P2, owl_error* ret_error, owl_error const* pass_through_error);

#endif // OWL_POLYNOMIAL_F32_H_INCLUDED
