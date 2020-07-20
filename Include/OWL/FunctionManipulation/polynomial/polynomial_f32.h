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
owl_Polynomial_f32* owl_Polynomial_f32_Create(unsigned int nb_coeff, owl_error* ret_error, owl_error const* pass_through_error);

//Destroy a polynomial
void owl_Polynomial_f32_Destroy(owl_Polynomial_f32* P);

//Copy
owl_Polynomial_f32* owl_Polynomial_f32_Copy(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, owl_error* ret_error, owl_error const* pass_through_error);

//P = 0
owl_Polynomial_f32* owl_Polynomial_f32_ZeroCoeff(owl_Polynomial_f32* P, owl_error* ret_error, owl_error const* pass_through_error);

//Adjust number of coefficients
owl_Polynomial_f32* owl_Polynomial_f32_AdjustNbCoeff(owl_Polynomial_f32* P, unsigned int new_nb_coeff, owl_error* ret_error, owl_error const* pass_through_error);

//Cuts null high degree coefficients
owl_Polynomial_f32* owl_Polynomial_f32_Cut(owl_Polynomial_f32* P, owl_error* ret_error, owl_error const* pass_through_error);

//deg(P1)
//deg(0) = -1
int owl_Polynomial_f32_Deg(owl_Polynomial_f32 const* P1, owl_error* ret_error, owl_error const* pass_through_error);

//P(x)
float owl_Polynomial_f32_Evaluate(owl_Polynomial_f32 const* P1, float x, owl_error* ret_error, owl_error const* pass_through_error);

//P = a * P1
owl_Polynomial_f32* owl_Polynomial_f32_ScalarMul(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, float a, owl_error* ret_error, owl_error const* pass_through_error);

//P = P1 + a * P2
owl_Polynomial_f32* owl_Polynomial_f32_AddScalarMul(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, owl_Polynomial_f32 const* P2, float a, owl_error* ret_error, owl_error const* pass_through_error);

//P = P1 * P2
owl_Polynomial_f32* owl_Polynomial_f32_Mul(owl_Polynomial_f32* P, owl_Polynomial_f32 const* P1, owl_Polynomial_f32 const* P2, owl_error* ret_error, owl_error const* pass_through_error);

#endif // OWL_POLYNOMIAL_F32_H_INCLUDED
