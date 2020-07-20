#ifndef OWL_TAYLOR_SERIES_F32_H_INCLUDED
#define OWL_TAYLOR_SERIES_F32_H_INCLUDED

#include <OWL/owl.h>

typedef struct
{
    float x0;
    unsigned int order;

    unsigned int max_order;
    float* terms;

} owl_TaylorDevf32;

//
owl_TaylorDevf32* owl_TaylorDevf32_Create(unsigned int max_order, owl_error* ret_error, owl_error const* pass_through_error);

//
void owl_TaylorDevf32_Destroy(owl_TaylorDevf32* D);

//D(x)
float owl_TaylorDevf32_Evaluate(owl_TaylorDevf32 const* D, float x, owl_error* ret_error, owl_error const* pass_through_error);

//Df = 0
owl_TaylorDevf32* owl_TaylorDevf32_Zero(owl_TaylorDevf32* Df, float x0, unsigned int order, owl_error* ret_error, owl_error const* pass_through_error);

//Df = a * D1
owl_TaylorDevf32* owl_TaylorDevf32_ScalarMul(owl_TaylorDevf32* Df, owl_TaylorDevf32 const* D1, float a, owl_error* ret_error, owl_error const* pass_through_error);

//Df = D1 + a * D2
owl_TaylorDevf32* owl_TaylorDevf32_AddScalarMul(owl_TaylorDevf32* Df, owl_TaylorDevf32 const* D1, owl_TaylorDevf32 const* D2, float a, owl_error* ret_error, owl_error const* pass_through_error);

//Df = D1 * D2
owl_TaylorDevf32* owl_TaylorDevf32_Mul(owl_TaylorDevf32* Df, owl_TaylorDevf32 const* D1, owl_TaylorDevf32 const* D2, owl_error* ret_error, owl_error const* pass_through_error);

//Df = D1 o D2
owl_TaylorDevf32* owl_TaylorDevf32_Composition(owl_TaylorDevf32* Df, owl_TaylorDevf32 const* D1, owl_TaylorDevf32 const* D2, owl_error* ret_error, owl_error const* pass_through_error);


#endif // OWL_TAYLOR_SERIES_F32_H_INCLUDED
