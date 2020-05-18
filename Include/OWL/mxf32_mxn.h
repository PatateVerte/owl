#ifndef OWL_MXF32_MXN_H_INCLUDED
#define OWL_MXF32_MXN_H_INCLUDED

#include <OWL/owl.h>

#include <OWL/vnf32.h>

typedef struct
{
    unsigned int m;
    unsigned int n;

    owl_Vnf32** column;

} owl_Mxf32_mxn;

//Create owl_Mxf32_mxn
owl_Mxf32_mxn* owl_Mxf32_mxn_Create(unsigned int m, unsigned int n, owl_error* ret_error, owl_error const* pass_through_error);

//Destroy owl_Mxf32_mxn
void owl_Mxf32_mxn_Destroy(owl_Mxf32_mxn* M);

//A[i, j]
float owl_Mxf32_mxn_GetCoefficient(owl_Mxf32_mxn const* A, unsigned int i, unsigned int j, owl_error* ret_error, owl_error const* pass_through_error);

//M[i, j] = a
owl_Mxf32_mxn* owl_Mxf32_mxn_SetCoefficient(owl_Mxf32_mxn* M, unsigned int i, unsigned int j, float a, owl_error* ret_error, owl_error const* pass_through_error);

//M = A + B
owl_Mxf32_mxn* owl_Mxf32_mxn_Add(owl_Mxf32_mxn* M, owl_Mxf32_mxn const* A, owl_Mxf32_mxn const* B, owl_error* ret_error, owl_error const* pass_through_error);

//M = A - B
owl_Mxf32_mxn* owl_Mxf32_mxn_Sub(owl_Mxf32_mxn* M, owl_Mxf32_mxn const* A, owl_Mxf32_mxn const* B, owl_error* ret_error, owl_error const* pass_through_error);

//M = a * A
owl_Mxf32_mxn* owl_Mxf32_mxn_ScalarMul(owl_Mxf32_mxn* M, owl_Mxf32_mxn const* A, float a, owl_error* ret_error, owl_error const* pass_through_error);

//M = A + a * B
owl_Mxf32_mxn* owl_Mxf32_mxn_AddScalarMul(owl_Mxf32_mxn* M, owl_Mxf32_mxn const* A, owl_Mxf32_mxn const* B, float a, owl_error* ret_error, owl_error const* pass_through_error);

//M = A * B
owl_Mxf32_mxn* owl_Mxf32_mxn_Mul(owl_Mxf32_mxn* M, owl_Mxf32_mxn const* A, owl_Mxf32_mxn const* B, owl_error* ret_error, owl_error const* pass_through_error);


#endif // OWL_MXF32_MXN_H_INCLUDED
