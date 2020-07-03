#ifndef OWL_MXF32BL_MXN_H_INCLUDED
#define OWL_MXF32BL_MXN_H_INCLUDED

#include <OWL/owl.h>

#include <OWL/mxf32_4x4.h>

typedef struct
{
    unsigned int m;
    unsigned int n;

	unsigned int nb_m_blocks;
	unsigned int nb_n_blocks;
	size_t nb_blocks;

    owl_mxf32_4x4* blocks;

} owl_Mxf32bl_mxn;

//Create owl_Mxf32bl_mxn
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Create(unsigned int m, unsigned int n, owl_error* ret_error, owl_error const* pass_through_error);

//Destroy owl_Mxf32bl_mxn
void owl_Mxf32bl_mxn_Destroy(owl_Mxf32bl_mxn* M);

//M = 0
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Zero(owl_Mxf32bl_mxn* M, owl_error* ret_error, owl_error const* pass_through_error);

//M = "Diag_mn(a)"
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Diag(owl_Mxf32bl_mxn* M, float a, owl_error* ret_error, owl_error const* pass_through_error);

//A[i, j]
float owl_Mxf32bl_mxn_GetElement(owl_Mxf32bl_mxn const* A, unsigned int i, unsigned int j, owl_error* ret_error, owl_error const* pass_through_error);

//M[i, j] = a
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_SetElement(owl_Mxf32bl_mxn* M, unsigned int i, unsigned int j, float a, owl_error* ret_error, owl_error const* pass_through_error);

//M = A
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Copy(owl_Mxf32bl_mxn* M, owl_Mxf32bl_mxn const* A, owl_error* ret_error, owl_error const* pass_through_error);

//M = A + B
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Add(owl_Mxf32bl_mxn* M, owl_Mxf32bl_mxn const* A, owl_Mxf32bl_mxn const* B, owl_error* ret_error, owl_error const* pass_through_error);

//M = A - B
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Sub(owl_Mxf32bl_mxn* M, owl_Mxf32bl_mxn const* A, owl_Mxf32bl_mxn const* B, owl_error* ret_error, owl_error const* pass_through_error);

//M = a * A
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_ScalarMul(owl_Mxf32bl_mxn* M, owl_Mxf32bl_mxn const* A, float a, owl_error* ret_error, owl_error const* pass_through_error);

//M = A + a * B
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_AddScalarMul(owl_Mxf32bl_mxn* M, owl_Mxf32bl_mxn const* A, owl_Mxf32bl_mxn const* B, float a, owl_error* ret_error, owl_error const* pass_through_error);

//M = A * B
owl_Mxf32bl_mxn* owl_Mxf32bl_mxn_Mul(owl_Mxf32bl_mxn* M, owl_Mxf32bl_mxn const* A, owl_Mxf32bl_mxn const* B, owl_error* ret_error, owl_error const* pass_through_error);

//Tr(A)
float owl_Mxf32bl_mxn_Tr(owl_Mxf32bl_mxn const* A, owl_error* ret_error, owl_error const* pass_through_error);

#endif // OWL_MXF32BL_MXN_H_INCLUDED
