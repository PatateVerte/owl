#ifndef OWL_VNF32_H_INCLUDED
#define OWL_VNF32_H_INCLUDED

#include <OWL/owl.h>

#include <math.h>

//SSE / 128-bits AVX
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>

#define OWL_HARD_BLOCK_F32_LEN 4

typedef struct
{
    unsigned int n;
    unsigned int nb_hard_blocks;
    float* data;

} owl_Vnf32;

//Create Vnf32
owl_Vnf32* owl_Vnf32_Create(unsigned int n, owl_error* ret_error, owl_error const* pass_through_error);

//Destroy Vnf32
void owl_Vnf32_Destroy(owl_Vnf32* V);

//Vf = 0
owl_Vnf32* owl_Vnf32_Zero(owl_Vnf32* Vf, owl_error* ret_error, owl_error const* pass_through_error);

//Vf = V
owl_Vnf32* owl_Vnf32_Copy(owl_Vnf32* Vf, owl_Vnf32 const* V1, owl_error* ret_error, owl_error const* pass_through_error);

//V1[i]
float owl_Vnf32_GetComponent(owl_Vnf32 const* V1, unsigned int i, owl_error* ret_error, owl_error const* pass_through_error);

//Vf[i] = a
owl_Vnf32* owl_Vnf32_SetComponent(owl_Vnf32* Vf, unsigned int i, float a, owl_error* ret_error, owl_error const* pass_through_error);

//Vf = a * V1
owl_Vnf32* owl_Vnf32_ScalarMul(owl_Vnf32* Vf, owl_Vnf32 const* V1, float a, owl_error* ret_error, owl_error const* pass_through_error);

//Vf = V1 + V2
owl_Vnf32* owl_Vnf32_Add(owl_Vnf32* Vf, owl_Vnf32 const* V1, owl_Vnf32 const* V2, owl_error* ret_error, owl_error const* pass_through_error);

//Vf = V1 - V2
owl_Vnf32* owl_Vnf32_Sub(owl_Vnf32* Vf, owl_Vnf32 const* V1, owl_Vnf32 const* V2, owl_error* ret_error, owl_error const* pass_through_error);

//Vf = V + a * V2
owl_Vnf32* owl_Vnf32_AddScalarMul(owl_Vnf32* Vf, owl_Vnf32 const* V1, owl_Vnf32 const* V2, float a, owl_error* ret_error, owl_error const* pass_through_error);

//Vf = (V1[0] * V2[0], V1[1] * V2[1], ... , V1[n-1] * V2[n-1])
owl_Vnf32* owl_Vnf32_Mul(owl_Vnf32* Vf, owl_Vnf32 const* V1, owl_Vnf32 const* V2, float a, owl_error* ret_error, owl_error const* pass_through_error);

//Dot product
float owl_Vnf32_Dot(owl_Vnf32 const* V1, owl_Vnf32 const* V2, owl_error* ret_error, owl_error const* pass_through_error);

//Gram Schmidt process assuming the vectors from V are independent_list
//All vectors must have the same dimension n, and list_len must be less than or equal to n
owl_Vnf32** owl_Vnf32_GramSchmidt(owl_Vnf32** Vf_list, owl_Vnf32 const** V1_list, unsigned int list_len, owl_error* ret_error, owl_error const* pass_through_error);

#endif // OWL_VNF32_H_INCLUDED
