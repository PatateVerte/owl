#ifndef OWL_TENSORF32_H_INCLUDED
#define OWL_TENSORF32_H_INCLUDED

#include <OWL/owl.h>

#include <OWL/Optimized3d/vector/v3f32.h>

typedef struct
{
    unsigned int stored_dim;    //>0
    unsigned int dim;   //dim can be changed manually as long as it stays <= stored_dim and > 0

    float* data;

} owl_Tensor3df32;

//Creates a new Tensor3d and sets all data to 0
owl_Tensor3df32* owl_Tensor3df32_Create(unsigned int stored_dim);

//Destroy a Tensor3d
void owl_Tensor3df32_Destroy(owl_Tensor3df32* T);

/*
    BASIC OPERATIONS
*/

//All the following functions assume that dimensions are correct
//T's (destination) dimension will always be trusted above all
//For that reason, they are potentially unsafe

//T = 0
owl_Tensor3df32* owl_Tensor3df32_Zero(owl_Tensor3df32* T);

//T = T1 + T2
owl_Tensor3df32* owl_Tensor3df32_Add(owl_Tensor3df32* T, owl_Tensor3df32 const* T1, owl_Tensor3df32 const* T2);

//T = T1 - T2
owl_Tensor3df32* owl_Tensor3df32_Sub(owl_Tensor3df32* T, owl_Tensor3df32 const* T1, owl_Tensor3df32 const* T2);

//T = a * T1
owl_Tensor3df32* owl_Tensor3df32_ScalarMul(owl_Tensor3df32* T, owl_Tensor3df32 const* T1, float a);

//T = T1
owl_Tensor3df32* owl_Tensor3df32_Copy(owl_Tensor3df32* T, owl_Tensor3df32 const* T1);

//T = T1 + a*T2
owl_Tensor3df32* owl_Tensor3df32_AddScalarMul(owl_Tensor3df32* T, owl_Tensor3df32 const* T1, owl_Tensor3df32 const* T2, float a);

//T[i0,i1,...]
//len(i_list) must be dim
float owl_Tensor3df32_GetComp(owl_Tensor3df32 const* T, unsigned int* i_list);

//T[i0,i1,...] = f
owl_Tensor3df32* owl_Tensor3df32_SetComp(owl_Tensor3df32* T, unsigned int* i_list, float f);

/*
    CONTRACTION
*/

//If T==T1, remember to update T's dimension (=dim(T1) - 1 > 0) before calling, because the function will always assume that dim(T1)=dim(T)+1
owl_Tensor3df32* owl_Tensor3df32_Contraction(owl_Tensor3df32* T, owl_Tensor3df32 const* T1, unsigned int contraction_index, owl_v3f32 v);

//T = T1 with the first vector fixed
//Mort optimized than the general contraction
owl_Tensor3df32* owl_Tensor3df32_LeftContraction(owl_Tensor3df32* T, owl_Tensor3df32 const* T1, owl_v3f32 v0);

#endif // OWL_TENSORF32_H_INCLUDED
