#ifndef OWL_TENSORF32_H_INCLUDED
#define OWL_TENSORF32_H_INCLUDED

#include <OWL/owl.h>

#include <OWL/Optimized3d/vector/v8f32.h>

//8 Tensors in one
typedef struct
{
    unsigned int vector_dim;

    unsigned int stored_dim;
    unsigned int dim;

    float* data;

} owl_Tensorp8f32;

//Creates a new Tensorp8 and sets all data to 0
owl_Tensorp8f32* owl_Tensorp8f32_Create(unsigned int vector_dim, unsigned int stored_dim);

//Destroy a Tensorp8
void owl_Tensorp8f32_Destroy(owl_Tensorp8f32* T);

/*
    BASIC OPERATIONS
*/

//All the following functions assume that dimensions are correct
//T's (destination) dimension will always be trusted above all
//For that reason, they are potentially unsafe

//T = 0
owl_Tensorp8f32* owl_Tensorp8f32_Zero(owl_Tensorp8f32* T);

//T = T1 + T2
owl_Tensorp8f32* owl_Tensorp8f32_Add(owl_Tensorp8f32* T, owl_Tensorp8f32 const* T1, owl_Tensorp8f32 const* T2);

//T = T1 - T2
owl_Tensorp8f32* owl_Tensorp8f32_Sub(owl_Tensorp8f32* T, owl_Tensorp8f32 const* T1, owl_Tensorp8f32 const* T2);

//T = a * T1
owl_Tensorp8f32* owl_Tensorp8f32_ScalarMul(owl_Tensorp8f32* T, owl_Tensorp8f32 const* T1, float a);

//T = T1
owl_Tensorp8f32* owl_Tensorp8f32_Copy(owl_Tensorp8f32* T, owl_Tensorp8f32 const* T1);

//T = T1 + a*T2
owl_Tensorp8f32* owl_Tensorp8f32_AddScalarMul(owl_Tensorp8f32* T, owl_Tensorp8f32 const* T1, owl_Tensorp8f32 const* T2, float a);

//T[i0,i1,...]
//len(index_list) must be == dim
owl_v8f32 owl_Tensorp8f32_GetComp(owl_Tensorp8f32 const* T, unsigned int const* index_list);

//T[i0,i1,...] = f
owl_Tensorp8f32* OWL_VECTORCALL owl_Tensorp8f32_SetComp(owl_Tensorp8f32* T, unsigned int const* index_list, owl_v8f32 f);

/*
    CONTRACTION
*/

//len(v)==vector_dim
//If T==T1, remember to update T's dimension (=dim(T1) - 1 > 0) before calling, because the function will always assume that dim(T1)=dim(T)+1
owl_Tensorp8f32* owl_Tensorp8f32_Contraction(owl_Tensorp8f32* T, owl_Tensorp8f32 const* T1, unsigned int contraction_index, owl_v8f32 const* v);

#endif // OWL_TENSORF32_H_INCLUDED
