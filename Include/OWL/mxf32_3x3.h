#ifndef OWL_MXF32_3X3_H_INCLUDED
#define OWL_MXF32_3X3_H_INCLUDED

#include <OWL/owl.h>
#include <OWL/v3f32.h>

//SSE / 128-bits AVX
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>

#include <string.h>

typedef struct
{
    owl_v3f32 column[3];

} owl_mxf32_3x3;

//M = 0
//Return M
static inline owl_mxf32_3x3* owl_mxf32_3x3_zero(owl_mxf32_3x3* M)
{
	for(int j = 0 ; j < 3 ; j++)
    {
        M->column[j] = owl_v3f32_zero();
    }

    return M;
}

//
owl_mxf32_3x3* owl_mxf32_3x3_diag(owl_mxf32_3x3* M, float diag_val);

//
float owl_mxf32_3x3_get_element(owl_mxf32_3x3* M, unsigned int i, unsigned int j);

//
owl_mxf32_3x3* owl_mxf32_3x3_set_element(owl_mxf32_3x3* M, float value, unsigned int i, unsigned int j);

//M = A + B
owl_mxf32_3x3* owl_mxf32_3x3_add(owl_mxf32_3x3* M, owl_mxf32_3x3 const* A, owl_mxf32_3x3 const* B);

//M = A - B
owl_mxf32_3x3* owl_mxf32_3x3_sub(owl_mxf32_3x3* M, owl_mxf32_3x3 const* A, owl_mxf32_3x3 const* B);

//M = a * A
owl_mxf32_3x3* owl_mxf32_3x3_scalar_mul(owl_mxf32_3x3* M, owl_mxf32_3x3 const* A, float a);

//M = tA
owl_mxf32_3x3* owl_mxf32_3x3_transp(owl_mxf32_3x3* M, owl_mxf32_3x3 const* A);

//Transform vector with A
owl_v3f32 owl_mxf32_3x3_transform_v3f32(owl_mxf32_3x3 const* A, owl_v3f32 v);


//M = A * B
owl_mxf32_3x3* owl_mxf32_3x3_mul(owl_mxf32_3x3* M, owl_mxf32_3x3 const* A, owl_mxf32_3x3 const* B);

//M = A
static inline owl_mxf32_3x3* owl_mxf32_3x3_copy(owl_mxf32_3x3* M, owl_mxf32_3x3 const* A)
{
    for(int j = 0 ; j < 3 ; j++)
    {
        M->column[j] = A->column[j];
    }
    return M;
}

//norm2(A)
float owl_mxf32_3x3_norm2(owl_mxf32_3x3 const* A);

//norminf(A)
float owl_mxf32_3x3_norminf(owl_mxf32_3x3 const* A);

//Tr(A)
static inline float owl_mxf32_3x3_trace(owl_mxf32_3x3 const* A)
{
    return owl_v3f32_unsafe_get_component(A->column[0], 0) +
            owl_v3f32_unsafe_get_component(A->column[1], 1) +
            owl_v3f32_unsafe_get_component(A->column[2], 2);
}

//det(A)
static inline float owl_mxf32_3x3_det(owl_mxf32_3x3 const* A)
{
    return owl_v3f32_triple(A->column[0], A->column[1], A->column[2]);
}

//M = A ^(-1)
owl_mxf32_3x3* owl_mxf32_3x3_Inv(owl_mxf32_3x3* M, owl_mxf32_3x3 const* A);

//Return dominant eigenvalue
float owl_mxf32_3x3_sym_dominant_eigenvalue(owl_v3f32* eigenvector_ptr, owl_mxf32_3x3 const* A);

//A = P * D * tP with A symmetric and D=diag(eigenvalue_list)
//Parameter P is optional
//P is an orthogonal direct base of R^3
//Return eigenvalue_list
float* owl_mxf32_3x3_sym_diagonalize(float* eigenvalue_list, owl_mxf32_3x3* P, owl_mxf32_3x3 const* A);

#endif // OWL_MXF32_3X3_H_INCLUDED
