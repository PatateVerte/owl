#ifndef OWL_MXF32_4X4_H_INCLUDED
#define OWL_MXF32_4X4_H_INCLUDED

#include <OWL/owl.h>

//SSE / 128-bits AVX
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>

#include <string.h>

typedef struct
{
    float data[16];

} __attribute__( (aligned(16)) ) owl_mxf32_4x4;

//M = 0
//Return M
static inline owl_mxf32_4x4* owl_mxf32_4x4_zero(owl_mxf32_4x4* M);

//
owl_mxf32_4x4* owl_mxf32_4x4_diag(owl_mxf32_4x4* M, float diag_val);

//
float owl_mxf32_4x4_get_element(owl_mxf32_4x4* M, int i, int j);

//
owl_mxf32_4x4* owl_mxf32_4x4_set_element(owl_mxf32_4x4* M, float value, int i, int j);

//M = A + B
owl_mxf32_4x4* owl_mxf32_4x4_add(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, owl_mxf32_4x4 const* B);

//M = A - B
owl_mxf32_4x4* owl_mxf32_4x4_sub(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, owl_mxf32_4x4 const* B);

//M = a * A
owl_mxf32_4x4* owl_mxf32_4x4_scalar_mul(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, float a);

//M = tA
owl_mxf32_4x4* owl_mxf32_4x4_transp(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A);


//M = A * B
owl_mxf32_4x4* owl_mxf32_4x4_mul(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, owl_mxf32_4x4 const* B);

//M = A
owl_mxf32_4x4* owl_mxf32_4x4_copy(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A);

//norm2(A)
float owl_mxf32_4x4_norm2(owl_mxf32_4x4 const* A);

//norminf(A)
float owl_mxf32_4x4_norminf(owl_mxf32_4x4 const* A);

//Tr(A)
float owl_mxf32_4x4_trace(owl_mxf32_4x4 const* A);

#endif // OWL_MXF32_4X4_H_INCLUDED
