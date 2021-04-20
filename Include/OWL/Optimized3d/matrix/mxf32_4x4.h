#ifndef OWL_MXF32_4X4_H_INCLUDED
#define OWL_MXF32_4X4_H_INCLUDED

#include <OWL/owl.h>

//SSE / 128-bits AVX
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>

#include <string.h>

typedef struct OWL_ALIGN16
{
    float data[16];

} owl_mxf32_4x4;

//M = 0
//Return M
owl_mxf32_4x4* OWL_DLL_EXPORT owl_mxf32_4x4_zero(owl_mxf32_4x4* M);

//
owl_mxf32_4x4* OWL_DLL_EXPORT owl_mxf32_4x4_diag(owl_mxf32_4x4* M, float diag_val);

//
static inline float owl_mxf32_4x4_get_element(owl_mxf32_4x4* M, unsigned int i, unsigned int j)
{
    return M->data[ (j % 4) * 4 + (i % 4) ];
}

//
static inline owl_mxf32_4x4* owl_mxf32_4x4_set_element(owl_mxf32_4x4* M, float value, unsigned int i, unsigned int j)
{
    M->data[ (j % 4) * 4 + (i % 4) ] = value;
    return M;
}

//M = A
owl_mxf32_4x4* OWL_DLL_EXPORT owl_mxf32_4x4_copy(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A);

//M = A + B
owl_mxf32_4x4* OWL_DLL_EXPORT owl_mxf32_4x4_add(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, owl_mxf32_4x4 const* B);

//M = A - B
owl_mxf32_4x4* OWL_DLL_EXPORT owl_mxf32_4x4_sub(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, owl_mxf32_4x4 const* B);

//M = a * A
owl_mxf32_4x4* OWL_DLL_EXPORT owl_mxf32_4x4_scalar_mul(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, float a);

//M = A + a * B
owl_mxf32_4x4* OWL_DLL_EXPORT owl_mxf32_4x4_add_scalar_mul(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, owl_mxf32_4x4 const* B, float a);

//M = tA
owl_mxf32_4x4* OWL_DLL_EXPORT owl_mxf32_4x4_transp(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A);


//M = A * B
owl_mxf32_4x4* OWL_DLL_EXPORT owl_mxf32_4x4_mul(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, owl_mxf32_4x4 const* B);

//M = A
owl_mxf32_4x4* OWL_DLL_EXPORT owl_mxf32_4x4_copy(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A);

//norm2(A)
float OWL_DLL_EXPORT owl_mxf32_4x4_norm2(owl_mxf32_4x4 const* A);

//norminf(A)
float OWL_DLL_EXPORT owl_mxf32_4x4_norminf(owl_mxf32_4x4 const* A);

//Tr(A)
float OWL_DLL_EXPORT owl_mxf32_4x4_trace(owl_mxf32_4x4 const* A);

#endif // OWL_MXF32_4X4_H_INCLUDED
