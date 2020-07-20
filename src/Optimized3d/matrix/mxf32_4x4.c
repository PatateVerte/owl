#include <OWL/Optimized3d/matrix/mxf32_4x4.h>

owl_mxf32_4x4* owl_mxf32_4x4_zero(owl_mxf32_4x4* M)
{
    for(int j = 0 ; j < 4 ; j++)
    {
        _mm_storeu_ps(M->data + 4*j, _mm_setzero_ps());
    }

    return M;
}

//
//
//
owl_mxf32_4x4* owl_mxf32_4x4_diag(owl_mxf32_4x4* M, float diag_val)
{
    __m128 const diag_val_v = _mm_set_ss(diag_val);

    _mm_storeu_ps(M->data + 0, _mm_insert_ps(diag_val_v, diag_val_v, 0b00001110) );
    _mm_storeu_ps(M->data + 4, _mm_insert_ps(diag_val_v, diag_val_v, 0b00011101) );
    _mm_storeu_ps(M->data + 8, _mm_insert_ps(diag_val_v, diag_val_v, 0b00101011) );
    _mm_storeu_ps(M->data + 12, _mm_insert_ps(diag_val_v, diag_val_v, 0b00110111) );

    return M;
}

//M = A
//
//
owl_mxf32_4x4* owl_mxf32_4x4_copy(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A)
{
    for(int j = 0 ; j < 4 ; j++)
    {
        _mm_storeu_ps(
                        M->data + 4*j,
                        _mm_loadu_ps(A->data + 4*j)
                      );
    }

    return M;
}

//M = A + B
//
//
owl_mxf32_4x4* owl_mxf32_4x4_add(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, owl_mxf32_4x4 const* B)
{
    for(int j = 0 ; j < 4 ; j++)
    {
        __m128 tmp = _mm_add_ps(
                                    _mm_loadu_ps(A->data + 4*j),
                                    _mm_loadu_ps(B->data + 4*j)
                                );
        _mm_storeu_ps(M->data + 4*j, tmp);
    }

    return M;
}

//M = A - B
//
//
owl_mxf32_4x4* owl_mxf32_4x4_sub(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, owl_mxf32_4x4 const* B)
{
    for(int j = 0 ; j < 4 ; j++)
    {
        __m128 tmp = _mm_sub_ps(
                                    _mm_loadu_ps(A->data + 4*j),
                                    _mm_loadu_ps(B->data + 4*j)
                                );
        _mm_storeu_ps(M->data + 4*j, tmp);
    }

    return M;
}

//M = a * A
//
//
owl_mxf32_4x4* owl_mxf32_4x4_scalar_mul(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, float a)
{
    __m128 broadcast = _mm_set1_ps(a);

    for(int j = 0 ; j < 4 ; j++)
    {
        __m128 tmp = _mm_mul_ps(
                                    _mm_loadu_ps(A->data + 4*j),
                                    broadcast
                                );
        _mm_storeu_ps(M->data + 4*j, tmp);
    }

    return M;
}

//M = A + a * B
//
//
owl_mxf32_4x4* owl_mxf32_4x4_add_scalar_mul(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, owl_mxf32_4x4 const* B, float a)
{
    __m128 broadcast = _mm_set1_ps(a);

    for(int j = 0 ; j < 4 ; j++)
    {
        __m128 tmp = _mm_add_ps(
                                    _mm_loadu_ps(A->data + 4*j),
                                    _mm_mul_ps(
                                                _mm_loadu_ps(B->data + 4*j),
                                               broadcast
                                               )
                                );
        _mm_storeu_ps(M->data + 4*j, tmp);
    }

    return M;
}

//M = tA
//
//
owl_mxf32_4x4* owl_mxf32_4x4_transp(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A)
{
    for(int j = 0 ; j < 4 ; j++)
    {
        M->data[4*j + j] = A->data[4*j + j];

        for(int i = 0 ; i < j ; i++)
        {
            float a1 = A->data[4*j + i];
            float a2 = A->data[4*i + j];

            M->data[4*j + i] = a2;
            M->data[4*i + j] = a1;
        }
    }

    return M;
}


//M = A * B
//
//
owl_mxf32_4x4* owl_mxf32_4x4_mul(owl_mxf32_4x4* M, owl_mxf32_4x4 const* A, owl_mxf32_4x4 const* B)
{
    owl_mxf32_4x4 T;
    owl_mxf32_4x4* P = (M == A) ? &T : M;

    for(int j = 0 ; j < 4 ; j++)
    {
        __m128 tmp = _mm_setzero_ps();

        for(int i = 0 ; i < 4 ; i++)
        {
            __m128 broadcast = _mm_set1_ps(B->data[4*j + i]);

            tmp = _mm_add_ps(
                                tmp,
                                _mm_mul_ps(
                                            broadcast,
                                            _mm_loadu_ps(A->data + 4*i)
                                           )
                             );
        }

        _mm_storeu_ps(P->data + 4*j, tmp);
    }

    if(P != M)
    {
        return owl_mxf32_4x4_copy(M, &T);
    }
    else
    {
        return M;
    }
}

//norm2(A)
//
//
float owl_mxf32_4x4_norm2(owl_mxf32_4x4 const* A)
{
    __m128 buffer = _mm_setzero_ps();

    for(int j = 0 ; j < 4 ; j++)
    {
        buffer = _mm_add_ps(
                                buffer,
                                _mm_mul_ps(
                                            _mm_loadu_ps(A->data + 4*j),
                                            _mm_loadu_ps(A->data + 4*j)
                                           )
                            );
    }

    return _mm_cvtss_f32( _mm_dp_ps(buffer, _mm_set1_ps(1.0), 0b11110001) );
}

//norminf(A)
//
//
float owl_mxf32_4x4_norminf(owl_mxf32_4x4 const* A)
{
    __m128 buffer = _mm_setzero_ps();

    for(int j = 0 ; j < 4 ; j++)
    {
        __m128 tmp = _mm_loadu_ps(A->data + 4*j);
        tmp = _mm_max_ps(tmp, _mm_sub_ps(_mm_setzero_ps(), tmp));
        buffer = _mm_max_ps(buffer, tmp);
    }

    __m128 maxv = _mm_max_ps(buffer, _mm_shuffle_ps(buffer, buffer, 0b11110101));
    maxv = _mm_max_ss(maxv, _mm_shuffle_ps(maxv, maxv, 0b11100110));
    return _mm_cvtss_f32(maxv);
}

//Tr(A)
//
//
float owl_mxf32_4x4_trace(owl_mxf32_4x4 const* A)
{
    float trace = 0.0;
    for(int j = 0 ; j < 4 ; j++)
    {
        trace += A->data[4*j + j];
    }

    return trace;
}
