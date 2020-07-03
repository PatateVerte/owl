#include <OWL/vnf32.h>

#include <stdlib.h>
#include <malloc.h>

//Create Vnf32
//
//
owl_Vnf32* owl_Vnf32_Create(unsigned int n, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;
    owl_Vnf32* V = NULL;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        V = malloc(sizeof(*V));
        if(V != NULL)
        {
            V->n = n;
            V->nb_hard_blocks = (n % OWL_HARD_BLOCK_F32_LEN == 0) ? (n / OWL_HARD_BLOCK_F32_LEN) : (n / OWL_HARD_BLOCK_F32_LEN) + 1;

            V->data = _aligned_malloc(OWL_HARD_BLOCK_F32_LEN * V->nb_hard_blocks * sizeof(*(V->data)), 16);
            if(V->data != NULL)
            {
                owl_Vnf32_Zero(V, &error, &error);
            }
            else
            {
                error = OWL_MEMORY_ERROR;
            }
        }
        else
        {
            error = OWL_MEMORY_ERROR;
        }
    }

    if(error != OWL_SUCCESS)
    {
        owl_Vnf32_Destroy(V);
        V = NULL;
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return V;
}

//Destroy Vnf32
//
//
void owl_Vnf32_Destroy(owl_Vnf32* V)
{
    if(V != NULL)
    {
        _aligned_free(V->data);
        free(V);
    }
}

//Vf = 0
//
//
owl_Vnf32* owl_Vnf32_Zero(owl_Vnf32* Vf, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        for(unsigned int b = 0 ; b < Vf->nb_hard_blocks ; b++)
        {
            unsigned int i = b * OWL_HARD_BLOCK_F32_LEN;
            _mm_storeu_ps(Vf->data + i, _mm_setzero_ps());
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Vf;
}

//Vf = V1
//
//
owl_Vnf32* owl_Vnf32_Copy(owl_Vnf32* Vf, owl_Vnf32 const* V1, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(Vf->n == V1->n)
        {
            for(unsigned int b = 0 ; b < Vf->nb_hard_blocks ; b++)
            {
                unsigned int i = b * OWL_HARD_BLOCK_F32_LEN;
                _mm_storeu_ps(
                                Vf->data + i,
                                _mm_loadu_ps(V1->data + i)
                              );
            }
        }
        else
        {
            error = OWL_DIMENSION_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Vf;
}

//V1[i]
//
//
float owl_Vnf32_GetComponent(owl_Vnf32 const* V1, unsigned int i, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    float a = 0.0;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(i < V1->n)
        {
            a = V1->data[i];
        }
        else
        {
            error = OWL_OUT_OF_BOUND_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return a;
}

//V1[i] = a
//
//
owl_Vnf32* owl_Vnf32_SetComponent(owl_Vnf32* Vf, unsigned int i, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(i < Vf->n)
        {
            Vf->data[i] = a;
        }
        else
        {
            error = OWL_OUT_OF_BOUND_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Vf;
}


//Vf = a * V1
//
//
owl_Vnf32* owl_Vnf32_ScalarMul(owl_Vnf32* Vf, owl_Vnf32 const* V1, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(Vf->n == V1->n)
        {
            __m128 broadcast = _mm_set1_ps(a);

            for(unsigned int b = 0 ; b < Vf->nb_hard_blocks ; b++)
            {
                unsigned int i = b * OWL_HARD_BLOCK_F32_LEN;
                _mm_storeu_ps(
                                Vf->data + i,
                                _mm_mul_ps(
                                            _mm_loadu_ps(V1->data + i),
                                           broadcast
                                           )
                              );
            }
        }
        else
        {
            error = OWL_DIMENSION_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Vf;
}

//Vf = V1 + V2
//
//
owl_Vnf32* owl_Vnf32_Add(owl_Vnf32* Vf, owl_Vnf32 const* V1, owl_Vnf32 const* V2, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(Vf->n == V1->n && Vf->n == V2->n)
        {
            for(unsigned int b = 0 ; b < Vf->nb_hard_blocks ; b++)
            {
                unsigned int i = b * OWL_HARD_BLOCK_F32_LEN;
                _mm_storeu_ps(
                                Vf->data + i,
                                _mm_add_ps(
                                            _mm_loadu_ps(V1->data + i),
                                            _mm_loadu_ps(V2->data + i)
                                           )
                              );
            }
        }
        else
        {
            error = OWL_DIMENSION_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Vf;
}

//Vf = V1 - V2
//
//
owl_Vnf32* owl_Vnf32_Sub(owl_Vnf32* Vf, owl_Vnf32 const* V1, owl_Vnf32 const* V2, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(Vf->n == V1->n && Vf->n == V2->n)
        {
            for(unsigned int b = 0 ; b < Vf->nb_hard_blocks ; b++)
            {
                unsigned int i = b * OWL_HARD_BLOCK_F32_LEN;
                _mm_storeu_ps(
                                Vf->data + i,
                                _mm_sub_ps(
                                            _mm_loadu_ps(V1->data + i),
                                            _mm_loadu_ps(V2->data + i)
                                           )
                              );
            }
        }
        else
        {
            error = OWL_DIMENSION_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Vf;
}

//Vf = V1 + a * V2
//
//
owl_Vnf32* owl_Vnf32_AddScalarMul(owl_Vnf32* Vf, owl_Vnf32 const* V1, owl_Vnf32 const* V2, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(Vf->n == V1->n && Vf->n == V2->n)
        {
            __m128 broadcast = _mm_set1_ps(a);

            for(unsigned int b = 0 ; b < Vf->nb_hard_blocks ; b++)
            {
                unsigned int i = b * OWL_HARD_BLOCK_F32_LEN;
                _mm_storeu_ps(
                                Vf->data + i,
                                _mm_add_ps(
                                            _mm_loadu_ps(V1->data + i),
                                            _mm_mul_ps(
                                                        _mm_loadu_ps(V2->data + i),
                                                        broadcast
                                                       )
                                           )
                              );
            }
        }
        else
        {
            error = OWL_DIMENSION_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Vf;
}

//Vf = (V1[0] * V2[0], V1[1] * V2[1], ... , V1[n-1] * V2[n-1])
//
//
owl_Vnf32* owl_Vnf32_Mul(owl_Vnf32* Vf, owl_Vnf32 const* V1, owl_Vnf32 const* V2, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(Vf->n == V1->n && Vf->n == V2->n)
        {
            for(unsigned int b = 0 ; b < Vf->nb_hard_blocks ; b++)
            {
                unsigned int i = b * OWL_HARD_BLOCK_F32_LEN;
                _mm_storeu_ps(
                                Vf->data + i,
                                _mm_mul_ps(
                                            _mm_loadu_ps(V1->data + i),
                                            _mm_loadu_ps(V2->data + i)
                                           )
                              );
            }
        }
        else
        {
            error = OWL_DIMENSION_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Vf;
}

//Dot product
//
//
float owl_Vnf32_Dot(owl_Vnf32 const* V1, owl_Vnf32 const* V2, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;
    float dot = 0.0;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        if(V1->n == V2->n)
        {
            __m128 acc = _mm_setzero_ps();

            for(unsigned int b = 0 ; b < V1->nb_hard_blocks ; b++)
            {
                unsigned int i = b * OWL_HARD_BLOCK_F32_LEN;
                acc = _mm_add_ps(
                                    acc,
                                    _mm_mul_ps(
                                                _mm_loadu_ps(V1->data + i),
                                                _mm_loadu_ps(V2->data + i)
                                               )
                                 );
            }

            dot = _mm_cvtss_f32( _mm_dp_ps(acc, _mm_set1_ps(1.0), 0b11110001) );
        }
        else
        {
            error = OWL_DIMENSION_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return dot;
}

