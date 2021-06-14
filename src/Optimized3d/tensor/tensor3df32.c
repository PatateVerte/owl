#include <OWL/Optimized3d/tensor/tensor3df32.h>

static inline size_t pow3(unsigned int n)
{
    size_t res = 1;
    size_t p = 3;

    for(unsigned int n_ = n ; n_ != 0 ; n_ >>= 1u)
    {
        if((n_ & 1u) != 0)
        {
            res *= p;
        }

        p *= p;
    }

    return res;
}

//Creates a new Tensor3d and sets all data to 0
//
//
owl_Tensor3df32* owl_Tensor3df32_Create(unsigned int stored_dim)
{
    owl_Tensor3df32* T = malloc(sizeof(*T));

    if(T != NULL)
    {
        T->stored_dim = stored_dim;
        T->dim = stored_dim;
        T->data = NULL;

        if(stored_dim > 0)
        {
            T->data = owl_v3f32_array_alloc(pow3(stored_dim - 1));
            if(T->data != NULL)
            {
                owl_Tensor3df32_Zero(T);
            }
            else
            {
                owl_Tensor3df32_Destroy(T);
                T = NULL;
            }
        }
        else
        {
            owl_Tensor3df32_Destroy(T);
            T = NULL;
        }
    }

    return T;
}

//Destroy a Tensor3d
//
//
void owl_Tensor3df32_Destroy(owl_Tensor3df32* T)
{
    if(T != NULL)
    {
        owl_v3f32_array_free(T->data);
        free(T);
    }
}

//T = 0
//
//
owl_Tensor3df32* owl_Tensor3df32_Zero(owl_Tensor3df32* T)
{
    owl_v3f32 zero = owl_v3f32_zero();
    size_t n = pow3(T->dim - 1);

    for(size_t k = 0 ; k < n ; k++)
    {
        owl_v3f32_store4(T->data + 4*k, zero);
    }

    return T;
}

//T = T1 + T2
//
//
owl_Tensor3df32* owl_Tensor3df32_Add(owl_Tensor3df32* T, owl_Tensor3df32 const* T1, owl_Tensor3df32 const* T2)
{
    size_t n = pow3(T->dim - 1);

    for(size_t k = 0 ; k < n ; k++)
    {
        owl_v3f32_store4(
                            T->data + 4*k,
                            owl_v3f32_add(
                                            owl_v3f32_load4(T1->data + 4*k),
                                            owl_v3f32_load4(T2->data + 4*k)
                                          )
                         );
    }

    return T;
}

//T = T1 - T2
//
//
owl_Tensor3df32* owl_Tensor3df32_Sub(owl_Tensor3df32* T, owl_Tensor3df32 const* T1, owl_Tensor3df32 const* T2)
{
    size_t n = pow3(T->dim - 1);

    for(size_t k = 0 ; k < n ; k++)
    {
        owl_v3f32_store4(
                            T->data + 4*k,
                            owl_v3f32_sub(
                                            owl_v3f32_load4(T1->data + 4*k),
                                            owl_v3f32_load4(T2->data + 4*k)
                                          )
                         );
    }

    return T;
}

//T = a * T1
//
//
owl_Tensor3df32* owl_Tensor3df32_ScalarMul(owl_Tensor3df32* T, owl_Tensor3df32 const* T1, float a)
{
    size_t n = pow3(T->dim - 1);

    for(size_t k = 0 ; k < n ; k++)
    {
        owl_v3f32_store4(
                            T->data + 4*k,
                            owl_v3f32_scalar_mul(owl_v3f32_load4(T1->data + 4*k), a)
                         );
    }

    return T;
}

//T = T1
//
//
owl_Tensor3df32* owl_Tensor3df32_Copy(owl_Tensor3df32* T, owl_Tensor3df32 const* T1)
{
    size_t n = pow3(T->dim - 1);

    for(size_t k = 0 ; k < n ; k++)
    {
        owl_v3f32_store4(
                            T->data + 4*k,
                            owl_v3f32_load4(T1->data + 4*k)
                         );
    }

    return T;
}

//T = T1 + a*T2
//
//
owl_Tensor3df32* owl_Tensor3df32_AddScalarMul(owl_Tensor3df32* T, owl_Tensor3df32 const* T1, owl_Tensor3df32 const* T2, float a)
{
    size_t n = pow3(T->dim - 1);

    for(size_t k = 0 ; k < n ; k++)
    {
        owl_v3f32_store4(
                            T->data + 4*k,
                            owl_v3f32_add_scalar_mul(
                                                        owl_v3f32_load4(T1->data + 4*k),
                                                        owl_v3f32_load4(T2->data + 4*k),
                                                        a
                                                    )
                         );
    }

    return T;
}

//T[i0,i1,...]
//
//
float owl_Tensor3df32_GetComp(owl_Tensor3df32 const* T, unsigned int* i_list)
{
    size_t k = 0;
    for(unsigned int m = 0 ; m < T->dim - 1 ; m++)
    {
        k = 3*k + i_list[m];
    }
    k = 4*k + i_list[T->dim - 1];

    return T->data[k];
}

//T[i0,i1,...] = f
//
//
owl_Tensor3df32* owl_Tensor3df32_SetComp(owl_Tensor3df32* T, unsigned int* i_list, float f)
{
    size_t k = 0;
    for(unsigned int m = 0 ; m < T->dim - 1 ; m++)
    {
        k = 3*k + i_list[m];
    }
    k = 4*k + i_list[T->dim - 1];

    T->data[k] = f;
    return T;
}

//
//
//
owl_Tensor3df32* owl_Tensor3df32_Contraction(owl_Tensor3df32* T, owl_Tensor3df32 const* T1, unsigned int contraction_index, owl_v3f32 v)
{
    unsigned int dim_T = T->dim;

    size_t n = pow3(dim_T - 1);

    size_t c_pow = pow3(dim_T - contraction_index);

    float OWL_ALIGN16 v_coords[4];
    owl_v3f32_store4(v_coords, v);

    for(size_t k = 0 ; k < n ; k++)
    {
        size_t low_order = k % c_pow;
        size_t high_order = k - low_order;

        float f = 0.0f;
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            size_t low_order_i = low_order + c_pow*i;
            f += T1->data[(low_order_i%3) + 4*(low_order_i/3) + 4*high_order] * v_coords[i];
        }
        T->data[(k%3) + 4*(k/3)] = f;
    }

    return T;
}

//
//
//
owl_Tensor3df32* owl_Tensor3df32_LeftContraction(owl_Tensor3df32* T, owl_Tensor3df32 const* T1, owl_v3f32 v0)
{
    unsigned int dim_T = T->dim;
    size_t h_pow = pow3(dim_T - 1);

    float const v0_0 = owl_v3f32_unsafe_get_component(v0, 0);
    float const v0_1 = owl_v3f32_unsafe_get_component(v0, 1);
    float const v0_2 = owl_v3f32_unsafe_get_component(v0, 2);

    for(size_t k = 0 ; k < h_pow ; k++)
    {
        owl_v3f32 s = owl_v3f32_scalar_mul(
                                            owl_v3f32_load4(T1->data + 4*k),
                                            v0_0
                                           );
        s = owl_v3f32_add_scalar_mul(
                                        s,
                                        owl_v3f32_load4(T1->data + 4*(k + h_pow)),
                                        v0_1
                                     );
        s = owl_v3f32_add_scalar_mul(
                                        s,
                                        owl_v3f32_load4(T1->data + 4*(k + 2*h_pow)),
                                        v0_2
                                     );
        owl_v3f32_store4(T->data + 4*k, s);
    }

    return T;
}
