#include <OWL/Optimized3d/tensor/tensorf32.h>

//Return a^n
static size_t powint(unsigned int a, unsigned int n)
{
    size_t res = 1;
    size_t p = a;

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

//
static size_t tensor_index(owl_Tensorp8f32 const* T, unsigned int const* index_list)
{
    size_t k = 0;
    for(unsigned int m = 0 ; m < T->dim ; m++)
    {
        k = (size_t)T->vector_dim * k + (size_t)index_list[m];
    }

    return k;
}

//Creates a new Tensorp8 and sets all data to 0
//
//
owl_Tensorp8f32* owl_Tensorp8f32_Create(unsigned int vector_dim, unsigned int stored_dim)
{
    owl_Tensorp8f32* T = malloc(sizeof(*T));

    if(T != NULL)
    {
        T->vector_dim = vector_dim;
        T->stored_dim = stored_dim;
        T->dim = stored_dim;
        T->data = NULL;

        if(vector_dim > 0 && stored_dim > 0)
        {
            T->data = owl_v8f32_array_alloc(powint(vector_dim, stored_dim));
            if(T->data != NULL)
            {
                owl_Tensorp8f32_Zero(T);
            }
            else
            {
                owl_Tensorp8f32_Destroy(T);
                T = NULL;
            }
        }
        else
        {
            owl_Tensorp8f32_Destroy(T);
            T = NULL;
        }
    }

    return T;
}

//Destroy a Tensorp8
//
//
void owl_Tensorp8f32_Destroy(owl_Tensorp8f32* T)
{
    if(T != NULL)
    {
        owl_v8f32_array_free(T->data);
        free(T);
    }
}

//T = 0
//
//
owl_Tensorp8f32* owl_Tensorp8f32_Zero(owl_Tensorp8f32* T)
{
    owl_v8f32 zero = owl_v8f32_zero();
    size_t n = powint(T->vector_dim, T->dim);

    for(size_t k = 0 ; k < n ; k++)
    {
        owl_v8f32_store8(T->data + 8*k, zero);
    }

    return T;
}

//T = T1 + T2
//
//
owl_Tensorp8f32* owl_Tensorp8f32_Add(owl_Tensorp8f32* T, owl_Tensorp8f32 const* T1, owl_Tensorp8f32 const* T2)
{
    size_t n = powint(T->vector_dim, T->dim);

    for(size_t k = 0 ; k < n ; k++)
    {
        owl_v8f32_store8(
                            T->data + 8*k,
                            owl_v8f32_add(
                                            owl_v8f32_load8(T1->data + 8*k),
                                            owl_v8f32_load8(T2->data + 8*k)
                                          )
                         );
    }

    return T;
}

//T = T1 - T2
//
//
owl_Tensorp8f32* owl_Tensorp8f32_Sub(owl_Tensorp8f32* T, owl_Tensorp8f32 const* T1, owl_Tensorp8f32 const* T2)
{
    size_t n = powint(T->vector_dim, T->dim);

    for(size_t k = 0 ; k < n ; k++)
    {
        owl_v8f32_store8(
                            T->data + 8*k,
                            owl_v8f32_sub(
                                            owl_v8f32_load8(T1->data + 8*k),
                                            owl_v8f32_load8(T2->data + 8*k)
                                          )
                         );
    }

    return T;
}

//T = a * T1
//
//
owl_Tensorp8f32* owl_Tensorp8f32_ScalarMul(owl_Tensorp8f32* T, owl_Tensorp8f32 const* T1, float a)
{
    size_t n = powint(T->vector_dim, T->dim);

    for(size_t k = 0 ; k < n ; k++)
    {
        owl_v8f32_store8(
                            T->data + 8*k,
                            owl_v8f32_scalar_mul(owl_v8f32_load8(T1->data + 8*k), a)
                         );
    }

    return T;
}

//T = T1
//
//
owl_Tensorp8f32* owl_Tensorp8f32_Copy(owl_Tensorp8f32* T, owl_Tensorp8f32 const* T1)
{
    size_t n = powint(T->vector_dim, T->dim);

    for(size_t k = 0 ; k < n ; k++)
    {
        owl_v8f32_store8(
                            T->data + 8*k,
                            owl_v8f32_load8(T1->data + 8*k)
                         );
    }

    return T;
}

//T = T1 + a*T2
//
//
owl_Tensorp8f32* owl_Tensorp8f32_AddScalarMul(owl_Tensorp8f32* T, owl_Tensorp8f32 const* T1, owl_Tensorp8f32 const* T2, float a)
{
    size_t n = powint(T->vector_dim, T->dim);

    for(size_t k = 0 ; k < n ; k++)
    {
        owl_v8f32_store8(
                            T->data + 8*k,
                            owl_v8f32_add_scalar_mul(
                                                        owl_v8f32_load8(T1->data + 8*k),
                                                        owl_v8f32_load8(T2->data + 8*k),
                                                        a
                                                     )
                         );
    }

    return T;
}

//T[i0,i1,...]
//
//
owl_v8f32 owl_Tensorp8f32_GetComp(owl_Tensorp8f32 const* T, unsigned int const* index_list)
{
    return owl_v8f32_load8(T->data + 8 * tensor_index(T, index_list));
}

//T[i0,i1,...] = f
//
//
owl_Tensorp8f32* OWL_VECTORCALL owl_Tensorp8f32_SetComp(owl_Tensorp8f32* T, unsigned int const* index_list, owl_v8f32 f)
{
    owl_v8f32_store8(T->data + 8 * tensor_index(T, index_list), f);
    return T;
}

//
//
//
owl_Tensorp8f32* owl_Tensorp8f32_Contraction(owl_Tensorp8f32* T, owl_Tensorp8f32 const* T1, unsigned int contraction_index, owl_v8f32 const* v)
{
    unsigned int dim_T = T->dim;

    size_t n = powint(T->vector_dim, dim_T);

    size_t c_pow = powint(T->vector_dim, dim_T - contraction_index);

    for(size_t k = 0 ; k < n ; k++)
    {
        size_t low_order = k % c_pow;
        size_t high_order = k - low_order;

        owl_v8f32 f = owl_v8f32_zero();

        for(unsigned int m = 0 ; m < T->vector_dim ; m++)
        {
            size_t low_order_i = low_order + c_pow * m;
            f = owl_v8f32_addmul(
                                    f,
                                    owl_v8f32_load8(T1->data + 8*(low_order_i + c_pow * T->vector_dim * high_order)),
                                    v[m]
                                 );
        }

        owl_v8f32_store8(T->data + 8*k, f);
    }

    return T;
}
