#include <OWL/taylor_series_f32.h>

#include <stdlib.h>
#include <math.h>

//
//
//
owl_TaylorDevf32* owl_TaylorDevf32_Create(unsigned int max_order, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;
    owl_TaylorDevf32* D = NULL;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        D = malloc(sizeof(*D));
        if(D != NULL)
        {
            D->max_order = max_order;
            D->order = 0;
            D->x0 = NAN;

            D->terms = malloc(((size_t)max_order + 1) * sizeof(*D->terms));
            if(D->terms == NULL)
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
        owl_TaylorDevf32_Destroy(D);
        D = NULL;
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return D;
}

//
//
//
void owl_TaylorDevf32_Destroy(owl_TaylorDevf32* D)
{
    if(D != NULL)
    {
        free(D->terms);
        free(D);
    }
}

//D(x)
//
//
float owl_TaylorDevf32_Evaluate(owl_TaylorDevf32 const* D, float x, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    float Dx = 0.0;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        for(unsigned int i = 0 ; i <= D->order ; i++)
        {
            Dx = (x - D->x0) * Dx + D->terms[D->order - i];
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Dx;
}

//Df = 0
//
//
owl_TaylorDevf32* owl_TaylorDevf32_Zero(owl_TaylorDevf32* Df, float x0, unsigned int order, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        unsigned int final_order = Df->max_order;
        if(final_order > order)
        {
            final_order = order;
        }

        Df->order = final_order;
        Df->x0 = x0;

        for(unsigned int i = 0 ; i <= final_order ; i++)
        {
            Df->terms[i] = 0.0;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Df;
}

//Df = a * D1
//
//
owl_TaylorDevf32* owl_TaylorDevf32_ScalarMul(owl_TaylorDevf32* Df, owl_TaylorDevf32 const* D1, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        unsigned int final_order = Df->max_order;
        if(final_order > D1->order)
        {
            final_order = D1->order;
        }

        Df->order = final_order;
        Df->x0 = D1->x0;

        for(unsigned int i = 0 ; i <= final_order ; i++)
        {
            Df->terms[i] = a * D1->terms[i];
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Df;
}

//Df = D1 + a * D2
//
//
owl_TaylorDevf32* owl_TaylorDevf32_AddScalarMul(owl_TaylorDevf32* Df, owl_TaylorDevf32 const* D1, owl_TaylorDevf32 const* D2, float a, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        unsigned int final_order = Df->max_order;
        if(final_order > D1->order)
        {
            final_order = D1->order;
        }
        if(final_order > D2->order)
        {
            final_order = D2->order;
        }

        if(D1->x0 == D2->x0)
        {
            Df->order = final_order;
            Df->x0 = D1->x0;

            for(unsigned int i = 0 ; i <= final_order ; i++)
            {
                Df->terms[i] = D1->terms[i] + a * D2->terms[i];
            }
        }
        else
        {
            error = OWL_TAYLOR_ERROR;
        }
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Df;
}

//Df = D1 * D2
//
//
owl_TaylorDevf32* owl_TaylorDevf32_Mul(owl_TaylorDevf32* Df, owl_TaylorDevf32 const* D1, owl_TaylorDevf32 const* D2, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        unsigned int final_order = Df->max_order;
        if(final_order > D1->order)
        {
            final_order = D1->order;
        }
        if(final_order > D2->order)
        {
            final_order = D2->order;
        }

        owl_TaylorDevf32* D_write = NULL;
        owl_TaylorDevf32* D_buffer = NULL;
        owl_TaylorDevf32 const* D_read1 = D1;
        owl_TaylorDevf32 const* D_read2 = D2;
        if(D2 != Df)
        {
            D_write = Df;
        }
        else if(D1 != Df)
        {
            D_write = Df;
            D_read1 = D2;
            D_read2 = D1;
        }
        else
        {
            D_buffer = owl_TaylorDevf32_Create(final_order, &error, &error);
            D_write = D_buffer;
        }

        if(error == OWL_SUCCESS)
        {
            if(D1->x0 == D2->x0)
            {
                D_write->order = final_order;
                D_write->x0 = D1->x0;

                for(unsigned int i = 0 ; i <= final_order ; i++)
                {
                    float a = 0.0;

                    for(unsigned int i1 = 0 ; i1 <= final_order - i ; i1++)
                    {
                        a += D_read1->terms[i1] * D_read2->terms[final_order - i - i1];
                    }

                    D_write->terms[final_order - i] = a;
                }
            }
            else
            {
                error = OWL_TAYLOR_ERROR;
            }
        }

        if(D_write != Df)
        {
            owl_TaylorDevf32_ScalarMul(Df, D_write, 1.0, &error, &error);
        }

        owl_TaylorDevf32_Destroy(D_buffer);
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Df;
}

//
//
//
owl_TaylorDevf32* owl_TaylorDevf32_Composition(owl_TaylorDevf32* Df, owl_TaylorDevf32 const* D1, owl_TaylorDevf32 const* D2, owl_error* ret_error, owl_error const* pass_through_error)
{
    owl_error error = OWL_SUCCESS;

    OWL_PASS_THROUGH_ERROR_VERIFICATION(error, pass_through_error)
    {
        unsigned int final_order = Df->max_order;
        if(final_order > D1->order)
        {
            final_order = D1->order;
        }
        if(final_order > D2->order)
        {
            final_order = D2->order;
        }

        owl_TaylorDevf32* D_write = NULL;
        owl_TaylorDevf32* D_buffer = NULL;
        if(D1 != Df && D2 != Df)
        {
            D_write = Df;
        }
        else
        {
            D_buffer = owl_TaylorDevf32_Create(final_order, &error, &error);
            D_write = D_buffer;
        }

        owl_TaylorDevf32* D2_pow = owl_TaylorDevf32_Create(final_order, &error, &error);

        if(error == OWL_SUCCESS)
        {
            if(D1->x0 == D2->terms[0])
            {
                D_write->order = final_order;
                D_write->x0 = D2->x0;

                owl_TaylorDevf32_ScalarMul(D2_pow, D2, 0.0, &error, &error);
                D2_pow->terms[0] = 1.0;

                owl_TaylorDevf32_ScalarMul(D_write, D2_pow, D1->terms[0], &error, &error);

                for(unsigned i1 = 1 ; i1 <= final_order && error == OWL_SUCCESS ; i1++)
                {
                    owl_TaylorDevf32_Mul(D2_pow, D2_pow, D2, &error, &error);
                    owl_TaylorDevf32_AddScalarMul(D_write, D_write, D2_pow, D1->terms[i1], &error, &error);
                }
            }
            else
            {
                error = OWL_TAYLOR_ERROR;
            }
        }

        if(Df != D_write)
        {
            owl_TaylorDevf32_ScalarMul(Df, D_write, 1.0, &error, &error);
        }
        owl_TaylorDevf32_Destroy(D_buffer);
        owl_TaylorDevf32_Destroy(D2_pow);
    }

    if(ret_error != NULL)
    {
        *ret_error = error;
    }

    return Df;
}
