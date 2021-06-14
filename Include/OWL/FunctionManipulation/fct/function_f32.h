#ifndef OWL_FUNCTION_F32_H_INCLUDED
#define OWL_FUNCTION_F32_H_INCLUDED

#include <OWL/owl.h>

typedef struct
{
    //float evaluate(float x)
    float (*evaluate)(float);

    //float derivative(float x)
    float (*derivative)(float);

} owl_derivable_fct_f32;

//
OWL_DLL_EXPORT extern owl_derivable_fct_f32 owl_id_der;
OWL_DLL_EXPORT extern owl_derivable_fct_f32 owl_square_der;
OWL_DLL_EXPORT extern owl_derivable_fct_f32 owl_sqrt_der;
//
OWL_DLL_EXPORT extern owl_derivable_fct_f32 owl_inv_der;
OWL_DLL_EXPORT extern owl_derivable_fct_f32 owl_log_der;
OWL_DLL_EXPORT extern owl_derivable_fct_f32 owl_exp_der;
//
OWL_DLL_EXPORT extern owl_derivable_fct_f32 owl_cos_der;
OWL_DLL_EXPORT extern owl_derivable_fct_f32 owl_sin_der;
OWL_DLL_EXPORT extern owl_derivable_fct_f32 owl_tan_der;

#endif // OWL_FUNCTION_F32_H_INCLUDED
