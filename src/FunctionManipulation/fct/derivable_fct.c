#include <OWL/FunctionManipulation/fct/function_f32.h>

#include <math.h>

//Id, Square
//
//
OWL_DLL_EXPORT float owl_id_eval(float x);
OWL_DLL_EXPORT float owl_id_deriv(float x);
OWL_DLL_EXPORT float owl_square_eval(float x);
OWL_DLL_EXPORT float owl_square_deriv(float x);
OWL_DLL_EXPORT float owl_sqrt_deriv(float x);

OWL_DLL_EXPORT float owl_id_eval(float x)
{
    return x;
}

OWL_DLL_EXPORT float owl_id_deriv(float x)
{
    (void)x;
    return 1.0f;
}

OWL_DLL_EXPORT float owl_square_eval(float x)
{
    return x*x;
}

OWL_DLL_EXPORT float owl_square_deriv(float x)
{
    return 2.0f * x;
}

OWL_DLL_EXPORT float owl_sqrt_deriv(float x)
{
    return 0.5f / sqrtf(x);
}

OWL_DLL_EXPORT owl_derivable_fct_f32 owl_id_der =
{
    .evaluate = owl_id_eval,
    .derivative = owl_id_deriv
};

OWL_DLL_EXPORT owl_derivable_fct_f32 owl_square_der =
{
    .evaluate = owl_square_eval,
    .derivative = owl_square_deriv
};

OWL_DLL_EXPORT owl_derivable_fct_f32 owl_sqrt_der =
{
    .evaluate = sqrtf,
    .derivative = owl_sqrt_deriv
};

//Inv, Exp, Log
//
//
OWL_DLL_EXPORT float owl_inv_eval(float x);
OWL_DLL_EXPORT float owl_inv_deriv(float x);
OWL_DLL_EXPORT float owl_log_deriv(float x);

OWL_DLL_EXPORT float owl_inv_eval(float x)
{
    return 1.0f / x;
}

OWL_DLL_EXPORT float owl_inv_deriv(float x)
{
    return -1.0f / (x*x);
}

OWL_DLL_EXPORT float owl_log_deriv(float x)
{
    if(x >= 0.0f)
    {
        return 1.0f / x;
    }
    else
    {
        return NAN;
    }
}

OWL_DLL_EXPORT owl_derivable_fct_f32 owl_inv_der =
{
    .evaluate = owl_inv_eval,
    .derivative = owl_inv_deriv
};

OWL_DLL_EXPORT owl_derivable_fct_f32 owl_log_der =
{
    .evaluate = logf,
    .derivative = owl_log_deriv
};

OWL_DLL_EXPORT owl_derivable_fct_f32 owl_exp_der =
{
    .evaluate = expf,
    .derivative = expf
};

//Cos, Sin, Tan
//
//
OWL_DLL_EXPORT float owl_cos_deriv(float x);
OWL_DLL_EXPORT float owl_tan_deriv(float x);

OWL_DLL_EXPORT float owl_cos_deriv(float x)
{
    return -sinf(x);
}

OWL_DLL_EXPORT float owl_tan_deriv(float x)
{
    float t = tanf(x);
    return 1.0f + t*t;
}

OWL_DLL_EXPORT owl_derivable_fct_f32 owl_cos_der =
{
    .evaluate = cosf,
    .derivative = owl_cos_deriv
};

OWL_DLL_EXPORT owl_derivable_fct_f32 owl_sin_der =
{
    .evaluate = sinf,
    .derivative = cosf
};

OWL_DLL_EXPORT owl_derivable_fct_f32 owl_tan_der =
{
    .evaluate = tanf,
    .derivative = owl_tan_deriv
};
