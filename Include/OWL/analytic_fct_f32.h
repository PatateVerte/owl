#ifndef OWL_ANALYTIC_FCT_F32_H_INCLUDED
#define OWL_ANALYTIC_FCT_F32_H_INCLUDED

typedef struct
{
    //float evaluate(float x)
    float (*evaluate)(float);

    //

} analytic_fct_f32;

#endif // OWL_ANALYTIC_FCT_F32_H_INCLUDED
