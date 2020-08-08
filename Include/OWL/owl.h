#ifndef OWL_H_INCLUDED
#define OWL_H_INCLUDED

#define OWL_ALIGN16 __attribute__( (aligned(16)) )
#define OWL_VECTORCALL __vectorcall

#define OWL_SQRT2	1.41421356237309504880
#define OWL_SQRT3   1.73205080756887729352
#define OWL_PI		3.14159265358979323846

typedef enum
{
    OWL_SUCCESS,

    OWL_MEMORY_ERROR,

    OWL_DIMENSION_ERROR,
    OWL_OUT_OF_BOUND_ERROR,
    OWL_TAYLOR_ERROR

} owl_error;

#define OWL_PASS_THROUGH_ERROR_VERIFICATION(err, pass_through_err) \
    if((pass_through_err) != NULL && *(pass_through_error) != OWL_SUCCESS) \
    { \
        err = *(pass_through_error); \
    } \
    else

#endif // OWL_H_INCLUDED
