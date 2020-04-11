#ifndef OWL_H_INCLUDED
#define OWL_H_INCLUDED

#define OWL_ALIGN16 __attribute__( (aligned(16)) )

typedef enum
{
    OWL_SUCCESS,

    OWL_MEMORY_ERROR,

    OWL_TAYLOR_ERROR

} owl_error;

#define OWL_PASS_THROUGH_ERROR_VERIFICATION(err, pass_through_err) \
    if((pass_through_err) != NULL && *(pass_through_error) != OWL_SUCCESS) \
    { \
        err = *(pass_through_error); \
    } \
    else

#endif // OWL_H_INCLUDED