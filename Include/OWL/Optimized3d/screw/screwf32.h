#ifndef OWL_SCREWF32_H_INCLUDED
#define OWL_SCREWF32_H_INCLUDED

#include <OWL/Optimized3d/vector/v3f32.h>

typedef struct
{
    owl_v3f32 P;
    owl_v3f32 R;
    owl_v3f32 M;

} owl_screwf32;

//
static inline owl_v3f32 owl_screwf32_moment(owl_screwf32 const* T1, owl_v3f32 A)
{
    return owl_v3f32_add(
                            T1->M,
                            owl_v3f32_cross(owl_v3f32_sub(T1->P, A), T1->R)
                         );
}

//(T1)_A
static inline owl_screwf32* owl_screwf32_change_point(owl_screwf32* Tf, owl_screwf32 const* T1, owl_v3f32 A)
{
    Tf->M = owl_screwf32_moment(T1, A);
    Tf->P = A;
    Tf->R = T1->R;

    return Tf;
}

//(T1 + T2)_(T1->P)
static inline owl_screwf32* owl_screwf32_add(owl_screwf32* Tf, owl_screwf32 const* T1, owl_screwf32 const* T2)
{
    Tf->M = owl_v3f32_add(T1->M, owl_screwf32_moment(T2, T1->P));
    Tf->P = T1->P;
    Tf->R = owl_v3f32_add(T1->R, T2->R);

    return Tf;
}

//(T1 + T2)_A
static inline owl_screwf32* owl_screwf32_add_with_point(owl_screwf32* Tf, owl_screwf32 const* T1, owl_screwf32 const* T2, owl_v3f32 A)
{
    Tf->M = owl_v3f32_add(owl_screwf32_moment(T1, A), owl_screwf32_moment(T2, A));
    Tf->P = T1->P;
    Tf->R = owl_v3f32_add(T1->R, T2->R);

    return Tf;
}

//
static inline float owl_screwf32_scalar_product(owl_screwf32 const* T1, owl_screwf32 const* T2)
{
    return owl_v3f32_dot(T1->R, owl_screwf32_moment(T2, T1->P)) + owl_v3f32_dot(T2->R, T1->M);
}

#endif // OWL_SCREWF32_H_INCLUDED
