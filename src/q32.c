#include <OWL/q32.h>

#include <math.h>

//q1 * q2
//
//
owl_q32 owl_q32_mul(owl_q32 q1, owl_q32 q2)
{
    //a
    __m128 qr = _mm_mul_ps(
                                _mm_shuffle_ps(q1, q1, 0b00000000),
                                _mm_shuffle_ps(q2, q2, 0b10110100)
                             );
    //c
    qr = _mm_addsub_ps(
                        qr,
                        _mm_mul_ps(
                                    _mm_shuffle_ps(q1, q1, 0b10101010),
                                    _mm_shuffle_ps(q2, q2, 0b00011110)
                                   )
                       );
    //d
    qr = _mm_addsub_ps(
                        _mm_shuffle_ps(qr, qr, 0b10011100),
                        _mm_mul_ps(
                                    _mm_shuffle_ps(q1, q1, 0b11111111),
                                    _mm_shuffle_ps(q2, q2, 0b00100111)
                                   )
                      );
    //b
    qr = _mm_addsub_ps(
                        _mm_shuffle_ps(qr, qr, 0b11011000),
                        _mm_mul_ps(
                                    _mm_shuffle_ps(q1, q1, 0b01010101),
                                    _mm_shuffle_ps(q2, q2, 0b10110001)
                                   )
                       );

    return qr;
}

//||v|| = 1
//Rotation of alpha around the unitary vector v
//
owl_q32 owl_q32_rotation(owl_v3f32 v, float alpha)
{
    float cos, sin;
    sincosf(0.5 * alpha, &sin, &cos);

    owl_v3f32 Im = owl_v3f32_scalar_mul(v, sin);
    Im = _mm_shuffle_ps(Im, Im, 0b10010011);

    return _mm_insert_ps(
                            Im,
                            _mm_set_ss(cos),
                            0b00000000
                        );
}

//|q| = 1
//(q) * u * (q^-1)
//
owl_v3f32 owl_q32_transform_v3f32(owl_q32 q, owl_v3f32 u)
{
    float a = owl_q32_Ref(q);
    owl_v3f32 v = owl_q32_Imv(q);

    owl_v3f32 rot = owl_v3f32_scalar_mul(owl_v3f32_cross(v, u), a);
    rot = owl_v3f32_add_scalar_mul(
                                        rot,
                                        v,
                                        owl_v3f32_dot(v, u)
                                     );
    rot = owl_v3f32_add_scalar_mul(
                                        rot,
                                        u,
                                        a*a - 0.5
                                     );

    return owl_v3f32_add(rot, rot);

}
