#include <OWL/Optimized3d/quaternion/q32.h>

#include <math.h>

//q1 * q2
//
//
owl_q32 OWL_VECTORCALL owl_q32_mul(owl_q32 q1, owl_q32 q2)
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
owl_q32 OWL_DLL_EXPORT OWL_VECTORCALL owl_q32_from_rotation(owl_v3f32 v, float alpha)
{
    float cos = cosf(0.5f * alpha);
    float sin = sinf(0.5f * alpha);

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
owl_v3f32 OWL_DLL_EXPORT OWL_VECTORCALL owl_q32_transform_v3f32(owl_q32 q, owl_v3f32 u)
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
                                        a*a - 0.5f
                                     );

    return owl_v3f32_add(rot, rot);

}

//Transform a rotation matrix to a quaternion
//
//
owl_q32 OWL_DLL_EXPORT owl_q32_from_rotation_matrix(owl_mxf32_3x3 const* O)
{
    //r = 0
    __m128 r = _mm_setzero_ps();

    //r33
    r = _mm_addsub_ps(r, _mm_shuffle_ps(O->column[2], O->column[2], 0b10101010));

    //r22
    r = _mm_shuffle_ps(r, r, 0b10110100);
    r = _mm_addsub_ps(r, _mm_shuffle_ps(O->column[1], O->column[1], 0b01010101));

    //r11
    r = _mm_shuffle_ps(r, r, 0b00110110);
    r = _mm_addsub_ps(r, _mm_shuffle_ps(O->column[0], O->column[0], 0b00000000));

    //Reordering
    r = _mm_shuffle_ps(r, r, 0b00101101);

    //Sign mask
    __m128 sign = _mm_setzero_ps();

    //Result 2
    __m128 p2;
    __m128 tmp;
    //r32 +- r23
    tmp = _mm_addsub_ps(
                            _mm_shuffle_ps(O->column[1], O->column[1], 0b10101010),
                            _mm_shuffle_ps(O->column[2], O->column[2], 0b01010101)
                       );
    sign = _mm_insert_ps(sign, tmp, 0b01010000);
    tmp = _mm_shuffle_ps(tmp, tmp, 0b01010000);
    p2 = _mm_mul_ps(tmp, tmp);
    //r13 +- r31
    tmp = _mm_addsub_ps(
                        _mm_shuffle_ps(O->column[2], O->column[2], 0b00000000),
                        _mm_shuffle_ps(O->column[0], O->column[0], 0b10101010)
                     );
    sign = _mm_insert_ps(sign, tmp, 0b01100000);
    tmp = _mm_shuffle_ps(tmp, tmp, 0b01000100);
    p2 = _mm_add_ps(p2, _mm_mul_ps(tmp, tmp));
    //r21 +- r12
    tmp = _mm_addsub_ps(
                            _mm_shuffle_ps(O->column[0], O->column[0], 0b01010101),
                            _mm_shuffle_ps(O->column[1], O->column[1], 0b00000000)
                        );
    sign = _mm_insert_ps(sign, tmp, 0b01110000);
    tmp = _mm_shuffle_ps(tmp, tmp, 0b00010100);
    p2 = _mm_add_ps(p2, _mm_mul_ps(tmp, tmp));
    //
    p2 = _mm_div_ps(
                        p2,
                        _mm_sub_ps(_mm_set1_ps(3.0), r)
                    );

    //Result 1
    __m128 p1 = _mm_add_ps(_mm_set1_ps(1.0), r);

    //Final result
    __m128 mask = _mm_cmpeq_ps(r, _mm_setzero_ps());
    __m128 p = _mm_or_ps(_mm_and_ps(mask, p1), _mm_andnot_ps(mask, p2));
    p = _mm_mul_ps(
                    _mm_set1_ps(0.5),
                    _mm_sqrt_ps(p)
                   );

    __m128 p_plus = p;
    __m128 p_minus = _mm_sub_ps(_mm_setzero_ps(), p);
    sign = _mm_cmpeq_ps(sign, _mm_setzero_ps());
    p = _mm_or_ps(_mm_and_ps(sign, p_plus), _mm_andnot_ps(sign, p_minus));

    return p;
}

owl_q32 OWL_DLL_EXPORT owl_q32_ext_mul(owl_q32 q1, owl_q32 q2) { return owl_q32_mul(q1, q2); }
owl_q32 OWL_DLL_EXPORT owl_q32_ext_from_rotation(owl_v3f32 v, float alpha) { return owl_q32_from_rotation(v, alpha); }
owl_v3f32 OWL_DLL_EXPORT owl_q32_ext_transform_v3f32(owl_q32 q, owl_v3f32 u) {return owl_q32_transform_v3f32(q, u); }
