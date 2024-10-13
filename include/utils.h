#include <iostream>
#include <vector>
#include <random>
#include "icicle/runtime.h"
#include "icicle/api/bn254.h"
#include "TwistedG2.h"
#include "Fr.h"
#include "Fp.h"
#include "EC.h"

#ifndef UTILS_H
#define UTILS_H

using namespace bn254;

// typedef ComplexExtensionField<fq_config, Field<fq_config>> ext_field;

storage<8> shift_left(const storage<8>& input);


Fr convert_scalar_to_fr(scalar_t& scalar);

scalar_t convert_fr_to_scalar(const Fr& fr);

Fp convert_point_field_to_fp(point_field_t& elem);

point_field_t convert_fp_to_point_field(const Fp& fp);

g1 convert_affine_to_g1(affine_t& p);

affine_t convert_g1_to_affine(g1& p);

Fp2 convert_g2_point_field_to_fp2(g2_point_field_t& elem);

TwistedG2 convert_g2_affine_to_twistedg2(g2_affine_t& p);

g2 convert_g2_affine_to_g2(g2_affine_t& p);

g2_affine_t convert_twistedg2_to_g2_affine(TwistedG2& p);





#endif