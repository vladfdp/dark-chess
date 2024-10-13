#include "utils.h"

#include <iostream>
#include <vector>
#include <random>
#include "icicle/runtime.h"
#include "icicle/api/bn254.h"
#include "TwistedG2.h"
#include "Fr.h"
#include "Fp.h"
#include "Fp2.h"
#include "EC.h"

using namespace bn254;

storage<8> shift_left(const storage<8>& input) {
    storage<8> result;
    uint32_t carry = 0;

    for (int i = 7; i >= 0; --i) {
        // Shift the current limb left by 1
        uint32_t new_limb = (input.limbs[i] << 1);


        // If this is not the last limb, XOR the carry into the next limb
        if (i) {
          // Calculate the new carry
        new_limb ^= input.limbs[i-1] >> 31;
        }

        result.limbs[i] = new_limb;
    }

    return result;
}


Fr convert_scalar_to_fr(scalar_t& scalar) {
    Fr result;
    uint32_t scalar_limbs[8];
    std::copy(scalar.export_limbs(), scalar.export_limbs() + 8, scalar_limbs);

    // Convert from 8 32-bit limbs (little-endian) to 4 64-bit limbs (big-endian)
    result.num.u0 = (static_cast<uint64_t>(scalar_limbs[1]) << 32) | scalar_limbs[0];
    result.num.u1 = (static_cast<uint64_t>(scalar_limbs[3]) << 32) | scalar_limbs[2];
    result.num.u2 = (static_cast<uint64_t>(scalar_limbs[5]) << 32) | scalar_limbs[4];
    result.num.u3 = (static_cast<uint64_t>(scalar_limbs[7]) << 32) | scalar_limbs[6];

    return result;
}

scalar_t convert_fr_to_scalar(const Fr& fr) {
    scalar_t result;

    // Convert from 4 64-bit limbs (big-endian) to 8 32-bit limbs (little-endian)
    result.limbs_storage.limbs[0] = static_cast<uint32_t>(fr.num.u0);
    result.limbs_storage.limbs[1] = static_cast<uint32_t>(fr.num.u0 >> 32);
    result.limbs_storage.limbs[2] = static_cast<uint32_t>(fr.num.u1);
    result.limbs_storage.limbs[3] = static_cast<uint32_t>(fr.num.u1 >> 32);
    result.limbs_storage.limbs[4] = static_cast<uint32_t>(fr.num.u2);
    result.limbs_storage.limbs[5] = static_cast<uint32_t>(fr.num.u2 >> 32);
    result.limbs_storage.limbs[6] = static_cast<uint32_t>(fr.num.u3);
    result.limbs_storage.limbs[7] = static_cast<uint32_t>(fr.num.u3 >> 32);

    return result;

}



Fp convert_point_field_to_fp(point_field_t& elem){
    Fp result;
    uint32_t scalar_limbs[8];
    std::copy(elem.export_limbs(), elem.export_limbs() + 8, scalar_limbs);

    // Convert from 8 32-bit limbs (little-endian) to 4 64-bit limbs (big-endian)
    result.num.u0 = (static_cast<uint64_t>(scalar_limbs[1]) << 32) | scalar_limbs[0];
    result.num.u1 = (static_cast<uint64_t>(scalar_limbs[3]) << 32) | scalar_limbs[2];
    result.num.u2 = (static_cast<uint64_t>(scalar_limbs[5]) << 32) | scalar_limbs[4];
    result.num.u3 = (static_cast<uint64_t>(scalar_limbs[7]) << 32) | scalar_limbs[6];

    return result;
}

point_field_t convert_fp_to_point_field(const Fp& fp) {
    point_field_t result;

    // Convert from 4 64-bit limbs (big-endian) to 8 32-bit limbs (little-endian)
    result.limbs_storage.limbs[0] = static_cast<uint32_t>(fp.num.u0);
    result.limbs_storage.limbs[1] = static_cast<uint32_t>(fp.num.u0 >> 32);
    result.limbs_storage.limbs[2] = static_cast<uint32_t>(fp.num.u1);
    result.limbs_storage.limbs[3] = static_cast<uint32_t>(fp.num.u1 >> 32);
    result.limbs_storage.limbs[4] = static_cast<uint32_t>(fp.num.u2);
    result.limbs_storage.limbs[5] = static_cast<uint32_t>(fp.num.u2 >> 32);
    result.limbs_storage.limbs[6] = static_cast<uint32_t>(fp.num.u3);
    result.limbs_storage.limbs[7] = static_cast<uint32_t>(fp.num.u3 >> 32);

    return result;

}


g1 convert_affine_to_g1(affine_t& p){
  g1 result;
  result.x = convert_point_field_to_fp(p.x);
  result.y = convert_point_field_to_fp(p.y);
  
  return result;
}


Fp2 convert_g2_point_field_to_fp2(g2_point_field_t& elem){
    Fp2 result;
    result.x0 = convert_point_field_to_fp(elem.real);
    result.x1 = convert_point_field_to_fp(elem.imaginary);
    return result;
}

g2_point_field_t convert_fp2_to_g2_point_field(Fp2& elem){
    g2_point_field_t result;
    result.real = convert_fp_to_point_field(elem.x0);
    result.imaginary = convert_fp_to_point_field(elem.x1);
    return result;
}

TwistedG2 convert_g2_affine_to_twistedg2(g2_affine_t& p){
    TwistedG2 result;
    result.x = convert_g2_point_field_to_fp2(p.x);
    result.y = convert_g2_point_field_to_fp2(p.y);
    return result;
}

g2 convert_g2_affine_to_g2(g2_affine_t& p){
    TwistedG2 result = convert_g2_affine_to_twistedg2(p);
    return g2_untwist(result);
}


g2_affine_t convert_twistedg2_to_g2_affine(TwistedG2& p){
    g2_affine_t result;
    result.x = convert_fp2_to_g2_point_field(p.x);
    result.y = convert_fp2_to_g2_point_field(p.y);
    return result;
}