#include <iostream>
#include <vector>
#include <random>
#include "icicle/runtime.h"
#include "icicle/api/bn254.h"



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