#include <array>
#include <iostream>
#include <vector>
#include <random>
#include "icicle/runtime.h"
#include "icicle/api/bn254.h"
#include "utils.h"
#include "icicle/polynomials/polynomials.h"
#include "Chess.h"
#include "icicle/ntt.h"
#include <botan/aes.h>
#include <botan/auto_rng.h>
#include <botan/cipher_mode.h>
#include <botan/hex.h>
#include <botan/hash.h>
#include "BiVOT.h"


BiVOT::BoardBiVOT(Polynomial_t board, scalar_t root, projective_t enemyBoard){
    //TODO: compute the different pairings for each piece value
    scalar_t index = scalar_t::one();
    for (int row = 0; row < size; ++row) {
        for (int col = 0; col < size; ++col) {
            //int index = row * size + col;
            scalar_t eval = poly(index);
            //TODO: KZG to generate proof
            //TODO: concatenate eval and the proof

             // Create a random number generator
            Botan::AutoSeeded_RNG rng;

            // Generate a random 128-bit (16-byte) key
            std::vector<uint8_t> key = rng.random_vec(16);

            // Create an initialization vector (IV)
            std::vector<uint8_t> iv = rng.random_vec(16);
                    
            //TODO: Encrypt eval and proof with the key

            //SquareBiVOT(cypher, key, board, enemyBoard)
            
            index = index * root;
        }
    }
}