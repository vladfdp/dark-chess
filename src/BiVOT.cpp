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
//#include <botan/aes.h>
#include <botan/auto_rng.h>
#include <botan/cipher_mode.h>
#include <botan/hex.h>
#include <botan/hash.h>
#include "BiVOT.h"
#include "utils.h"
#include "TwistedG2.h"
#include "Fr.h"
#include "Fp.h"
#include "Fp2.h"
#include "Fp12.h"
#include "EC.h"
#include "Pairing.h"


static std::unique_ptr<scalar_t[]> generate_pows(scalar_t tau, uint32_t size)
{
  auto vec = std::make_unique<scalar_t[]>(size);
  vec[0] = scalar_t::one();
  for (size_t i = 1; i < size; ++i) {
    vec[i] = vec[i - 1] * tau;
  }
  return std::move(vec);
}

std::unique_ptr<affine_t[]> generate_SRS(uint32_t size, scalar_t secret_scalar)
{
  auto gen = projective_t::generator();
  auto pows_of_tau = generate_pows(secret_scalar, size);
  auto SRS = std::make_unique<affine_t[]>(size);
  for (size_t i = 0; i < size; ++i) {
    SRS[i] = projective_t::to_affine(pows_of_tau[i] * gen);
  }
  return std::move(SRS);
}

affine_t commit(Polynomial_t& Q, affine_t* SRS){
    MSMConfig config = default_msm_config();
    auto [viewQ, sizeQ] = Q.get_coefficients_view();

    projective_t commitment;
    msm(viewQ.get(), SRS, sizeQ, config, &commitment);
    affine_t afComm = projective_t::to_affine(commitment);
    return afComm;
}


//create proof that board(index) = eval
affine_t KZGProofGen(Polynomial_t& P, scalar_t eval, scalar_t index, affine_t* SRS){
    auto vec = std::make_unique<scalar_t[]>(2);
    vec[0] = scalar_t::neg(index);
    vec[1] = scalar_t::one();

    P.sub_monomial_inplace(eval);
    Polynomial_t divident = Polynomial_t::from_coefficients(vec.get(), 2);
    Polynomial_t Q = P / divident;
    P.add_monomial_inplace(eval); //get it back to normal

    return commit(Q, SRS);
}

bool KZGProofVerif(affine_t commitment, affine_t proof, g2_affine_t g2_srs, scalar_t eval, scalar_t index){

    projective_t scaled_generator_g1 = eval * projective_t::generator();
    
    // Calculate the difference between the commitment and the scaled G1 generator
    affine_t commitment_affine_diff = projective_t::to_affine(projective_t::from_affine(commitment) - scaled_generator_g1);
    g1 commitmentDifference = convert_affine_to_g1(commitment_affine_diff);

    // Compute the left-hand side of the pairing equation
    g2_affine_t g2_gen = g2_projective_t::to_affine(g2_projective_t::generator());
    Fp12 lhs = Tate_pairing(commitmentDifference, convert_g2_affine_to_g2(g2_gen));

    // Compute index * G2 generator
    g2_projective_t scaled_generator_g2 = index * g2_projective_t::generator();
    
    // Calculate the difference between the G2 SRS and the scaled G2 generator
    g2_affine_t g2_affine_diff = g2_projective_t::to_affine(g2_projective_t::from_affine(g2_srs) - scaled_generator_g2);
    g2 g2_offset = convert_g2_affine_to_g2(g2_affine_diff);

    // Compute the right-hand side of the pairing equation
    Fp12 rhs = Tate_pairing(convert_affine_to_g1(proof), g2_offset);

    // Return true if both sides are equal, proving the KZG proof
    return Fp12_equal(lhs, rhs);
}


//fills the board with (piece, proof)
BoardProofs MakeBoardProofs(Polynomial_t& board, scalar_t root, affine_t* srs){

    BoardProofs result;
    scalar_t index = scalar_t::one();
    for(int i = 0; i < TOTAL_SQUARES; ++i){
        scalar_t eval = board(index);
        affine_t proof = KZGProofGen(board, eval, index, srs);
        result.board[i] = {eval.limbs_storage.limbs[0], proof};
        index = index * root;
    }
    return result;
}


// BoardBiVOT::BoardBiVOT(Polynomial_t& board, projective_t& enemyBoard, scalar_t root, affine_t* srs){
//     //TODO: compute the different pairings for each piece value
//     BoardProofs board_proofs MakeBoardProofs(board, root, srs);

//     scalar_t index = scalar_t::one();
//     for (int row = 0; row < BOARD_SIZE; ++row) {
//         for (int col = 0; col < BOARD_SIZE; ++col) {
//             //int index = row * size + col;
//             scalar_t eval = board(index);
//             //TODO: KZG to generate proof
//             //TODO: make eval and the proof

//              // Create a random number generator
//             Botan::AutoSeeded_RNG rng;

//             // // Generate a random 128-bit (16-byte) key
//             // std::vector<uint8_t> key = rng.random_vec(16);

//             // // Create an initialization vector (IV)
//             // std::vector<uint8_t> iv = rng.random_vec(16);
                    
//             //TODO: Encrypt eval and proof with the key

//             //SquareBiVOT(cypher, key, board, enemyBoard)
            
//             index = index * root;
//         }
//     }
// }