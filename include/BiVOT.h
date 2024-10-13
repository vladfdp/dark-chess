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
#include "EC.h"

constexpr size_t FIXED_CIPHER_LENGTH = 32;
constexpr size_t FIXED_SQUARE_LENGTH = 32 + (32 * 12);

struct EncryptedMessage {
    int attackingSquare;
    std::array<uint8_t, FIXED_CIPHER_LENGTH> message;
};

struct EncryptedSquare {
    std::array<uint8_t, FIXED_SQUARE_LENGTH> square;
};

struct DecryptedSquare {
    uint32_t piece;
    affine_t proof;
};


struct AttackVectors{
    std::vector<EncryptedMessage> av;
};

struct SquareBiVOT{
    EncryptedSquare square;
    AttackVectors King; //TODO: add other pieces
};

struct BoardProofs{
    std::array<DecryptedSquare,TOTAL_SQUARES> board;
};

struct BoardBiVOT{
    std::array<SquareBiVOT,TOTAL_SQUARES> Pieces;
    //std::array<SBiVOTVariant,TOTAL_SQUARES> Positions;


    // BoardBiVOT(Polynomial_t board, scalar_t root, projective_t enemyBoard);

};


static std::unique_ptr<scalar_t[]> generate_pows(scalar_t tau, uint32_t size);

std::unique_ptr<affine_t[]> generate_SRS(uint32_t size, scalar_t secret_scalar);

affine_t commit(Polynomial_t& P, affine_t* SRS);


affine_t KZGProofGen(Polynomial_t& P, scalar_t eval, scalar_t index, affine_t* SRS);

bool KZGProofVerif(affine_t commitment, affine_t proof, g2_affine_t g2_srs, scalar_t eval, scalar_t index);

BoardProofs MakeBoardProofs(Polynomial_t& board, scalar_t root, affine_t* srs);