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


constexpr size_t FIXED_CIPHER_LENGTH = 32;
constexpr size_t FIXED_SQUARE_LENGTH = 32 + (32 * 12);

struct EncryptedMessage {
    int attackingSquare;
    std::array<uint8_t, FIXED_CIPHER_LENGTH> message;
};

struct EncryptedSquare {
    std::array<uint8_t, FIXED_SQUARE_LENGTH> square;
};

// struct DecryptedSquare {

// }

struct AttackVectors{
    std::vector<EncryptedMessage> av;
};

struct SquareBiVOT{
    EncryptedMessage square;
    AttackVectors King; //TODO: add other pieces
};

struct BoardBiVOT{
    std::array<SquareBiVOT,64> Pieces;
    //std::array<SBiVOTVariant,64> Positions;


    BoardBiVOT(Polynomial_t board, scalar_t root, projective_t enemyBoard);

};