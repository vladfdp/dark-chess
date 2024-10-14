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

constexpr size_t FIXED_CIPHER_LENGTH = 16;
constexpr size_t FIXED_SQUARE_LENGTH = 80; //Decrypted square is 72 bytes but encryption works by 16 bytes block so we round to multiple of 16

struct EncryptedKeys {
    uint64_t attackingSquare;
    std::array<uint8_t, FIXED_CIPHER_LENGTH> message;
};



struct DecryptedSquare {
    uint64_t piece;
    affine_t proof;
};

struct EncryptedSquare {
    std::array<uint8_t, FIXED_SQUARE_LENGTH> square;
};

struct AttackVectors{
    std::vector<EncryptedKeys> av;
};


struct SquareWE{    //TODO: add other pieces
    Fp12 empty;
    Fp12 king;
    g2_affine_t ct;
};

struct BoardWE{
    std::array<SquareWE, TOTAL_SQUARES> board;
};

struct BoardProofs{
    std::array<DecryptedSquare,TOTAL_SQUARES> board;
};

struct SquareBiVOT{
    EncryptedSquare square;
    AttackVectors King; //TODO: add other pieces
};

struct BoardBiVOT{
    std::array<SquareBiVOT,TOTAL_SQUARES> Pieces;
    //std::array<SBiVOTVariant,TOTAL_SQUARES> Positions;
    std::array<g2_affine_t,TOTAL_SQUARES> cts;


    BoardBiVOT(Polynomial_t& board, affine_t enemyBoard, scalar_t root, affine_t* srs, g2_affine_t g2_srs);

};


static std::unique_ptr<scalar_t[]> generate_pows(scalar_t tau, uint32_t size);

std::unique_ptr<affine_t[]> generate_SRS(uint32_t size, scalar_t secret_scalar);

affine_t commit(Polynomial_t& P, affine_t* SRS);

BoardWE MakeBoardWE(affine_t enemyBoard, g2_affine_t g2_srs, scalar_t root);

SquareWE MakeSquareWE(affine_t enemyBoard, g2_affine_t g2_srs, int index);

affine_t KZGProofGen(Polynomial_t& P, scalar_t eval, scalar_t index, affine_t* SRS);

bool KZGProofVerif(affine_t commitment, affine_t proof, g2_affine_t g2_srs, scalar_t eval, scalar_t index);

BoardProofs MakeBoardProofs(Polynomial_t& board, scalar_t root, affine_t* srs);

ChessBoard DecryptBoardBiVOT(BoardBiVOT& enemyBiVOT, Polynomial_t& poly_board, affine_t enemyBoardCom, scalar_t root, affine_t* g1_srs, g2_affine_t g2_srs);