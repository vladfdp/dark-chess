#include <array>
#include "icicle/runtime.h"
#include "icicle/api/bn254.h"
#include "utils.h"
#include "icicle/polynomials/polynomials.h"
#include "icicle/ntt.h"

using namespace bn254;
#ifndef CHESS_H
#define CHESS_H

const int BOARD_SIZE = 4;
const int TOTAL_SQUARES = BOARD_SIZE*BOARD_SIZE;

typedef Polynomial<scalar_t> Polynomial_t;

enum PieceType {
    Empty = 0,
    King = 1,
    Queen = 2,
    Rook = 3,
    Bishop = 4,
    Knight = 5,
    Pawn = 6,
    Unknown = 7
};

struct ChessBoard {
    std::array<std::array<int, BOARD_SIZE>, BOARD_SIZE> board;

    void initializeWBoard();
    void initializeBBoard();
    void visualizeBoard();
    static char getPieceChar(int piece);
    std::array<scalar_t, TOTAL_SQUARES> toBoardArray() const;
    Polynomial_t toPoly();
    static ChessBoard fromPoly(const Polynomial_t &poly, scalar_t root);

};



#endif