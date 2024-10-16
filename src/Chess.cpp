#include "Chess.h"
#include <array>
#include <iostream>
#include "icicle/runtime.h"
#include "icicle/api/bn254.h"
#include "utils.h"
#include "icicle/polynomials/polynomials.h"
#include "icicle/ntt.h"

using namespace bn254;

typedef Polynomial<scalar_t> Polynomial_t;

void ChessBoard::initializeWBoard() {   //TODO: Make a version with random input
    // Initialize all squares to Empty
    for (int i = 0; i < BOARD_SIZE; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {
            board[i][j] = Empty;
        }
    }

    // Set up white pieces
    board[0][0] = board[0][7] = Rook;
    // board[0][1] = board[0][6] = Knight;
    board[0][2] = board[0][5] = Bishop;
    board[0][3] = Queen;
    board[0][4] = King;
    // for (int i = 0; i < BOARD_SIZE; ++i) {
    //     board[1][i] = Pawn;
    // }

}

void ChessBoard::initializeBBoard() {   //just the king for the moment
    // Initialize all squares to Empty
    for (int i = 0; i < BOARD_SIZE; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {
            board[i][j] = Empty;
        }
    }

    // Set up black pieces
    board[7][0] = board[7][7] = Rook;
    // board[7][1] = board[7][6] = Knight;
    board[7][2] = board[7][5] = Bishop;
    board[7][3] = Queen;
    board[7][4] = King;
    // board[0][0] = Knight;
    // for (int i = 0; i < BOARD_SIZE; ++i) {
    //     board[6][i] = Pawn;
    // }
}

void ChessBoard::visualizeBoard(){
        std::cout << "  a b c d e f g h" << std::endl;
        for (int i = BOARD_SIZE - 1; i >= 0; --i) {
            std::cout << i + 1 << " ";
            for (int j = 0; j < BOARD_SIZE; ++j) {
                char piece = getPieceChar(board[i][j]);
                std::cout << piece << " ";
            }
            std::cout << i + 1 << std::endl;
        }
        std::cout << "  a b c d e f g h" << std::endl;
}

char ChessBoard::getPieceChar(int piece) {
        char pieceChar;
        switch (static_cast<PieceType>(piece)) {
            case King: pieceChar = 'K'; break;
            case Queen: pieceChar = 'Q'; break;
            case Rook: pieceChar = 'R'; break;
            case Bishop: pieceChar = 'B'; break;
            case Knight: pieceChar = 'N'; break;
            case Pawn: pieceChar = 'P'; break;
            case Unknown: pieceChar = 'X'; break;
            default: return '.';
        }
        return pieceChar;
}


std::array<scalar_t, TOTAL_SQUARES> ChessBoard::toBoardArray() const {
        std::array<scalar_t, TOTAL_SQUARES> boardArray;
        int index = 0;

        for (const auto& row : board) {
            for (int piece : row) {
                boardArray[index++] = scalar_t::from(piece);
            }
        }

        return boardArray;
}

Polynomial_t ChessBoard::toPoly(){
    //std::array<scalar_t, TOTAL_SQUARES> boardArray;
    auto boardArray = std::make_unique<scalar_t[]>(TOTAL_SQUARES);
    int index = 0;

    for (const auto& row : board) {
        for (int piece : row) {
            boardArray[index++] = scalar_t::from(piece);
        }
    }

    return Polynomial_t::from_rou_evaluations(boardArray.get(), 64);
}

Polynomial_t ChessBoard::toPosPoly(){
    //std::array<scalar_t, TOTAL_SQUARES> boardArray;
    auto boardArray = std::make_unique<scalar_t[]>(TOTAL_SQUARES);
    int index = 0;

    for (const auto& row : board) {
        for (int piece : row) {
            if(piece){
            boardArray[index++] = scalar_t::one();
            }else{
            boardArray[index++] = scalar_t::zero();
            }
        }
    }

    return Polynomial_t::from_rou_evaluations(boardArray.get(), 64);
}

ChessBoard ChessBoard::fromPoly(const Polynomial_t &poly, scalar_t root){  //TODO: use ntt instead
    scalar_t i = scalar_t::one();
    std::array<std::array<int, BOARD_SIZE>, BOARD_SIZE> board;
    for(int j = 0; j < TOTAL_SQUARES; ++j){
        board[j / BOARD_SIZE][j % BOARD_SIZE] = poly(i).export_limbs()[0];
        i = i * root;
    }
    return {board};
}