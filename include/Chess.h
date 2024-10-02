#ifndef CHESS_H
#define CHESS_H

const int BOARD_SIZE = 8;

enum PieceType {
    Empty = 0,
    King = 1,
    Queen = 2,
    Rook = 3,
    Bishop = 4,
    Knight = 5,
    Pawn = 6
};

struct ChessBoard {
    std::array<std::array<int, BOARD_SIZE>, BOARD_SIZE> board;

    void initializeWBoard();
    void initializeBBoard();
    void visualizeBoard();
    char getPieceChar(field piece, bool isBlack);
};