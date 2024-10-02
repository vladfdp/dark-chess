#include "Chess.h"


void ChessBoard::initializeWBoard() {
    // Initialize all squares to Empty
    for (int i = 0; i < BOARD_SIZE; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {
            board[i][j] = Empty;
        }
    }

    // Set up white pieces
    board[0][0] = board[0][7] = Rook;
    board[0][1] = board[0][6] = Knight;
    board[0][2] = board[0][5] = Bishop;
    board[0][3] = Queen;
    board[0][4] = King;
    for (int i = 0; i < BOARD_SIZE; ++i) {
        board[1][i] = Pawn;
    }

}

void ChessBoard::initializeBBoard() {
    // Initialize all squares to Empty
    for (int i = 0; i < BOARD_SIZE; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {
            board[i][j] = Empty;
        }
    }

    // Set up black pieces
    board[7][0] = board[7][7] = Rook;
    board[7][1] = board[7][6] = Knight;
    board[7][2] = board[7][5] = Bishop;
    board[7][3] = Queen;
    board[7][4] = King;
    for (int i = 0; i < BOARD_SIZE; ++i) {
        board[6][i] = Pawn;
    }
}

void ChessBoard::visualizeBoard() const {
        std::cout << "  a b c d e f g h" << std::endl;
        for (int i = BOARD_SIZE - 1; i >= 0; --i) {
            std::cout << i + 1 << " ";
            for (int j = 0; j < BOARD_SIZE; ++j) {
                char piece = getPieceChar(board[i][j], i >= 4);
                std::cout << piece << " ";
            }
            std::cout << i + 1 << std::endl;
        }
        std::cout << "  a b c d e f g h" << std::endl;
}

char ChessBoard::getPieceChar(field piece, bool isBlack) const {
        char pieceChar;
        switch (static_cast<PieceType>(piece)) {
            case King: pieceChar = 'K'; break;
            case Queen: pieceChar = 'Q'; break;
            case Rook: pieceChar = 'R'; break;
            case Bishop: pieceChar = 'B'; break;
            case Knight: pieceChar = 'N'; break;
            case Pawn: pieceChar = 'P'; break;
            default: return '.';
        }
        return isBlack ? std::tolower(pieceChar) : pieceChar;
}
