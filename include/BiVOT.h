#include <array>
#include <variant>
#include <vector>
#include <cstdint>


constexpr size_t FIXED_CIPHER_LENGTH = 64;
constexpr size_t FIXED_SQUARE_LENGTH = 64 + (64 * 12);

struct EncryptedMessage {
    std::array<uint8_t, FIXED_CIPHER_LENGTH> message;
};

struct EncryptedSquare {
    std::array<uint8_t, FIXED_SQUARE_LENGTH> square;
};

template<size_t NumAttack>
struct AttackVectors{
    std::array<EncryptedMessage, NumAttack> av;
};

template<size_t NumKA>
struct SquareBiVOT{
    EncryptedMessage square;
    AttackVectors<NumKA> King; //TODO: add other pieces
};

using SBiVOTVariant = std::variant<SquareBiVOT<3>, SquareBiVOT<5>, SquareBiVOT<8>>;

struct BoardBiVOT{
    std::array<SBiVOTVariant,64> Pieces;
    std::array<SBiVOTVariant,64> Positions;
};