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
#include <botan/auto_rng.h>
#include <botan/rng.h>
#include <botan/cipher_mode.h>
#include <botan/hex.h>
#include <botan/hash.h>
#include "Fp12.h"
using namespace bn254;

typedef Polynomial<scalar_t> Polynomial_t;


int main(){

     // Load installed backends
  icicle_load_backend_from_env_or_default();

  // trying to choose CUDA if available, or fallback to CPU otherwise (default device)
  const bool is_cuda_device_available = (eIcicleError::SUCCESS == icicle_is_device_available("CUDA"));
  if (is_cuda_device_available) {
    Device device = {"CUDA", 0};             // GPU-0
    ICICLE_CHECK(icicle_set_device(device)); // ICICLE_CHECK asserts that the api call returns eIcicleError::SUCCESS
  } // else we stay on CPU backend

    const int MAX_NTT_LOG_SIZE = 6;
    const scalar_t basic_root = scalar_t::omega(MAX_NTT_LOG_SIZE);
    ntt_init_domain(basic_root, default_ntt_init_domain_config());




    ChessBoard whiteboard;
    ChessBoard blackboard;
    whiteboard.initializeWBoard();
    blackboard.initializeBBoard();
    whiteboard.visualizeBoard();
    blackboard.visualizeBoard();

    Polynomial_t wb_poly = whiteboard.toPoly();
    Polynomial_t bb_poly = blackboard.toPoly();

    ChessBoard bbfrompoly = ChessBoard::fromPoly(bb_poly, basic_root);
    bbfrompoly.visualizeBoard();


    printf("\n\n\n\n___________________n___________________n___________________n___________________\n");

    scalar_t key = scalar_t::rand_host();
    //scalar_t message = scalar_t::rand_host();
    Fp12 message = Fp12_one();
    //std::cout << "message:" << message << std::endl;
    Botan::AutoSeeded_RNG rng;

    const auto hash = Botan::HashFunction::create_or_throw("SHA-256");
    hash->update(reinterpret_cast<const uint8_t*>(&key), sizeof(scalar_t));

    auto hashed_key = hash->final();
    std::vector<uint8_t> final_key(hashed_key.begin(), hashed_key.begin() + 16);


    const auto enc = Botan::Cipher_Mode::create_or_throw("AES-128/CBC/PKCS7", Botan::Cipher_Dir::Encryption);
    enc->set_key(final_key);

    Botan::secure_vector<uint8_t> iv = rng.random_vec(enc->default_nonce_length());


    Botan::secure_vector<uint8_t> pt(reinterpret_cast<const uint8_t*>(&message), 
                                      reinterpret_cast<const uint8_t*>(&message) + sizeof(scalar_t));

    enc->start(iv);
    enc->finish(pt);

   std::cout << enc->name() << " with iv " << Botan::hex_encode(iv) << " " << Botan::hex_encode(pt) << '\n';



    const auto dec = Botan::Cipher_Mode::create_or_throw("AES-128/CBC/PKCS7", Botan::Cipher_Dir::Decryption);
    dec->set_key(final_key);

    // Decrypt
    dec->start(iv);
    Botan::secure_vector<uint8_t> decrypted_pt = pt;
    dec->finish(decrypted_pt);

    // Cast back to scalar_t
    Fp12 decrypted_message;
    // assert(decrypted_pt.size() >= sizeof(scalar_t));
    std::memcpy(&decrypted_message, decrypted_pt.data(), sizeof(Fp12));

    std::cout << "Decrypted message: " ; 
    print_Fp12(decrypted_message);
    std::cout << std::endl;
































    // auto coeff = std::make_unique<scalar_t[]>(64);
    // ChessBoard board;
    // board.initializeWBoard();
    // board.visualizeBoard();

    // auto arr = board.toBoardArray();
    // for(int i = 0; i < 64; ++i){
    //     coeff[i] = arr[i];
    // }

    // Polynomial_t PolyBoard = board.toPoly();
    // scalar_t root = basic_root;

    // for(int i = 0; i < 64; ++i){
    //     std::cout << PolyBoard(root) << std::endl;
    //     root = root * basic_root;
    // }




























//     auto P = projective_t::rand_host();

//   auto Q = projective_t::rand_host();

//    const auto& modulus = fp_config::modulus;

//     // Convert the modulus to an array for easier access
//     // Print the modulus
//     std::cout << "BN254 scalar field modulus:" << std::endl;
//     for (int i = 7; i >= 0; --i) {
//         printf("%x ", modulus.limbs[i]);
//     }
//     std::cout << std::endl;


//     auto shifted = shift_left(modulus);

//     std::cout << "BN254 scalar field modulus:" << std::endl;
//     for (int i = 7; i >= 0; --i) {
//         printf("%x ", shifted.limbs[i]);
//     }
//     std::cout << std::endl;


    return 0;
}