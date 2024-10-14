#include <array>
#include <iostream>
#include <vector>
#include <random>
#include "icicle/runtime.h"
#include "icicle/api/bn254.h"
#include "utils.h"
#include "icicle/polynomials/polynomials.h"
#include "Chess.h"
#include "BiVOT.h"
#include "icicle/ntt.h"
#include <botan/auto_rng.h>
#include <botan/rng.h>
#include <botan/cipher_mode.h>
#include <botan/hex.h>
#include <botan/hash.h>
extern "C" {
#include "Fp12.h"
#include "Fp.h"
#include "EC.h"
#include "TwistedG2.h"
#include "Pairing.h"
}

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


    scalar_t secret_scalar = scalar_t::rand_host();
    auto g1_srs = generate_SRS(100, secret_scalar);
    auto g2_srs = g2_projective_t::to_affine(secret_scalar*g2_projective_t::generator());



    ChessBoard whiteboard;
    whiteboard.initializeWBoard();
    whiteboard.visualizeBoard();

    Polynomial_t wb_poly = whiteboard.toPoly();
    affine_t wb_com = commit(wb_poly, g1_srs.get());

    ChessBoard blackboard;
    blackboard.initializeBBoard();
    blackboard.visualizeBoard();

    Polynomial_t bb_poly = blackboard.toPoly();
    affine_t bb_com = commit(bb_poly, g1_srs.get());



    BoardBiVOT wb_bivot(wb_poly, bb_com, basic_root, g1_srs.get(), g2_srs);

    ChessBoard blackboard2 = DecryptBoardBiVOT(wb_bivot, bb_poly, wb_com, basic_root, g1_srs.get(), g2_srs);

    blackboard2.visualizeBoard();













    // BoardWE wb_board_we = MakeBoardWE(com, g2_srs, basic_root);

    // scalar_t index = scalar_t::one();
    // scalar_t eval = wb_poly(index);

    // affine_t proof = KZGProofGen(wb_poly, eval, index, g1_srs.get());
    // std::cout << eval << std::endl;

    // g1 g1_proof = convert_affine_to_g1(proof);

    
    // scalar_t r = scalar_t::rand_host();
    // projective_t proj_com = projective_t::from_affine(com);
    // g2_projective_t proj_g2_srs = g2_projective_t:: from_affine(g2_srs);

    // g2_projective_t gen = g2_projective_t::generator();
    // g2_affine_t aff_gen = g2_projective_t::to_affine(gen);
    // g2 g2_gen = convert_g2_affine_to_g2(aff_gen);
    // //empty
    // affine_t empty_dif = projective_t::to_affine( r * (proj_com - projective_t::generator()));
    // Fp12 key = Tate_pairing(convert_affine_to_g1(empty_dif), g2_gen);

    // g2_projective_t alpha = index * gen;
    // g2_projective_t ct = r * (proj_g2_srs - alpha);
    // g2_affine_t aff_ct = g2_projective_t::to_affine(ct);
    // g2 g2_ct = convert_g2_affine_to_g2(aff_ct);

    // Fp12 key2 = Tate_pairing(g1_proof, g2_ct);





    // g2 ct = convert_g2_affine_to_g2(wb_board_we.board[0].ct);

    // Fp12 key = Tate_pairing(g1_proof, ct);

    // Fp12 key2 = wb_board_we.board[0].king;


    // print_Fp12(key);
    // printf("\n\n");
    // print_Fp12(key2);
   






















    // BoardProofs board_proof = MakeBoardProofs(wb_poly, basic_root, g1_srs.get());
    // scalar_t index = scalar_t::one();

    // for(int i = 0; i < 64; ++i){
    //     scalar_t piece = scalar_t::from(board_proof.board[i].piece);
    //     affine_t proof = board_proof.board[i].proof;
    //     bool ok = KZGProofVerif(com, proof, g2_srs, piece, index);
    //     printf("piece is :%d, and verif gives:%d\n", board_proof.board[i].piece, ok);

    //     index = index * basic_root;
    // }




































    // auto coeff = std::make_unique<scalar_t[]>(64);
    // for(int i = 0; i < 64; ++i){
    //     coeff[i] = scalar_t::zero();
    // }
    // coeff[2] = scalar_t::one();
    // scalar_t index = basic_root*basic_root;
    // Polynomial_t poly = Polynomial_t::from_rou_evaluations(coeff.get(), 64);

    // affine_t com = commit(poly, g1_srs.get());

    // affine_t proof = KZGProofGen(poly, scalar_t::one(), index, g1_srs.get());

    // bool ok = KZGProofVerif(com, proof, g2_srs, scalar_t::one(), index);

    // printf("ok is :%d\n", ok);
    // printf("should be:%d\n", poly(index) == scalar_t::one());


























    // ChessBoard whiteboard;
    // ChessBoard blackboard;
    // whiteboard.initializeWBoard();
    // blackboard.initializeBBoard();
    // whiteboard.visualizeBoard();
    // blackboard.visualizeBoard();

    // Polynomial_t wb_poly = whiteboard.toPoly();
    // Polynomial_t bb_poly = blackboard.toPoly();

    // ChessBoard bbfrompoly = ChessBoard::fromPoly(bb_poly, basic_root);
    // bbfrompoly.visualizeBoard();




    // scalar_t tf = scalar_t::from(34);
    // g2_projective_t P = g2_projective_t::rand_host();
    // g2_affine_t P2 = g2_projective_t::to_affine(tf*P);
    // g2_affine_t P3 = g2_projective_t::to_affine(P);
    // std::cout << P2 << std::endl;

    // TwistedG2 Q = TwistedG2_mul_by_int(convert_g2_affine_to_twistedg2(P3) ,{0,0,0,34});

    // std::cout << convert_twistedg2_to_g2_affine(Q) << std::endl;













































    // printf("sould be equal %d and %d", sizeof(EncryptedSquare), sizeof(DecryptedSquare));

    // printf("\n\n\n\n___________________n___________________n___________________n___________________\n");

//     scalar_t key = scalar_t::rand_host();
//     //scalar_t message = scalar_t::rand_host();
//     Fp12 message = Fp12_one();
//     //std::cout << "message:" << message << std::endl;
//     Botan::AutoSeeded_RNG rng;

//     const auto hash = Botan::HashFunction::create_or_throw("SHA-256");
//     hash->update(reinterpret_cast<const uint8_t*>(&key), sizeof(scalar_t));

//     auto hashed_key = hash->final();
//     std::vector<uint8_t> final_key(hashed_key.begin(), hashed_key.begin() + 16);


//     const auto enc = Botan::Cipher_Mode::create_or_throw("AES-128/CBC/PKCS7", Botan::Cipher_Dir::Encryption);
//     enc->set_key(final_key);

//     Botan::secure_vector<uint8_t> iv = rng.random_vec(enc->default_nonce_length());


//     Botan::secure_vector<uint8_t> pt(reinterpret_cast<const uint8_t*>(&message), 
//                                       reinterpret_cast<const uint8_t*>(&message) + sizeof(Fp12));

//     enc->start(iv);
//     enc->finish(pt);

//    std::cout << enc->name() << " with iv " << Botan::hex_encode(iv) << " " << Botan::hex_encode(pt) << '\n';



//     const auto dec = Botan::Cipher_Mode::create_or_throw("AES-128/CBC/PKCS7", Botan::Cipher_Dir::Decryption);
//     dec->set_key(final_key);

//     // Decrypt
//     dec->start(iv);
//     Botan::secure_vector<uint8_t> decrypted_pt = pt;
//     dec->finish(decrypted_pt);

//     // Cast back to scalar_t
//     Fp12 decrypted_message;
//     // assert(decrypted_pt.size() >= sizeof(scalar_t));
//     std::memcpy(&decrypted_message, decrypted_pt.data(), sizeof(Fp12));

//     std::cout << "Decrypted message: " ; 
//     print_Fp12(decrypted_message);
//     std::cout << std::endl;
































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