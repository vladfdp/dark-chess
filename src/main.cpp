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
    //auto ntt_config = ntt::default_ntt_config<scalar_t>();
    const scalar_t basic_root = scalar_t::omega(MAX_NTT_LOG_SIZE);
    ntt_init_domain(basic_root, default_ntt_init_domain_config());

    auto coeff = std::make_unique<scalar_t[]>(64);
    ChessBoard board;
    board.initializeWBoard();
    board.visualizeBoard();

    auto arr = board.toBoardArray();
    for(int i = 0; i < 64; ++i){
        coeff[i] = arr[i];
    }

    Polynomial_t PolyBoard = Polynomial_t::from_rou_evaluations(coeff.get(), 64);
    scalar_t root = basic_root;

    for(int i = 0; i < 64; ++i){
        std::cout << PolyBoard(root) << std::endl;
        root = root * basic_root;
    }




























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