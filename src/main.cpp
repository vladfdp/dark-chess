#include <array>
#include <iostream>
#include <vector>
#include <random>
#include "icicle/runtime.h"
#include "icicle/api/bn254.h"
#include "utils.h"

using namespace bn254;




int main(){

     // Load installed backends
  icicle_load_backend_from_env_or_default();

  // trying to choose CUDA if available, or fallback to CPU otherwise (default device)
  const bool is_cuda_device_available = (eIcicleError::SUCCESS == icicle_is_device_available("CUDA"));
  if (is_cuda_device_available) {
    Device device = {"CUDA", 0};             // GPU-0
    ICICLE_CHECK(icicle_set_device(device)); // ICICLE_CHECK asserts that the api call returns eIcicleError::SUCCESS
  } // else we stay on CPU backend



    auto P = projective_t::rand_host();

  auto Q = projective_t::rand_host();

   const auto& modulus = fp_config::modulus;

    // Convert the modulus to an array for easier access
    // Print the modulus
    std::cout << "BN254 scalar field modulus:" << std::endl;
    for (int i = 7; i >= 0; --i) {
        printf("%x ", modulus.limbs[i]);
    }
    std::cout << std::endl;


    auto shifted = shift_left(modulus);

    std::cout << "BN254 scalar field modulus:" << std::endl;
    for (int i = 7; i >= 0; --i) {
        printf("%x ", shifted.limbs[i]);
    }
    std::cout << std::endl;


    return 0;
}