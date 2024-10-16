#include <array>
#include <iostream>
#include <vector>
#include <random>
#include <cassert>
#include "icicle/runtime.h"
#include "icicle/api/bn254.h"
#include "utils.h"
#include "icicle/polynomials/polynomials.h"
#include "Chess.h"
#include "icicle/ntt.h"
//#include <botan/aes.h>
#include <botan/auto_rng.h>
#include <botan/cipher_mode.h>
#include <botan/hex.h>
#include <botan/hash.h>
#include "BiVOT.h"
#include "utils.h"
#include "TwistedG2.h"
#include "Fr.h"
#include "Fp.h"
#include "Fp2.h"
#include "Fp12.h"
#include "EC.h"
#include "Pairing.h"


static std::unique_ptr<scalar_t[]> generate_pows(scalar_t tau, uint32_t size)
{
  auto vec = std::make_unique<scalar_t[]>(size);
  vec[0] = scalar_t::one();
  for (size_t i = 1; i < size; ++i) {
    vec[i] = vec[i - 1] * tau;
  }
  return std::move(vec);
}

std::unique_ptr<affine_t[]> generate_SRS(uint32_t size, scalar_t secret_scalar)
{
  auto gen = projective_t::generator();
  auto pows_of_tau = generate_pows(secret_scalar, size);
  auto SRS = std::make_unique<affine_t[]>(size);
  for (size_t i = 0; i < size; ++i) {
    SRS[i] = projective_t::to_affine(pows_of_tau[i] * gen);
  }
  return std::move(SRS);
}

affine_t commit(Polynomial_t& Q, affine_t* SRS){
    MSMConfig config = default_msm_config();
    auto [viewQ, sizeQ] = Q.get_coefficients_view();

    projective_t commitment;
    msm(viewQ.get(), SRS, sizeQ, config, &commitment);
    affine_t afComm = projective_t::to_affine(commitment);
    return afComm;
}


//create proof that board(index) = eval
affine_t KZGProofGen(Polynomial_t& P, scalar_t eval, scalar_t index, affine_t* SRS){
    auto vec = std::make_unique<scalar_t[]>(2);
    vec[0] = scalar_t::neg(index);
    vec[1] = scalar_t::one();

    P.sub_monomial_inplace(eval);
    Polynomial_t divident = Polynomial_t::from_coefficients(vec.get(), 2);
    Polynomial_t Q = P / divident;
    P.add_monomial_inplace(eval); //get it back to normal

    return commit(Q, SRS);
}

bool KZGProofVerif(affine_t commitment, affine_t proof, g2_affine_t g2_srs, scalar_t eval, scalar_t index){

    projective_t scaled_generator_g1 = eval * projective_t::generator();
    
    // Calculate the difference between the commitment and the scaled G1 generator
    affine_t commitment_affine_diff = projective_t::to_affine(projective_t::from_affine(commitment) - scaled_generator_g1);
    g1 commitmentDifference = convert_affine_to_g1(commitment_affine_diff);

    // Compute the left-hand side of the pairing equation
    g2_affine_t g2_gen = g2_projective_t::to_affine(g2_projective_t::generator());
    Fp12 lhs = Tate_pairing(commitmentDifference, convert_g2_affine_to_g2(g2_gen));

    // Compute index * G2 generator
    g2_projective_t scaled_generator_g2 = index * g2_projective_t::generator();
    
    // Calculate the difference between the G2 SRS and the scaled G2 generator
    g2_affine_t g2_affine_diff = g2_projective_t::to_affine(g2_projective_t::from_affine(g2_srs) - scaled_generator_g2);
    g2 g2_offset = convert_g2_affine_to_g2(g2_affine_diff);

    // Compute the right-hand side of the pairing equation
    Fp12 rhs = Tate_pairing(convert_affine_to_g1(proof), g2_offset);

    // Return true if both sides are equal, proving the KZG proof
    return Fp12_equal(lhs, rhs);
}


//fills the board with (piece, proof)
BoardProofs MakeBoardProofs(Polynomial_t& board, scalar_t root, affine_t* srs){

    BoardProofs result;
    scalar_t index = scalar_t::one();
    for(int i = 0; i < TOTAL_SQUARES; ++i){
        scalar_t eval = board(index);
        if( eval != scalar_t::one() && eval != scalar_t::zero()) printf("problem with MakeBoardProof\n\n\n");
        affine_t proof = KZGProofGen(board, eval, index, srs);
        result.board[i] = {eval.limbs_storage.limbs[0], proof};
        index = index * root;
    }
    return result;
}

SquareWE MakeSquareWE(affine_t enemyBoard, g2_affine_t g2_srs, scalar_t index, bool debug){//TODO: add other pieces
    SquareWE result;
    scalar_t r = scalar_t::rand_host();
    projective_t proj_enemyBoard = projective_t::from_affine(enemyBoard);
    g2_projective_t proj_g2_srs = g2_projective_t:: from_affine(g2_srs);

    g2_projective_t gen = g2_projective_t::generator();
    g2_affine_t aff_gen = g2_projective_t::to_affine(gen);
    g2 g2_gen = convert_g2_affine_to_g2(aff_gen);
    //empty
    affine_t empty_dif = projective_t::to_affine(r * proj_enemyBoard);
    result.empty = Tate_pairing(convert_affine_to_g1(empty_dif), g2_gen);

    //king
    affine_t king_diff = projective_t::to_affine( r * (proj_enemyBoard - projective_t::generator()) );
    result.king = Tate_pairing(convert_affine_to_g1(king_diff), g2_gen);

    //queen
    affine_t queen_diff = projective_t::to_affine( r * (proj_enemyBoard - (scalar_t::from(2) * projective_t::generator())) );
    result.queen = Tate_pairing(convert_affine_to_g1(queen_diff), g2_gen);


    //encode position
    g2_projective_t alpha = index * gen;
    g2_projective_t ct = r * (proj_g2_srs - alpha);
    result.ct = g2_projective_t::to_affine(ct);

    if(debug){
        print_Fp12(result.king);
        printf("<-this is from MakeSquareWE\n\n");
        std::cout << "\nalso eval should be 1" << "index:" << index << "ct" << ct << "\n";
    }

    return result;
}

BoardWE MakeBoardWE(affine_t enemyBoard, g2_affine_t g2_srs, scalar_t root){
    BoardWE result;
    scalar_t index = scalar_t::one();
    bool debug;
    for(int i = 0; i < TOTAL_SQUARES; i++){
        debug = false;
        if(i == 11){debug = true;}
        result.board[i] = MakeSquareWE(enemyBoard, g2_srs, index, debug);
        index = index*root;
    }
    return result;
}


PosBoardWE MakePosBoardWE(affine_t enemyPosBooard, g2_affine_t g2_srs, scalar_t root){
    PosBoardWE result;
    scalar_t index = scalar_t::one();
    for(int i = 0; i < TOTAL_SQUARES; i++){
         scalar_t r = scalar_t::rand_host();
        projective_t proj_enemyPosBooard = projective_t::from_affine(enemyPosBooard);
        g2_projective_t proj_g2_srs = g2_projective_t:: from_affine(g2_srs);

        g2_projective_t gen = g2_projective_t::generator();
        g2_affine_t aff_gen = g2_projective_t::to_affine(gen);
        g2 g2_gen = convert_g2_affine_to_g2(aff_gen);

        //non-empty
        affine_t piece_diff = projective_t::to_affine( r * (proj_enemyPosBooard - projective_t::generator()) );
        result.keys[i] = Tate_pairing(convert_affine_to_g1(piece_diff), g2_gen);

        //encode position
        g2_projective_t alpha = index * gen;
        g2_projective_t ct = r * (proj_g2_srs - alpha);
        result.cts[i] = g2_projective_t::to_affine(ct);

        index = index*root;      
    }
    return result;
}

std::vector<uint8_t> GetKey(Fp12 key){
    const auto hash = Botan::HashFunction::create_or_throw("SHA-256");
    hash->update(reinterpret_cast<const uint8_t*>(&key), sizeof(Fp12));

    auto hashed_key = hash->final();
    std::vector<uint8_t> final_key(hashed_key.begin(), hashed_key.begin() + 16);
    return final_key;
}

std::vector<uint8_t> GetKey_final(Fp12 key, affine_t proof){
    const auto hash = Botan::HashFunction::create_or_throw("SHA-256");
    hash->update(reinterpret_cast<const uint8_t*>(&key), sizeof(Fp12));
    hash->update(reinterpret_cast<const uint8_t*>(&proof), sizeof(affine_t));

    auto hashed_key = hash->final();
    std::vector<uint8_t> final_key(hashed_key.begin(), hashed_key.begin() + 16);
    return final_key;
}

AttackVectors MakeKingAttackVectors(int row, int col, BoardProofs& pos_board_proofs, BoardWE& board_we, Botan::secure_vector<uint8_t> key){

    AttackVectors KingAttackVectors;

    const int dx[] = {-1, -1, -1, 0, 0, 1, 1, 1};   //directions the king can go in
    const int dy[] = {-1, 0, 1, -1, 1, -1, 0, 1};
    
    for (int i = 0; i < 8; ++i) {
        int newRow = row + dx[i];
        int newCol = col + dy[i];
        
        // Check if the new position is within the board
        if (newRow >= 0 && newRow < BOARD_SIZE && newCol >= 0 && newCol < BOARD_SIZE) {
                        
            
            
            
            
            EncryptedKeys enc_keys;
            uint64_t index = newRow * BOARD_SIZE + newCol;
            enc_keys.attackingSquare = index;
            Fp12 av_key = board_we.board[index].king;
            DecryptedSquare pos_proof = pos_board_proofs.board[index];
            affine_t random_point = projective_t::to_affine(projective_t::rand_host());
            affine_t av_empty_proof = pos_proof.piece ? random_point : pos_proof.proof;

            //hash to get the key
            std::vector<uint8_t> final_key = GetKey_final(av_key, av_empty_proof);

            const auto enc = Botan::Cipher_Mode::create_or_throw("AES-128/CBC/NoPadding", Botan::Cipher_Dir::Encryption);
            enc->set_key(final_key);

            Botan::secure_vector<uint8_t> iv(16, 0);


            Botan::secure_vector<uint8_t> pt(reinterpret_cast<const uint8_t*>(key.data()), 
                                            reinterpret_cast<const uint8_t*>(key.data()) + 16);

            

            enc->start(iv);
            enc->finish(pt);

            if((row == 1)&&(col == 2)&& (index == 11)){
                print_Fp12(av_key);
                printf("\n\nfor attacking square:%llu\n",index);
                std::cout << "we have msg:" << Botan::hex_encode(key) << "with key:" << Botan::hex_encode(final_key) <<"\n";
                std::cout << "that gives cipher " << Botan::hex_encode(pt) << " Cipher size after encryption: " << pt.size() << std::endl;
                std::cout << "IV:" << Botan::hex_encode(iv) << "\n";
                //print_Fp12(final_key)
            }

            std::memcpy(&enc_keys.message, pt.data(), 16);

            KingAttackVectors.av.push_back(enc_keys);
        }
    }
    return KingAttackVectors;
}

AttackVectors MakeQueenAttackVectors(int row, int col, BoardProofs& board_proofs, BoardProofs& pos_board_proofs, BoardWE& board_we, Botan::secure_vector<uint8_t> key){

    AttackVectors QueenAttackVectors;

    const int dx[] = {-1, -1, -1, 0, 0, 1, 1, 1};   //directions the queen can go in
    const int dy[] = {-1, 0, 1, -1, 1, -1, 0, 1};
    
    for (int i = 0; i < 8; ++i) {
        Fp12 empty_inbetween_squares = Fp12_one();
        projective_t inbetween_proofs = projective_t::zero();
        for( int step = 1; step < BOARD_SIZE; ++step){
            int newRow = row + (step*dx[i]);
            int newCol = col + (step*dy[i]);
            // Check if the new position is within the board
            if (newRow >= 0 && newRow < BOARD_SIZE && newCol >= 0 && newCol < BOARD_SIZE) {
                            
                EncryptedKeys enc_keys;
                uint64_t index = newRow * BOARD_SIZE + newCol;
                enc_keys.attackingSquare = index;
                Fp12 av_key = Fp12_mul(board_we.board[index].queen, empty_inbetween_squares);

                DecryptedSquare pos_proof = pos_board_proofs.board[index];
                affine_t random_point = projective_t::to_affine(projective_t::rand_host());
                affine_t av_empty_proof = pos_proof.piece ? random_point : pos_proof.proof; //if square is not empty we just use a random point
                affine_t proof_key = projective_t::to_affine( inbetween_proofs + av_empty_proof);

                //hash to get the key
                std::vector<uint8_t> final_key = GetKey_final(av_key, proof_key);

                const auto enc = Botan::Cipher_Mode::create_or_throw("AES-128/CBC/NoPadding", Botan::Cipher_Dir::Encryption);
                enc->set_key(final_key);

                Botan::secure_vector<uint8_t> iv(16, 0);


                Botan::secure_vector<uint8_t> pt(reinterpret_cast<const uint8_t*>(key.data()), 
                                                reinterpret_cast<const uint8_t*>(key.data()) + 16);

                

                enc->start(iv);
                enc->finish(pt);

                std::memcpy(&enc_keys.message, pt.data(), 16);

                QueenAttackVectors.av.push_back(enc_keys);
                
                empty_inbetween_squares = Fp12_mul(empty_inbetween_squares, board_we.board[index].empty);
                DecryptedSquare board_proof = board_proofs.board[index];
                affine_t random_point2 = projective_t::to_affine(projective_t::rand_host());
                affine_t empty_proof = board_proof.piece ? random_point2 : board_proof.proof; //if square is not empty we just use a random point
                inbetween_proofs = inbetween_proofs + empty_proof;
            }
        }
    }
    return QueenAttackVectors;
}


BoardBiVOT::BoardBiVOT(Polynomial_t& board, Polynomial_t& pos_board, affine_t enemyBoardCom, affine_t enemyPosBoardCom, scalar_t root, affine_t* srs, g2_affine_t g2_srs){
    BoardProofs board_proofs = MakeBoardProofs(board, root, srs);
    BoardWE board_we = MakeBoardWE(enemyBoardCom, g2_srs, root);

    BoardProofs pos_board_proofs = MakeBoardProofs(pos_board, root, srs);
    PosBoardWE pos_board_we = MakePosBoardWE(enemyPosBoardCom, g2_srs, root);

    scalar_t index = scalar_t::one();
    for (int row = 0; row < BOARD_SIZE; ++row) {
        for (int col = 0; col < BOARD_SIZE; ++col) {

            uint64_t i = row * BOARD_SIZE + col;
            const auto enc = Botan::Cipher_Mode::create_or_throw("AES-128/CBC/PKCS7", Botan::Cipher_Dir::Encryption);
            // // Create an initialization vector (IV)
            Botan::secure_vector<uint8_t> iv(16, 0);  // 128-bit IV (16 bytes)
            std::memcpy(iv.data(), &i, sizeof(uint64_t));  // Copy counter into the IV

            //_____ Fill Positions board BiVOT ______________________________

            
            DecryptedSquare pos_square = pos_board_proofs.board[i];
            std::vector<uint8_t> pos_key = GetKey(pos_board_we.keys[i]);
            enc->set_key(pos_key);

            Botan::secure_vector<uint8_t> pt2(reinterpret_cast<const uint8_t*>(&pos_square), 
                                      reinterpret_cast<const uint8_t*>(&pos_square) + sizeof(DecryptedSquare));

            enc->start(iv);
            enc->finish(pt2);

            EncryptedSquare pos_enc_square;
            std::memcpy(&pos_enc_square.square, pt2.data(), FIXED_SQUARE_LENGTH);

            Positions[i] = pos_enc_square;
            PositionsCts[i] = pos_board_we.cts[i];

           



            //_____ Fill Pieces board BiVOT ______________________________

            SquareBiVOT square_bivot;

            scalar_t eval = board(index);
            DecryptedSquare square = board_proofs.board[i];

             // Create a random number generator
            Botan::AutoSeeded_RNG rng;

            // // Generate a random 128-bit (16-byte) key
            Botan::secure_vector<uint8_t> key = rng.random_vec(16);

            
                    
            //Encrypt eval and proof with the key

            enc->set_key(key);

            Botan::secure_vector<uint8_t> pt(reinterpret_cast<const uint8_t*>(&square), 
                                      reinterpret_cast<const uint8_t*>(&square) + sizeof(DecryptedSquare));

            enc->start(iv);
            enc->finish(pt);

             if(row == 1 && col == 2){
                printf("i:%d\n", i);
                std::cout << "we have key:" << Botan::hex_encode(key) <<"\n";
                std::cout << "that gives cipher " << Botan::hex_encode(pt) << " Cipher size after encryption: " << pt.size() << std::endl;
                std::cout << "IV:" << Botan::hex_encode(iv) << "\n";
            }


            EncryptedSquare enc_square;
            std::memcpy(&enc_square.square, pt.data(), FIXED_SQUARE_LENGTH);
            square_bivot.square = enc_square;

            square_bivot.King = MakeKingAttackVectors(row, col, pos_board_proofs, board_we, key);    //TODO: add other pieces
            square_bivot.Queen = MakeQueenAttackVectors(row, col, board_proofs, pos_board_proofs, board_we, key);
            
            Pieces[i] = square_bivot;
            PiecesCts[i] = board_we.board[i].ct;


            index = index * root;
        }
    }
}

std::array<Fp12, TOTAL_SQUARES> DecryptWE(BoardBiVOT& enemyBiVOT, Polynomial_t& poly_board, scalar_t root, affine_t* g1_srs){
    std::array<Fp12, TOTAL_SQUARES> result;
    scalar_t index = scalar_t::one();
    for(int i = 0; i < TOTAL_SQUARES; ++i){
        scalar_t eval = poly_board(index);
        affine_t proof = KZGProofGen(poly_board, eval, index, g1_srs);
        g2 g2_ct = convert_g2_affine_to_g2(enemyBiVOT.PiecesCts[i]);
        result[i] = Tate_pairing(convert_affine_to_g1(proof), g2_ct);

        if(i==11){
            printf("\n\nthis is from decrytpWE");
            print_Fp12(result[i]);
            std::cout << "\nalso eval:" << eval << " index: " << index << " ct: " << enemyBiVOT.PiecesCts[i] << "\n";
        }

        index = index * root;
    }
    return result;
}

std::array<affine_t, TOTAL_SQUARES> DecryptPosProofs(ChessBoard& board, BoardBiVOT& enemyBiVOT, affine_t enemyPosBoardCom, Polynomial_t& pos_poly_board, scalar_t root, affine_t* g1_srs, g2_affine_t g2_srs){
    std::array<affine_t, TOTAL_SQUARES> result;
    scalar_t index = scalar_t::one();
    for(uint64_t i = 0; i < TOTAL_SQUARES; ++i){
        if(pos_poly_board(index) == scalar_t::one()){
        
            
            affine_t proof = KZGProofGen(pos_poly_board, scalar_t::one(), index, g1_srs);
            g2 g2_ct = convert_g2_affine_to_g2(enemyBiVOT.PositionsCts[i]);
            Fp12 key = Tate_pairing(convert_affine_to_g1(proof), g2_ct);
            std::vector<uint8_t> final_key = GetKey(key);

            const auto dec = Botan::Cipher_Mode::create_or_throw("AES-128/CBC/PKCS7", Botan::Cipher_Dir::Decryption);
            dec->set_key(final_key);

            Botan::secure_vector<uint8_t> iv(16, 0);
            std::memcpy(iv.data(), &i, sizeof(uint64_t));

            EncryptedSquare square = enemyBiVOT.Positions[i];

            Botan::secure_vector<uint8_t> cipher(reinterpret_cast<const uint8_t*>(&square), 
                                      reinterpret_cast<const uint8_t*>(&square) + FIXED_SQUARE_LENGTH);


            // std::cout << "we have key:" << Botan::hex_encode(final_key) <<"\n";
            //     std::cout << "with cipher " << Botan::hex_encode(cipher) << " Cipher size before encryption: " << cipher.size() << std::endl;
            //     std::cout << "IV:" << Botan::hex_encode(iv) << "\n";

            dec->start(iv);
            dec->finish(cipher);

            DecryptedSquare decrypted_square;
            std::memcpy(&decrypted_square, cipher.data(), sizeof(DecryptedSquare)); //we don't need to copy the full 80 bytes as the last 8 are padding

            result[i] = decrypted_square.proof;

            if(decrypted_square.piece){
                printf("eaten piece at position i:%d", i);
                board.board[i/BOARD_SIZE][i%BOARD_SIZE] = Unknown;
            }
            
                
            


            assert(KZGProofVerif(enemyPosBoardCom, decrypted_square.proof, g2_srs, scalar_t::from(decrypted_square.piece), index));

            
        }
        index = index * root;
    }
    return result;
}

void DecryptKingBiVOT(int row, int col, ChessBoard& dec_board, BoardBiVOT& enemyBiVOT, affine_t enemyBoardCom, std::array<Fp12, TOTAL_SQUARES> keys, g2_affine_t g2_srs, scalar_t root, affine_t proof_not_taken){

    const int dx[] = {-1, -1, -1, 0, 0, 1, 1, 1};   //directions the king can go in
    const int dy[] = {-1, 0, 1, -1, 1, -1, 0, 1};
    uint64_t index = row * BOARD_SIZE + col;
    
    for (int i = 0; i < 8; ++i) {
        int newRow = row + dx[i];
        int newCol = col + dy[i];
        
        // Check if the new position is within the board
        if (newRow >= 0 && newRow < BOARD_SIZE && newCol >= 0 && newCol < BOARD_SIZE && (dec_board.board[newRow][newCol] == Unknown)) {
            EncryptedKeys enc_keys;
            uint64_t newIndex = newRow * BOARD_SIZE + newCol;

            Botan::secure_vector<uint8_t> cipher(16);

            for(EncryptedKeys ek : enemyBiVOT.Pieces[newIndex].King.av){
                if (ek.attackingSquare == index){
                    std::memcpy(cipher.data(), &ek.message, 16);
                }
            }
            std::cout << "Cipher size before decryption: " << cipher.size() << std::endl;


            Fp12 key = keys[index];
            print_Fp12(key);

            //hash to get the key
            std::vector<uint8_t> hashed_key = GetKey_final(key, proof_not_taken);

            printf("\nfor newIndex :%llu and attacking square:%llu you have :", newIndex, index);
            std::cout << "\ncipher:" << Botan::hex_encode(cipher) << "\n key: " << Botan::hex_encode(hashed_key) << '\n';


            const auto dec = Botan::Cipher_Mode::create_or_throw("AES-128/CBC/NoPadding", Botan::Cipher_Dir::Decryption);
            dec->set_key(hashed_key);

            Botan::secure_vector<uint8_t> iv(16, 0);

            dec->start(iv);
            dec->finish(cipher);

            std::cout << "decrytpion gives the key" << Botan::hex_encode(cipher) << "   using IV:" << Botan::hex_encode(iv) << "\n\n";

            const auto dec2 = Botan::Cipher_Mode::create_or_throw("AES-128/CBC/PKCS7", Botan::Cipher_Dir::Decryption);

            dec2->set_key(cipher);   //use the decrypted message as key to open the (eval,proof) cipher

            std::memcpy(iv.data(), &newIndex, sizeof(uint64_t)); 

            dec2->start(iv);


            Botan::secure_vector<uint8_t> pt(reinterpret_cast<const uint8_t*>(&enemyBiVOT.Pieces[newIndex].square.square), 
                                      reinterpret_cast<const uint8_t*>(&enemyBiVOT.Pieces[newIndex].square.square) + FIXED_SQUARE_LENGTH);

            std::cout << "key is: " << Botan::hex_encode(cipher) << ", IV is: " << Botan::hex_encode(iv) << "\nAnd cipher is: " << Botan::hex_encode(pt) << " of size: " << pt.size() << "\n";
            printf("index:%d, newIndex:%d \n", index, newIndex);
            dec2->finish(pt);
            DecryptedSquare decrypted_square;
            std::memcpy(&decrypted_square, pt.data(), sizeof(DecryptedSquare)); //we don't need to copy the full 80 bytes as the last 8 are padding


            scalar_t rou_index = scalar_t::one();
            for(int k = 0; k < newIndex; ++k){
                rou_index = rou_index*root;
            }

            assert(KZGProofVerif(enemyBoardCom, decrypted_square.proof, g2_srs, scalar_t::from(decrypted_square.piece), rou_index));

            dec_board.board[newRow][newCol] = decrypted_square.piece;
        }
    }

}

void DecryptQueenBiVOT(int row, int col, ChessBoard& dec_board, BoardBiVOT& enemyBiVOT, affine_t enemyBoardCom, std::array<Fp12, TOTAL_SQUARES> keys, g2_affine_t g2_srs, scalar_t root, affine_t proof_not_taken){

    const int dx[] = {-1, -1, -1, 0, 0, 1, 1, 1};   //directions the queen can go in
    const int dy[] = {-1, 0, 1, -1, 1, -1, 0, 1};
    uint64_t index = row * BOARD_SIZE + col;
    
    for (int i = 0; i < 8; ++i) {
        Fp12 empty_inbetween_squares = Fp12_one();
        projective_t aggregate_empty_proofs = projective_t::from_affine(proof_not_taken);
        for( int step = 1; step < BOARD_SIZE; ++step){
            int newRow = row + (step*dx[i]);
            int newCol = col + (step*dy[i]);
            // Check if the new position is within the board
            if (newRow >= 0 && newRow < BOARD_SIZE && newCol >= 0 && newCol < BOARD_SIZE) {
                EncryptedKeys enc_keys;
                uint64_t newIndex = newRow * BOARD_SIZE + newCol;

                Botan::secure_vector<uint8_t> cipher(16);

                for(EncryptedKeys ek : enemyBiVOT.Pieces[newIndex].Queen.av){
                    if (ek.attackingSquare == index){
                        std::memcpy(cipher.data(), &ek.message, 16);
                    }
                }
                Fp12 key = Fp12_mul(keys[index],empty_inbetween_squares);
                

                //hash to get the key
                std::vector<uint8_t> hashed_key = GetKey_final(key, projective_t::to_affine(aggregate_empty_proofs));

        //         printf("for newIndex :%llu and attacking square:%llu you have :", newIndex, index);
        //         std::cout << "\ncipher:" << Botan::hex_encode(cipher) << "\n key: " << Botan::hex_encode(hashed_key) << '\n';


                const auto dec = Botan::Cipher_Mode::create_or_throw("AES-128/CBC/NoPadding", Botan::Cipher_Dir::Decryption);
                dec->set_key(hashed_key);

                Botan::secure_vector<uint8_t> iv(16, 0);

                dec->start(iv);
                dec->finish(cipher);

        //         std::cout << "decrytpion gives the key" << Botan::hex_encode(cipher) << "   using IV:" << Botan::hex_encode(iv) << "\n\n";

                const auto dec2 = Botan::Cipher_Mode::create_or_throw("AES-128/CBC/PKCS7", Botan::Cipher_Dir::Decryption);

                dec2->set_key(cipher);   //use the decrypted message as key to open the (eval,proof) cipher

                std::memcpy(iv.data(), &newIndex, sizeof(uint64_t)); 

                dec2->start(iv);


                Botan::secure_vector<uint8_t> pt(reinterpret_cast<const uint8_t*>(&enemyBiVOT.Pieces[newIndex].square.square), 
                                        reinterpret_cast<const uint8_t*>(&enemyBiVOT.Pieces[newIndex].square.square) + FIXED_SQUARE_LENGTH);

                //std::cout << "key is: " << Botan::hex_encode(cipher) << ", IV is: " << Botan::hex_encode(iv) << "\nAnd crebzervserbzbzebrebzyzrybzipher is: " << Botan::hex_encode(pt) << " of size: " << pt.size() << "\n";
                //printf("index:%d, newIndex:%d, et step:%d", index, newIndex, step);
                dec2->finish(pt);

                DecryptedSquare decrypted_square;
                std::memcpy(&decrypted_square, pt.data(), sizeof(DecryptedSquare)); //we don't need to copy the full 80 bytes as the last 8 are padding


                scalar_t rou_index = scalar_t::one();
                for(int k = 0; k < newIndex; ++k){
                    rou_index = rou_index*root;
                }

                assert(KZGProofVerif(enemyBoardCom, decrypted_square.proof, g2_srs, scalar_t::from(decrypted_square.piece), rou_index));

                if (dec_board.board[newRow][newCol] == Unknown) dec_board.board[newRow][newCol] = decrypted_square.piece;

                if (dec_board.board[newRow][newCol] == Empty){
                    empty_inbetween_squares = Fp12_mul(empty_inbetween_squares, keys[newIndex]);
                    aggregate_empty_proofs = aggregate_empty_proofs + decrypted_square.proof;
                }else{
                    break;
                }
            }
        }
    }

}

ChessBoard DecryptBoardBiVOT(BoardBiVOT& enemyBiVOT, Polynomial_t& poly_board, Polynomial_t& pos_poly_board, affine_t enemyBoardCom, affine_t enemyPosBoardCom, scalar_t root, affine_t* g1_srs, g2_affine_t g2_srs){
    ChessBoard dec_board;
    ChessBoard board = ChessBoard::fromPoly(poly_board, root);
    std::array<Fp12, TOTAL_SQUARES> keys = DecryptWE(enemyBiVOT, poly_board, root, g1_srs);

    printf("\n ________________________________________________________\n\n");

    for (int i = 0; i < BOARD_SIZE; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {
            dec_board.board[i][j] = Unknown;
            if(board.board[i][j]) dec_board.board[i][j] = board.board[i][j];
        }
    }
    dec_board.visualizeBoard();
    std::array<affine_t, TOTAL_SQUARES> pos_proofs = DecryptPosProofs(dec_board, enemyBiVOT, enemyPosBoardCom, pos_poly_board, root, g1_srs, g2_srs);
    dec_board.visualizeBoard();
    //TODO: get rid of eaten pieces

    for (int i = 0; i < BOARD_SIZE; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {

            int index = i * BOARD_SIZE + j;
            affine_t proof_not_taken = pos_proofs[index];
            
            switch(dec_board.board[i][j]){         //TODO: add other pieces
                case PieceType::King:
                DecryptKingBiVOT(i, j, dec_board, enemyBiVOT, enemyBoardCom, keys, g2_srs, root, proof_not_taken);
                break;
                case PieceType::Queen:
                DecryptQueenBiVOT(i, j, dec_board, enemyBiVOT, enemyBoardCom, keys, g2_srs, root, proof_not_taken);
                break;
                default:
                break;
            }
        
        }
    }
    return dec_board;

}