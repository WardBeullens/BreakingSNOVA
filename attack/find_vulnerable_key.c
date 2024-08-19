#include <stdio.h>

#define NO_MT4B
#include "../deriv_params.h"
#include "../snova.h"
#include "../util/util.h"
#include "../api.h"
#include <time.h>

#if rank == 2
    #define PRINT_EVERY 250000
#else
    #define PRINT_EVERY 5
#endif

int main(int argc, char **argv) {

    if(argc != 3){
        printf("Supply start and number of keys!\n");
        return 0;
    }

    size_t start = atoi(argv[1]);
    size_t number_of_keys = atoi(argv[2]);

    printf("start: %ld, number_of_keys: %ld\n", start, number_of_keys);

    snova_init();
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t entropy_input[48];
    for (int i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }

    randombytes_init(entropy_input, NULL, 256);

    for (size_t i = 0; i < start; i++)
    {
        randombytes(pk, seed_length);
    }
    

    int lowest_ranks[sq_rank+1] = {0};
    size_t tries = 0;
    int vulnerable_key_found = 0;
    time_t seconds = time(NULL);
    while (tries < number_of_keys)
    {
        tries++;
        crypto_sign_keypair(pk, sk);

        int lr = lowest_rank(pk);
        lowest_ranks[lr] ++;

        if(lr <= ATTACK_RANK && vulnerable_key_found == 0){
            vulnerable_key_found = 1;
            print_public_key(pk);
            printf("first key found after %ld attempts\n", tries);
        }

        if (tries % PRINT_EVERY == 0){
            printf("tries: %ld, tries/second = %f \n", tries, (tries+0.0)/(time(NULL)+0.001-seconds));
            printf("lowest ranks: \n");
            for (size_t r = 0; r < sq_rank+1; r++)
            {
                printf("%2ld: %d, %f \n", r, lowest_ranks[r], (lowest_ranks[r]+0.0)/tries);
            }
        }
    }

    printf("tries: %ld, tries/second = %f \n", tries, (tries+0.0)/(time(NULL)+0.001-seconds));
    printf("lowest ranks: \n");
    for (size_t r = 0; r < sq_rank+1; r++)
    {
        printf("%2ld: %d, %f \n", r, lowest_ranks[r], (lowest_ranks[r]+0.0)/tries);
    }

    return 0;
}
