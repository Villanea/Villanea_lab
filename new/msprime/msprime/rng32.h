#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <Python.h>
#define NN 624
#define MM 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

typedef struct
  {
    unsigned long mt[NN];
    int mti;
  }
mt_state_t;


void init_by_array(unsigned long long init_key[], int key_length, void *vstate)
{
    int i, j, k;
    //init_genrand(19650218UL, state);
    i=1; j=0;
    mt_state_t *state = (mt_state_t *) vstate;
    k = (NN>key_length ? NN : key_length);
    for (; k; k--) {
        state->mt[i] = (state->mt[i] ^ ((state->mt[i-1] ^ (state->mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        state->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=NN) { state->mt[0] = state->mt[NN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=NN-1; k; k--) {
        state->mt[i] = (state->mt[i] ^ ((state->mt[i-1] ^ (state->mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        state->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=NN) { state->mt[0] = state->mt[NN-1]; i=1; }
    }

    state->mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}
