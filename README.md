# Population Genetics
## Instruction
- Load blanca: ``` module load slurm/blanca```
- simulation code/network model code and submit to slurm: ``` sbatch script/[sim/ML]```


#
## Files
- **vigt**: 22 chr files using for simulation.
- **new**: Msprime program with accepting 64-bit random seed.
- **orginal**: Orginal msprime source code.
- **ouput**: Some outputs by running network model.
    - The way naming the files:
        - fix: fix f values value in model parameters
        - random: using msp seed generator generate parameters
        - preseed: pre-generate all the seeds
    ```
    [number of simulation/model]_[orginal msp or msp with new rng]_[parameter info or seed info].out
    ```
- **script**: Blanca script for running simulation and network model
    - [More Blanca script info](https://curc.readthedocs.io/en/latest/access/blanca.html)
- **data**: All the data files(matrix and statistics) that simulation generated.
    - **mat**: 20000 matrix data each model by new msprime and random model parameters
    - **mat1**: 4000 matrix data each model by old msprime and random model parameters, but fixed f values.
    - **mat2**: 4000 matrix data each model by old msprime and random model parameters.
    - **mat1**: 4000 matrix data each model by old msprime and random model parameters, but preseed.

    Data naming format follows:
```
symmetry_matrix_[model number]_[simulation number]_[node number]

symmetry_stat_[model number]_[simulation number]_[node number]
```
- **model.txt**: Orginal model parameter generator
- **ndmix9**: The simulation code
    - ndmix9_oldmsp and ndmix9_oldrand both generates 64x64 matrix
    - ndmix9_newmsp generates 101x101 matrix
- **network**: Machine Learning code
    - 64x64/101x101 is the size of matrix data
- **test**: both file has the main function of running simulation in parallel
#
## Random seeds via SHA256
First, I create a unique string for each simulation and parameter generating. The formating follows: 
```
[chr file number]_[simluation number]_[slurm job id]_[node id]_[model number]
```
Then, put into hash function to generate 64-bit random seed:
```
int.from_bytes(hashlib.sha256(pass_str.encode('utf-8')).digest()[:8], 'little')
```
#
## New Msprime VS old Msprime
- Disable 32-bit seed checking in **msprime/msprime/_msprimemodule.c**:
```
    // if (seed == 0 || seed >= (1ULL<<32)) {
    //     PyErr_Format(PyExc_ValueError,
    //         "%llu",seed);
    //     goto out;
    // }
```
- Adding init_by_array function to state vector. This is using 64-bit seed as the speical key to mix up with state vector. Mersenne Twister has 19937 bits of state space that it uses to iterate through the sequence of values it produces. If you initialize it with a 32 bit integer, you are restricting it to just 232 out of the 2^19937 possible starting points, and there are a massive number of sample trajectories that you will never see. The init_by_array() function allows you to specify more bits for the initial state, giving the potential to achieve any of the sampling trajectories which MT is capable of generating.
```
init_by_array(&self->seed, 1, self->rng->state);

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
```

