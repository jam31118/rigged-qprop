- [x] `MPI_init()`
- [x] get rank, size

### Retrieve wavefunction data from `current-wf.bin` 

- [x] specify which range of wf to read - depending on rank
- [x] open file collectively : `current-wf.bin`
- [x] read wf into memory - collectively through `MPI_File_read_at()`

### Evaulate propagators

- [x] evaluate propagator and store it to array for stacked propagators

### Propagate

- [x] open file collectively : `tsurffpsi.raw` and `tsurff-dpsidr.raw` in append write mode
- [x] iteration over wf that has been read (iter over `lm_index`):
  - [x] get propagator from stacked propagators
  - [ ] iteration over `time_index`
    - [x] store `tsurff`-related quantities: `phi_lm@R`, `dphi-drho@R` to an array, intended to be a cache-like entity
    - [x] wf propagation
    - [ ] (for every period1) save `current-wf.bin` as a backup or mid-product
    - [ ] (for every period2 - depending on the size of cache-like) save `tsurff-related quantities` into file (should be collective, but may be not)
- [ ] save remaining tsurff-quantities cache-array
- [x] save `current-wf.bin`
- [x] close file: `tsurffpsi.raw` and `tsurff-dpsidr.raw` 
- [x] close file: `current-wf.bin`

### Finalize MPI program

- [x] `MPI_Filnalize`