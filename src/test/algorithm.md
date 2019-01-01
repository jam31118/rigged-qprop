- [ ] `MPI_init()`
- [ ] get rank, size

### Retrieve wavefunction data from `current-wf.bin` 

- [ ] specify which range of wf to read - depending on rank
- [ ] open file collectively : `current-wf.bin`
- [ ] read wf into memory - collectively through `MPI_File_read_at()`
- [ ] close file collectively

### Propagate

- [ ] open (again) file collectively : `current-wf.bin`
- [ ] open file collectively : `tsurffpsi.raw` and `tsurff-dpsidr.raw` in append write mode
- [ ] iteration over wf that has been read (iter over `lm_index`):
  - [ ] construct propagator
  - [ ] iteration over `time_index`
    - [ ] wf propagation
    - [ ] store `tsurff`-related quantities: `phi_lm@R`, `dphi-drho@R` to an array, intended to be a cache-like entity
    - [ ] (for every period1) save `current-wf.bin` as a backup or mid-product
    - [ ] (for every period2 - depending on the size of cache-like) save `tsurff-related quantities` into file (should be collective, but may be not)
- [ ] save remaining tsurff-quantities cache-array
- [ ] save `current-wf.bin`
- [ ] close file: `tsurffpsi.raw` and `tsurff-dpsidr.raw` 
- [ ] close file: `current-wf.bin`

### Finalize MPI program

- [ ] `MPI_Filnalize`