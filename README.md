# sexyBWT
* `boost` is required
* Intel `tbb` is required
* google test is required for testing

## Hybrid
* [README.md](https://github.com/WilliamHsieh/sexyBWT/tree/hybrid)

## pSAIS
* This is the parallel version of SAIS.
* paper: https://doi.org/10.1007/s11227-018-2395-5

### Test
* All the testing is writen in google test.
```
make test
./test
```

### Run example
* The input file is in fasta format.
```
make psais
./psais < {fasta_file}
```


## radix_sort
* [README.md](https://github.com/WilliamHsieh/sexyBWT/tree/radix_sort)

## reference
* [slide](https://docs.google.com/presentation/d/1_wfaj8DifSW6FZVzTrDeEVucdw74WgY7prA6zdyCS1o/edit#slide=id.ge4f602b3f0_0_213)
* [api](https://hackmd.io/@williamhsieh/SyT9fFN9d)