# sexyBWT
* `boost` is required
* Intel `tbb` is required
* google test is required for testing

## Hybrid
* This is the hybrid version of the algorithm (pSAIS + radix_sort).

### Test
* All the testing is writen in google test.
```
make test
./test
```

### Run example
* The input file is in fasta format.
```
make hybrid
./hybrid < {fasta_file}
```

## pSAIS
* [README.md](https://github.com/WilliamHsieh/sexyBWT/blob/psais/README.md)

## radix_sort
* [README.md](https://github.com/WilliamHsieh/sexyBWT/blob/radix_sort/README.md)

## reference
* [slide](https://docs.google.com/presentation/d/1_wfaj8DifSW6FZVzTrDeEVucdw74WgY7prA6zdyCS1o/edit#slide=id.ge4f602b3f0_0_213)
* [api](https://hackmd.io/@williamhsieh/SyT9fFN9d)