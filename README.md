# sexyBWT
* `boost` is required
* Intel `tbb` is required
* google test is required for testing
* compiler with `C++20` support

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

### Result
* runtime of two dataset (drosophila and hs37d5)
* the result are in seconds

#### drosophila
![](https://i.imgur.com/buX5PaH.png)

#### hs37d5
![](https://i.imgur.com/7zNRN7I.png)

## pSAIS
* [README.md](https://github.com/WilliamHsieh/sexyBWT/tree/psais)

## radix_sort
* [README.md](https://github.com/WilliamHsieh/sexyBWT/tree/radix_sort)

## reference
* [slide](https://docs.google.com/presentation/d/1_wfaj8DifSW6FZVzTrDeEVucdw74WgY7prA6zdyCS1o/edit#slide=id.ge4f602b3f0_0_213)
* [api](https://hackmd.io/@williamhsieh/SyT9fFN9d)