CC = g++
CCFLAG = -Iinclude -std=c++20 -Wall -Wextra -Wshadow -fopenmp -pthread -Ofast

.PHONY: test

test:
	$(CC) $(CCFLAG) tests/main.cpp -lgtest -o test

example-fasta:
	$(CC) $(CCFLAG) example/fasta.cpp -o psais

clean:
	rm -f a.out test
