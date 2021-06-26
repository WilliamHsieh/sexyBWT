CC = g++
CCFLAG = -Iinclude -std=c++20 -Wall -Wextra -Wshadow -fopenmp -pthread -Ofast

.PHONY: test

test: tests/main.cpp
	$(CC) $(CCFLAG) $^ -lgtest -o $@

psais: example/fasta.cpp
	$(CC) $(CCFLAG) $^ -o $@

clean:
	rm -f a.out test psais
