CC = g++
CCFLAG = -Iinclude -std=c++20 -Wall -Wextra -Wshadow -fopenmp -pthread -Ofast -ltbb -g3

.PHONY: test psais

test: tests/main.cpp
	$(CC) $(CCFLAG) $^ -lgtest -o $@

psais: example/psais.cpp
	$(CC) $(CCFLAG) $^ -o $@

hybrid: example/hybrid.cpp
	$(CC) $(CCFLAG) $^ -o $@

radix_sort: example/radix_sort.cpp
	$(CC) $(CCFLAG) $^ -o $@

clean:
	rm -f a.out test psais hybrid radix_sort
