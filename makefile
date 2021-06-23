CC = g++
CCFLAG = -Iinclude -std=c++20 -Wall -Wextra -Wshadow -fopenmp -pthread -Ofast

test-psais:
	$(CC) $(CCFLAG) tests/psais/psais.cpp

test-ThreadPool:
	$(CC) $(CCFLAG) tests/psais/utility/thread_pool.cpp

test-parallel:
	$(CC) $(CCFLAG) tests/psais/utility/parallel.cpp

clean:
	rm -rf a.out
