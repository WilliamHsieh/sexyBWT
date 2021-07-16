#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1

#include <gtest/gtest.h>
#include "psais/psais.cpp"
#include "psais/utility/parallel.cpp"
#include "psais/utility/thread_pool.cpp"
#include "radix_sort/radix_sort.cpp"

int main(int argc, char **argv) {
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
