#include "psais/utility/thread_pool.hpp"

struct Uncopyable {
	Uncopyable() = default;
	Uncopyable(const Uncopyable&) = delete;

	int operator()(int n) {
		if (n <= 1) return n;
		return operator()(n - 1) + operator()(n - 2);
	}
};

TEST(utility, thread_pool) {
	auto pool = psais::utility::ThreadPool(5);

	// function / functor
	Uncopyable fib;
	auto res1 = pool.enqueue(std::ref(fib), 40);

	// lambda
	auto lam = [](int x, int y) { return x + y; };
	auto res2 = pool.enqueue(lam, 3, 5);

	ASSERT_EQ(res1.get(), Uncopyable{}(40));
	ASSERT_EQ(res2.get(), 8);
}
