#include <iostream>
#include <future>
#include <random>
#include <vector>
#include "psais/utility/thread_pool.hpp"

std::atomic<int> num_job1 = 0;
std::condition_variable cv;
std::mutex cv_m;

void print_1() {
	if (--num_job1 == 0) {
		cv.notify_all();
	}
}

struct print_2 {
	void operator()() {
		std::unique_lock<std::mutex> lock(cv_m);
		cv.wait(lock, []{ return num_job1 == 0; });
		lock.unlock();
	}
};

int fib(int n) {
	if (n <= 1) return n;
	return fib(n - 1) + fib(n - 2);
}

TEST(utility, thread_pool) {
	auto pool = psais::utility::ThreadPool(5);
	auto result = std::vector<std::future<void>>{};

	// send jobs
	for (int i = 0; i < 8; i++, num_job1++) {
		result.push_back(pool.enqueue(print_1));
	}
	for (int i = 0; i < 4; i++) {
		result.push_back(pool.enqueue(print_2{}));
	}

	// get result from future
	for (auto& r : result) {
		r.get();
	}
	ASSERT_EQ(num_job1, 0);

	// different return type(int)
	auto res = pool.enqueue(fib, 30);
	res.wait();

	// lambda
	auto lam = [](int x, int y) { return x + y; };
	auto res2 = pool.enqueue(lam, 3, 5);
	res2.wait();

	ASSERT_EQ(res.get(), fib(30));
	ASSERT_EQ(res2.get(), 8);
}
