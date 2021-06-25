#include <iostream>
#include <algorithm>
#include <random>
#include <cassert>
#include "psais/psais.hpp"

template<typename IndexType>
auto suffix_array(std::string_view s) {
	auto res = std::vector<IndexType>(s.size() + 1);
	std::iota(res.begin(), res.end(), 0);
	std::sort(res.begin(), res.end(), [s](IndexType a, IndexType b) {
		return s.substr(a) < s.substr(b);
	});
	return res;
}

int main(int argc, char**) {
	using IndexType = uint32_t;
	if (argc == 2) {
		std::cout << "pid: " << getpid() << std::endl;
		sleep(10);
	}

	auto len = 10000000;
	auto gen = std::mt19937(std::random_device{}());
	auto dis = std::uniform_int_distribution(0, 25);

	for (int i = 0; i < 10; i++) {
		std::cout << i + 1 << '\n';

		std::string str;
		for (int j = 0; j < len; j++) {
			str += 'a' + dis(gen);
		}

		auto ans = suffix_array<IndexType>(str);
		auto suf = psais::suffix_array<IndexType>(str);

		for (int j = 0; j <= len; j++) {
			assert(ans[j] == suf[j]);
		}
	}
}
