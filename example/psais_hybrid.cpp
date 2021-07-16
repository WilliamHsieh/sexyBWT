#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1

#include <bits/stdc++.h>
#include "psais/psais_hybrid.hpp"

using namespace std;

int main() {
	cin.tie(0) -> sync_with_stdio(0);

	string str, buf;
	getline(cin, buf);
	while (getline(cin, buf)) {
		str += buf;
	}

	auto result = psais::suffix_array_hybrid<uint32_t>(str, 64);
}
