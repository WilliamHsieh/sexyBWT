#include <bits/stdc++.h>
#include "psais/psais.hpp"

using namespace std;

int main() {
	cin.tie(0) -> sync_with_stdio(0);

	string str, buf;
	getline(cin, buf);
	while (getline(cin, buf)) {
		str += buf;
	}

	auto result = psais::suffix_array<uint32_t>(str);
}
