// RUN: marco --omc-bypass -c -O0 -c-wrappers -o %basename_t.mo %s.mo
// RUN: g++ %s -c -std=c++1z -I %marco_include_dirs -I %llvm_include_dirs -o %basename_t.o
// RUN: g++ %basename_t.o %basename_t.mo.o %runtime_lib $(llvm-config --ldflags --libs) -Wl,-R%libs -o %t
// RUN: %t | FileCheck %s

// CHECK-LABEL: results
// CHECK-NEXT: false
// CHECK-NEXT: false
// CHECK-NEXT: false
// CHECK-NEXT: true
// CHECK-NEXT: true

#include <array>
#include <iostream>
#include <llvm/ADT/STLExtras.h>

extern "C" bool __modelica_ciface_foo(long x, double y);

using namespace std;

int main() {
	array<long, 5> x = { 0, 10, 10, 11, -10 };
	array<double, 5> y = { 0.5, 10.5, 10, 10.5, -10.5 };

	cout << "results" << endl;

	for (const auto& [x, y] : llvm::zip(x, y))
		cout << boolalpha << __modelica_ciface_foo(x, y) << endl;

	return 0;
}