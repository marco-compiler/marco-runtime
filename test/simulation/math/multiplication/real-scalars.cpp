// RUN: marco --omc-bypass -c -O0 -c-wrappers -o %basename_t.mo %s.mo
// RUN: g++ %s -c -std=c++1z -I %marco_include_dirs -I %llvm_include_dirs -o %basename_t.o
// RUN: g++ %basename_t.o %basename_t.mo.o %runtime_lib $(llvm-config --ldflags --libs) -Wl,-R%libs -o %t
// RUN: %t | FileCheck %s

// CHECK-LABEL: results
// CHECK-NEXT: 0.000000e+00
// CHECK-NEXT: 1.925000e+01
// CHECK-NEXT: -1.925000e+01
// CHECK-NEXT: -1.925000e+01

#include <array>
#include <iostream>
#include <llvm/ADT/STLExtras.h>

extern "C" double __modelica_ciface_foo(double x, double y);

using namespace std;

int main() {
	array<double, 4> x = { 0, 3.5, 3.5, -5.5 };
	array<double, 4> y = { 10.5, 5.5, -5.5, 3.5 };

	cout << "results" << endl;

	for (const auto& [x, y] : llvm::zip(x, y))
		cout << scientific << __modelica_ciface_foo(x, y) << endl;

	return 0;
}