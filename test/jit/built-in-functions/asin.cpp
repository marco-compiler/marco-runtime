// RUN: marco %s.mo --O0 --emit-c-wrappers | llc -o %basename_t.mo.s
// RUN: clang++ %s -g -c -std=c++1z -I %runtime_h -I %llvm_include_dirs -o %basename_t.o
// RUN: clang++ %basename_t.o %basename_t.mo.s -g $(llvm-config --ldflags --libs) %runtime_lib -Wl,-R%libs/runtime -o %t
// RUN: %t | FileCheck %s

// CHECK-LABEL: results
// CHECK-NEXT: 1.570796e+00
// CHECK-NEXT: 1.047198e+00
// CHECK-NEXT: 7.853982e-01
// CHECK-NEXT: 0.000000e+00
// CHECK-NEXT: -7.853982e-01
// CHECK-NEXT: -1.047198e+00
// CHECK-NEXT: -1.570796e+00

#include <array>
#include <iostream>
#include <marco/runtime/ArrayDescriptor.h>

extern "C" double __modelica_ciface_foo(double x);

using namespace std;

int main() {
	array<double, 7> x = { 1, 0.866025403, 0.707106781, 0, -0.707106781, -0.866025403, -1 };
	cout << "results" << endl;

	for (const auto& value : x)
		cout << scientific << __modelica_ciface_foo(value) << endl;

	return 0;
}
