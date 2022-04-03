// RUN: marco --omc-bypass -c -O0 -c-wrappers -o %basename_t.mo %s.mod
// RUN: g++ %s -c -std=c++1z -I %marco_include_dirs -I %llvm_include_dirs -o %basename_t.o
// RUN: g++ %basename_t.o %basename_t.mo.o %runtime_lib $(llvm-config --ldflags --libs) -Wl,-R%libs -o %t
// RUN: %t | FileCheck %s

// CHECK-LABEL: results
// CHECK-NEXT{LITERAL}: [1, 2, 3]

#include "marco/Runtime/ArrayDescriptor.h"
#include <array>
#include <iostream>

extern "C" void __modelica_ciface_foo(ArrayDescriptor<long, 1>* y);

using namespace std;

int main() {
	array<long, 3> y = { 0, 0, 0 };
	ArrayDescriptor<long, 1> yDescriptor(y);

	__modelica_ciface_foo(&yDescriptor);

	cout << "results" << endl;
	cout << yDescriptor << endl;

	return 0;
}