// RUN: marco %s.mo --emit-c-wrappers | llc -o %basename_t.mo.s
// RUN: clang++ %s -g -c -std=c++1z -I %runtime_h -I %llvm_include_dirs -o %basename_t.o
// RUN: clang++ %basename_t.o %basename_t.mo.s -g $(llvm-config --ldflags --libs) -o %t
// RUN: %t | FileCheck %s

// CHECK-LABEL: results
// CHECK-NEXT: -10
// CHECK-NEXT: -23
// CHECK-NEXT: 57

#include <array>
#include <iostream>
#include <modelica/runtime/ArrayDescriptor.h>

extern "C" void __modelica_ciface_foo(
		ArrayDescriptor<double, 1>* y, ArrayDescriptor<double, 1>* x);

using namespace std;

int main() {
	array<double, 3> x = { 10, 23, -57 };
	ArrayDescriptor<double, 1> xDescriptor(x);

	array<double, 3> y = { 10, 23, -57 };
	ArrayDescriptor<double, 1> yDescriptor(y);

	__modelica_ciface_foo(&yDescriptor, &xDescriptor);

	cout << "results" << endl;

	for (const auto& value : yDescriptor)
		cout << value << endl;

	delete[] yDescriptor.getData();
	return 0;
}
