// RUN: marco %s.mo --emit-c-wrappers | llc -o %basename_t.mo.s
// RUN: clang++ %s -g -c -std=c++1z -I %runtime_h -I %llvm_include_dirs -o %basename_t.o
// RUN: clang++ %basename_t.o %basename_t.mo.s -g -L%libs/runtime $(llvm-config --ldflags --libs) -lruntime -Wl,-R%libs/runtime -o %t
// RUN: %t | FileCheck %s

// CHECK-LABEL: results
// CHECK-NEXT{LITERAL}: [2, 4, 6]

#include <array>
#include <iostream>
#include <modelica/runtime/ArrayDescriptor.h>

extern "C" void __modelica_ciface_foo(long x, ArrayDescriptor<long, 1>* y, ArrayDescriptor<long, 1>* z);

using namespace std;

int main() {
	long x = 2;

	array<long, 3> y = { 1, 2, 3 };
	ArrayDescriptor<long, 1> yDescriptor(y);

	array<long, 3> z = { 0, 0, 0 };
	ArrayDescriptor<long, 1> zDescriptor(z);

	__modelica_ciface_foo(x, &yDescriptor, &zDescriptor);

	cout << "results" << endl;
	cout << zDescriptor << endl;

	return 0;
}