// RUN: marco %s.mo --O0 --emit-c-wrappers | llc -o %basename_t.mo.s
// RUN: clang++ %s -g -c -std=c++1z -I %runtime_h -I %llvm_include_dirs -o %basename_t.o
// RUN: clang++ %basename_t.o %basename_t.mo.s -g -L%libs/runtime $(llvm-config --ldflags --libs) -lruntime -Wl,-R%libs/runtime -o %t
// RUN: %t | FileCheck %s

// CHECK-LABEL: results
// CHECK-NEXT{LITERAL}: [[1, 2, 3], [2, 5, 6], [3, 6, 9]]

#include <array>
#include <iostream>
#include <modelica/runtime/ArrayDescriptor.h>

extern "C" void __modelica_ciface_foo(ArrayDescriptor<long, 2>* y, ArrayDescriptor<long, 2>* x);

using namespace std;

int main() {
	array<long, 9> x = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	ArrayDescriptor<long, 2> xDescriptor(x.data(), { 3, 3 });

	ArrayDescriptor<long, 2> yDescriptor(nullptr, { 1, 1 });

	__modelica_ciface_foo(&yDescriptor, &xDescriptor);

	cout << "results" << endl;
	cout << yDescriptor << endl;
	free(yDescriptor.getData());

	return 0;
}
