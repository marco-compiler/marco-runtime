// RUN: marco %s.mo --emit-c-wrappers | llc -o %basename_t.mo.s
// RUN: clang++ %s -g -c -std=c++1z -I %runtime_h -I %llvm_include_dirs -o %basename_t.o
// RUN: clang++ %basename_t.o %basename_t.mo.s -g $(llvm-config --ldflags --libs) -o %t
// RUN: %t | FileCheck %s

// CHECK-LABEL: results
// CHECK-NEXT: 57

#include <array>
#include <iostream>

extern "C" long __modelica_ciface_foo();

using namespace std;

int main() {
	cout << "results" << endl;
	cout << __modelica_ciface_foo() << endl;

	return 0;
}