// RUN: marco %s --omc-bypass --emit-mlir -o %t
// RUN: cat %t | FileCheck %s

// CHECK-LABEL: @foo
// CHECK: modelica.mod
// CHECK-SAME: (!modelica.int, !modelica.int) -> !modelica.int

function foo
    input Integer x;
    input Integer y;
    output Integer z;
algorithm
    z := mod(x, y);
end foo;
