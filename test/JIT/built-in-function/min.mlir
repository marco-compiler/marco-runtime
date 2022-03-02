// RUN: modelica-opt %s                             \
// RUN:     --convert-modelica                      \
// RUN:     --convert-modelica-to-cfg               \
// RUN:     --convert-to-llvm                       \
// RUN: | mlir-opt                                  \
// RUN:      --convert-scf-to-std                   \
// RUN: | mlir-cpu-runner                           \
// RUN:     -e main -entry-point-result=void -O0    \
// RUN:     -shared-libs=%runtime_lib               \
// RUN: | FileCheck %s

// CHECK: 1.500000e+00
// CHECK-NEXT: 1.500000e+00
// CHECK-NEXT: -2.500000e+00
// CHECK-NEXT: -2.500000e+00
// CHECK-NEXT: -1.500000e+00
// CHECK-NEXT: -1.500000e+00

func @test_scalars() -> () {
    %size = constant 6 : index

    %y = modelica.alloca %size : index -> !modelica.array<stack, ?x!modelica.real>
    %x = modelica.alloca %size : index -> !modelica.array<stack, ?x!modelica.real>

    %c0 = constant 0 : index
    %x0 = modelica.constant #modelica.real<1.5> : !modelica.real
    %y0 = modelica.constant #modelica.real<2.5> : !modelica.real
    modelica.store %x[%c0], %x0 : !modelica.array<stack, ?x!modelica.real>
    modelica.store %y[%c0], %y0 : !modelica.array<stack, ?x!modelica.real>

    %c1 = constant 1 : index
    %x1 = modelica.constant #modelica.real<2.5> : !modelica.real
    %y1 = modelica.constant #modelica.real<1.5> : !modelica.real
    modelica.store %x[%c1], %x1 : !modelica.array<stack, ?x!modelica.real>
    modelica.store %y[%c1], %y1 : !modelica.array<stack, ?x!modelica.real>

    %c2 = constant 2 : index
    %x2 = modelica.constant #modelica.real<-1.5> : !modelica.real
    %y2 = modelica.constant #modelica.real<-2.5> : !modelica.real
    modelica.store %x[%c2], %x2 : !modelica.array<stack, ?x!modelica.real>
    modelica.store %y[%c2], %y2 : !modelica.array<stack, ?x!modelica.real>

    %c3 = constant 3 : index
    %x3 = modelica.constant #modelica.real<-2.5> : !modelica.real
    %y3 = modelica.constant #modelica.real<-1.5> : !modelica.real
    modelica.store %x[%c3], %x3 : !modelica.array<stack, ?x!modelica.real>
    modelica.store %y[%c3], %y3 : !modelica.array<stack, ?x!modelica.real>

    %c4 = constant 4 : index
    %x4 = modelica.constant #modelica.real<1.5> : !modelica.real
    %y4 = modelica.constant #modelica.real<-1.5> : !modelica.real
    modelica.store %x[%c4], %x4 : !modelica.array<stack, ?x!modelica.real>
    modelica.store %y[%c4], %y4 : !modelica.array<stack, ?x!modelica.real>

    %c5 = constant 5 : index
    %x5 = modelica.constant #modelica.real<-1.5> : !modelica.real
    %y5 = modelica.constant #modelica.real<1.5> : !modelica.real
    modelica.store %x[%c5], %x5 : !modelica.array<stack, ?x!modelica.real>
    modelica.store %y[%c5], %y5 : !modelica.array<stack, ?x!modelica.real>

    scf.for %i = %c0 to %size step %c1 {
      %xi = modelica.load %x[%i] : !modelica.array<stack, ?x!modelica.real>
      %yi = modelica.load %y[%i] : !modelica.array<stack, ?x!modelica.real>
      %result = modelica.min %xi, %yi : (!modelica.real, !modelica.real) -> !modelica.real
      modelica.print %result : !modelica.real
    }

    return
}

// CHECK-NEXT: -3.500000e+00

func @test_array() -> () {
    %size = constant 6 : index
    %array = modelica.alloca %size : index -> !modelica.array<stack, ?x!modelica.real>

    %c0 = constant 0 : index
    %0 = modelica.constant #modelica.real<-1.5> : !modelica.real
    modelica.store %array[%c0], %0 : !modelica.array<stack, ?x!modelica.real>

    %c1 = constant 1 : index
    %1 = modelica.constant #modelica.real<9.5> : !modelica.real
    modelica.store %array[%c1], %1 : !modelica.array<stack, ?x!modelica.real>

    %c2 = constant 2 : index
    %2 = modelica.constant #modelica.real<-3.5> : !modelica.real
    modelica.store %array[%c2], %2 : !modelica.array<stack, ?x!modelica.real>

    %c3 = constant 3 : index
    %3 = modelica.constant #modelica.real<0.0> : !modelica.real
    modelica.store %array[%c3], %3 : !modelica.array<stack, ?x!modelica.real>

    %c4 = constant 4 : index
    %4 = modelica.constant #modelica.real<4.5> : !modelica.real
    modelica.store %array[%c4], %4 : !modelica.array<stack, ?x!modelica.real>

    %c5 = constant 5 : index
    %5 = modelica.constant #modelica.real<3.5> : !modelica.real
    modelica.store %array[%c5], %5 : !modelica.array<stack, ?x!modelica.real>

    %result = modelica.min %array : !modelica.array<stack, ?x!modelica.real> -> !modelica.real
    modelica.print %result : !modelica.real

    return
}

func @main() -> () {
    call @test_scalars() : () -> ()
    call @test_array() : () -> ()
    return
}