// RUN: modelica-opt %s --split-input-file --canonicalize | FileCheck %s

// CHECK-LABEL: @test
// CHECK-NEXT: %[[VALUE:[a-zA-Z0-9]*]] = modelica.constant #modelica.real
// CHECK-SAME: -1.000000e+00
// CHECK-NEXT: return %[[VALUE]]

func @test() -> (!modelica.real) {
    %x = modelica.constant #modelica.real<-2.5>
    %result = modelica.sign %x : !modelica.real -> !modelica.real
    return %result : !modelica.real
}

// -----

// CHECK-LABEL: @test
// CHECK-NEXT: %[[VALUE:[a-zA-Z0-9]*]] = modelica.constant #modelica.real
// CHECK-SAME: 0.000000e+00
// CHECK-NEXT: return %[[VALUE]]

func @test() -> (!modelica.real) {
    %x = modelica.constant #modelica.real<0.0>
    %result = modelica.sign %x : !modelica.real -> !modelica.real
    return %result : !modelica.real
}

// -----

// CHECK-LABEL: @test
// CHECK-NEXT: %[[VALUE:[a-zA-Z0-9]*]] = modelica.constant #modelica.real
// CHECK-SAME: 1.000000e+00
// CHECK-NEXT: return %[[VALUE]]

func @test() -> (!modelica.real) {
    %x = modelica.constant #modelica.real<2.5>
    %result = modelica.sign %x : !modelica.real -> !modelica.real
    return %result : !modelica.real
}