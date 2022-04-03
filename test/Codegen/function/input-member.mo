// RUN: marco %s --omc-bypass --emit-modelica-dialect | FileCheck %s

// CHECK-LABEL: @inputBooleanScalar
// CHECK: modelica.member_create
// CHECK-SAME: name = "x"
// CHECK-SAME: !modelica.member<!modelica.bool, input>

function inputBooleanScalar
    input Boolean x;

algorithm
end inputBooleanScalar;


// CHECK-LABEL: @inputIntegerScalar
// CHECK: modelica.member_create
// CHECK-SAME: name = "x"
// CHECK-SAME: !modelica.member<!modelica.int, input>

function inputIntegerScalar
    input Integer x;

algorithm
end inputIntegerScalar;


// CHECK-LABEL: @inputRealScalar
// CHECK: modelica.member_create
// CHECK-SAME: name = "x"
// CHECK-SAME: !modelica.member<!modelica.real, input>

function inputRealScalar
    input Real x;

algorithm
end inputRealScalar;


// CHECK-LABEL: @inputBooleanStaticArray
// CHECK: modelica.member_create
// CHECK-SAME: name = "x"
// CHECK-SAME: !modelica.member<3x2x!modelica.bool, input>

function inputBooleanStaticArray
    input Boolean[3,2] x;

algorithm
end inputBooleanStaticArray;


// CHECK-LABEL: @inputBooleanDynamicArray
// CHECK: modelica.member_create
// CHECK-SAME: name = "x"
// CHECK-SAME: !modelica.member<?x?x!modelica.bool, input>

function inputBooleanDynamicArray
    input Boolean[:,:] x;

algorithm
end inputBooleanDynamicArray;


// CHECK-LABEL: @inputIntegerStaticArray
// CHECK: modelica.member_create
// CHECK-SAME: name = "x"
// CHECK-SAME: !modelica.member<3x2x!modelica.int, input>

function inputIntegerStaticArray
    input Integer[3,2] x;

algorithm
end inputIntegerStaticArray;


// CHECK-LABEL: @inputIntegerDynamicArray
// CHECK: modelica.member_create
// CHECK-SAME: name = "x"
// CHECK-SAME: !modelica.member<?x?x!modelica.int, input>

function inputIntegerDynamicArray
    input Integer[:,:] x;

algorithm
end inputIntegerDynamicArray;


// CHECK-LABEL: @inputRealStaticArray
// CHECK: modelica.member_create
// CHECK-SAME: name = "x"
// CHECK-SAME: !modelica.member<3x2x!modelica.real, input>

function inputRealStaticArray
    input Real[3,2] x;

algorithm
end inputRealStaticArray;


// CHECK-LABEL: @inputRealDynamicArray
// CHECK: modelica.member_create
// CHECK-SAME: name = "x"
// CHECK-SAME: !modelica.member<?x?x!modelica.real, input>

function inputRealDynamicArray
    input Real[:,:] x;

algorithm
end inputRealDynamicArray;