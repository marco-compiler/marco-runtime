// RUN: marco %s --omc-bypass --emit-modelica-dialect | FileCheck %s

// CHECK-LABEL: body
// CHECK-NEXT:  ^bb0(%[[arg0:[a-zA-Z0-9]*]]: !modelica.array<3x!modelica.real>):
// CHECK-NEXT:    modelica.equation {
// CHECK-NEXT:      %[[array:[a-zA-Z0-9]*]] = modelica.alloc  : !modelica.array<3x!modelica.real>
// CHECK-NEXT:      %[[index0:[a-zA-Z0-9]*]] = modelica.constant 0 : index
// CHECK-NEXT:      %[[subscription0:[a-zA-Z0-9]*]] = modelica.subscription %[[array]][%[[index0]]] : !modelica.array<!modelica.real>
// CHECK-NEXT:      %[[cst0:[a-zA-Z0-9]*]] = modelica.constant #modelica.int<1> : !modelica.int
// CHECK-NEXT:      modelica.assignment %[[subscription0]], %[[cst0]]
// CHECK-NEXT:      %[[index1:[a-zA-Z0-9]*]] = modelica.constant 1 : index
// CHECK-NEXT:      %[[subscription1:[a-zA-Z0-9]*]] = modelica.subscription %[[array]][%[[index1]]] : !modelica.array<!modelica.real>
// CHECK-NEXT:      %[[cst1:[a-zA-Z0-9]*]] = modelica.constant #modelica.real<2.000000e+00> : !modelica.real
// CHECK-NEXT:      modelica.assignment %[[subscription1]], %[[cst1]]
// CHECK-NEXT:      %[[index2:[a-zA-Z0-9]*]] = modelica.constant 2 : index
// CHECK-NEXT:      %[[subscription2:[a-zA-Z0-9]*]] = modelica.subscription %[[array]][%[[index2]]] : !modelica.array<!modelica.real>
// CHECK-NEXT:      %[[cst2:[a-zA-Z0-9]*]] = modelica.constant #modelica.int<3> : !modelica.int
// CHECK-NEXT:      modelica.assignment %[[subscription2]], %[[cst2]]
// CHECK-DAG:       %[[lhs:[a-zA-Z0-9]*]] = modelica.equation_side %[[arg0]] : tuple<!modelica.array<3x!modelica.real>>
// CHECK-DAG:       %[[rhs:[a-zA-Z0-9]*]] = modelica.equation_side %[[array]] : tuple<!modelica.array<3x!modelica.real>>
// CHECK-NEXT:      modelica.equation_sides %[[lhs]], %[[rhs]]
// CHECK-NEXT:    }

model MemberWithStartArray
	Real[3] x = {1, 2.0, 3};
equation
end MemberWithStartArray;