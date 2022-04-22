// RUN: modelica-opt %s --split-input-file --convert-modelica-to-cfg | FileCheck %s

// CHECK:       func @foo() {
// CHECK:           br ^[[while_condition:[a-zA-Z0-9]*]]
// CHECK:       ^[[while_condition]]:
// CHECK:           cond_br %{{[a-zA-Z0-9]*}}, ^[[while_body:[a-zA-Z0-9]*]], ^[[while_out:[a-zA-Z0-9]*]]
// CHECK-NEXT:  ^[[while_body:[a-zA-Z0-9]*]]:
// CHECK:           cond_br %{{[a-zA-Z0-9]*}}, ^[[if_then:[a-zA-Z0-9]*]], ^[[if_out:[a-zA-Z0-9]*]]
// CHECK-NEXT:  ^[[if_then]]:
// CHECK:           br ^[[while_out]]
// CHECK-NEXT:  ^[[if_out]]:
// CHECK:           br ^[[while_condition]]
// CHECK-NEXT:  ^[[while_out]]:
// CHECK:           br ^[[out:[a-zA-Z0-9]*]]
// CHECK-NEXT:  ^[[out]]:
// CHECK:           return
// CHECK-NEXT:  }

modelica.function @foo : () -> () {
    modelica.while {
        %0 = modelica.constant #modelica.bool<true>
        modelica.condition (%0 : !modelica.bool)
    } do {
        %0 = modelica.constant #modelica.bool<true>

        modelica.if (%0 : !modelica.bool) {
            modelica.break
        }
    }
}

// -----

// CHECK:       func @foo() {
// CHECK:           br ^[[while_1_condition:[a-zA-Z0-9]*]]
// CHECK:       ^[[while_1_condition]]:
// CHECK:           cond_br %{{[a-zA-Z0-9]*}}, ^[[while_1_body:[a-zA-Z0-9]*]], ^[[while_1_out:[a-zA-Z0-9]*]]
// CHECK-NEXT:  ^[[while_1_body]]:
// CHECK:           br ^[[while_2_condition:[a-zA-Z0-9]*]]
// CHECK:       ^[[while_2_condition]]:
// CHECK:           cond_br %{{[a-zA-Z0-9]*}}, ^[[while_2_body:[a-zA-Z0-9]*]], ^[[while_2_out:[a-zA-Z0-9]*]]
// CHECK-NEXT:  ^[[while_2_body]]:
// CHECK:           br ^[[while_2_out]]
// CHECK-NEXT:  ^[[while_2_out]]:
// CHECK:           cond_br %{{[a-zA-Z0-9]*}}, ^[[if_then:[a-zA-Z0-9]*]], ^[[if_out:[a-zA-Z0-9]*]]
// CHECK-NEXT:  ^[[if_then]]:
// CHECK:           br ^[[while_1_out]]
// CHECK-NEXT:  ^[[if_out]]:
// CHECK:           br ^[[while_1_condition]]
// CHECK-NEXT:  ^[[while_1_out]]:
// CHECK:           br ^[[out:[a-zA-Z0-9]*]]
// CHECK-NEXT:  ^[[out]]:
// CHECK:           return
// CHECK-NEXT:  }

modelica.function @foo : () -> () {
    modelica.while {
        %0 = modelica.constant #modelica.bool<true>
        modelica.condition (%0 : !modelica.bool)
    } do {
        modelica.while {
            %0 = modelica.constant #modelica.bool<true>
            modelica.condition (%0 : !modelica.bool)
        } do {
            modelica.break
        }

        %0 = modelica.constant #modelica.bool<true>

        modelica.if (%0 : !modelica.bool) {
            modelica.break
        }
    }
}

// -----

// CHECK:       func @foo() {
// CHECK:           br ^[[while_condition:[a-zA-Z0-9]*]]
// CHECK:       ^[[while_condition]]:
// CHECK:           cond_br %{{[a-zA-Z0-9]*}}, ^[[while_body:[a-zA-Z0-9]*]], ^[[while_out:[a-zA-Z0-9]*]]
// CHECK-NEXT:  ^[[while_body]]:
// CHECK:           br ^[[for_condition:[a-zA-Z0-9]*]]
// CHECK-NEXT:  ^[[for_condition]]:
// CHECK:           cond_br %{{[a-zA-Z0-9]*}}, ^[[for_body:[a-zA-Z0-9]*]], ^[[for_out:[a-zA-Z0-9]*]]
// CHECK-NEXT:  ^[[for_body]]:
// CHECK:           cond_br %{{[a-zA-Z0-9]*}}, ^[[if_then:[a-zA-Z0-9]*]], ^[[if_out:[a-zA-Z0-9]*]]
// CHECK-NEXT:  ^[[if_then]]:
// CHECK:           br ^[[for_out]]
// CHECK-NEXT:  ^[[if_out]]:
// CHECK:           br ^[[for_step:[a-zA-Z0-9]*]]
// CHECK-NEXT:  ^[[for_step]]:
// CHECK:           br ^[[for_condition]]
// CHECK-NEXT:  ^[[for_out]]:
// CHECK:           br ^[[while_condition]]
// CHECK-NEXT:  ^[[while_out]]:
// CHECK-NEXT:      br ^[[out:[a-zA-Z0-9]*]]
// CHECK-NEXT:  ^[[out]]:
// CHECK:           return
// CHECK-NEXT:  }

modelica.function @foo : () -> () {
    modelica.while {
        %3 = modelica.constant #modelica.bool<true>
        modelica.condition (%3 : !modelica.bool)
    } do {
        modelica.for condition {
            %3 = modelica.constant #modelica.bool<true>
            modelica.condition (%3 : !modelica.bool)
        } body {
            %3 = modelica.constant #modelica.bool<true>

            modelica.if (%3 : !modelica.bool) {
                modelica.break
            }

            modelica.yield
        } step {
            modelica.yield
        }
    }
}