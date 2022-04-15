// RUN: marco --omc-bypass --model=InductionUsage --end-time=1 --time-step=0.1 --solver=ida --ida-equidistant-time-grid -o %basename_t %s
// RUN: ./%basename_t | FileCheck %s

// CHECK: time;x[1];x[2];x[3];x[4];x[5]
// CHECK-NEXT: 0.000000000000;0.500000000000;0.500000000000;0.500000000000;0.500000000000;0.500000000000
// CHECK-NEXT: 0.100000000000;0.520000000000;0.580000000000;0.620000000000;0.660000000000;0.700000000000
// CHECK-NEXT: 0.200000000000;0.540000000000;0.660000000000;0.740000000000;0.820000000000;0.900000000000
// CHECK-NEXT: 0.300000000000;0.560000000000;0.740000000000;0.860000000000;0.980000000000;1.100000000000
// CHECK-NEXT: 0.400000000000;0.580000000000;0.820000000000;0.980000000000;1.140000000000;1.300000000000
// CHECK-NEXT: 0.500000000000;0.600000000000;0.900000000000;1.100000000000;1.300000000000;1.500000000000
// CHECK-NEXT: 0.600000000000;0.620000000000;0.980000000000;1.220000000000;1.460000000000;1.700000000000
// CHECK-NEXT: 0.700000000000;0.640000000000;1.060000000000;1.340000000000;1.620000000000;1.900000000000
// CHECK-NEXT: 0.800000000000;0.660000000000;1.140000000000;1.460000000000;1.780000000000;2.100000000000
// CHECK-NEXT: 0.900000000000;0.680000000000;1.220000000000;1.580000000000;1.940000000000;2.300000000000
// CHECK-NEXT: 1.000000000000;0.700000000000;1.300000000000;1.700000000000;2.100000000000;2.500000000000

model InductionUsage
	Real[5] x(start = 0.5);
equation
	5.0 * der(x[1]) = 1.0;

	for i in 2:5 loop
		5.0 * der(x[i]) = 2.0 * i;
	end for;
end InductionUsage;
