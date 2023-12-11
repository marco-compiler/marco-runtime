#ifndef MARCO_RUNTIME_MANGLING_H
#define MARCO_RUNTIME_MANGLING_H

#include <cstdint>

#define NUM_ARGS_H1(dummy, x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0, ...) x0
#define NUM_ARGS(...) NUM_ARGS_H1(dummy, ##__VA_ARGS__, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)

#define APPLY0(t, n, dummy)
#define APPLY1(t, n, a) t(n, a)
#define APPLY2(t, n, a, b) APPLY1(t, 1, a), t(n, b)
#define APPLY3(t, n, a, b, c) APPLY2(t, 2, a, b), t(n, c)
#define APPLY4(t, n, a, b, c, d) APPLY3(t, 3, a, b, c), t(n, d)
#define APPLY5(t, n, a, b, c, d, e) APPLY4(t, 4, a, b, c, d), t(n, e)
#define APPLY6(t, n, a, b, c, d, e, f) APPLY5(t, 5, a, b, c, d, e), t(n, f)
#define APPLY7(t, n, a, b, c, d, e, f, g) APPLY6(t, 6, a, b, c, d, e, f), t(n, g)
#define APPLY8(t, n, a, b, c, d, e, f, g, h) APPLY7(t, 7, a, b, c, d, e, f, g), t(n, h)
#define APPLY9(t, n, a, b, c, d, e, f, g, h, i) APPLY8(t, 8, a, b, c, d, e, f, g, h), t(n, i)
#define APPLY10(t, n, a, b, c, d, e, f, g, h, i, j) APPLY9(t, 9, a, b, c, d, e, f, g, h, i), t(n, j)
#define APPLY11(t, n, a, b, c, d, e, f, g, h, i, j, k) APPLY10(t, 9, a, b, c, d, e, f, g, h, i, j), t(n, k)

#define APPLY_ALL_H3(t, n, ...) APPLY##n(t, n, __VA_ARGS__)
#define APPLY_ALL_H2(t, n, ...) APPLY_ALL_H3(t, n, __VA_ARGS__)
#define APPLY_ALL(t, ...) APPLY_ALL_H2(t, NUM_ARGS(__VA_ARGS__), __VA_ARGS__)

#define CONCAT0(dummy)
#define CONCAT1(a) a
#define CONCAT2(a, b) a ##b
#define CONCAT3(a, b, c) a ##b ##c
#define CONCAT4(a, b, c, d) a ##b ##c ##d
#define CONCAT5(a, b, c, d, e) a ##b ##c ##d ##e
#define CONCAT6(a, b, c, d, e, f) a ##b ##c ##d ##e ##f
#define CONCAT7(a, b, c, d, e, f, g) a ##b ##c ##d ##e ##f ##g
#define CONCAT8(a, b, c, d, e, f, g, h) a ##b ##c ##d ##e ##f ##g ##h
#define CONCAT9(a, b, c, d, e, f, g, h, i) a ##b ##c ##d ##e ##f ##g ##h ##i
#define CONCAT10(a, b, c, d, e, f, g, h, i, j) a ##b ##c ##d ##e ##f ##g ##h ##i ##j
#define CONCAT11(a, b, c, d, e, f, g, h, i, j, k) a ##b ##c ##d ##e ##f ##g ##h ##i ##j ##k

#define CONCAT_ALL_H3(n, ...) CONCAT##n(__VA_ARGS__)
#define CONCAT_ALL_H2(n, ...) CONCAT_ALL_H3(n, __VA_ARGS__)
#define CONCAT_ALL(...) CONCAT_ALL_H2(NUM_ARGS(__VA_ARGS__), __VA_ARGS__)

#define void_CPP void
#define bool_CPP bool
#define int32_t_CPP int32_t
#define uint32_t_CPP uint32_t
#define int64_t_CPP int64_t
#define uint64_t_CPP uint64_t
#define float_CPP float
#define double_CPP double

#define _MANGLED
#define void_MANGLED _void
#define bool_MANGLED _i1
#define int32_t_MANGLED _i32
#define uint32_t_MANGLED _i32
#define int64_t_MANGLED _i64
#define uint64_t_MANGLED _i64
#define float_MANGLED _f32
#define double_MANGLED _f64

#define TYPE_CPP(n, type) type ##_CPP
#define TYPES_CPP(...) APPLY_ALL(TYPE_CPP, __VA_ARGS__)

#define ARRAY(type) ARRAY_ ##type

#define ARRAY_bool_CPP UnrankedMemRefType<bool_CPP>*
#define ARRAY_int32_t_CPP UnrankedMemRefType<int32_t_CPP>*
#define ARRAY_uint32_t_CPP UnrankedMemRefType<uint32_t_CPP>*
#define ARRAY_int64_t_CPP UnrankedMemRefType<int64_t_CPP>*
#define ARRAY_uint64_t_CPP UnrankedMemRefType<uint64_t_CPP>*
#define ARRAY_float_CPP UnrankedMemRefType<float_CPP>*
#define ARRAY_double_CPP UnrankedMemRefType<double_CPP>*

#define ARRAY_bool_MANGLED _ai1
#define ARRAY_int32_t_MANGLED _ai32
#define ARRAY_uint32_t_MANGLED _ai32
#define ARRAY_int64_t_MANGLED _ai64
#define ARRAY_uint64_t_MANGLED _ai64
#define ARRAY_float_MANGLED _af32
#define ARRAY_double_MANGLED _af64

#define PTR(type) PTR_ ##type

#define PTR_void_CPP void*
#define PTR_bool_CPP bool*
#define PTR_int32_t_CPP int32_t*
#define PTR_uint32_t_CPP uint32_t*
#define PTR_int64_t_CPP int64_t*
#define PTR_uint64_t_CPP uint64_t*
#define PTR_float_CPP float*
#define PTR_double_CPP double*

#define PTR_void_MANGLED _pvoid
#define PTR_bool_MANGLED _pi1
#define PTR_int32_t_MANGLED _pi32
#define PTR_uint32_t_MANGLED _pi32
#define PTR_int64_t_MANGLED _pi64
#define PTR_uint64_t_MANGLED _pi64
#define PTR_float_MANGLED _pf32
#define PTR_double_MANGLED _pf64

#define TYPE_MANGLED(n, type) type ##_MANGLED
#define TYPES_MANGLED(...) APPLY_ALL(TYPE_MANGLED, __VA_ARGS__)

#define MODELICA_PREFIX _M

#define NAME_MANGLED(name, ...) CONCAT_ALL(MODELICA_PREFIX, name, TYPES_MANGLED(__VA_ARGS__))

#define ARG_NAME(n, type) arg ##n
#define ARGS_NAMES(...) APPLY_ALL(ARG_NAME, __VA_ARGS__)

#define ARG_DECLARATION(n, type) type ARG_NAME(n, type)
#define ARGS_DECLARATIONS(...) APPLY_ALL(ARG_DECLARATION, __VA_ARGS__)

#define RUNTIME_FUNC_SIGNATURE(name, returnType, ...) \
	TYPES_CPP(returnType) NAME_MANGLED(name, returnType, __VA_ARGS__) (ARGS_DECLARATIONS(TYPES_CPP(__VA_ARGS__)))

#define RUNTIME_FUNC_DECL(name, returnType, ...) \
	extern "C" TYPES_CPP(returnType) NAME_MANGLED(name, returnType, __VA_ARGS__) (ARGS_DECLARATIONS(TYPES_CPP(__VA_ARGS__)));

#define RUNTIME_FUNC_DEF(name, returnType, ...) \
	TYPES_CPP(returnType) NAME_MANGLED(name, returnType, __VA_ARGS__) (ARGS_DECLARATIONS(TYPES_CPP(__VA_ARGS__))) \
	{ \
		return CONCAT_ALL(name, TYPES_MANGLED(returnType))(ARGS_NAMES(__VA_ARGS__)); \
	}

#endif	// MARCO_RUNTIME_MANGLING_H
