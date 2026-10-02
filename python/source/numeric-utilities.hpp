/***************************************************************************
 *            numeric-utilities.hpp
 *
 *  Copyright  2026  Pieter Collins
 *
 ****************************************************************************/

#ifndef ARIADNE_PYTHON_NUMERIC_UTILITIES_HPP
#define ARIADNE_PYTHON_NUMERIC_UTILITIES_HPP

#include "utilities.hpp"
#include "paradigm-utilities.hpp"
#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

using namespace PyBind11;

#define __py_div__ (PY_MAJOR_VERSION>=3) ? "__truediv__" : "__div__"
#define __py_rdiv__ (PY_MAJOR_VERSION>=3) ? "__rtruediv__" : "__rdiv__"

template<class... TS> struct Tag { };

template<class A, class RET=Return<decltype(+declval<A>())>>
ReturnType<RET> __pos__(const A& a) { return static_cast<ReturnType<RET>>(+a); }
template<class A, class RET=Return<decltype(-declval<A>())>>
ReturnType<RET> __neg__(const A& a) { return static_cast<ReturnType<RET>>(-a); }
template<class A1, class A2, class RET=Return<decltype(declval<A1>()+declval<A2>())>>
ReturnType<RET> __add__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1+a2); }
template<class A1, class A2, class RET=Return<decltype(declval<A1>()-declval<A2>())>>
ReturnType<RET> __sub__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1-a2); }
template<class A1, class A2, class RET=Return<decltype(declval<A1>()*declval<A2>())>>
ReturnType<RET> __mul__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1*a2); }
template<class A1, class A2, class RET=Return<decltype(declval<A1>()/declval<A2>())>>
ReturnType<RET> __div__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1/a2); }
template<class A1, class A2, class RET=Return<decltype(declval<A2>()+declval<A1>())>>
ReturnType<RET> __radd__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a2+a1); }
template<class A1, class A2, class RET=Return<decltype(declval<A2>()-declval<A1>())>>
ReturnType<RET> __rsub__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a2-a1); }
template<class A1, class A2, class RET=Return<decltype(declval<A2>()*declval<A1>())>>
ReturnType<RET> __rmul__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a2*a1); }
template<class A1, class A2, class RET=Return<decltype(declval<A2>()/declval<A1>())>>
ReturnType<RET> __rdiv__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a2/a1); }
template<class A1, class A2, class RET=Return<decltype(pow(declval<A1>(),declval<A2>()))>>
ReturnType<RET> __pow__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(pow(a1,a2)); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()==declval<A2>())>>
ReturnType<RET> __eq__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1==a2); }
template<class A1, class A2, class RET=Return<decltype(declval<A1>()!=declval<A2>())>>
ReturnType<RET> __ne__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1!=a2); }
template<class A1, class A2, class RET=Return<decltype(declval<A1>()>declval<A2>())>>
ReturnType<RET> __gt__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1>a2); }
template<class A1, class A2, class RET=Return<decltype(declval<A1>()<declval<A2>())>>
ReturnType<RET> __lt__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1<a2); }
template<class A1, class A2, class RET=Return<decltype(declval<A1>()>=declval<A2>())>>
ReturnType<RET> __ge__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1>=a2); }
template<class A1, class A2, class RET=Return<decltype(declval<A1>()<=declval<A2>())>>
ReturnType<RET> __le__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1<=a2); }

template<class... AS> auto _add_(AS const& ... as) -> decltype(add(as...)) { return add(as...); }
template<class... AS> auto _sub_(AS const& ... as) -> decltype(sub(as...)) { return sub(as...); }
template<class... AS> auto _mul_(AS const& ... as) -> decltype(mul(as...)) { return mul(as...); }
template<class... AS> auto _div_(AS const& ... as) -> decltype(div(as...)) { return div(as...); }
template<class... AS> auto _fma_(AS const& ... as) -> decltype(fma(as...)) { return fma(as...); }
template<class... AS> auto _nul_(AS const& ... as) -> decltype(nul(as...)) { return nul(as...); }
template<class... AS> auto _pos_(AS const& ... as) -> decltype(pos(as...)) { return pos(as...); }
template<class... AS> auto _neg_(AS const& ... as) -> decltype(neg(as...)) { return neg(as...); }
template<class... AS> auto _hlf_(AS const& ... as) -> decltype(hlf(as...)) { return hlf(as...); }
template<class... AS> auto _rec_(AS const& ... as) -> decltype(rec(as...)) { return rec(as...); }
template<class... AS> auto _sqr_(AS const& ... as) -> decltype(sqr(as...)) { return sqr(as...); }
template<class... AS> auto _pow_(AS const& ... as) -> decltype(pow(as...)) { return pow(as...); }
template<class... AS> auto _sqrt_(AS const& ... as) -> decltype(sqrt(as...)) { return sqrt(as...); }
template<class... AS> auto _exp_(AS const& ... as) -> decltype(exp(as...)) { return exp(as...); }
template<class... AS> auto _log_(AS const& ... as) -> decltype(log(as...)) { return log(as...); }
template<class... AS> auto _sin_(AS const& ... as) -> decltype(sin(as...)) { return sin(as...); }
template<class... AS> auto _cos_(AS const& ... as) -> decltype(cos(as...)) { return cos(as...); }
template<class... AS> auto _tan_(AS const& ... as) -> decltype(tan(as...)) { return tan(as...); }
template<class... AS> auto _asin_(AS const& ... as) -> decltype(asin(as...)) { return asin(as...); }
template<class... AS> auto _acos_(AS const& ... as) -> decltype(acos(as...)) { return acos(as...); }
template<class... AS> auto _atan_(AS const& ... as) -> decltype(atan(as...)) { return atan(as...); }
template<class... AS> auto _max_(AS const& ... as) -> decltype(max(as...)) { return max(as...); }
template<class... AS> auto _min_(AS const& ... as) -> decltype(min(as...)) { return min(as...); }
template<class... AS> auto _abs_(AS const& ... as) -> decltype(abs(as...)) { return abs(as...); }
template<class... AS> auto _sgn_(AS const& ... as) -> decltype(sgn(as...)) { return sgn(as...); }
template<class... AS> auto _mag_(AS const& ... as) -> decltype(mag(as...)) { return mag(as...); }
template<class... AS> auto _mig_(AS const& ... as) -> decltype(mig(as...)) { return mig(as...); }
template<class... AS> auto _arg_(AS const& ... as) -> decltype(arg(as...)) { return arg(as...); }
template<class... AS> auto _conj_(AS const& ... as) -> decltype(conj(as...)) { return conj(as...); }
template<class... AS> auto _dist_(AS const& ... as) -> decltype(dist(as...)) { return dist(as...); }
template<class A> auto _log2_(A const& a) -> decltype(log2(a)) { return log2(a); }
template<class... AS> auto _coarsening_(AS... as) -> decltype(coarsening(as...)) { return coarsening(as...); }
template<class... AS> auto _refinement_(AS... as) -> decltype(refinement(as...)) { return refinement(as...); }
template<class... AS> auto _refines_(AS... as) -> decltype(refines(as...)) { return refines(as...); }
template<class... AS> auto _inconsistent_(AS... as) -> decltype(inconsistent(as...)) { return inconsistent(as...); }

template<class X> pybind11::class_<X>& define_lattice(pybind11::module& module, pybind11::class_<X>& pyclass) {
    module.def("abs", &_abs_<X>); module.def("max", &_max_<X,X>); module.def("min", &_min_<X,X>); return pyclass;
}
template<class X, class Y> pybind11::class_<X>& define_mixed_lattice(pybind11::module& module, pybind11::class_<X>& pyclass, Tag<Y> = Tag<Y>()) {
    module.def("max", &_max_<X,Y>); module.def("max", &_max_<Y,X>); module.def("min", &_min_<X,Y>); module.def("min", &_min_<Y,X>); return pyclass;
}
template<class X> pybind11::class_<X>& define_comparisons(pybind11::module&, pybind11::class_<X>& pyclass) {
    pyclass.def("__eq__", &__eq__<X,X>, pybind11::is_operator()); pyclass.def("__ne__", &__ne__<X,X>, pybind11::is_operator());
    pyclass.def("__le__", &__le__<X,X>, pybind11::is_operator()); pyclass.def("__ge__", &__ge__<X,X>, pybind11::is_operator());
    pyclass.def("__lt__", &__lt__<X,X>, pybind11::is_operator()); pyclass.def("__gt__", &__gt__<X,X>, pybind11::is_operator()); return pyclass;
}
template<class X, class Y> pybind11::class_<X>& define_mixed_comparisons(pybind11::module&, pybind11::class_<X>& pyclass, Tag<Y> = Tag<Y>()) {
    pyclass.def("__eq__", &__eq__<X,Y>, pybind11::is_operator()); pyclass.def("__ne__", &__ne__<X,Y>, pybind11::is_operator());
    pyclass.def("__le__", &__le__<X,Y>, pybind11::is_operator()); pyclass.def("__ge__", &__ge__<X,Y>, pybind11::is_operator());
    pyclass.def("__lt__", &__lt__<X,Y>, pybind11::is_operator()); pyclass.def("__gt__", &__gt__<X,Y>, pybind11::is_operator());
    pyclass.def("__eq__", &__eq__<Y,X>, pybind11::is_operator()); pyclass.def("__ne__", &__ne__<Y,X>, pybind11::is_operator());
    pyclass.def("__le__", &__le__<Y,X>, pybind11::is_operator()); pyclass.def("__ge__", &__ge__<Y,X>, pybind11::is_operator());
    pyclass.def("__lt__", &__lt__<Y,X>, pybind11::is_operator()); pyclass.def("__gt__", &__gt__<Y,X>, pybind11::is_operator()); return pyclass;
}
template<class X> pybind11::class_<X>& define_arithmetic(pybind11::module&, pybind11::class_<X>& pyclass) {
    pyclass.def("__pos__", &__pos__<X>, pybind11::is_operator()); pyclass.def("__neg__", &__neg__<X>, pybind11::is_operator());
    pyclass.def("__add__", &__add__<X,X>, pybind11::is_operator()); pyclass.def("__sub__", &__sub__<X,X>, pybind11::is_operator());
    pyclass.def("__mul__", &__mul__<X,X>, pybind11::is_operator()); if constexpr(CanDivide<X,X>) pyclass.def(__py_div__, &__div__<X,X>, pybind11::is_operator());
    pyclass.def("__radd__", &__radd__<X,X>, pybind11::is_operator()); pyclass.def("__rsub__", &__rsub__<X,X>, pybind11::is_operator());
    pyclass.def("__rmul__", &__rmul__<X,X>, pybind11::is_operator()); if constexpr(CanDivide<X,X>) pyclass.def(__py_rdiv__, &__rdiv__<X,X>, pybind11::is_operator()); return pyclass;
}
template<class X> pybind11::class_<X>& define_self_arithmetic(pybind11::module& module, pybind11::class_<X>& pyclass) { return define_arithmetic(module,pyclass); }
template<class X, class Y> pybind11::class_<X>& define_mixed_arithmetic(pybind11::module&, pybind11::class_<X>& pyclass, Tag<Y> = Tag<Y>()) {
    pyclass.def("__add__", &__add__<X,Y>, pybind11::is_operator()); pyclass.def("__radd__", &__radd__<X,Y>, pybind11::is_operator());
    pyclass.def("__sub__", &__sub__<X,Y>, pybind11::is_operator()); pyclass.def("__rsub__", &__rsub__<X,Y>, pybind11::is_operator());
    pyclass.def("__mul__", &__mul__<X,Y>, pybind11::is_operator()); pyclass.def("__rmul__", &__rmul__<X,Y>, pybind11::is_operator());
    if constexpr(CanDivide<X,Y>) pyclass.def(__py_div__, &__div__<X,Y>, pybind11::is_operator());
    if constexpr(CanDivide<Y,X>) pyclass.def(__py_rdiv__, &__rdiv__<X,Y>, pybind11::is_operator()); return pyclass;
}
template<class X> pybind11::class_<X>& define_transcendental(pybind11::module& module, pybind11::class_<X>& pyclass) {
    module.def("nul", &_nul_<X>); module.def("pos", &_pos_<X>); module.def("neg", &_neg_<X>); module.def("sqr", &_sqr_<X>);
    module.def("hlf", &_hlf_<X>); module.def("rec", &_rec_<X>); module.def("pow", &_pow_<X,Int>); module.def("sqrt", &_sqrt_<X>);
    module.def("exp", &_exp_<X>); module.def("log", &_log_<X>); module.def("sin", &_sin_<X>); module.def("cos", &_cos_<X>);
    module.def("tan", &_tan_<X>); module.def("asin", &_asin_<X>); module.def("acos", &_acos_<X>); module.def("atan", &_atan_<X>); return pyclass;
}
template<class X> pybind11::class_<X>& define_elementary(pybind11::module& module, pybind11::class_<X>& pyclass) { define_arithmetic(module,pyclass); define_transcendental(module,pyclass); return pyclass; }
template<class X> pybind11::class_<X>& define_monotonic(pybind11::module& module, pybind11::class_<X>& pyclass) {
    using NX=decltype(-declval<X>()); pyclass.def("__pos__", &__pos__<X>, pybind11::is_operator()); pyclass.def("__neg__", &__neg__<X>, pybind11::is_operator());
    pyclass.def("__add__", &__add__<X,X>, pybind11::is_operator()); pyclass.def("__sub__", &__sub__<X,NX>, pybind11::is_operator());
    module.def("sqrt", &_sqrt_<X>); module.def("exp", &_exp_<X>); module.def("log", &_log_<X>); module.def("atan", &_atan_<X>); return pyclass;
}
template<class X, class Y> pybind11::class_<X>& define_mixed_monotonic(pybind11::module&, pybind11::class_<X>& pyclass, Tag<Y> = Tag<Y>()) {
    using NY=decltype(-declval<Y>()); pyclass.def("__add__", &__add__<X,Y>, pybind11::is_operator()); pyclass.def("__radd__", &__radd__<X,Y>, pybind11::is_operator());
    pyclass.def("__sub__", &__sub__<X,NY>, pybind11::is_operator()); pyclass.def("__rsub__", &__rsub__<X,NY>, pybind11::is_operator()); return pyclass;
}

} // namespace Ariadne

#endif
