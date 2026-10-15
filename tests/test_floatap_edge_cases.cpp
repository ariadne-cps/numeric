// FloatAP special values, rounding ties and object semantics.
#include "numeric/floatap.hpp"
#include <gmp.h>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
using K=Ariadne::FloatAPKernel;
using R=K::Rounding;
namespace {
struct Z{mpz_t z;explicit Z(long n){mpz_init_set_si(z,n);}~Z(){mpz_clear(z);}};
struct Q{mpq_t q;Q(){mpq_init(q);}~Q(){mpq_clear(q);}};
void check(bool x,char const* msg){if(!x)throw std::runtime_error(msg);}
template<class E,class F>void raises(F f,char const* m){bool ok=false;try{f();}catch(E const&){ok=true;}check(ok,m);}
K v(long a,long b=1,unsigned long p=16,R r=R::nearest){Z n(a),d(b);return K::from_ratio(n.z,d.z,0,p,r);}
void same(K const& a,K const& b,char const* m){
check(a.kind()==b.kind(),m);
if(a.is_nan())return;
check(a.negative()==b.negative(),m);
if(!a.is_finite())return;
Q x,y;a.to_rational(x.q);b.to_rational(y.q);check(mpq_cmp(x.q,y.q)==0,m);
}
void conversions(){
raises<std::invalid_argument>([]{K x(1);},"precision");
raises<std::domain_error>([]{v(1,0);},"denominator");
same(v(2,-3),v(-2,3),"signed denominator");
same(v(-2,-3),v(2,3),"both negative");
Z i(-129);same(K::from_integer(i.z,7,R::downward),v(-129,1,7,R::downward),"integer");
same(v(5,4,2),v(1,1,2),"nearest even low");
same(v(7,4,2),v(2,1,2),"nearest even high");
same(v(5,4,2,R::upward),v(3,2,2),"round up");
same(v(-5,4,2,R::downward),v(-3,2,2),"negative downward");
same(v(-5,4,2,R::upward),v(-1,1,2),"negative upward");
same(v(-5,4,2,R::toward_zero),v(-1,1,2),"negative toward zero");
raises<std::invalid_argument>([]{v(5,3,7,static_cast<R>(99));},"invalid rounding");
raises<std::overflow_error>([]{Z n(7),d(1);K::from_ratio(n.z,d.z,std::numeric_limits<K::Exponent>::max(),2,R::nearest);},"overflow");
auto x=v(-11,8,8);check(x.precision()==8&&x.exponent()<0&&mpz_sgn(x.significand())>0,"accessors");
K y=x;same(y,x,"copy");
K z=K::zero(32);z=x;same(z,x,"assignment");
K moved(std::move(y));same(moved,x,"move");
K other=v(3);z.swap(other);same(other,x,"swap");
K assigned=K::zero(32);assigned=std::move(moved);same(assigned,x,"move assign");
check(K::negate(K::nan(8)).is_nan(),"negate nan");
check(K::negate(K::zero(8)).negative(),"negate zero");
raises<std::domain_error>([]{Q q;K::infinity(16).to_rational(q.q);},"infinite rational");
raises<std::domain_error>([]{Q q;K::nan(16).to_rational(q.q);},"nan rational");
}
void specials(){
auto p=K::zero(16),n=K::zero(16,true),x=v(3),neg=v(-3);
auto inf=K::infinity(16),minus_inf=K::infinity(16,true),nan=K::nan(16);
check(p.is_finite()&&!p.is_inf()&&!p.is_nan(),"finite classification");
check(inf.is_inf()&&!inf.is_zero(),"inf classification");
check(nan.is_nan()&&!nan.is_finite(),"nan classification");
check(K::add(nan,x,16,R::nearest).is_nan(),"add lhs nan");
check(K::add(x,nan,16,R::nearest).is_nan(),"add rhs nan");
same(K::add(inf,x,16,R::nearest),inf,"add inf lhs");
same(K::add(x,minus_inf,16,R::nearest),minus_inf,"add inf rhs");
same(K::add(inf,inf,16,R::nearest),inf,"add inf inf");
check(K::add(inf,minus_inf,16,R::nearest).is_nan(),"infinity cancellation");
same(K::add(n,n,16,R::nearest),n,"negative zero");
same(K::add(p,n,16,R::downward),n,"downward zero");
same(K::add(x,neg,16,R::downward),n,"negative cancellation zero");
same(K::add(x,neg,16,R::upward),p,"positive cancellation zero");
check(K::mul(nan,x,16,R::nearest).is_nan(),"mul lhs nan");
check(K::mul(x,nan,16,R::nearest).is_nan(),"mul rhs nan");
check(K::mul(inf,p,16,R::nearest).is_nan(),"infinity times zero");
check(K::mul(p,inf,16,R::nearest).is_nan(),"zero times infinity");
same(K::mul(inf,neg,16,R::nearest),minus_inf,"mul inf");
same(K::mul(neg,inf,16,R::nearest),minus_inf,"mul inf rhs");
same(K::mul(n,x,16,R::nearest),n,"mul negative zero");
same(K::mul(x,n,16,R::nearest),n,"mul rhs zero");
check(K::div(nan,x,16,R::nearest).is_nan(),"div lhs nan");
check(K::div(x,nan,16,R::nearest).is_nan(),"div rhs nan");
check(K::div(p,n,16,R::nearest).is_nan(),"0/0");
check(K::div(inf,minus_inf,16,R::nearest).is_nan(),"inf/inf");
same(K::div(x,n,16,R::nearest),minus_inf,"finite / negative zero");
same(K::div(minus_inf,x,16,R::nearest),minus_inf,"infinity/finite");
same(K::div(n,x,16,R::nearest),n,"zero/finite");
same(K::div(x,minus_inf,16,R::nearest),n,"finite/infinity");
check(K::fma(nan,x,x,16,R::nearest).is_nan(),"fma lhs nan");
check(K::fma(x,nan,x,16,R::nearest).is_nan(),"fma rhs nan");
check(K::fma(x,x,nan,16,R::nearest).is_nan(),"fma addend nan");
check(K::fma(p,inf,x,16,R::nearest).is_nan(),"fma zero infinity");
check(K::fma(inf,p,x,16,R::nearest).is_nan(),"fma infinity zero");
check(K::fma(inf,x,minus_inf,16,R::nearest).is_nan(),"fma inf cancellation");
same(K::fma(inf,x,x,16,R::nearest),inf,"fma inf lhs");
same(K::fma(x,inf,x,16,R::nearest),inf,"fma inf rhs");
same(K::fma(x,x,minus_inf,16,R::nearest),minus_inf,"fma infinite addend");
same(K::fma(p,x,x,16,R::nearest),x,"fma zero lhs");
same(K::fma(x,p,x,16,R::nearest),x,"fma zero rhs");
same(K::fma(n,p,n,16,R::nearest),n,"fma signed zero");
}
}
int main(){try{conversions();specials();std::cout<<"FloatAP edge cases OK\n";return 0;}catch(std::exception const& e){std::cerr<<e.what()<<'\n';return 1;}}
