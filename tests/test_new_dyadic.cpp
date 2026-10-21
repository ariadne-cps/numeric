/***************************************************************************
 * test_new_dyadic.cpp — exact native Dyadic candidate regressions
 ****************************************************************************/
#include "numeric/new_dyadic.hpp"
#include <stdexcept>
#include <iostream>

using namespace Ariadne;

int main() {
    try {
        auto require=[](bool condition,char const* message) {
            if(!condition)throw std::runtime_error(message);
        };
        const NewDyadic zero(0),one(1),half(Integer(1),1u);
        const NewDyadic third_half(Integer(3),1u);
        require(zero.is_zero(),"zero");
        require(NewDyadic(Integer(12),3u)==third_half,"canonicalize factor of 2");
        require(NewDyadic(Integer(-12),3u)==-third_half,"negative canonicalization");
        require(half+half==one,"simple addition");
        require(half-half==zero,"exact subtraction");
        require(half*half==NewDyadic(Integer(1),2u),"exact multiplication");
        require(NewDyadic(3)/2==third_half,"divide by two");
        require(pow(half,3u)==NewDyadic(Integer(1),3u),"power");
        require(floor(NewDyadic(Integer(-5),1u))==Integer(-3),"negative floor");
        require(ceil(NewDyadic(Integer(-5),1u))==Integer(-2),"negative ceil");
        require(round(NewDyadic(Integer(-1),1u))==Integer(-1),"negative tie");
        require(round(NewDyadic(Integer(1),1u))==Integer(1),"positive tie");
        require(round(NewDyadic(Integer(1),2u))==Integer(0),"below half");
        require(round(NewDyadic(Integer(-3),2u))==Integer(-1),"above half");
        require(NewDyadic::inf()>one,"positive infinity ordering");
        require(NewDyadic::inf(-1)<one,"negative infinity ordering");
        require(NewDyadic::inf(-1)<NewDyadic::inf(),"infinity ordering");
        require(NewDyadic::nan()!=NewDyadic::nan(),"NaN unordered");
        require(NewDyadic::nan().kind()==NewDyadic::Kind::NaN,"NaN kind");
        for(unsigned bits:{53u,128u,256u,1024u,2048u,4096u}) {
            const auto high=shift_left(BigUnsigned(1u),bits-1u);
            const Integer a=Integer::from_magnitude(add(high,BigUnsigned(13u)),true);
            const Integer b=Integer::from_magnitude(add(high,BigUnsigned(31u)),false);
            const NewDyadic x(a,bits/5u+1u),y(b,bits/5u+2u);
            const Rational rx=x.to_rational(),ry=y.to_rational();
            require((x+y).to_rational()==rx+ry,"Rational addition oracle");
            require((x-y).to_rational()==rx-ry,"Rational subtraction oracle");
            require((x*y).to_rational()==rx*ry,"Rational multiplication oracle");
            require((x<y)==(rx<ry),"Rational comparison oracle");
            require(floor(x)==floor(rx),"Rational floor oracle");
            require(ceil(x)==ceil(rx),"Rational ceil oracle");
            require(round(x)==round(rx),"Rational round oracle");
            require(x.get_d()==rx.get_d(),"Rational binary64 oracle");
        }
        std::cout<<"NewDyadic exact regression PASS\n";
        return 0;
    } catch(std::exception const& e) {
        std::cerr<<"NewDyadic exact regression FAILED: "<<e.what()<<'\n';
        return 1;
    }
}
