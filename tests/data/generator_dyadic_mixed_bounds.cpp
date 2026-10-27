// SPDX-License-Identifier: MIT
// Bounded exact rational interval arithmetic and three-valued comparisons.
#include "oracle_exact_integer.hpp"
#include <numeric>
namespace {
using I=std::int64_t;
struct Q {
    I n,d;
    Q(I p,I q=1):n(p),d(q) {
        if (d==0) throw std::invalid_argument("zero denominator");
        if (d<0) { n=-n; d=-d; }
        const auto g=std::gcd(n,d); n/=g; d/=g;
    }
};
bool less(Q a,Q b) { return a.n*b.d<b.n*a.d; }
Q add(Q a,Q b) { return {a.n*b.d+b.n*a.d,a.d*b.d}; }
Q sub(Q a,Q b) { return {a.n*b.d-b.n*a.d,a.d*b.d}; }
Q mul(Q a,Q b) { return {a.n*b.n,a.d*b.d}; }
Q divide(Q a,Q b) { return {a.n*b.d,a.d*b.n}; }
Q minimum(Q a,Q b) { return less(a,b)?a:b; }
Q maximum(Q a,Q b) { return less(a,b)?b:a; }
struct Interval { Q l,u; };
bool contains_zero(Interval a) { return !less(a.u,Q(0)) && !less(Q(0),a.l); }
Interval arithmetic(std::string const& op,Interval a,Interval b) {
    if (op=="add") return {add(a.l,b.l),add(a.u,b.u)};
    if (op=="sub") return {sub(a.l,b.u),sub(a.u,b.l)};
    if (op=="max") return {maximum(a.l,b.l),maximum(a.u,b.u)};
    if (op=="min") return {minimum(a.l,b.l),minimum(a.u,b.u)};
    std::array<Q,4> corners={op=="mul"?mul(a.l,b.l):divide(a.l,b.l),
        op=="mul"?mul(a.l,b.u):divide(a.l,b.u),
        op=="mul"?mul(a.u,b.l):divide(a.u,b.l),
        op=="mul"?mul(a.u,b.u):divide(a.u,b.u)};
    Q lower=corners.front(),upper=lower;
    for (Q x:corners) { lower=minimum(lower,x); upper=maximum(upper,x); }
    return {lower,upper};
}
int compare(std::string const& op,Interval a,Interval b) {
    const bool equal=!less(a.l,a.u) && !less(b.l,b.u) && !less(a.l,b.l) && !less(b.l,a.l);
    const bool disjoint=less(a.u,b.l)||less(b.u,a.l);
    if (op=="eq") return equal?1:disjoint?-1:0;
    if (op=="ne") return equal?-1:disjoint?1:0;
    if (op=="gt") return compare("lt",b,a);
    if (op=="ge") return compare("le",b,a);
    if (op=="lt") return less(a.u,b.l)?1:!less(a.l,b.u)?-1:0;
    return !less(b.l,a.u)?1:less(b.u,a.l)?-1:0;
}
std::pair<std::string,std::string> fixtures() {
    const std::array<std::pair<I,I>,5> dyadic={{{-3,-1},{-1,3},{1,3},{0,0},{2,2}}};
    const std::array<std::pair<I,I>,5> rational={{{-2,-1},{-1,2},{1,2},{3,3},{6,9}}};
    std::string values="# Independent exact mixed interval arithmetic.\n# id op reverse dyadic_low_n dyadic_high_n dyadic_den rational_low_n rational_high_n rational_den expected_low_n expected_low_d expected_high_n expected_high_d\n";
    std::string logic="# Independent three-valued mixed interval comparisons (1 true, -1 false, 0 indeterminate).\n# id op reverse dyadic_low_n dyadic_high_n dyadic_den rational_low_n rational_high_n rational_den expected\n";
    unsigned vi=0u,li=0u;
    for (auto [dl,du]:dyadic) for (auto [ql,qu]:rational) for (bool reverse:{false,true}) {
        const Interval d{{dl,2},{du,2}},q{{ql,3},{qu,3}};
        const auto a=reverse?q:d,b=reverse?d:q;
        const auto key=[&](unsigned id,std::string const& op) {
            return std::vector<std::string>{std::to_string(id),op,reverse?"1":"0",
                std::to_string(dl),std::to_string(du),"2",std::to_string(ql),std::to_string(qu),"3"};
        };
        for (std::string const op:{"add","sub","mul","div","max","min"}) {
            if (op=="div" && contains_zero(b)) continue;
            const auto expected=arithmetic(op,a,b);
            auto fields=key(vi++,op);
            for (I n:{expected.l.n,expected.l.d,expected.u.n,expected.u.d}) fields.push_back(std::to_string(n));
            values+=independent_oracle::row(fields);
        }
        for (std::string const op:{"eq","ne","lt","le","gt","ge"}) {
            auto fields=key(li++,op); fields.push_back(std::to_string(compare(op,a,b)));
            logic+=independent_oracle::row(fields);
        }
    }
    return {values,logic};
}
std::string powers() {
    std::string result="# Independent interval powers on dyadic endpoints.\n# id low_n high_n power expected_low_n expected_low_d expected_high_n expected_high_d\n";
    unsigned id=0u;
    for (auto [ln,un]:std::array<std::pair<I,I>,5>{{{-3,-1},{-1,3},{1,3},{0,0},{2,2}}}) {
        for (unsigned degree:{0u,1u,2u,3u,4u}) {
            auto power=[&](I n) { Q q(1); for (unsigned i=0u;i<degree;++i) q=mul(q,Q(n,2)); return q; };
            const Q left=power(ln),right=power(un);
            const Q lower=degree!=0u && degree%2u==0u && ln<0 && un>0?Q(0):minimum(left,right);
            const Q upper=maximum(left,right);
            result+=independent_oracle::row({std::to_string(id++),std::to_string(ln),std::to_string(un),
                std::to_string(degree),std::to_string(lower.n),std::to_string(lower.d),
                std::to_string(upper.n),std::to_string(upper.d)});
        }
    }
    return result;
}

}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        const auto [values,logic]=fixtures();
        independent_oracle::output(dir,"dyadic_mixed_bounds.tsv",values,write);
        independent_oracle::output(dir,"dyadic_mixed_comparisons.tsv",logic,write);
        independent_oracle::output(dir,"dyadic_bounds_power.tsv",powers(),write);
    } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
