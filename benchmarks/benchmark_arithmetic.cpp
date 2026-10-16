/***************************************************************************
 *            benchmark_arithmetic.cpp
 *
 *  Copyright 2008--17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <vector>

#include "utility/stopwatch.hpp"
#include "numeric/rounding.hpp"
#include "logging/logging.hpp"
#include "logging/command_line.hpp"
#include "benchmark_command_line.hpp"
#include <stdexcept>
#include <string>

#include <gmp.h>

using namespace Ariadne;

// Machine epsilon, approximately 2.2e-16;
const double eps=1./(1<<26)/(1<<26);


double rndm() {
    double w=static_cast<double>(1u<<16)*(1u<<15);
    Nat r1=2u*static_cast<Nat>(std::rand())/2u;
    Nat r2=static_cast<Nat>(std::rand());
    Nat r3=static_cast<Nat>(std::rand());
    double r=((static_cast<double>(r3)/w+static_cast<double>(r2))/w+static_cast<double>(r1))/(1u<<28);
    //std::cerr<<"r="<<r<<" r1,2,3="<<r1<<","<<r2<<","<<r3<<"\n";
    return r;
}

Void dot_lu_rat(double& l, double& u, SizeType n, const double* x, const double* y);
Void dot_md_rat(mpq_t m, SizeType n, const double* x, const double* y);
Void dot_lu_std(double& l, double& u, SizeType n, const double* x, const double* y);
Void dot_lu_opp(double& l, double& u, SizeType n, const double* x, const double* y);
Void dot_lu_ivl(double& l, double& u, SizeType n, const double* x, const double* y);
Void dot_mr_std(double& m, double& r, SizeType n, const double* x, const double* y);
Void dot_mr_mid(double& m, double& r, SizeType n, const double* x, const double* y);
Void dot_mr_csy(double& m, double& r, SizeType n, const double* x, const double* y);

Void add_mr_std(double& e, SizeType n, double* r, const double* x, const double* y);
Void add_mr_buf(double& e, SizeType n, [[maybe_unused]] double* r, const double* x, const double* y);
Void add_mr_csy(double& e, SizeType n, double* r, const double* x, const double* y);

Void scal_mr_std(double& e, SizeType n, double* r, const double* x, const double& c);
Void scal_mr_csy(double& e, SizeType n, double* r, const double* x, const double& c);

double rational_error(double x, mpq_srcptr exact) {
    mpq_t error;
    mpq_init(error);
    mpq_set_d(error,x);
    mpq_sub(error,error,exact);
    if(mpq_sgn(error)<0) {
        mpq_neg(error,error);
    }
    const double result=mpq_get_d(error);
    mpq_clear(error);
    return result;
}

Bool rational_contains(double x, double radius, mpq_srcptr exact) {
    mpq_t error;
    mpq_t bound;
    mpq_inits(error,bound,nullptr);

    mpq_set_d(error,x);
    mpq_sub(error,error,exact);
    if(mpq_sgn(error)<0) {
        mpq_neg(error,error);
    }
    mpq_set_d(bound,radius);

    const Bool result=(mpq_cmp(error,bound)<=0);
    mpq_clears(error,bound,nullptr);
    return result;
}


Void benchmark_dot(SizeType n, SizeType nn) {
    std::vector<double> x(n);
    std::vector<double> y(n);

    for(SizeType i=0; i!=n; ++i) {
        x[i]=rndm();
        y[i]=rndm();
    }

    double m=0,r=0,l=0,u=0;

    Stopwatch<Microseconds> sw; double t=0;

    double ql=0; double qu=0;
    dot_lu_rat(ql,qu,n,x.data(),y.data());
    LOGGING_PRINTLN("lu_rat:        e="<<(qu-ql)/2<<" l="<<ql<<" u="<<qu)
    assert(ql<=qu);

    mpq_t qm;
    mpq_init(qm);
    mpq_set_ui(qm,0u,1u);
    dot_md_rat(qm,n,x.data(),y.data());

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        l=0; u=0;
        dot_lu_std(l,u,n,x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    LOGGING_PRINTLN("lu_std: t(us)="<<std::setprecision(5)<<t
             <<std::setprecision(20)<<" e="<<(u-l)/2<<" l="<<l<<" u="<<u)
    assert(l<=ql && qu<=u);

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        l=0; u=0;
        dot_lu_opp(l,u,n,x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    LOGGING_PRINTLN("lu_opp: t(us)="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" e="<<(u-l)/2<<" l="<<l<<" u="<<u)
    assert(l<=ql && qu<=u);

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        l=0; u=0;
        dot_lu_ivl(l,u,n,x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    LOGGING_PRINTLN("lu_ivl: t(us)="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" e="<<(u-l)/2<<" l="<<l<<" u="<<u)
    assert(l<=ql && qu<=u);



    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        m=0; r=0;
        dot_mr_std(m,r,n,x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    LOGGING_PRINTLN("mr_std: t(us)="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r<<" m="<<m<<" e="<<rational_error(m,qm))
    assert(rational_contains(m,r,qm));

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        m=0; r=0;
        dot_mr_mid(m,r,n,x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    LOGGING_PRINTLN("mr_mid: t(us)="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r<<" m="<<m<<" e="<<rational_error(m,qm))
    assert(rational_contains(m,r,qm));

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        m=0; r=0;
        dot_mr_csy(m,r,n,x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    LOGGING_PRINTLN("mr_csy: t(us)="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r<<" m="<<m<<" e="<<rational_error(m,qm))
    assert(rational_contains(m,r,qm));

    mpq_clear(qm);
}


Void benchmark_add(SizeType n, SizeType nn) {
    std::vector<double> x(n);
    std::vector<double> y(n);
    std::vector<double> z(n);

    for(SizeType i=0; i!=n; ++i) {
        x[i]=rndm();
        y[i]=rndm();
    }

    double r=0;

    Stopwatch<Microseconds> sw; double t=0;

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        r=0;
        add_mr_std(r,n,z.data(),x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    LOGGING_PRINTLN("add_mr_std: t="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r)

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        r=0;
        add_mr_buf(r,n,z.data(),x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    LOGGING_PRINTLN("add_mr_buf: t="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r)

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        r=0;
        add_mr_csy(r,n,z.data(),x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    LOGGING_PRINTLN("add_mr_csy: t="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r)
}


 Void benchmark_scal(SizeType n, SizeType nn) {
    std::vector<double> x(n);
    std::vector<double> z(n);

    double c=rndm();
    for(SizeType i=0; i!=n; ++i) {
        x[i]=rndm();
    }


    double r=0;

    Stopwatch<Microseconds> sw; double t=0;

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        r=0;
        scal_mr_std(r,n,z.data(),x.data(),c);
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    LOGGING_PRINTLN("scal_mr_std: t="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r)

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        r=0;
        scal_mr_csy(r,n,z.data(),x.data(),c);
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    LOGGING_PRINTLN("scal_mr_csy: t="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r)
}



struct ArithmeticBenchmarkSelection {
    SizeType vector_exponent=10u;
    SizeType repetition_exponent=12u;
    SizeType runs=1u;
};

List<ArgumentParser> make_arithmetic_argument_parsers(ArithmeticBenchmarkSelection& selection) {
    List<ArgumentParser> parsers;
    parsers.append(make_valued_argument_parser(
        "x","vector-exponent","Log2 of vector size (nonnegative)",
        [&selection](String const& value) -> VoidFunction {
            const SizeType exponent=BenchmarkCommandLine::parse_unsigned<SizeType>(value);
            if(exponent>=std::numeric_limits<SizeType>::digits) {
                throw std::invalid_argument("vector exponent out of range");
            }
            return [&selection,exponent] { selection.vector_exponent=exponent; };
        }));

    parsers.append(make_valued_argument_parser(
        "e","repetition-exponent","Log2 of repetition count (nonnegative)",
        [&selection](String const& value) -> VoidFunction {
            const SizeType exponent=BenchmarkCommandLine::parse_unsigned<SizeType>(value);
            if(exponent>=std::numeric_limits<SizeType>::digits) {
                throw std::invalid_argument("repetition exponent out of range");
            }
            return [&selection,exponent] { selection.repetition_exponent=exponent; };
        }));

    parsers.append(make_valued_argument_parser(
        "r","runs","Number of benchmark runs (positive)",
        [&selection](String const& value) -> VoidFunction {
            const SizeType runs=BenchmarkCommandLine::parse_unsigned<SizeType>(value);
            if(runs==0u) throw std::invalid_argument("runs must be positive");
            return [&selection,runs] { selection.runs=runs; };
        }));
    return parsers;
}

Int main(Int argc, const char* argv[]) {
    BenchmarkCommandLine::initialize_logging();
    try {
        ArithmeticBenchmarkSelection selection;
        CommandLineInterface command_line;
        command_line.add(logging_command_line());
        command_line.add(CommandLineContribution(
            "numeric-arithmetic-benchmark",
            make_arithmetic_argument_parsers(selection)));

        const auto arguments=BenchmarkCommandLine::normalize_legacy_arguments(
            argc,argv,{"vector-exponent","repetition-exponent","runs"},
            {"--vector-exponent","-x","--repetition-exponent","-e","--runs","-r",
             "--scheduler","-s","--theme","-t","--verbosity","-v"});
        if (!command_line.acquire(arguments)) {
            return BenchmarkCommandLine::requests_help(arguments) ? 0 : 1;
        }

        const SizeType vector_exponent=selection.vector_exponent;
        const SizeType repetition_exponent=selection.repetition_exponent;
        const SizeType runs=selection.runs;

        const SizeType vector_size=SizeType{1} << vector_exponent;
        const SizeType repetitions=SizeType{1} << repetition_exponent;
        LOGGING_PRINTLN(std::string(78u,'='))
        LOGGING_PRINTLN("Arithmetic benchmark")
        LOGGING_PRINTLN("Vector size: " << vector_size
                        << "  |  Repetitions: " << repetitions
                        << "  |  Runs: " << runs)
        LOGGING_PRINTLN(std::string(78u,'='))

        for(SizeType run=1u; run<=runs; ++run) {
            LOGGING_PRINTLN(" ")
            LOGGING_PRINTLN("Run " << run << "/" << runs)
            LOGGING_PRINTLN(std::string(78u,'-'))

            Stopwatch<Milliseconds> sw;
            LOGGING_PRINTLN("Scalar operations")
            benchmark_scal(vector_size,repetitions);
            LOGGING_PRINTLN(std::string(78u,'-'))
            LOGGING_PRINTLN("Vector addition")
            benchmark_add(vector_size,repetitions);
            LOGGING_PRINTLN(std::string(78u,'-'))
            LOGGING_PRINTLN("Dot products")
            benchmark_dot(vector_size,repetitions);

            sw.click();
            LOGGING_PRINTLN(std::string(78u,'-'))
            LOGGING_PRINTLN("Run total: " << sw.duration().count() << " ms")
            LOGGING_PRINTLN(std::string(78u,'='))
        }
        return 0;
    } catch(std::exception const& e) {
        LOGGING_PRINTLN("benchmark_arithmetic: " << e.what())
        return 1;
    }
}

Void dot_lu_rat(double& l, double& u, SizeType n, const double* x, const double* y) {
    set_builtin_rounding_to_nearest();
    assert(l==u);

    mpq_t exact;
    mpq_t xq;
    mpq_t yq;
    mpq_t product;
    mpq_t candidate;
    mpq_inits(exact,xq,yq,product,candidate,nullptr);

    mpq_set_d(exact,l);
    for(SizeType i=0; i!=n; ++i) {
        mpq_set_d(xq,x[i]);
        mpq_set_d(yq,y[i]);
        mpq_mul(product,xq,yq);
        mpq_add(exact,exact,product);
    }

    const double d=mpq_get_d(exact);
    const double e=std::abs(d)*(eps/16);
    assert(e>0);

    l=d;
    Int m=0;
    mpq_set_d(candidate,l);
    while(mpq_cmp(candidate,exact)<0) {
        ++m;
        l=d+static_cast<double>(m)*e;
        mpq_set_d(candidate,l);
    }
    while(mpq_cmp(candidate,exact)>0) {
        --m;
        l=d+static_cast<double>(m)*e;
        mpq_set_d(candidate,l);
    }
    ++m;
    u=d+static_cast<double>(m)*e;

    mpq_set_d(candidate,l);
    assert(mpq_cmp(candidate,exact)<=0);
    mpq_set_d(candidate,u);
    assert(mpq_cmp(candidate,exact)>=0);

    mpq_clears(exact,xq,yq,product,candidate,nullptr);
}


Void dot_md_rat(mpq_t m, SizeType n, const double* x, const double* y) {
    set_builtin_rounding_to_nearest();

    mpq_t xq;
    mpq_t yq;
    mpq_t product;
    mpq_inits(xq,yq,product,nullptr);

    for(SizeType i=0; i!=n; ++i) {
        mpq_set_d(xq,x[i]);
        mpq_set_d(yq,y[i]);
        mpq_mul(product,xq,yq);
        mpq_add(m,m,product);
    }

    mpq_clears(xq,yq,product,nullptr);
}

Void dot_lu_ivl(double& l, double& u, SizeType n, const double* x, const double* y) {
    for(SizeType i=0; i!=n; ++i) {
        set_builtin_rounding_upward();
        u=u+x[i]*y[i];
        set_builtin_rounding_downward();
        l=l+x[i]*y[i];
    }
    set_builtin_rounding_to_nearest();
}

Void dot_lu_std(double& l, double& u, SizeType n, const double* x, const double* y) {
    set_builtin_rounding_upward();
    for(SizeType i=0; i!=n; ++i) {
        u=u+x[i]*y[i];
        //std::cerr<<" i="<<i<<" u="<<u<<"\n";
    }
    set_builtin_rounding_downward();
    for(SizeType i=0; i!=n; ++i) {
        l=l+x[i]*y[i];
        //std::cerr<<" i="<<i<<" l="<<l<<"\n";
    }
    set_builtin_rounding_to_nearest();
}

Void dot_lu_opp(double& l, double& u, SizeType n, const double* x, const double* y) {
    set_builtin_rounding_upward();
    volatile double uu=u;
    volatile double ll=-l;
    volatile double t;
    for(SizeType i=0; i!=n; ++i) {
        uu=uu+x[i]*y[i];
        t=-x[i];
        t=t*y[i];
        ll=ll+t;
        //l+=t*y[i];
        //std::cerr<<" i="<<i<<" l="<<(-l)<<" u="<<u<<"\n";
    }
    u=uu;
    l=-ll;
    set_builtin_rounding_to_nearest();
}

Void dot_mr_std(double& m, double& r, SizeType n, const double* x, const double* y) {
    set_builtin_rounding_downward();
    volatile double l=m-r;
    for(SizeType i=0; i!=n; ++i) {
        l=l+x[i]*y[i];
    }
    set_builtin_rounding_upward();
    volatile double u=m+r;
    for(SizeType i=0; i!=n; ++i) {
        u=u+x[i]*y[i];
    }
    set_builtin_rounding_to_nearest();
    m=(u+l)/2;
    set_builtin_rounding_upward();
    r=std::max(u-m,m-l);
    set_builtin_rounding_to_nearest();
}

Void dot_mr_mid(double& m, double& r, SizeType n, const double* x, const double* y) {
    volatile double a=m;
    for(SizeType i=0; i!=n; ++i) {
        a=a+x[i]*y[i];
    }
    set_builtin_rounding_downward();
    volatile double l=m-r;
    for(SizeType i=0; i!=n; ++i) {
        l=l+x[i]*y[i];
    }
    set_builtin_rounding_upward();
    volatile double u=m+r;
    for(SizeType i=0; i!=n; ++i) {
        u=u+x[i]*y[i];
    }
    m=a;
    r=r+std::max(u-m,m-l);
    set_builtin_rounding_to_nearest();
}

Void dot_mr_csy(double& m, double& r, SizeType n, const double* x, const double* y) {
    double e=0;
    for(SizeType i=0; i!=n; ++i) {
        double p=x[i]*y[i];
        m+=p;
        e+=std::abs(p)*(eps/2);
        e+=std::abs(m)*(eps/2);
    }
    r+=e;
}




Void add_mr_std(double& e, SizeType n, double* r, const double* x, const double* y)
{
    for(SizeType i=0; i!=n; ++i) {
        volatile double a=x[i]+y[i];
        set_builtin_rounding_downward();
        volatile double l=x[i]+y[i];
        set_builtin_rounding_upward();
        volatile double u=x[i]+y[i];
        r[i]=a;
        e+=(u-l)/2;
        set_builtin_rounding_to_nearest();
    }
}

Void add_mr_buf(double& e, SizeType n, [[maybe_unused]] double* r, const double* x, const double* y)
{
    std::vector<double> z(n);
    for(SizeType i=0; i!=n; ++i) {
        z[i]=x[i]+y[i];
    }
    set_builtin_rounding_upward();
    for(SizeType i=0; i!=n; ++i) {
        volatile double u=x[i]+y[i];
        volatile double t=-x[i];
        volatile double ml=t-y[i];
        e+=(u+ml)/2;
    }
    set_builtin_rounding_to_nearest();
}

Void add_mr_csy(double& e, SizeType n, double* r, const double* x, const double* y)
{
    double d=0;
    for(SizeType i=0; i!=n; ++i) {
        r[i]=x[i]+y[i];
        d+=std::abs(r[i]);
    }
    d=d*(1+eps*n/2)*eps/2;
    e=e+d;
    e=e*(1+eps/2);
}


Void scal_mr_std(double& e, SizeType n, double* r, const double* x, const double& c)
{
    for(SizeType i=0; i!=n; ++i) {
        volatile double a=x[i]*c;
        set_builtin_rounding_downward();
        volatile double l=x[i]*c;
        set_builtin_rounding_upward();
        volatile double u=x[i]*c;
        r[i]=a;
        e+=(u-l)/2;
        set_builtin_rounding_to_nearest();
    }
}

Void scal_mr_csy(double& e, SizeType n, double* r, const double* x, const double& c)
{
    double d=0;
    for(SizeType i=0; i!=n; ++i) {
        r[i]=x[i]*c;
        d+=std::abs(r[i]);
    }
    d=d*(1+eps*n/2)*eps/2;
    e=e+d;
    e=e*(1+eps/2);
}
