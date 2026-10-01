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

#include <gmpxx.h>

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
Void dot_md_rat(mpq_class& m, SizeType n, const double* x, const double* y);
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
    std::cout<<"lu_rat:        e="<<(qu-ql)/2<<" l="<<ql<<" u="<<qu<<std::endl;
    assert(ql<=qu);

    mpq_class qm=0;
    dot_md_rat(qm,n,x.data(),y.data());

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        l=0; u=0;
        dot_lu_std(l,u,n,x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    std::cout<<"lu_std: t(us)="<<std::setprecision(5)<<t
             <<std::setprecision(20)<<" e="<<(u-l)/2<<" l="<<l<<" u="<<u<<std::endl;
    assert(l<=ql && qu<=u);

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        l=0; u=0;
        dot_lu_opp(l,u,n,x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    std::cout<<"lu_opp: t(us)="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" e="<<(u-l)/2<<" l="<<l<<" u="<<u<<std::endl;
    assert(l<=ql && qu<=u);

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        l=0; u=0;
        dot_lu_ivl(l,u,n,x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    std::cout<<"lu_ivl: t(us)="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" e="<<(u-l)/2<<" l="<<l<<" u="<<u<<std::endl;
    assert(l<=ql && qu<=u);



    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        m=0; r=0;
        dot_mr_std(m,r,n,x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    std::cout<<"mr_std: t(us)="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r<<" m="<<m<<" e="<<mpq_class(abs(mpq_class(m)-qm)).get_d()<<std::endl;
    assert(abs(mpq_class(m)-qm)<=r);

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        m=0; r=0;
        dot_mr_mid(m,r,n,x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    std::cout<<"mr_mid: t(us)="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r<<" m="<<m<<" e="<<mpq_class(abs(mpq_class(m)-qm)).get_d()<<std::endl;
    assert(abs(mpq_class(m)-qm)<=r);

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        m=0; r=0;
        dot_mr_csy(m,r,n,x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    std::cout<<"mr_csy: t(us)="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r<<" m="<<m<<" e="<<mpq_class(abs(mpq_class(m)-qm)).get_d()<<std::endl;
    assert(abs(mpq_class(m)-qm)<=r);
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
    std::cout<<"add_mr_std: t="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r<<std::endl;

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        r=0;
        add_mr_buf(r,n,z.data(),x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    std::cout<<"add_mr_buf: t="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r<<std::endl;

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        r=0;
        add_mr_csy(r,n,z.data(),x.data(),y.data());
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    std::cout<<"add_mr_csy: t="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r<<std::endl;
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
    std::cout<<"scal_mr_std: t="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r<<std::endl;

    sw.restart();
    for(SizeType i=0; i!=nn; ++i) {
        r=0;
        scal_mr_csy(r,n,z.data(),x.data(),c);
    }
    sw.click();
    t=static_cast<double>(sw.duration().count());
    std::cout<<"scal_mr_csy: t="<<std::setprecision(5)<<t<<std::setprecision(20)
             <<" r="<<r<<std::endl;
}


Int main(Int argc, const char* argv[]) {
    std::cout << std::setprecision(20);
    std::cerr << std::setprecision(20);

    SizeType vector_size = (SizeType{1} << 10u);
    SizeType repetitions = (SizeType{1} << 12u);
    SizeType runs = 1u;

    if (argc > 1) {
        const auto exponent = std::strtoul(argv[1], nullptr, 10);
        if (exponent >= static_cast<unsigned long>(std::numeric_limits<SizeType>::digits)) {
            std::cerr << "Vector-size exponent is too large" << std::endl;
            return 1;
        }
        vector_size = (SizeType{1} << exponent);
    }
    if (argc > 2) {
        const auto exponent = std::strtoul(argv[2], nullptr, 10);
        if (exponent >= static_cast<unsigned long>(std::numeric_limits<SizeType>::digits)) {
            std::cerr << "Repetition exponent is too large" << std::endl;
            return 1;
        }
        repetitions = (SizeType{1} << exponent);
    }
    if (argc > 3) {
        runs = static_cast<SizeType>(std::strtoul(argv[3], nullptr, 10));
    }

    std::cout << "\nVector size=" << vector_size
              << "\nRepetitions=" << repetitions
              << "\nRuns=" << runs << "\n" << std::endl;

    for (SizeType run = 1u; run <= runs; ++run) {
        std::cout << "=== Run " << run << "/" << runs << " ===" << std::endl;

        Stopwatch<Milliseconds> sw;

        benchmark_scal(vector_size, repetitions);
        benchmark_add(vector_size, repetitions);
        benchmark_dot(vector_size, repetitions);

        sw.click();
        std::cout << "run_total: t(ms)=" << sw.duration().count() << "\n" << std::endl;
    }

    return 0;
}


Void dot_lu_rat(double& l, double& u, SizeType n, const double* x, const double* y) {
    set_builtin_rounding_to_nearest();
    assert(l==u);
    mpq_class r=l;
    for(SizeType i=0; i!=n; ++i) {
        r+=mpq_class(x[i])*mpq_class(y[i]);
    }
    //mpf_class::set_precision(128);
    //std::cerr<<"r="<<r<<"~"<<mpf_class(r)<<"\n";
    mpq_class a(r.get_d());
    double d=a.get_d();
    assert(mpq_class(d)==a);
    assert(d==a);
    double e=std::abs(d)*(eps/16);
    assert(e>0);
    l=d;
    Int m=0;
    while(l<r) {
        ++m;
        l=d+m*e;
    }
    while(l>r) {
        --m;
        l=d+m*e;
    }
    ++m;
    u=d+m*e;

    assert(mpq_class(l)<=r);
    assert(mpq_class(u)>=r);
}


Void dot_md_rat(mpq_class& m, SizeType n, const double* x, const double* y) {
    set_builtin_rounding_to_nearest();
    for(SizeType i=0; i!=n; ++i) {
        m+=mpq_class(x[i])*mpq_class(y[i]);
    }
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
