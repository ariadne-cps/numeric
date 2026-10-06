/***************************************************************************
 *            test_operators.cpp
 *
 *  This file is part of Ariadne.
 ****************************************************************************/

#include <array>
#include <sstream>

#include "numeric/operators.hpp"
#include "numeric/integer.hpp"
#include "numeric/rational.hpp"
#include "numeric/dyadic.hpp"
#include "utility/test.hpp"

using namespace Ariadne;

namespace {

struct CodeName {
    OperatorCode code;
    const char* text;
};

struct CodeKind {
    OperatorCode code;
    OperatorKind kind;
};

template<class OP>
void test_metadata(OP op, OperatorCode code, OperatorKind kind) {
    Operator generic(op);
    ARIADNE_TEST_EQUALS(generic.code(),code);
    ARIADNE_TEST_EQUALS(generic.kind(),kind);
}

} // namespace

class TestOperators {
  public:
    void test() {
        ARIADNE_TEST_CALL(test_kinds());
        ARIADNE_TEST_CALL(test_names_and_symbols());
        ARIADNE_TEST_CALL(test_objects());
        ARIADNE_TEST_CALL(test_simple_operations());
        ARIADNE_TEST_CALL(test_inverses());
    }

  private:
    void test_kinds() {
        const std::array kinds={
            CodeKind{OperatorCode::CNST,OperatorKind::NULLARY},
            CodeKind{OperatorCode::VAR,OperatorKind::VARIABLE},
            CodeKind{OperatorCode::IND,OperatorKind::COORDINATE},
            CodeKind{OperatorCode::GET,OperatorKind::COORDINATE},
            CodeKind{OperatorCode::VEC,OperatorKind::NULLARY},
            CodeKind{OperatorCode::ADD,OperatorKind::BINARY},
            CodeKind{OperatorCode::SUB,OperatorKind::BINARY},
            CodeKind{OperatorCode::MUL,OperatorKind::BINARY},
            CodeKind{OperatorCode::DIV,OperatorKind::BINARY},
            CodeKind{OperatorCode::FMA,OperatorKind::TERNARY},
            CodeKind{OperatorCode::MAX,OperatorKind::BINARY},
            CodeKind{OperatorCode::MIN,OperatorKind::BINARY},
            CodeKind{OperatorCode::AND,OperatorKind::BINARY},
            CodeKind{OperatorCode::OR,OperatorKind::BINARY},
            CodeKind{OperatorCode::XOR,OperatorKind::BINARY},
            CodeKind{OperatorCode::IMPL,OperatorKind::BINARY},
            CodeKind{OperatorCode::SADD,OperatorKind::SCALAR},
            CodeKind{OperatorCode::SSUB,OperatorKind::SCALAR},
            CodeKind{OperatorCode::SMUL,OperatorKind::SCALAR},
            CodeKind{OperatorCode::SDIV,OperatorKind::SCALAR},
            CodeKind{OperatorCode::NUL,OperatorKind::UNARY},
            CodeKind{OperatorCode::POS,OperatorKind::UNARY},
            CodeKind{OperatorCode::NEG,OperatorKind::UNARY},
            CodeKind{OperatorCode::HLF,OperatorKind::UNARY},
            CodeKind{OperatorCode::REC,OperatorKind::UNARY},
            CodeKind{OperatorCode::SQR,OperatorKind::UNARY},
            CodeKind{OperatorCode::SQRT,OperatorKind::UNARY},
            CodeKind{OperatorCode::EXP,OperatorKind::UNARY},
            CodeKind{OperatorCode::LOG,OperatorKind::UNARY},
            CodeKind{OperatorCode::SIN,OperatorKind::UNARY},
            CodeKind{OperatorCode::COS,OperatorKind::UNARY},
            CodeKind{OperatorCode::TAN,OperatorKind::UNARY},
            CodeKind{OperatorCode::TANH,OperatorKind::UNARY},
            CodeKind{OperatorCode::ASIN,OperatorKind::UNARY},
            CodeKind{OperatorCode::ACOS,OperatorKind::UNARY},
            CodeKind{OperatorCode::ATAN,OperatorKind::UNARY},
            CodeKind{OperatorCode::ABS,OperatorKind::UNARY},
            CodeKind{OperatorCode::ITOR,OperatorKind::UNARY},
            CodeKind{OperatorCode::NOT,OperatorKind::UNARY},
            CodeKind{OperatorCode::POW,OperatorKind::GRADED},
            CodeKind{OperatorCode::ROOT,OperatorKind::GRADED},
            CodeKind{OperatorCode::SGN,OperatorKind::PREDICATE},
            CodeKind{OperatorCode::SUBS,OperatorKind::PREDICATE},
            CodeKind{OperatorCode::DISJ,OperatorKind::PREDICATE},
            CodeKind{OperatorCode::EQ,OperatorKind::COMPARISON},
            CodeKind{OperatorCode::NEQ,OperatorKind::COMPARISON},
            CodeKind{OperatorCode::GEQ,OperatorKind::COMPARISON},
            CodeKind{OperatorCode::LEQ,OperatorKind::COMPARISON},
            CodeKind{OperatorCode::GT,OperatorKind::COMPARISON},
            CodeKind{OperatorCode::LT,OperatorKind::COMPARISON}
        };

        for(auto const& entry:kinds) {
            Operator inferred(entry.code);
            ARIADNE_TEST_EQUALS(inferred.code(),entry.code);
            ARIADNE_TEST_EQUALS(inferred.kind(),entry.kind);

            Operator explicit_kind(entry.code,entry.kind);
            ARIADNE_TEST_EQUALS(static_cast<OperatorCode>(explicit_kind),entry.code);
            ARIADNE_TEST_EQUALS(explicit_kind.kind(),entry.kind);
        }

        const std::array kind_names={
            std::pair{OperatorKind::VARIABLE,"VARIABLE"},
            std::pair{OperatorKind::COORDINATE,"COORDINATE"},
            std::pair{OperatorKind::NULLARY,"NULLARY"},
            std::pair{OperatorKind::UNARY,"UNARY"},
            std::pair{OperatorKind::BINARY,"BINARY"},
            std::pair{OperatorKind::TERNARY,"TERNARY"},
            std::pair{OperatorKind::GRADED,"GRADED"},
            std::pair{OperatorKind::SCALAR,"SCALAR"},
            std::pair{OperatorKind::PREDICATE,"PREDICATE"},
            std::pair{OperatorKind::COMPARISON,"COMPARISON"}
        };
        for(auto const& entry:kind_names) {
            std::ostringstream stream;
            stream << entry.first;
            ARIADNE_TEST_EQUALS(stream.str(),std::string(entry.second));
        }
        std::ostringstream invalid_kind;
        invalid_kind << static_cast<OperatorKind>(127);
        ARIADNE_TEST_EQUALS(invalid_kind.str(),std::string("UNKNOWN"));

        ARIADNE_TEST_FAIL(static_cast<void>(Operator(static_cast<OperatorCode>(127))));
    }

    void test_names_and_symbols() {
        const std::array names={
            CodeName{OperatorCode::CNST,"cnst"}, CodeName{OperatorCode::VAR,"var"},
            CodeName{OperatorCode::IND,"ind"}, CodeName{OperatorCode::ADD,"add"},
            CodeName{OperatorCode::SUB,"sub"}, CodeName{OperatorCode::MUL,"mul"},
            CodeName{OperatorCode::DIV,"div"}, CodeName{OperatorCode::FMA,"fma"},
            CodeName{OperatorCode::POW,"pow"}, CodeName{OperatorCode::ROOT,"root"},
            CodeName{OperatorCode::NUL,"nul"}, CodeName{OperatorCode::POS,"pos"},
            CodeName{OperatorCode::NEG,"neg"}, CodeName{OperatorCode::HLF,"hlf"},
            CodeName{OperatorCode::REC,"rec"}, CodeName{OperatorCode::SQR,"sqr"},
            CodeName{OperatorCode::SQRT,"sqrt"}, CodeName{OperatorCode::EXP,"exp"},
            CodeName{OperatorCode::LOG,"log"}, CodeName{OperatorCode::SIN,"sin"},
            CodeName{OperatorCode::COS,"cos"}, CodeName{OperatorCode::TAN,"tan"},
            CodeName{OperatorCode::ASIN,"asin"}, CodeName{OperatorCode::ACOS,"acos"},
            CodeName{OperatorCode::ATAN,"atan"}, CodeName{OperatorCode::SADD,"sadd"},
            CodeName{OperatorCode::SSUB,"ssub"}, CodeName{OperatorCode::SMUL,"smul"},
            CodeName{OperatorCode::SDIV,"sdiv"}, CodeName{OperatorCode::SFMA,"sfma"},
            CodeName{OperatorCode::ABS,"abs"}, CodeName{OperatorCode::MAX,"max"},
            CodeName{OperatorCode::MIN,"min"}, CodeName{OperatorCode::NOT,"not"},
            CodeName{OperatorCode::AND,"and"}, CodeName{OperatorCode::OR,"or"},
            CodeName{OperatorCode::XOR,"xor"}, CodeName{OperatorCode::IMPL,"impl"},
            CodeName{OperatorCode::ITOR,"itor"}, CodeName{OperatorCode::GET,"get"},
            CodeName{OperatorCode::VEC,"vec"}, CodeName{OperatorCode::PUSH,"push"},
            CodeName{OperatorCode::PULL,"pull"}, CodeName{OperatorCode::TANH,"tanh"},
            CodeName{OperatorCode::EQ,"eq"}, CodeName{OperatorCode::NEQ,"neq"},
            CodeName{OperatorCode::GEQ,"geq"}, CodeName{OperatorCode::LEQ,"leq"},
            CodeName{OperatorCode::GT,"gt"}, CodeName{OperatorCode::LT,"lt"},
            CodeName{OperatorCode::SGN,"sgn"}, CodeName{OperatorCode::SUBS,"subs"},
            CodeName{OperatorCode::DISJ,"disj"}
        };
        for(auto const& entry:names) {
            ARIADNE_TEST_EQUALS(std::string(name(entry.code)),std::string(entry.text));
            std::ostringstream stream;
            stream << entry.code;
            ARIADNE_TEST_EQUALS(stream.str(),std::string(entry.text));
        }
        ARIADNE_TEST_EQUALS(std::string(name(static_cast<OperatorCode>(127))),std::string("UNKNOWN"));

        const std::array symbols={
            CodeName{OperatorCode::POS,"+"}, CodeName{OperatorCode::NEG,"-"},
            CodeName{OperatorCode::ADD,"+"}, CodeName{OperatorCode::SUB,"-"},
            CodeName{OperatorCode::MUL,"*"}, CodeName{OperatorCode::DIV,"/"},
            CodeName{OperatorCode::POW,"^"}, CodeName{OperatorCode::SADD,"+"},
            CodeName{OperatorCode::SMUL,"*"}, CodeName{OperatorCode::NOT,"!"},
            CodeName{OperatorCode::AND,"&"}, CodeName{OperatorCode::OR,"|"},
            CodeName{OperatorCode::EQ,"=="}, CodeName{OperatorCode::NEQ,"!="},
            CodeName{OperatorCode::LEQ,"<="}, CodeName{OperatorCode::GEQ,">="},
            CodeName{OperatorCode::LT,"<"}, CodeName{OperatorCode::GT,">"}
        };
        for(auto const& entry:symbols) {
            ARIADNE_TEST_EQUALS(std::string(symbol(entry.code)),std::string(entry.text));
        }
        ARIADNE_TEST_EQUALS(std::string(symbol(OperatorCode::CNST)),std::string("???"));
    }

    void test_objects() {
        test_metadata(Gtr{},OperatorCode::GT,OperatorKind::COMPARISON);
        test_metadata(Less{},OperatorCode::LT,OperatorKind::COMPARISON);
        test_metadata(Geq{},OperatorCode::GEQ,OperatorKind::COMPARISON);
        test_metadata(Leq{},OperatorCode::LEQ,OperatorKind::COMPARISON);
        test_metadata(Equal{},OperatorCode::EQ,OperatorKind::COMPARISON);
        test_metadata(Unequal{},OperatorCode::NEQ,OperatorKind::COMPARISON);
        test_metadata(XOrOp{},OperatorCode::XOR,OperatorKind::BINARY);
        test_metadata(AndOp{},OperatorCode::AND,OperatorKind::BINARY);
        test_metadata(OrOp{},OperatorCode::OR,OperatorKind::BINARY);
        test_metadata(NotOp{},OperatorCode::NOT,OperatorKind::UNARY);
        test_metadata(Cnst{},OperatorCode::CNST,OperatorKind::NULLARY);
        test_metadata(Ind{},OperatorCode::IND,OperatorKind::COORDINATE);
        test_metadata(Var{},OperatorCode::VAR,OperatorKind::VARIABLE);
        test_metadata(Add{},OperatorCode::ADD,OperatorKind::BINARY);
        test_metadata(Sub{},OperatorCode::SUB,OperatorKind::BINARY);
        test_metadata(Mul{},OperatorCode::MUL,OperatorKind::BINARY);
        test_metadata(Div{},OperatorCode::DIV,OperatorKind::BINARY);
        test_metadata(Pow{},OperatorCode::POW,OperatorKind::GRADED);
        test_metadata(Root{},OperatorCode::ROOT,OperatorKind::GRADED);
        test_metadata(Fma{},OperatorCode::FMA,OperatorKind::TERNARY);
        test_metadata(Nul{},OperatorCode::NUL,OperatorKind::UNARY);
        test_metadata(Pos{},OperatorCode::POS,OperatorKind::UNARY);
        test_metadata(Neg{},OperatorCode::NEG,OperatorKind::UNARY);
        test_metadata(Hlf{},OperatorCode::HLF,OperatorKind::UNARY);
        test_metadata(Rec{},OperatorCode::REC,OperatorKind::UNARY);
        test_metadata(Sqr{},OperatorCode::SQR,OperatorKind::UNARY);
        test_metadata(Sqrt{},OperatorCode::SQRT,OperatorKind::UNARY);
        test_metadata(Exp{},OperatorCode::EXP,OperatorKind::UNARY);
        test_metadata(Log{},OperatorCode::LOG,OperatorKind::UNARY);
        test_metadata(Sin{},OperatorCode::SIN,OperatorKind::UNARY);
        test_metadata(Cos{},OperatorCode::COS,OperatorKind::UNARY);
        test_metadata(Tan{},OperatorCode::TAN,OperatorKind::UNARY);
        test_metadata(Tanh{},OperatorCode::TANH,OperatorKind::UNARY);
        test_metadata(Asin{},OperatorCode::ASIN,OperatorKind::UNARY);
        test_metadata(Acos{},OperatorCode::ACOS,OperatorKind::UNARY);
        test_metadata(Atan{},OperatorCode::ATAN,OperatorKind::UNARY);
        test_metadata(Max{},OperatorCode::MAX,OperatorKind::BINARY);
        test_metadata(Min{},OperatorCode::MIN,OperatorKind::BINARY);
        test_metadata(Abs{},OperatorCode::ABS,OperatorKind::UNARY);
        test_metadata(Sgn{},OperatorCode::SGN,OperatorKind::PREDICATE);

        std::ostringstream comparisons;
        comparisons << Less{} << Gtr{} << Leq{} << Geq{} << Equal{} << Unequal{};
        ARIADNE_TEST_EQUALS(comparisons.str(),std::string("<><=>===!="));

        std::ostringstream logical;
        logical << AndOp{} << OrOp{} << NotOp{};
        ARIADNE_TEST_EQUALS(logical.str(),std::string("&&||!"));

        std::ostringstream builtin;
        builtin << Plus{} << Minus{} << Times{} << Divides{};
        ARIADNE_TEST_EQUALS(builtin.str(),std::string("+-*/"));

        std::ostringstream arithmetic;
        arithmetic << Add{} << Sub{} << Mul{} << Div{} << Pow{} << Root{} << Fma{}
                   << Nul{} << Pos{} << Neg{} << Hlf{} << Rec{} << Sqr{} << Sqrt{}
                   << Exp{} << Log{} << Sin{} << Cos{} << Tan{} << Tanh{} << Asin{}
                   << Acos{} << Atan{} << Max{} << Min{} << Abs{};
        ARIADNE_TEST_ASSERT(not arithmetic.str().empty());
    }

    void test_simple_operations() {
        Integer z1(1);
        Integer z2(2);
        Integer z3(3);

        ARIADNE_TEST_ASSERT(Less{}(z1,z2));
        ARIADNE_TEST_ASSERT(Gtr{}(z2,z1));
        ARIADNE_TEST_ASSERT(Leq{}(z1,z1));
        ARIADNE_TEST_ASSERT(Geq{}(z2,z1));
        ARIADNE_TEST_ASSERT(Equal{}(z2,z2));
        ARIADNE_TEST_ASSERT(Unequal{}(z1,z2));

        ARIADNE_TEST_ASSERT(AndOp{}(true,true));
        ARIADNE_TEST_ASSERT(OrOp{}(false,true));
        ARIADNE_TEST_ASSERT(XOrOp{}(true,false));
        ARIADNE_TEST_ASSERT(NotOp{}(false));

        ARIADNE_TEST_EQUALS(Plus{}(2),2);
        ARIADNE_TEST_EQUALS(Plus{}(2,3),5);
        ARIADNE_TEST_EQUALS(Minus{}(2),-2);
        ARIADNE_TEST_EQUALS(Minus{}(5,3),2);
        ARIADNE_TEST_EQUALS(Times{}(2,3),6);
        ARIADNE_TEST_EQUALS(Divides{}(6,3),2);

        ARIADNE_TEST_EQUALS(Add{}(z2,z3),Integer(5));
        ARIADNE_TEST_EQUALS(Sub{}(z3,z2),Integer(1));
        ARIADNE_TEST_EQUALS(Mul{}(z2,z3),Integer(6));
        ARIADNE_TEST_EQUALS(Pow{}(z2,3u),Integer(8));
        ARIADNE_TEST_EQUALS(Nul{}(z3),Integer(0));
        ARIADNE_TEST_EQUALS(Pos{}(z3),Integer(3));
        ARIADNE_TEST_EQUALS(Neg{}(z3),Integer(-3));
        ARIADNE_TEST_EQUALS(Sqr{}(z3),Natural(9u));
        ARIADNE_TEST_EQUALS(Max{}(z2,z3),Integer(3));
        ARIADNE_TEST_EQUALS(Min{}(z2,z3),Integer(2));
        ARIADNE_TEST_EQUALS(Abs{}(Integer(-3)),Natural(3u));
        ARIADNE_TEST_EQUALS(Sgn{}(Integer(-3)),Sign::NEGATIVE);

        ARIADNE_TEST_EQUALS(Hlf{}(Integer(4)),Dyadic(2));
        ARIADNE_TEST_EQUALS(Rec{}(Integer(4)),Rational(1,4));

        BinaryComparisonOperator less(Less{});
        ARIADNE_TEST_ASSERT(less(z1,z2));
        BinaryRingOperator add(Add{});
        ARIADNE_TEST_EQUALS(add(z2,z3),Integer(5));
        GradedRingOperator power(Pow{});
        ARIADNE_TEST_EQUALS(power(z2,3u),Integer(8));
    }

    void test_inverses() {
        ARIADNE_TEST_EQUALS(inverse(Pos{}).code(),OperatorCode::POS);
        ARIADNE_TEST_EQUALS(inverse(Neg{}).code(),OperatorCode::NEG);
        ARIADNE_TEST_EQUALS(inverse(Rec{}).code(),OperatorCode::REC);
        ARIADNE_TEST_EQUALS(inverse(Sqr{}).code(),OperatorCode::SQRT);
        ARIADNE_TEST_EQUALS(inverse(Sqrt{}).code(),OperatorCode::SQR);
        ARIADNE_TEST_EQUALS(inverse(Exp{}).code(),OperatorCode::LOG);
        ARIADNE_TEST_EQUALS(inverse(Log{}).code(),OperatorCode::EXP);
        ARIADNE_TEST_EQUALS(inverse(Tan{}).code(),OperatorCode::ATAN);
        ARIADNE_TEST_EQUALS(inverse(Atan{}).code(),OperatorCode::TAN);

        UnaryTranscendentalOperator exponential(Exp{});
        UnaryTranscendentalOperator logarithm(Log{});
        UnaryTranscendentalOperator sine(Sin{});
        ARIADNE_TEST_ASSERT(are_inverses(exponential,logarithm));
        ARIADNE_TEST_ASSERT(are_inverses(logarithm,exponential));
        ARIADNE_TEST_ASSERT(not are_inverses(sine,logarithm));
    }
};

int main() {
    ARIADNE_TEST_CLASS(Operators,TestOperators());
    return ARIADNE_TEST_FAILURES;
}
