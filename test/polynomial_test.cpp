#define CATCH_CONFIG_MAIN 
#include <catch/catch.hpp>

#include <sisl/sisl.hpp>
#include <sisl/polynomial.hpp>

TEST_CASE("Constant constructor test", "polycon"){
	using namespace sisl;

	ConstantTerm<3> p(100.);
	ConstantTerm<3> q(p);


	REQUIRE(p.evaluate({0.1, 0.1, 0.1}) == 100.);
	REQUIRE(p.evaluate({0.1, 0.1, 0.1}) == q.evaluate({0.1, 0.1, 0.1}));
}


TEST_CASE("linear constructor test", "polylin"){
	using namespace sisl;

	VariableTerm<3> p(2);
	VariableTerm<3> q(p);

	REQUIRE(p.evaluate({0.1, 0.1, 0.1}) == 0.1);
	REQUIRE(p.evaluate({0.1, 0.1, 0}) == 0);
	REQUIRE(p.evaluate({0.1, 0.1, 0.1}) == q.evaluate({0.1, 0.1, 0.1}));
}

TEST_CASE("polynomial constructor test", "polyconst"){
	using namespace sisl;

	Polynomial<3> p;
	Polynomial<3> q(p);

	REQUIRE(p.evaluate({0.1, 0.1, 0}) == 0);
	REQUIRE(p.evaluate({0.1, 0.1, 0.1}) == q.evaluate({0.1, 0.1, 0.1}));
}

TEST_CASE("polynomial operator test", "polyconst"){
	using namespace sisl;

	Polynomial<3> p = 1. + PVar<3>(1);
	Polynomial<3> q = 2. * PVar<3>(0);

	REQUIRE(p.evaluate({1, 0, 1}) == 1);
	REQUIRE(p.evaluate({1, 1, 1}) == 2);

	REQUIRE(q.evaluate({1, 0, 1}) == 2);
	REQUIRE(q.evaluate({0, 1, 1}) == 0);

	REQUIRE((p+q).evaluate({1, 1, 0}) == 4);
	REQUIRE((q*p).evaluate({1, 1, 0}) == 4);

}