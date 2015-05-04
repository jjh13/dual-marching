#define CATCH_CONFIG_MAIN 
#include <catch/catch.hpp>
#include <sisl/sisl.hpp>

#include <iostream>

TEST_CASE("Matrix Mult Test","mmtest"){
	using namespace sisl;
	vector3<float> v(1,2,3);

	matrix4x4<float> m;
	m.mat[0][0] = 1;
	m.mat[1][0] = 1;
	m.mat[2][0] = 1;
	m.mat[3][0] = 0;
	m.mat[0][1] = 2;
	m.mat[1][1] = 2;
	m.mat[2][1] = 2;
	m.mat[3][1] = 0;
	m.mat[0][2] = 3;
	m.mat[1][2] = 3;
	m.mat[2][2] = 3;
	m.mat[3][2] = 0;
	m.mat[0][3] = 0;
	m.mat[1][3] = 0;
	m.mat[2][3] = 0;
	m.mat[3][3] = 1;

	v = m * v;
	REQUIRE(v.i == 6);
	REQUIRE(v.j == 12); 
	REQUIRE(v.k == 18);
}

TEST_CASE("Index Out of Bounds",""){
}

TEST_CASE("FFTW Alloc",""){
}
