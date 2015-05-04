#define CATCH_CONFIG_MAIN 
#include <catch/catch.hpp>
#include <sisl/array.hpp>
#include <sisl/sparse_array.hpp>
#include <sisl/fftwalloc.hpp>

#include <iostream>

TEST_CASE("Array Write Test","arraywrite"){
	sisl::array3<float> a(200,200,200);
	REQUIRE(a.size() == 200*200*200);

	printf("AW Test\n");

	float val = 100;
	bool valsWereWrong = false;

	for(int i = 0; i < 200; i++)
		for(int j = 0; j < 200; j++)
			for(int k = 0; k < 200; k++){
				a(i,j,k) = val++;
			}
	
	val = 100;
	for(int i = 0; i < 200; i++)
		for(int j = 0; j < 200; j++)
			for(int k = 0; k < 200; k++){
				valsWereWrong |= (a(i,j,k) != val++);
			}

	REQUIRE(valsWereWrong == false);
}

TEST_CASE("4D Array Write Test","arraywrite"){
	sisl::array4<float> a(10,10,10,10);
	REQUIRE(a.size() == 200*200*200);

	printf("AW Test\n");

	float val = 100;
	bool valsWereWrong = false;

	for(int i = 0; i < 10; i++)
		for(int j = 0; j < 10; j++)
			for(int k = 0; k < 10; k++)
				for(int l = 0; l < 10; l++){
					a(i,j,k,l) = val++;
				}
	
	val = 100;
	for(int i = 0; i < 10; i++)
		for(int j = 0; j < 10; j++)
			for(int k = 0; k < 10; k++)
				for(int l = 0; l < 10; l++){
					valsWereWrong |= (a(i,j,k,l) != val++);
				}

	REQUIRE(valsWereWrong == false);
}
TEST_CASE("2d Array Write Test","arraywrite2d"){
	sisl::array2<float> a(200,200);
	REQUIRE(a.size() == 200*200);
	printf("2D Test\n");
	
	float val = 100;
	bool valsWereWrong = false;

	for(int i = 0; i < 200; i++)
		for(int j = 0; j < 200; j++){
				a(i,j) = val++;
			}
	
	val = 100;
	for(int i = 0; i < 200; i++)
		for(int j = 0; j < 200; j++){
			valsWereWrong |= (a(i,j) != val++);
		}

	REQUIRE(valsWereWrong == false);
}
TEST_CASE("Index Out of Bounds",""){
	sisl::array3<float> a(200,200,200);
	printf("IOB Test\n");

	bool exceptionThrew = false;
	try{
		a(150,150,200) = 10;
	}catch(...){
		exceptionThrew = true;
	}
	REQUIRE(exceptionThrew == true);

	exceptionThrew = false;
	try{
		a(150,150,-1) = 10;
	}catch(...){
		exceptionThrew = true;
	}
	REQUIRE(exceptionThrew == true);

	exceptionThrew = false;
	try{
		a(200,150,150) = 10;
	}catch(...){
		exceptionThrew = true;
	}
	REQUIRE(exceptionThrew == true);

	exceptionThrew = false;
	try{
		a(-1,150,150) = 10;
	}catch(...){
		exceptionThrew = true;
	}
	REQUIRE(exceptionThrew == true);

	exceptionThrew = false;
	try{
		a(150,200,150) = 10;
	}catch(...){
		exceptionThrew = true;
	}
	REQUIRE(exceptionThrew == true);

	exceptionThrew = false;
	try{
		a(150,-1,150) = 10;
	}catch(...){
		exceptionThrew = true;
	}
	REQUIRE(exceptionThrew == true);
}

TEST_CASE("FFTW Alloc",""){
		printf("FFTW Test\n");
	sisl::array3<float, _fftwalloc<float> > a(200,200,200);


	bool exceptionThrew = false;
	try{
		a(150,150,200) = 10;
	}catch(...){
		exceptionThrew = true;
	}
	REQUIRE(exceptionThrew == true);

	exceptionThrew = false;
	try{
		a(150,150,-1) = 10;
	}catch(...){
		exceptionThrew = true;
	}
	REQUIRE(exceptionThrew == true);

	exceptionThrew = false;
	try{
		a(200,150,150) = 10;
	}catch(...){
		exceptionThrew = true;
	}
	REQUIRE(exceptionThrew == true);

	exceptionThrew = false;
	try{
		a(-1,150,150) = 10;
	}catch(...){
		exceptionThrew = true;
	}
	REQUIRE(exceptionThrew == true);

	exceptionThrew = false;
	try{
		a(150,200,150) = 10;
	}catch(...){
		exceptionThrew = true;
	}
	REQUIRE(exceptionThrew == true);

	exceptionThrew = false;
	try{
		a(150,-1,150) = 10;
	}catch(...){
		exceptionThrew = true;
	}
	REQUIRE(exceptionThrew == true);
}

TEST_CASE("Sparse map",""){
	sisl::sparse_array3<float> a(200,200,200,0.f);
	float val = 100;
	bool valsWereWrong = false;

	REQUIRE(a.hasValue(10,10,10) == false);

	for(int i = 0; i < 200; i++)
		for(int j = 0; j < 200; j++)
			for(int k = 0; k < 200; k++){
				a(i,j,k) = val++;
			}
	
	val = 100;
	for(int i = 0; i < 200; i++)
		for(int j = 0; j < 200; j++)
			for(int k = 0; k < 200; k++){
				valsWereWrong |= (a(i,j,k) != val++);
			}

	REQUIRE(valsWereWrong == false);
	REQUIRE(a.hasValue(10,10,10) == true);

	a(10,10,10) = 100;
	REQUIRE(a(10,10,10) == 100);
}