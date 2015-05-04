#define CATCH_CONFIG_MAIN 
#include <catch/catch.hpp>

#include <sisl/utility/ply_writer.hpp>

#include <iostream>

using namespace sisl;
using namespace sisl::utility;

TEST_CASE("Ply binary writer test", "plytest"){
		ply_writer<double> ply;

		ply.addTriangle(sisl::vertex3<double>(1,1,1,1,1,1),
						sisl::vertex3<double>(1,1,-1,1,1,-1),
						sisl::vertex3<double>(1,-1,-1,1,-1,-1));

		ply.addTriangle(sisl::vertex3<double>(1,1,1,1,1,1),
						sisl::vertex3<double>(1,-1,1,1,-1,1),
						sisl::vertex3<double>(1,-1,-1,1,-1,-1));

		ply.addTriangle(sisl::vertex3<double>(1,1,1,1,1,1),
						sisl::vertex3<double>(-1,1,1,-1,1,1),
						sisl::vertex3<double>(1,-1,1,1,-1,1));

		ply.addTriangle(
						sisl::vertex3<double>(-1,-1,1,-1,-1,1),
						sisl::vertex3<double>(-1,1,1,-1,1,1),
						sisl::vertex3<double>(1,-1,1,1,-1,1));

		ply.addTriangle(sisl::vertex3<double>(1,1,1,1,1,1),
						sisl::vertex3<double>(1,1,-1,1,1,-1),
						sisl::vertex3<double>(-1,1,1,-1,1,1));
		
		ply.addTriangle(sisl::vertex3<double>(-1,1,-1,-1,1,-1),
						sisl::vertex3<double>(-1,1,1,-1,1,1),
						sisl::vertex3<double>(1,1,-1,1,1,-1));

		ply.addTriangle(sisl::vertex3<double>(-1,1,1,1,1,1),
						sisl::vertex3<double>(-1,1,-1,1,1,-1),
						sisl::vertex3<double>(-1,-1,-1,1,-1,-1));

		ply.addTriangle(sisl::vertex3<double>(-1,1,1,1,1,1),
						sisl::vertex3<double>(-1,-1,1,1,-1,1),
						sisl::vertex3<double>(-1,-1,-1,1,-1,-1));

		ply.addTriangle(sisl::vertex3<double>(1,1,-1,1,1,1),
						sisl::vertex3<double>(-1,1,-1,-1,1,1),
						sisl::vertex3<double>(1,-1,-1,1,-1,1));

		ply.addTriangle(
						sisl::vertex3<double>(-1,-1,-1,-1,-1,1),
						sisl::vertex3<double>(-1,1,-1,-1,1,1),
						sisl::vertex3<double>(1,-1,-1,1,-1,1));

		ply.addTriangle(sisl::vertex3<double>(1,-1,1,1,1,1),
						sisl::vertex3<double>(1,-1,-1,1,1,-1),
						sisl::vertex3<double>(-1,-1,1,-1,1,1));
		
		ply.addTriangle(sisl::vertex3<double>(-1,-1,-1,-1,1,-1),
						sisl::vertex3<double>(-1,-1,1,-1,1,1),
						sisl::vertex3<double>(1,-1,-1,1,1,-1));

		if(ply.countVertices() != 8) {
			printf("The test model should have 8 vertices, not %d \n", ply.countVertices());
		//	FAIL();
		}

		if(ply.countFaces() != 12) {
			printf("The test model should have 12 faces, not %d \n", ply.countFaces());
		//	FAIL();
		}

		if(!ply.writePly("ply_tests_1.ply")) {
			printf("Couldn't write the PLY file!\n");
		//	FAIL();
		}

	}
