// FEM1D-Test.cpp : Defines the entry point for the console application.
//
//Include files
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "FEM1.h"
#include "writeSolutions.h"

//The main program, using the FEM class
int main(int argc, char *argv[]) {
	try {
		deallog.depth_console(0);

		//Specify the basis function order: 1, 2, or 3
		unsigned int order = 1;

		//Specify the subproblem: 1 or 2
		unsigned int problem = 1;

		std::istringstream ss(argv[1]);
		if (argc > 0)
			ss >> order;
		if (argc > 1)
			ss >> problem;

		FEM<1> problemObject(order, problem);

		//Define the number of elements as an input to "generate_mesh"
		problemObject.generate_mesh(10); //e.g. a 10 element mesh
		problemObject.setup_system();
		problemObject.assemble_system();
		problemObject.solve();
		double l2norm = problemObject.l2norm_of_error();
		std::cout << l2norm << std::endl;

		//write output file in vtk format for visualization
		problemObject.output_results();

		//write solutions to h5 file
		char tag[21];
		snprintf(tag, 21, "CA1_Order%d_Problem%d", order, problem);
		writeSolutionsToFileCA1(problemObject.D, l2norm, tag);
	}
	catch (std::exception &exc) {
		std::cerr << std::endl << std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		std::cerr << "Exception on processing: " << std::endl
			<< exc.what() << std::endl
			<< "Aborting!" << std::endl
			<< "----------------------------------------------------"
			<< std::endl;

		return 1;
	}
	catch (...) {
		std::cerr << std::endl << std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		std::cerr << "Unknown exception!" << std::endl
			<< "Aborting!" << std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		return 1;
	}

	return 0;
}
