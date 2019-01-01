/*This is a template file for use with 1D finite elements.
  The portions of the code you need to fill in are marked with the comment "//EDIT".

  Do not change the name of any existing functions, but feel free
  to create additional functions, variables, and constants.
  It uses the deal.II FEM library.*/

//Include files
//Data structures and solvers
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/numerics/data_out.h>
//Mesh related classes
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
//Finite element implementation classes
#include <deal.II/fe/fe_system.h>
//Standard C++ libraries
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

#include "FEM1DGauSrc.h"

using namespace dealii;

class FEM1DWave
{
	static const unsigned int N = 1000;
	static constexpr double t_max = 1, T = t_max / N, theta = 0.5;
	unsigned int prob, n;
	double		 L;

	//Class objects
	Triangulation<1>   triangulation; //mesh
	FESystem<1>        fe; 	      //FE element
	DoFHandler<1>      dof_handler;   //Connectivity matrices
	QGauss<1>		   basis;	//basis function
	MappingQ1<1, 1>    mapping;
	ConstraintMatrix   constraints;
	FEM1DGauSrc		   src;

	public:
		//Class functions
		FEM1DWave (unsigned int order,unsigned int problem); // Class constructor 
		~FEM1DWave(); //Class destructor

		//Solution steps
		void generate_mesh(unsigned int numberOfElements);
		void define_boundary_conds();
		void setup_system();
		void assemble_system();
		void solve_p();
		void solve_v();
		void run();
		void output_results();

		//Function to calculate the l2 norm of the error in the finite element sol'n vs. the exact solution
		double l2norm_of_error();
		double u_exact(double x);
				    
		//Data structures
		SparsityPattern       sparsity_pattern; //Sparse matrix pattern
		SparseMatrix<double>  K, M, Au;		 //Global stiffness (sparse) matrix
		Vector<double>        P1, P2, V, F, RHS;  //Global vectors - Solutions (P1,P2), Derivatives (V), Global force vector (F), and Global right-hand side (RHS)
		std::vector<double>   nodeLocation;	 //Vector of the x-coordinate of nodes by global dof number
		//std::map<unsigned int,double> boundary_values;	//Map of dirichlet boundary conditions

		//solution name array
		std::vector<std::string> nodal_solution_names;
		std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};
