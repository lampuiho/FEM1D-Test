#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/base/logstream.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include "FEM1DWave.h"

// Class constructor for a vector field
FEM1DWave::FEM1DWave(unsigned int order, unsigned int problem)
	:
	fe(FE_Q<1>(order), 1),
	dof_handler(triangulation),
	basis(QGauss<1>(order + 1)),
	src(FEM1DGauSrc(Point<1>(0.5), 0.05, 0.05, 0.1 / T)),
	n(0)
{
	if (problem == 1 || problem == 2) {
		prob = problem;
	}
	else {
		std::cout << "Error: problem number should be 1 or 2.\n";
		exit(0);
	}

	//Nodal Solution names - this is for writing the output file
	nodal_solution_names.push_back("p");
	nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
}

//Class destructor
FEM1DWave::~FEM1DWave() {
	dof_handler.clear();
}

//Define the problem domain and generate the mesh
void FEM1DWave::generate_mesh(unsigned int numberOfElements) {

	//Define the limits of your domain
	L = 1; //EDIT
	double x_min = 0.;
	double x_max = L;

	Point<1, double> min(x_min), max(x_max);
	std::vector<unsigned int> meshDimensions(1, numberOfElements);
	GridGenerator::subdivided_hyper_rectangle(triangulation, meshDimensions, min, max);
}

//Setup data structures (sparse matrix, vectors)
void FEM1DWave::setup_system() {
	//Let deal.II organize degrees of freedom
	dof_handler.distribute_dofs(fe);
	
	//Enter the global node x-coordinates into the vector "nodeLocation"
	std::vector< Point<1, double> > dof_coords(dof_handler.n_dofs());
	nodeLocation.resize(dof_handler.n_dofs());
	DoFTools::map_dofs_to_support_points<1, 1>(mapping, dof_handler, dof_coords);
	for (unsigned int i = 0; i < dof_coords.size(); i++) {
		nodeLocation[i] = dof_coords[i][0];
	}

	//Specify boundary condtions (call the function)
	define_boundary_conds();

	//Define the size of the global matrices and vectors
	sparsity_pattern.reinit(dof_handler.n_dofs(), dof_handler.n_dofs(),
		dof_handler.max_couplings_between_dofs());
	DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
	sparsity_pattern.compress();
	K.reinit(sparsity_pattern);
	M.reinit(sparsity_pattern);
	Au.reinit(sparsity_pattern);
	F.reinit(dof_handler.n_dofs());
	P1.reinit(dof_handler.n_dofs());
	P2.reinit(dof_handler.n_dofs());
	V.reinit(dof_handler.n_dofs());
	RHS.reinit(dof_handler.n_dofs());
	constraints.close();
	
	//Just some notes...
	std::cout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
	std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;
}

//Form elmental vectors and matrices and assemble to the global vector (F) and matrix (K)
void FEM1DWave::assemble_system() {
	MatrixCreator::create_mass_matrix(dof_handler, basis, M);
	MatrixCreator::create_laplace_matrix(dof_handler, basis,K);
	Au = M;
	Au.add(T*T*theta*theta, K);

	//zero initial conditions
	P1 = 0;
	P2 = 0;
	V = 0;
}

void FEM1DWave::solve_p()
{
	SolverControl solver_control(1000, 1e-8*RHS.l2_norm());
	SolverCG<>    cg(solver_control);

	cg.solve(Au, P1, RHS, PreconditionIdentity());
	std::cout << "   u-equation: " << solver_control.last_step() << " CG iterations." << std::endl;
}

void FEM1DWave::solve_v()
{
	SolverControl           solver_control(1000, 1e-8*F.l2_norm());
	SolverCG<>              cg(solver_control);
	cg.solve(M, P2, RHS, PreconditionIdentity());
	V += P2;
	std::cout << "   v-equation: " << solver_control.last_step() << " CG iterations." << std::endl;
}

void FEM1DWave::run() {
	Vector<double> tmp(P1.size());
	double t;

	for (; n <= N;)
	{
		++n;
		t = n * T;

		std::cout << "Time step " << n << " at t=" << t	<< std::endl;
		
		P2 = P1;
		tmp = P2;
		tmp *= T * T*theta*(1 - theta);
		K.vmult_add(RHS, tmp);
		tmp = V;
		tmp.sadd(T, P2);
		M.vmult_add(RHS, tmp);

		tmp = F;
		tmp *= (1 - theta);
		src.set_time(t);
		VectorTools::create_right_hand_side(mapping, dof_handler, basis, src, F);
		tmp.add(theta, F);
		RHS.add(T * T * theta, tmp);

		//if dirichlet boundary required... update matrix here
		solve_p();

		RHS = tmp;
		RHS *= T;
		P2 *= -T * (1 - theta);
		P2.add(-T * theta,P1);
		K.vmult_add(RHS, P2);
		solve_v();

		//write output file in vtk format for visualization
		output_results();

		std::cout << "   Total energy: " << (M.matrix_norm_square(V) + K.matrix_norm_square(P1)) / 2 << std::endl;
	}
}

//Specify the Dirichlet boundary conditions
void FEM1DWave::define_boundary_conds() {
	//const unsigned int totalNodes = dof_handler.n_dofs(); //Total number of nodes

	//Identify dirichlet boundary nodes and specify their values.
	//This function is called from within "setup_system"

	/*The vector "nodeLocation" gives the x-coordinate in the real domain of each node,
	  organized by the global node number.*/

	  /*This loops through all nodes in the system and checks to see if they are
		at one of the boundaries. If at a Dirichlet boundary, it stores the node number
		and the applied displacement value in the std::map "boundary_values". Deal.II
		will use this information later to apply Dirichlet boundary conditions.
		Neumann boundary conditions are applied when constructing Flocal in "assembly"*/
		/*
		for(unsigned int globalNode=0; globalNode<totalNodes; globalNode++){
		  if(nodeLocation[globalNode] == 0){
			  boundary_values[globalNode] = g1;
		  }
		  if(nodeLocation[globalNode] == L){
			  if(prob == 1){
			  boundary_values[globalNode] = g2;
			  }
		  }
		}
		*/
}

//Output results
void FEM1DWave::output_results() {
	const std::string tag = "./sln/wave_o" + Utilities::int_to_string(basis.size() - 1) + "_step_" + Utilities::int_to_string(n) + ".vtk";
	//Write results to VTK file
	std::ofstream output1(tag);
	DataOut<1> data_out;
	data_out.attach_dof_handler(dof_handler);

	//Add nodal DOF data
	data_out.add_data_vector(P1, nodal_solution_names, DataOut<1>::type_dof_data, nodal_data_component_interpretation);
	data_out.build_patches();
	data_out.write_vtk(output1);
	output1.close();
}

double FEM1DWave::u_exact(double x) {
	double result = 0;
	//diriclet bc
	if (prob == 2)
	{
	}
	else
	{
	}
	return result;
}

double FEM1DWave::l2norm_of_error() {
	/*
	double l2norm = 0.;

	//Find the l2 norm of the error between the finite element sol'n and the exact sol'n
	const unsigned int   			dofs_per_elem = fe.dofs_per_cell; //This gives you dofs per element
	std::vector<unsigned int> local_dof_indices(dofs_per_elem);
	Vector<double> u_ex; u_ex.reinit(dof_handler.n_dofs());
	double u_h, x, h_e;

	//loop over elements  
	for (unsigned int q = 0; q < dof_handler.n_dofs(); q++)
		u_ex[q] = u_exact(nodeLocation[q]);
	u_ex.print(std::cout);
	for (auto i = nodeLocation.begin(); i != nodeLocation.end(); ++i)
		std::cout << *i << ' ';
	std::cout << std::endl;

	typename DoFHandler<1>::active_cell_iterator elem = dof_handler.begin_active(),
		endc = dof_handler.end();
	for (; elem != endc; ++elem) {
		//Retrieve the effective "connectivity matrix" for this element
		elem->get_dof_indices(local_dof_indices);

		//Find the element length
		h_e = nodeLocation[local_dof_indices[1]] - nodeLocation[local_dof_indices[0]];

		double sum = 0.;
		for (unsigned int q = 0; q < quadRule; q++) {
			x = 0.;
			u_h = 0.;
			//Find the values of x and u_h (the finite element solution) at the quadrature points
			for (unsigned int B = 0; B < dofs_per_elem; B++) {
				u_h += D[local_dof_indices[B]] * basis_function(B, quad_points[q]);
				x += nodeLocation[local_dof_indices[B]] * basis_function(B, quad_points[q]);
			}
			//EDIT - Find the l2-norm of the error through numerical integration.
			//This includes evaluating the exact solution at the quadrature points
			sum += pow(u_exact(x) - u_h, 2) * quad_weight[q];
		}
		sum *= h_e / 2;
		l2norm += sum;
	}
	return sqrt(l2norm);
	*/
	return 0.0;
}
