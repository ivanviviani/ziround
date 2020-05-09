/**
 * @file main.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

/**
 * @brief Main.
 * 
 * @param argc Number of command line arguments (includes the program name).
 * @param argv Command line arguments (strings).
 * @return Exit status.
 */
int main(int argc, char** argv) {

	// Local variables
	instance inst;       /**< General support structure. */
	int rounded = 1;     /**< Support flag. */
	FILE* output = NULL; /**< Output file pointer. */
	
	// Initialize instance
	init_inst(&inst);
	
	// Parse command line arguments
	parse_cmd(argc, argv, &inst);

	// Set up CPLEX environment
	setup_CPLEX_env(&inst);                   print_verbose(100, "[INFO][OK]: CPLEX env setup.\n");

	// Read MIP problem from input file (.mps)
	read_MIP_problem(&inst, inst.input_file); print_verbose(100, "[INFO][OK]: MIP problem read.\n");

	// Remember integer/binary variables of the original MIP
	save_integer_variables(&inst);            print_verbose(100, "[INFO][OK]: Integer variables saved.\n");
	
	// Solve continuous relaxation of the MIP problem
	solve_continuous_relaxation(&inst);       print_verbose(100, "[INFO][OK]: Continuous relaxation solved.\n");

	// Print objective value of continuous relaxation
	print_verbose(100, "[INFO]: Continuous relaxation objective value: %.10g.\n", inst.objval);

	// Populate the instance with problem data retrieved from the CPLEX lp
	populate_inst(&inst);                     print_verbose(100, "[INFO][OK]: Problem data read.\n");

	// Print problem info to file
	if (VERBOSE >= 200) print_problem_info(&inst, 1, 1);
	if (VERBOSE >= 150) {
		output = fopen("output.txt", "a");
		fprintf(output, "\n[INFO]: Continuous relaxation objective value: %.10g.\n", inst.objval);
		fprintf(output, "[INFO]: Solution of continuous relaxation: \n");
		for (int j = 0; j < inst.ncols; j++) fprintf(output, "%f ", inst.x[j]);
		fclose(output);
	}

	print_verbose(100, "\n**************************************************************************\n");
	print_verbose(100, "[INFO]: ... Starting ZI-Round ...\n");
	print_verbose(100, "**************************************************************************\n");
//******************************************* ZI-ROUND *******************************************
	zi_round(&inst);
//************************************************************************************************
	print_verbose(100, "**************************************************************************\n");
	print_verbose(100, "[INFO][OK]: ZI-Round terminated. Solution fractionality: %f.\n", sol_fractionality(inst.x, inst.int_var, inst.ncols));
	print_verbose(100, "**************************************************************************\n");

	// Print objective value of candidate rounded solution
	print_verbose(100, "[INFO]: Candidate objective value: %f\n", inst.objval);

	// Print candidate rounded solution and its objective value to file
	if (VERBOSE >= 150) {
		output = fopen("output.txt", "a");
		fprintf(output, "\n[INFO]: Candidate objective value: %f\n", inst.objval);
		fprintf(output, "[INFO]: Candidate rounded x: ");
		for (int j = 0; j < inst.ncols; j++) fprintf(output, "%f ", inst.x[j]);
		fprintf(stdout, "\n");
		fclose(output);
	}

	// Check variable bounds and constraints for the candidate rounded solution
	print_verbose(100, "[INFO]: ... verifying variable bounds and constraints ...\n");
	check_bounds(&inst, inst.x);
	check_constraints(&inst, inst.x);
	print_verbose(100, "[INFO][OK]: Variable bounds and constraints satisfied.\n");

	// Check whether all integer/binary variables of the original MIP have been rounded
	print_verbose(100, "[INFO]: ... verifying whether all integer variables of the original MIP have been rounded ...\n");
	check_rounding(&inst);
	print_verbose(100, "[INFO][OK]: All integer/binary variables of the MIP have been rounded.\n");
	
	// Print objective value of rounded solution
	fprintf(stdout, "[FINAL RESULT]: Objective value of rounded solution: %f\n\n", inst.objval);

	// Print rounded solution to file
	if (VERBOSE >= 150) {
		output = fopen("output.txt", "a");
		fprintf(output, "[INFO][FINAL ROUNDED SOLUTION]: Rounded x: \n");
		for (int j = 0; j < inst.ncols; j++) fprintf(output, "%f ", inst.x[j]);
		fclose(output);
	}

	// Free
	free_inst(&inst);

	return EXIT_SUCCESS;
}