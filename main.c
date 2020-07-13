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

	instance inst;       /**< General support structure. */
	FILE* output = NULL; /**< Output file pointer. */
	
	// Initialize instance
	init_inst(&inst);
	
	// Parse command line arguments
	parse_cmd(argc, argv, &inst);

	// Set up CPLEX environment
	setup_CPLEX_env(&inst);

	// Read MIP problem from input file (.mps)
	read_MIP_problem(&inst, inst.input_file);

	// Remember integer/binary variables of the original MIP
	save_integer_variables(&inst);

	if (VERBOSE >= 201) print_problem_info(&inst, 0, 1);
	
	// Solve continuous relaxation of the MIP problem
	solve_continuous_relaxation(&inst);

	// Print objective value of continuous relaxation
	print_verbose(100, "[INFO]: Continuous relaxation objective value: %.10g.\n", inst.objval);

	// Populate the instance with problem data retrieved from the CPLEX lp
	populate_inst(&inst);

	// Print problem info to file
	if (VERBOSE >= 201) print_problem_info(&inst, 1, 1);
	if (VERBOSE >= 201) {
		output = fopen("output.txt", "a");
		fprintf(output, "\n[INFO]: Continuous relaxation objective value: %.10g.\n", inst.objval);
		fprintf(output, "[INFO]: Solution of continuous relaxation: \n");
		for (int j = 0; j < inst.ncols; j++) fprintf(output, "%f ", inst.x[j]);
		fclose(output);
	}

	print_verbose(10, "[INFO]: ... Starting ZI-Round ...\n");
//******************************************* ZI-ROUND *******************************************
	zi_round(&inst);
//************************************************************************************************
	print_verbose(10, "[INFO]: ZI-Round terminated.\n");
	print_verbose(10, "[INFO]: Solution fractionality: %f.\n", sol_fractionality(inst.x, inst.int_var, inst.ncols));
	print_verbose(100, "[INFO]: Candidate objective value: %f\n", inst.objval);

	// Print candidate rounded solution and its objective value to file
	if (VERBOSE >= 201) {
		output = fopen("output.txt", "a");
		fprintf(output, "\n[INFO]: Candidate objective value: %f\n", inst.objval);
		fprintf(output, "[INFO]: Candidate rounded x: ");
		for (int j = 0; j < inst.ncols; j++) fprintf(output, "%f ", inst.x[j]);
		fprintf(stdout, "\n");
		fclose(output);
	}

	// Check variable bounds and constraints for the candidate rounded solution
	check_bounds(inst.x, inst.lb, inst.ub, inst.ncols);
	check_constraints(inst.x, inst.ncols, inst.nrows, inst.nzcnt, inst.rmatbeg, 
					  inst.rmatind, inst.rmatval, inst.sense, inst.rhs);

	// Check whether all integer/binary variables of the original MIP have been rounded
	if (check_rounding(inst.x, inst.ncols, inst.int_var, inst.vartype)) {
		print_verbose(10, "[INFO]: All integer/binary variables of the MIP have been rounded.\n");
		fprintf(stdout, "[INFO]: Objective value of rounded solution: %f\n\n", inst.objval);
	}
	else {
		print_verbose(10, "[INFO]: Failed to round all integer/binary variables of the MIP ...\n");
	}

	// Print rounded solution to file
	if (VERBOSE >= 201) {
		output = fopen("output.txt", "a");
		fprintf(output, "[INFO][FINAL ROUNDED SOLUTION]: Rounded x: \n");
		for (int j = 0; j < inst.ncols; j++) fprintf(output, "%f ", inst.x[j]);
		fclose(output);
	}

	// Free instance
	free_inst(&inst);

	return EXIT_SUCCESS;
}