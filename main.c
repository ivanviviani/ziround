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
	instance inst; /**< General support structure. */
	int rounded;   /**< Support flag. */
	FILE* output;  /**< Output file pointer. */

	// Initialize
	rounded = 1;
	output = NULL;
	
	// Initialize instance
	init_inst(&inst);
	
	// Parse command line arguments
	parse_cmd(argc, argv, &inst);

	// Set up CPLEX environment
	setup_CPLEX_env(&inst);					  
	fprintf(stdout, "[INFO][OK]: CPLEX env setup.\n");

	// Read MIP problem from input file (.mps)
	read_MIP_problem(&inst, inst.input_file); 
	fprintf(stdout, "[INFO][OK]: MIP problem read.\n");

	// Read problem sizes
	read_problem_sizes(&inst);
	fprintf(stdout, "[INFO][OK]: Problem sizes read.\n");

	// Remember integer/binary variables of the original MIP
	save_integer_variables(&inst);
	fprintf(stdout, "[INFO][OK]: Integer variables saved.\n");
	
	// Solve continuous relaxation of the MIP problem
	solve_continuous_relaxation(&inst);
	fprintf(stdout, "[INFO][OK]: Continuous relaxation solved.\n");
	fprintf(stdout, "[INFO]: Continuous relaxation objective value: %.10g.\n", inst.objval);

	// Populate the instance with problem data retrieved from the CPLEX lp
	populate_inst(&inst);

	if (VERBOSE >= 200) print_problem_info(&inst, 1, 1);
	if (VERBOSE >= 150) {
		output = fopen("output.txt", "a");
		fprintf(output, "\n[INFO]: Continuous relaxation objective value: %.10g.\n", inst.objval);
		fprintf(output, "[INFO]: Solution of continuous relaxation: \n");
		for (int j = 0; j < inst.ncols; j++) fprintf(output, "%f ", inst.x[j]);
		fclose(output);
	}

	fprintf(stdout, "[INFO][OK]: Problem data read.\n[INFO]: ... Starting ZI-Round ...\n");
	
//******************************************* ZI-ROUND *******************************************
	zi_round(&inst);
//************************************************************************************************

	fprintf(stdout, "[INFO][OK]: ZI-Round terminated.\n");
	fprintf(stdout, "[INFO]: Candidate objective value: %f\n", inst.objval);

	if (VERBOSE >= 150) {
		output = fopen("output.txt", "a");
		fprintf(output, "\n[INFO]: Candidate objective value: %f\n", inst.objval);
		fprintf(output, "[INFO]: Candidate rounded x: ");
		for (int j = 0; j < inst.ncols; j++) fprintf(output, "%f ", inst.x[j]);
		fprintf(stdout, "\n");
		fclose(output);
	}

	// Verify variable bounds and constraints for the rounded candidate solution
	fprintf(stdout, "[INFO]: ... verifying variable bounds and constraints ...\n");
	check_bounds(&inst, inst.x);
	check_constraints(&inst, inst.x);
	fprintf(stdout, "[INFO][OK]: variable bounds and constraints satisfied.\n");

	// Verify whether all integer/binary variables of the original MIP have been rounded
	fprintf(stdout, "[INFO]: ... verifying whether all integer variables of the original MIP have been rounded ...\n");
	for (int j = 0; j < inst.ncols; j++) {
		if (inst.int_var[j] && is_fractional(inst.x[j])) {
			fprintf(stdout, "[INFO]: Integer variable x_%d = %f has not been rounded!\n", j + 1, inst.x[j]);
			rounded = 0;
		}
	}

	// Print final result (fail/success)
	if (!rounded) fprintf(stdout, "[INFO][FINAL RESULT]: ... Failed to round all integer/binary variables of the MIP ...\n"); 
	else fprintf(stdout, "[INFO][FINAL RESULT]: All integer/binary variables of the MIP have been rounded!\n");

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