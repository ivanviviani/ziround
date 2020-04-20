#include "ziround.h"

/**
 * @brief Main.
 * 
 * @param argc Number of command line arguments (includes the program name).
 * @param argv Command line arguments (strings).
 * @return Exit status.
 */
int main(int argc, char** argv) {

	int status = 0;
	instance inst;
	int rounded = 1;
	FILE* output = fopen("output.txt", "w");
	
	initialize_instance(&inst);
	
	parse_command_line(argc, argv, &inst);
	
	status = setup_CPLEX_env(&inst); if (status) { fprintf(stderr, "[ERR]: Error inside setup_CPLEX_env.\n"); goto TERMINATE; }
	fprintf(stdout, "[INFO][OK]: CPLEX env setup.\n");
	
	status = read_MIP_problem(&inst, inst.input_file); if (status) { fprintf(stderr, "[ERR]: Error inside read_problem.\n"); goto TERMINATE; }
	fprintf(stdout, "[INFO][OK]: MIP problem read.\n");
	
	read_problem_sizes(&inst);
	fprintf(stdout, "[INFO][OK]: Problem sizes read.\n");

	status = save_integer_variables(&inst); if (status) { fprintf(stderr, "[ERR]: Error inside save_integer_variables.\n"); goto TERMINATE; }
	fprintf(stdout, "[INFO][OK]: Integer variables saved.\n");
	
	status = solve_continuous_relaxation(&inst);
	if (status) { fprintf(stderr, "[ERR]: Error inside solve_continuous_relaxation.\n"); goto TERMINATE; }
	fprintf(stdout, "[INFO][OK]: Continuous relaxation solved.\n");

	status = read_problem_data(&inst);
	if (status) { fprintf(stderr, "[ERR]: Error inside read_problem_data.\n"); goto TERMINATE; }

	print_problem_info(&inst, 1, 1); // flags solution available (1) and print problem to file (1)

	fprintf(output, "[INFO]: Solution of continuous relaxation: ");
	for (int j = 0; j < inst.ncols; j++) fprintf(output, "%f ", inst.x[j]);
	fprintf(stdout, "\n[INFO]: Continuous relaxation objective value: %.10g.\n", inst.objval);
	fprintf(output, "\n[INFO]: Continuous relaxation objective value: %.10g.\n", inst.objval);
	fprintf(stdout, "[INFO][OK]: Problem data read.\n[INFO]: ... Starting ZI-Round ...\n");
	
	//******************************************* ZI-ROUND *******************************************************************************************
	status = zi_round(&inst);
	if (status) { fprintf(stderr, "[ERR]: Error inside zi_round.\n"); goto TERMINATE; }
	//************************************************************************************************************************************************

	// Print candidate solution found
	fprintf(output, "[INFO]: Candidate rounded x: ");
	for (int j = 0; j < inst.ncols; j++) fprintf(output, "%f ", inst.x[j]);
	fprintf(stdout, "\n");
	fprintf(stdout, "[INFO]: Candidate objective value: %f\n", inst.objval);
	fprintf(output, "\n\n\n\n[INFO]: Candidate objective value: %f\n", inst.objval);

	// Verify variable bounds and constraints
	fprintf(stdout, "[INFO]: ... verifying variable bounds and constraints ...\n");
	status = check_bounds(&inst, inst.x);
	if (status) { fprintf(stderr, "[ERR]: Error inside check_bounds.\n"); goto TERMINATE; }
	status = check_constraints(&inst, inst.x);
	if (status) { fprintf(stderr, "[ERR]: Error inside check_constraints.\n"); goto TERMINATE; }

	// Verify whether all integer variables of the original MIP have been rounded
	fprintf(stdout, "[INFO]: ... verifying whether all integer variables of the original MIP have been rounded ...\n");
	for (int j = 0; j < inst.ncols; j++) {
		if (inst.int_var[j] && is_fractional(inst.x[j])) {
			fprintf(stdout, "[INFO]: Integer variable x_%d = %f has not been rounded!\n", j + 1, inst.x[j]);
			rounded = 0;
		}
	}

	// Print final result
	if (!rounded) { fprintf(stdout, "[INFO][FINAL RESULT]: All the integer variables have NOT been rounded! :(\n"); goto TERMINATE; }
	fprintf(stdout, "[INFO][FINAL ROUNDED SOLUTION]: Rounded x in output.txt");
	fprintf(output, "[INFO][FINAL ROUNDED SOLUTION]: Rounded x: ");
	for (int j = 0; j < inst.ncols; j++) fprintf(output, "%f ", inst.x[j]);
	fprintf(stdout, "\n\n");

	/*
		Prova altri esempi a mano e sicuramente la MIPLIB 2003
	*/

TERMINATE:
	fclose(output);
	status = free_instance(&inst); if (status) fprintf(stderr, "[ERR]: Error inside free_instance.\n");
	return status;
}