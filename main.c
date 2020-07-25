/**
 * @file main.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"
#include "dirent.h"

/**
 * @brief Main.
 * 
 * @param argc Number of command line arguments (includes the program name).
 * @param argv Command line arguments (strings).
 * @return Exit status.
 */
int main(int argc, char** argv) {

	instance inst; /**< General support structure. */
	
	init_inst(&inst);
	
	parse_cmd(argc, argv, &inst);

	strcmp(inst.input_file, "NULL")   ? test_instance(&inst) :
	strcmp(inst.input_folder, "NULL") ? test_folder(&inst) :
	print_error("Input file or folder required! See help.\n");

	free_inst(&inst);
	
	return EXIT_SUCCESS;
}

void test_instance(instance* inst) {

	LARGE_INTEGER lpfreq, lpstart, lpend; /**< Variables for measuring execution time for solving the initial continuous relaxation. */
	LARGE_INTEGER zifreq, zistart, ziend; /**< Variables for measuring execution time of ZI-Round. */
	LONG64 lp_solve_exec_time = 0;        /**< Execution time (in milliseconds) for solving the initial continuous relaxation. */
	LONG64 ziround_exec_time = 0;         /**< Execution time (in milliseconds) of ZI-Round. */
	int numrounds = 0;                    /**< Number of rounds (outer loops) of ZI-Round. */
		
	setup_CPLEX_env(inst);
	read_MIP_problem(inst, inst->input_file);
	save_integer_variables(inst);
	
	// Measure execution time (in milliseconds) for solving the continuous relaxation (time limit of 5 minutes)
	QueryPerformanceFrequency(&lpfreq);
	QueryPerformanceCounter(&lpstart);
	solve_continuous_relaxation(inst);
	QueryPerformanceCounter(&lpend);
	lp_solve_exec_time = (lpend.QuadPart - lpstart.QuadPart) * 1000 / lpfreq.QuadPart;

	populate_inst(inst);

	// Measure execution time (in milliseconds) of ZI-Round
	QueryPerformanceFrequency(&zifreq);
	QueryPerformanceCounter(&zistart);
	zi_round(inst, &numrounds);
	QueryPerformanceCounter(&ziend);
	ziround_exec_time = (ziend.QuadPart - zistart.QuadPart) * 1000 / zifreq.QuadPart;

	print_verbose(10, "[INFO]: ZI-Round terminated. #Rounds: %d\n", numrounds);
	print_verbose(10, "[INFO]: LP solve execution time (in milliseconds): %lld ms\n", lp_solve_exec_time);
	print_verbose(10, "[INFO]: ZI-Round execution time (in milliseconds): %lld ms\n", ziround_exec_time);
	print_verbose(10, "[INFO]: Sum of LP solve + ZI-Round execution time (in milliseconds): %lld ms\n", lp_solve_exec_time + ziround_exec_time);
	assert(equals(inst->solfrac, sol_fractionality(inst->x, inst->int_var, inst->ncols)));
	assert(equals(inst->objval, dot_product(inst->obj, inst->x, inst->ncols)));
	print_verbose(10, "[INFO]: Solution fractionality: %f\n", inst->solfrac);
	print_verbose(20, "[INFO]: Candidate objective value: %f\n", inst->objval);
	
	check_bounds(inst->x, inst->lb, inst->ub, inst->ncols);
	check_constraints(inst->x, inst->ncols, inst->nrows, inst->nzcnt, inst->rmatbeg,
		inst->rmatind, inst->rmatval, inst->sense, inst->rhs);

	if (VERBOSE >= 10) plot(inst);

	if (check_rounding(inst->x, inst->ncols, inst->int_var, inst->vartype)) {
		print_verbose(10, "[INFO]: All integer/binary variables of the MIP have been rounded.\n");
		fprintf(stdout, "[INFO]: Objective value of rounded solution: %f\n\n", inst->objval);
	}
	else print_verbose(10, "[INFO]: Failed to round all integer/binary variables of the MIP ...\n");
}

void test_folder(instance* inst) {

	LARGE_INTEGER lpfreq, lpstart, lpend; /**< Variables for measuring execution time for solving the initial continuous relaxation. */
	LARGE_INTEGER zifreq, zistart, ziend; /**< Variables for measuring execution time of ZI-Round. */
	LONG64 lp_solve_exec_time = 0;        /**< Execution time (in milliseconds) for solving the initial continuous relaxation. */
	LONG64 ziround_exec_time = 0;         /**< Execution time (in milliseconds) of ZI-Round. */
	int numrounds = 0;                    /**< Number of rounds (outer loops) of ZI-Round. */
	char* input_folder_name = NULL;       /**< Input folder name. */
	char* output_path = NULL;             /**< Test results output path. */

	// Allocate input folder name, test results output path and test-bed aggregate measures output path
	input_folder_name = (char*)calloc(30, sizeof(char));
	output_path = (char*)calloc(60, sizeof(char));
	if (input_folder_name == NULL || output_path == NULL) print_error("[test_folder]: Failed to allocate output related structures.\n");

	// Set file/folder names
	sprintf(input_folder_name, inst->input_folder);
	sprintf(output_path, "test_results_nogap.csv");

	// Initialize the directory and the directory element that represents a single file  
	DIR* dir = opendir(input_folder_name); if (dir == NULL) print_error("[test_folder]: Failed to open directory %s.\n", input_folder_name);
	struct dirent* direlem;

	// Print file headers
	FILE* output = fopen(output_path, "w");
	fprintf(output, "Instance;Seed;Cost;Fractionality;Rounds;LPtime(ms);ZItime(ms);SumLPZI(ms)\n");
	fclose(output);

	// Scan files
	while ((direlem = readdir(dir)) != NULL) {

		// Check whether the input filename (direlem->d_name) is a .mps file
		if (strlen(direlem->d_name) < 5 || strstr(direlem->d_name, ".mps") == NULL) continue;

		// Create a new instance (clone)
		instance test_inst;
		init_inst(&test_inst);
		sprintf(test_inst.input_file, "%s/%s", ((strcmp(input_folder_name, "NULL")) ? input_folder_name : "."), direlem->d_name);
		test_inst.singletons = inst->singletons;
		test_inst.shift_nonfracvars = inst->shift_nonfracvars;
		test_inst.timelimit = inst->timelimit;
		test_inst.rseed = inst->rseed;

		print_verbose(10, "TEST INSTANCE ------------------------------------------------------------------\n");
		print_verbose(10, "[] Instance name: %s\n", test_inst.input_file);
		print_verbose(10, "[] Use singletons: %d\n", test_inst.singletons);
		print_verbose(10, "[] Shift non-fractional integer variables: %d\n", test_inst.shift_nonfracvars);
		print_verbose(10, "[] Random seed: %d\n", test_inst.rseed);
		print_verbose(10, "--------------------------------------------------------------------------------\n");

		setup_CPLEX_env(&test_inst);
		read_MIP_problem(&test_inst, test_inst.input_file);
		save_integer_variables(&test_inst);

		// Measure execution time (in milliseconds) for solving the continuous relaxation (time limit of 5 minutes)
		QueryPerformanceFrequency(&lpfreq);
		QueryPerformanceCounter(&lpstart);
		solve_continuous_relaxation(&test_inst);
		QueryPerformanceCounter(&lpend);
		lp_solve_exec_time = (lpend.QuadPart - lpstart.QuadPart) * 1000 / lpfreq.QuadPart;

		populate_inst(&test_inst);

		// Measure execution time (in milliseconds) of ZI-Round
		QueryPerformanceFrequency(&zifreq);
		QueryPerformanceCounter(&zistart);
		zi_round(&test_inst, &numrounds);
		QueryPerformanceCounter(&ziend);
		ziround_exec_time = (ziend.QuadPart - zistart.QuadPart) * 1000 / zifreq.QuadPart;

		assert(equals(test_inst.solfrac, sol_fractionality(test_inst.x, test_inst.int_var, test_inst.ncols)));
		assert(equals(test_inst.objval, dot_product(test_inst.obj, test_inst.x, test_inst.ncols)));
		check_bounds(test_inst.x, test_inst.lb, test_inst.ub, test_inst.ncols);
		check_constraints(test_inst.x, test_inst.ncols, test_inst.nrows, test_inst.nzcnt, test_inst.rmatbeg,
						  test_inst.rmatind, test_inst.rmatval, test_inst.sense, test_inst.rhs);

		// if (VERBOSE >= 10) plot(&test_inst);

		// Print test results to file
		output = fopen(output_path, "a");
		fprintf(output, "%s;%d;%f;%f;%d;%lld;%lld;%lld\n", 
			strtok(direlem->d_name, "."), test_inst.rseed, test_inst.objval, test_inst.solfrac, numrounds, lp_solve_exec_time, ziround_exec_time, lp_solve_exec_time + ziround_exec_time);
		fclose(output);

		print_verbose(10, "TEST RESULT --------------------------------------------------------------------\n");
		print_verbose(10, "[] Solution cost: %.2f\n", test_inst.objval);
		print_verbose(10, "[] Solution fractionality: %.2f\n", test_inst.solfrac);
		print_verbose(10, "[] LP solve execution time (in milliseconds): %lld ms\n", lp_solve_exec_time);
		print_verbose(10, "[] ZI-Round execution time (in milliseconds): %lld ms\n", ziround_exec_time);
		print_verbose(10, "[] Sum of LP solve + ZI-Round execution time (in milliseconds): %lld ms\n", lp_solve_exec_time + ziround_exec_time);
		print_verbose(10, "--------------------------------------------------------------------------------\n\n\n");

		free_inst(&test_inst);
	} // end while
	
	free_all(2, output_path, input_folder_name);
	closedir(dir);
}