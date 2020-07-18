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
	
	// Initialize instance
	init_inst(&inst);
	
	// Parse command line arguments
	parse_cmd(argc, argv, &inst);

	// [] Test a single instance
	if (strcmp(inst.input_file, "NULL")) {

		test_instance(&inst);
	}
	// [] Test all instances in a folder
	else if (strcmp(inst.input_folder, "NULL")) {

		test_folder(&inst);
	}
	// [] More parameters are required
	else {
		print_error("Input file or folder required! See help.\n");
	}

	// Free instance
	free_inst(&inst);
	
	return EXIT_SUCCESS;
}

void test_instance(instance* inst) {

	FILE* output = NULL; /**< Output file pointer. */
	int numrounds = 0;   /**< Number of rounds (outer loops) of ZI-Round. */

	// Set up CPLEX environment
	setup_CPLEX_env(inst);

	// Read MIP problem from input file (.mps)
	read_MIP_problem(inst, inst->input_file);

	// Remember integer/binary variables of the original MIP
	save_integer_variables(inst);

	// [DEBUG ONLY]: Print problem info to file
	if (VERBOSE >= 201) print_problem_info(inst, 0, 1);

	// Solve continuous relaxation of the MIP problem
	solve_continuous_relaxation(inst);

	// Print objective value of continuous relaxation
	print_verbose(100, "[INFO]: Continuous relaxation objective value: %.10g.\n", inst->objval);

	// Populate the instance with problem data retrieved from the CPLEX lp
	populate_inst(inst);

	// Print problem info to file
	if (VERBOSE >= 201) print_problem_info(inst, 1, 1);
	if (VERBOSE >= 201) {
		output = fopen("output.txt", "a");
		fprintf(output, "\n[INFO]: Continuous relaxation objective value: %.10g.\n", inst->objval);
		fprintf(output, "[INFO]: Solution of continuous relaxation: \n");
		for (int j = 0; j < inst->ncols; j++) fprintf(output, "%f ", inst->x[j]);
		fclose(output);
	}

	print_verbose(10, "[INFO]: ... Starting ZI-Round ...\n");
	//******************************************* ZI-ROUND *******************************************
	zi_round(inst, &numrounds);
	//************************************************************************************************
	print_verbose(10, "[INFO]: ZI-Round terminated. #Rounds: %d\n", numrounds);
	print_verbose(10, "[INFO]: Solution fractionality: %f.\n", sol_fractionality(inst->x, inst->int_var, inst->ncols));
	print_verbose(100, "[INFO]: Candidate objective value: %f\n", inst->objval);

	// [DEBUG ONLY]: Print candidate rounded solution and its objective value to file
	if (VERBOSE >= 201) {
		output = fopen("output.txt", "a");
		fprintf(output, "\n[INFO]: Candidate objective value: %f\n", inst->objval);
		fprintf(output, "[INFO]: Candidate rounded x: ");
		for (int j = 0; j < inst->ncols; j++) fprintf(output, "%f ", inst->x[j]);
		fprintf(stdout, "\n");
		fclose(output);
	}

	// Check variable bounds and constraints for the candidate rounded solution
	check_bounds(inst->x, inst->lb, inst->ub, inst->ncols);
	check_constraints(inst->x, inst->ncols, inst->nrows, inst->nzcnt, inst->rmatbeg,
		inst->rmatind, inst->rmatval, inst->sense, inst->rhs);

	// Plot trackers
	if (VERBOSE >= 10) plot(inst);

	// Check whether all integer/binary variables of the original MIP have been rounded
	if (check_rounding(inst->x, inst->ncols, inst->int_var, inst->vartype)) {
		print_verbose(10, "[INFO]: All integer/binary variables of the MIP have been rounded.\n");
		fprintf(stdout, "[INFO]: Objective value of rounded solution: %f\n\n", inst->objval);
	}
	else {
		print_verbose(10, "[INFO]: Failed to round all integer/binary variables of the MIP ...\n");
	}

	// [DEBUG ONLY]: Print rounded solution to file
	if (VERBOSE >= 201) {
		output = fopen("output.txt", "a");
		fprintf(output, "[INFO][FINAL ROUNDED SOLUTION]: Rounded x: \n");
		for (int j = 0; j < inst->ncols; j++) fprintf(output, "%f ", inst->x[j]);
		fclose(output);
	}
}

void test_folder(instance* inst) {

	time_t start, exec_time;   /**< Execution time variables. */
	double solfrac = LONG_MIN; /**< Solution fractionality. */
	int numrounds = 0;         /**< Number of rounds (outer loops) of ZI-Round. */

	// Allocate input folder and output path
	char* input_folder_name = (char*)calloc(30, sizeof(char));
	char* output_path = (char*)calloc(60, sizeof(char)); if (input_folder_name == NULL || output_path == NULL) print_error("[test_folder]: Failed to allocate input filename, folder or output path.\n");
	sprintf(input_folder_name, "%s", inst->input_folder);

	// Set output file name
	sprintf(output_path, "ziround_test_results.csv");

	// Initialize the directory and the directory element that represents a single file  
	DIR* dir = opendir(input_folder_name); if (dir == NULL) print_error("[test_folder]: Failed to open directory %s.\n", input_folder_name);
	struct dirent* direlem;

	// Print file header
	FILE* output = fopen(output_path, "w");
	fprintf(output, "Instance;Cost;Fractionality;#Rounds;Time\n");
	fclose(output);

	// Scan files
	while ((direlem = readdir(dir)) != NULL) {

		// Check whether the input filename (direlem->d_name) is a .mps file
		if (strlen(direlem->d_name) < 5 || strstr(direlem->d_name, ".mps") == NULL) continue;

		// Create a new instance (clone)
		instance test_inst;
		init_inst(&test_inst);
		sprintf(test_inst.input_file, "%s/%s", ((strcmp(input_folder_name, "NULL")) ? input_folder_name : "."), direlem->d_name);
		test_inst.extension = inst->extension;
		test_inst.timelimit = inst->timelimit;
		test_inst.rseed = inst->rseed;

		print_verbose(10, "TEST INSTANCE ------------------------------------------------------------------\n");
		print_verbose(10, "[] Instance name: %s\n", test_inst.input_file);
		print_verbose(10, "[] Extension flag: %d\n", test_inst.extension);
		print_verbose(10, "[] Random seed: %d\n", test_inst.rseed);
		print_verbose(10, "--------------------------------------------------------------------------------\n");

		// Test
		setup_CPLEX_env(&test_inst);
		read_MIP_problem(&test_inst, test_inst.input_file);
		save_integer_variables(&test_inst);
		solve_continuous_relaxation(&test_inst);
		populate_inst(&test_inst);
		//******************************************* ZI-ROUND *******************************************
		start = time(NULL);
		zi_round(&test_inst, &numrounds);
		exec_time = time(NULL) - start;
		//************************************************************************************************
		check_bounds(test_inst.x, test_inst.lb, test_inst.ub, test_inst.ncols);
		check_constraints(test_inst.x, test_inst.ncols, test_inst.nrows, test_inst.nzcnt, test_inst.rmatbeg,
						  test_inst.rmatind, test_inst.rmatval, test_inst.sense, test_inst.rhs);

		// Compute final solution fractionality
		solfrac = sol_fractionality(test_inst.x, test_inst.int_var, test_inst.ncols);

		// Plot trackers
		// if (VERBOSE >= 10) plot(&test_inst);

		// Print results to file
		output = fopen(output_path, "a");
		fprintf(output, "%s;%f;%f;%d;%d\n", strtok(direlem->d_name, "."), test_inst.objval, solfrac, numrounds, (long)exec_time);
		fclose(output);

		// Print CPLEX solution cost and execution time
		print_verbose(10, "TEST RESULT --------------------------------------------------------------------\n");
		print_verbose(10, "[] Solution cost: %f\n", test_inst.objval);
		print_verbose(10, "[] Solution fractionality: %f\n", solfrac);
		print_verbose(10, "[] Execution time: %d s\n", exec_time);
		print_verbose(10, "--------------------------------------------------------------------------------\n\n\n");

		// Free test instance
		free_inst(&test_inst);
	} // end while

	// Free
	free_all(2, output_path, input_folder_name);
	// Close directory
	closedir(dir);
}

void plot(instance* inst) {

	// Plot trackers of solution fractionality, cost or number of rounded variables
	if (PLOT_SOL_FRAC && PLOT_NUM_VARS_TO_ROUND) {
		char** label = (char**)calloc(3, sizeof(char*));
		label[0] = (char*)calloc(6, sizeof(char));
		label[1] = (char*)calloc(30, sizeof(char));
		label[2] = (char*)calloc(30, sizeof(char));
		char** name = (char**)calloc(2, sizeof(char*));
		name[0] = (char*)calloc(20, sizeof(char));
		name[1] = (char*)calloc(20, sizeof(char));
		char* temp = (char*)calloc(20, sizeof(char));
		sprintf(label[0], "Round");
		sprintf(label[1], "Solution Fractionality (SF)");
		sprintf(label[2], "#Variables to Round (#VR)");
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name[0], "%s (SF)", strtok(NULL, "/."));
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name[1], "%s (#VR)", strtok(NULL, "/."));
		assert(inst->size_frac == inst->size_toround);
		plot_tracker_pair(inst->tracker_sol_frac, inst->tracker_toround, name, label, inst->size_frac, NULL);
		free_all(8, label[0], label[1], label[2], label, name[0], name[1], name, temp);
	}
	else if (PLOT_SOL_FRAC) {
		char** label = (char**)calloc(2, sizeof(char*));
		label[0] = (char*)calloc(6, sizeof(char));
		label[1] = (char*)calloc(30, sizeof(char));
		char* name = (char*)calloc(20, sizeof(char));
		char* temp = (char*)calloc(20, sizeof(char));
		sprintf(label[0], "Round");
		sprintf(label[1], "Solution Fractionality");
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name, strtok(NULL, "/."));
		plot_tracker(inst->tracker_sol_frac, name, label, inst->size_frac, NULL);
		free_all(5, label[0], label[1], label, name, temp);
	}
	else if (PLOT_NUM_VARS_TO_ROUND) {
		char** label = (char**)calloc(2, sizeof(char*));
		label[0] = (char*)calloc(6, sizeof(char));
		label[1] = (char*)calloc(30, sizeof(char));
		char* name = (char*)calloc(20, sizeof(char));
		char* temp = (char*)calloc(20, sizeof(char));
		sprintf(label[0], "Round");
		sprintf(label[1], "#Variables to Round");
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name, strtok(NULL, "/."));
		plot_tracker(inst->tracker_toround, name, label, inst->size_toround, NULL);
		free_all(5, label[0], label[1], label, name, temp);
	}
	
	if (PLOT_SOL_COST) {
		char** label = (char**)calloc(2, sizeof(char*));
		label[0] = (char*)calloc(6, sizeof(char));
		label[1] = (char*)calloc(30, sizeof(char));
		char* name = (char*)calloc(20, sizeof(char));
		char* temp = (char*)calloc(20, sizeof(char));
		sprintf(label[0], "Round");
		sprintf(label[1], "Solution Cost");
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name, strtok(NULL, "/."));
		plot_tracker(inst->tracker_sol_cost, name, label, inst->size_cost, NULL);
		free_all(5, label[0], label[1], label, name, temp);
	}
}