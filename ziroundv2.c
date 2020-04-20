/**
 * @file ziroundv2.c
 * @author Ivan Viviani
 * @brief main file for testing the ZI-Round MIP rounding heuristic.
 */
/*  README.md
	ZI-Round MIP Rounding Heuristic (version 2)

	Usage: ziroundv2 --input-mip [filename.mps] --obj-sense [max/min]

	Variables:

		Primal feasible point x = [x1, x2, ..., xj, ...]
		Upper bounds of all xj
		Lower bounds of all xj
		Slacks si
		Constraint matrix coefficients aij (both by rows and by columns)
		Constraints senses {'L','E','G'}
		Objective coefficients
		Problem type {CPX_MAX,CPX_MIN}
		Threshold epsilon = 0.00001

	Formulas:

		Fractionality of xj: ZI(xj) = min{xj - floor(xj), ceil(xj) - xj}
		Integer infeasibility of x[]: ZI(x) = sum{ZI(xj)}
		For 'L' (<=) constraints: (si non-negative)
			ub_j = min_i{si/aij : aij > 0}
			lb_j = min_i{-si/aij : aij < 0}
		For 'G' (>=) constraints: (si non-positive)
			ub_j = min_i{si/aij : aij < 0}
			lb_j = min_i{-si/aij : aij > 0}
		UBj = min{ub_j, ub(xj) - xj}
		LBj = min{lb_j, xj - lb(xj)}

	ZI-Round Algorithm (version 2):
		Input: Primal feasible point x

		Repeat
		  Loop For each integer variable xj
			If xj non-fractional
			  Then
				Calculate UBj,LBj,1
				If cj > 0 and LBj = 1 or cj < 0 and UB = 1
				  Then Update xj to improve objective and Update slacks
			If xj fractional
			  Then
				Calculate UBj,LBj,epsilon
				If ZI(xj + UBj) = ZI(xj - LBj) and ZI(xj + UBj) < ZI(xj)
				  Then Update xj to improve objective and Update slacks
				Else If ZI(xj + UBj) < ZI(xj - LBj) and ZI(xj + UBj) < ZI(xj)
				  Then Set xj = xj + UBj and Update slacks
				Else If ZI(xj - LBj) < ZI(xj + UBj) and ZI(xj - LBj) < ZI(xj)
				  Then Set xj = xj - LBj and Update slacks
		Until No more updates found

		Detail: stop calculating UBj,LBj if they both fall below epsilon
*/
#include <cplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * @brief Verbosity level.
 */
#define VERBOSE 5

/**
 * @brief Tolerance for non-integer numbers as considered by CPLEX.
 */
#define TOLERANCE 1e-6

/**
 * @brief Threshold used in the ZI-Round heuristic.
 */
#define EPSILON 1e-5

/**
 * @brief Structure holding local and global information about a
 * 	      problem instance, parameters included.
 */
typedef struct {
	// Variables
	int cur_numrows; 		/**< Number of rows of the coefficients matrix. */
	int cur_numcols; 		/**< Number of variables, also columns of the coefficients matrix. */
	double* x; 				/**< Current problem solution. Will be rounded. */
	double* obj; 			/**< Objective function coefficients. */
	double* lb; 			/**< Variable lower bounds. */
	double* ub; 			/**< Variable upper bounds. */
	double* slack; 			/**< Row (constraint) slacks, defined as right hand side minus row activity. */
	double objval; 			/**< Current objective value (for current problem solution). */
	int objsen; 			/**< Objective function sense, CPX_MIN (default) or CPX_MAX (specified from command line). */
	char* mip_ctype; 		/**< Variable types (before converting MIP to LP), integer/binary or continuous. */
	int* int_var; 			/**< Flags array that keeps track of integer/binary (value 1) and continuous (value 0) variables. */
	// Coefficient matrix
	int* rmatbeg; 			/**< Begin row indices of non-zero coefficients for rmatind and rmatval. */
	int* rmatind; 			/**< Column indices of non-zero coefficients. */
	double* rmatval; 		/**< Non-zero coefficients (row major). */
	int* cmatbeg; 			/**< Begin column indices of non-zero coefficients for cmatind and cmatval. */
	int* cmatind; 			/**< Row indices of non-zero coefficients. */
	double* cmatval; 		/**< Non-zero coefficients (column major). */
	char* sense; 			/**< Constraint (row) senses, 'L' (<=) or 'G' (>=) or 'E' (=). */
	double* rhs; 			/**< Constraint right hand sides (rhs). */
	// ZI-Round variables
	double* UB; 			/**< Maximum variable up-shifts. */
	double* LB; 			/**< Maximum variable down-shifts. */
	double ZI; 				/**< Fractionality of a variable (used in function zi_round). */
	double ZIplus; 			/**< Fractionality of a shifted up variable (used in function zi_round). */
	double ZIminus; 		/**< Fractionality of a shifted down variable (used in fucntion zi_round). */
	double obj_plusUBj; 	/**< Objective value for a shifted up variable (used in function zi_round). */
	double obj_minusLBj; 	/**< Objective value for a shifted down variable (used in function zi_round). */
	// double* x_prev; 		/**< Solution before a single rounding (up/down shifting). */
	// double* x_updated; 	/**< Support copy needed to calculate obj_plusUBj and obj_minusLBj. */
	int updated; 			/**< Flag set to 1 whenever a solution has been rounded (in multiple rounds). */
	// Parameters
	int status; 			/**< Error status flag set to 1 whenever any error occurs. */
	CPXENVptr env; 			/**< CPLEX environment pointer. */
	CPXLPptr lp; 			/**< CPLEX lp pointer. */
	int solnstat; 			/**< Solution status according to CPLEX. */
	int solnmethod; 		/**< Solution method according to CPLEX. */
	int solntype; 			/**< Solution type according to CPLEX. */
	int nzcnt; 				/**< Number of non-zero coefficients. */
	int surplus; 			/**< CPLEX parameter. */
	char* input_file; 		/**< Input filename (mps format, specified from command line). */
} instance;

/**
 * @brief Initialize the appropriate fields of the instance.
 * 
 * @param inst Pointer to the instance.
 */
void initialize_instance(instance* inst);

/**
 * @brief Parse the arguments from the command line and populate the
 *        appropriate fields of the already initialized instance.
 * 
 * @param argc Number of command line arguments (includes the program name).
 * @param argv Program arguments (strings).
 * @param inst Pointer to the instance.
 */
void parse_command_line(int argc, char** argv, instance* inst);

/**
 * @brief Setup the CPLEX environment for the problem represented by the instance.
 * 		  Also turn on CPLEX screen output.
 * 
 * @param inst Pointer to the already initialized and populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int setup_CPLEX_env(instance* inst);

/**
 * @brief Create the lp into the CPLEX env of the instance and populate the lp
 * 		  with problem data read from a file (mps format).
 * 
 * @param inst Pointer to the already populated instance.
 * @param filename Name of the input file (mps format).
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int read_MIP_problem(instance* inst, char* filename);

/**
 * @brief Read the current number of rows and columns of the coefficients matrix
 * 	 	  and populate the corresponding fields of the instance.
 * 
 * @param inst Pointer to the already populated instance.
 */
void read_problem_sizes(instance* inst);

/**
 * @brief Read and save variable types from the MIP and save the integer/binary ones in a
 * 		  flags array. Both arrays are fields of the instance.
 * 
 * @param inst Pointer to the already populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int save_integer_variables(instance* inst);

/**
 * @brief Change the problem type from MIP to LP and solve its continuous relaxation.
 * 
 * @param inst Pointer to the already populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int solve_continuous_relaxation(instance* inst);

/**
 * @brief Read the continuous relaxation solution from the CPLEX lp
 *        and populate the corresponding array field of the instance.
 * 
 * @param inst Pointer to the already populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int read_solution(instance* inst);

/**
 * @brief Read the variable bounds from the CPLEX lp
 * 		  and populate the corresponding array fields of the instance.
 * 
 * @param inst Pointer to the already populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int read_variable_bounds(instance* inst);

/**
 * @brief Read the objective value from the CPLEX lp
 * 		  and populate the corresponding field of the instance.
 * 
 * @param inst Pointer to the already populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int read_objective_value(instance* inst);

/**
 * @brief Read the objective function coefficients from the CPLEX lp
 * 		  and populate the corresponding array field of the instance.
 * 
 * @param inst Pointer to the already populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int read_objective_coefficients(instance* inst);

/**
 * @brief Read the non-zero constraint coefficients both by row and column
 * 		  from the CPLEX lp, according to the data structures used by CPLEX,
 * 	      and populate the corresponding array fields of the instance.
 * 
 * @param inst Pointer to the already populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int read_constraints_coefficients(instance* inst);

/**
 * @brief Read the constraint senses from the CPLEX lp,
 * 		  and populate the corresponding array field of the instance.
 * 
 * @param inst Pointer to the already populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int read_constraints_senses(instance* inst);

/**
 * @brief Read the constraint right hand sides from the CPLEX lp,
 * 		  and populate the corresponding array field of the instance.
 * 
 * @param inst Pointer to the already populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int read_constraints_right_hand_sides(instance* inst);

/**
 * @brief Read the constraint (row) slacks for the continuous relaxation
 * 		  solution from the CPLEX lp, and populate the corresponding
 * 		  array field of the instance. They will then be updated during
 * 		  the rounding phase of zi_round.
 * 
 * @param inst Pointer to the already populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int read_row_slacks(instance* inst);

/**
 * @brief Read the problem data from the CPLEX lp using all the read functions,
 * 		  and populate the corresponding fields of the instance.
 * 
 * @param inst Pointer to the already populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int read_problem_data(instance* inst);

/**
 * @brief Use the ZI-Round heuristic to round the continuous relaxation solution
 * 	      until an integer feasible solution is found or the procedure terminates
 * 		  unsuccessfully.
 * 
 * @param inst Pointer to the already populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int zi_round(instance* inst);

/**
 * @brief Update variable xj to improve objective, according to the values
 * of UBj and LBj calculated before. Also update the row slacks.
 * 
 * @param inst Pointer to the already populated instance.
 * @param j Column (variable xj) index.
 * @param is_fractional Flag that indicates whether xj is fractional or integer.
 */
void update_xj_to_improve_objective(instance* inst, int j, int is_fractional);

/**
 * @brief Calculate the activity of the row of index i, for the solution x.
 * 
 * @param inst Pointer to the already populated instance.
 * @param i Row index.
 * @param x Solution to be used for calculating the row activity.
 * @return The row activity.
 */
double row_activity(instance* inst, int i, double* x);

/**
 * @brief Update the slack array field of the instance (incrementally) 
 *        after a rounding of the current solution (x) field of the instance.
 *
 * @details Whenever it is called, only one variable xj has been updated.
 * 
 * @param inst Pointer to the already populated instance.
 * @param j Index of the (only) variable just updated.
 * @param signed_delta Signed delta xj.
 */
void update_slacks(instance* inst, int j, double signed_delta);

/**
 * @brief Update the objective value (incrementally) for the current updated solution x.
 *
 * @details Whenever it is called, only one variable xj has been updated.
 * 
 * @param inst Pointer to the already populated instance.
 * @param j Index of the (only) variable just updated.
 * @param signed_delta Signed delta xj.
 */
void update_objective_value(instance* inst, int j, double signed_delta);

/**
 * @brief 
 * 
 * @param inst 
 * @param j 
 * @param epsilon 
 */
void calculate_UBjLBj(instance* inst, int j, const double epsilon);

/**
 * @brief Check that all the variable bounds are satisfied, for the
 * 		  given solution x.
 * 
 * @param inst Pointer to the already populated instance.
 * @param x Solution to be used for evaluating bounds satisfiability.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int check_bounds(instance* inst, double* x);

/**
 * @brief Check that all the constraints are satisfied, for the
 * 		  given solution x.
 * 
 * @param inst Pointer to the already populated instance.
 * @param x Solution to be used for evaluating constraints satisfiability.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int check_constraints(instance* inst, double* x);

/**
 * @brief Deallocate all the fields of the instance.
 * 
 * @param inst Pointer to the already populated instance.
 * @return Error status (1 if any error occured, 0 otherwise).
 */
int free_instance(instance* inst);

/**
 * @brief [DEBUG] Print problem information to standard output or file.
 * 
 * @param inst Pointer to the already populated instance.
 * @param sol_available Flag set to 1 if the instance contains a solution, 0 otherwise.
 * @param to_file Flag set to 1 if the problem info is to be printed to file, 0 otherwise.
 */
void print_problem_info(instance* inst, int sol_available, int to_file);

/**
 * @brief Deallocate a resource (previously allocated memory area).
 * 
 * @param ptr Pointer to the resource to deallocate.
 */
void free_and_null(char** ptr);

/**
 * @brief Calculate the fractionality of the given value.
 * 
 * @param xj Value.
 * @return Fractionality of the value.
 */
double fractionality(double xj);

/**
 * @brief Check the integrality of the given value.
 * 
 * @param num Value.
 * @return 1 if the value is fractional, 0 otherwise.
 */
int is_fractional(double num);

/**
 * @brief Calculate the dot product of two arrays (of the same length).
 * 
 * @param coef First array.
 * @param var_value Second array.
 * @param len Length of both arrays.
 * @return Dot product of the two arrays.
 */
double dot_product(double* coef, double* var_value, int len);

/**
 * @brief Check whether two arrays (of the same length) are equal, 
 * 		  according to an integer threshold (see macros).
 * 
 * @param prev First array.
 * @param new Second array.
 * @param len Length of both arrays
 * @return 1 if the two arrays are different, 0 otherwise.
 */
int different_arr(double* prev, double* new, int len);

/**
 * @brief Clone an array into a copy (both already allocated).
 * 
 * @param arr The array to be cloned.
 * @param clo Clone copy of the given array.
 * @param len Length of the array.
 */
void clone_array(double* arr, double* clo, int len);

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
	for (int j = 0; j < inst.cur_numcols; j++) fprintf(output, "%f ", inst.x[j]);
	fprintf(stdout, "\n[INFO]: Continuous relaxation objective value: %.10g.\n", inst.objval);
	fprintf(output, "\n[INFO]: Continuous relaxation objective value: %.10g.\n", inst.objval);
	fprintf(stdout, "[INFO][OK]: Problem data read.\n[INFO]: ... Starting ZI-Round ...\n");
	
	//******************************************* ZI-ROUND *******************************************************************************************
	status = zi_round(&inst);
	if (status) { fprintf(stderr, "[ERR]: Error inside zi_round.\n"); goto TERMINATE; }
	//************************************************************************************************************************************************

	// Print candidate solution found
	fprintf(output, "[INFO]: Candidate rounded x: ");
	for (int j = 0; j < inst.cur_numcols; j++) fprintf(output, "%f ", inst.x[j]);
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
	for (int j = 0; j < inst.cur_numcols; j++) {
		if (inst.int_var[j] && is_fractional(inst.x[j])) {
			fprintf(stdout, "[INFO]: Integer variable x_%d = %f has not been rounded!\n", j + 1, inst.x[j]);
			rounded = 0;
		}
	}

	// Print final result
	if (!rounded) { fprintf(stdout, "[INFO][FINAL RESULT]: All the integer variables have NOT been rounded! :(\n"); goto TERMINATE; }
	fprintf(stdout, "[INFO][FINAL ROUNDED SOLUTION]: Rounded x in output.txt");
	fprintf(output, "[INFO][FINAL ROUNDED SOLUTION]: Rounded x: ");
	for (int j = 0; j < inst.cur_numcols; j++) fprintf(output, "%f ", inst.x[j]);
	fprintf(stdout, "\n\n");

	/*
		Prova altri esempi a mano e sicuramente la MIPLIB 2003
	*/

TERMINATE:
	fclose(output);
	status = free_instance(&inst); if (status) fprintf(stderr, "[ERR]: Error inside free_instance.\n");
	return status;
}

void initialize_instance(instance* inst) {
	inst->x = NULL; inst->obj = NULL; inst->lb = NULL; inst->ub = NULL;
	inst->slack = NULL; inst->objsen = CPX_MIN; inst->mip_ctype = NULL; 
	inst->int_var = NULL; inst->rmatbeg = NULL; inst->rmatind = NULL; 
	inst->rmatval = NULL; inst->cmatbeg = NULL; inst->cmatind = NULL; 
	inst->cmatval = NULL; inst->sense = NULL; inst->rhs = NULL; inst->UB = NULL;
	inst->LB = NULL; inst->updated = 0; // inst->x_prev = NULL; inst->x_updated = NULL; 
	inst->status = 0; inst->env = NULL; inst->lp = NULL; inst->input_file = NULL;
}

void parse_command_line(int argc, char** argv, instance* inst) {
	// Parse command line arguments
	int help = (argc < 2) ? 1 : 0;
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--input-mip") == 0) { inst->input_file = argv[++i]; continue; }
		// if (strcmp(argv[i], "--obj-sense") == 0) { inst->objsen = (strcmp(argv[++i], "min") == 0) ? CPX_MIN : CPX_MAX; continue; }
		if (strcmp(argv[i], "--help") == 0) { help = 1; continue; }
		if (strcmp(argv[i], "-help") == 0) { help = 1; continue; }
		help = 1;
	}
	// Print chosen parameters
	if (help || VERBOSE > 10) {
		fprintf(stdout, "\n\n**** Parameters ***************************\n");
		fprintf(stdout, "**   --input-mip %s\n", inst->input_file);
		// fprintf(stdout, "**   --obj-sense %s\n", ((inst->objsen == CPX_MIN) ? "min" : "max"));
		fprintf(stdout, "*******************************************\n");
	}
	if (help) exit(EXIT_FAILURE);
}

int setup_CPLEX_env(instance* inst) {
	inst->env = CPXopenCPLEX(&(inst->status));
	if (inst->env == NULL) {
		char errmsg[CPXMESSAGEBUFSIZE];
		fprintf(stderr, "[ERR][setup_CPLEX_env]: Could not open CPLEX environment.\n");
		CPXgeterrorstring(inst->env, inst->status, errmsg);
		fprintf(stderr, "[ERR][setup_CPLEX_env]: %s", errmsg); return inst->status;
	}
	// Set CPLEX parameters
	inst->status = CPXsetintparam(inst->env, CPXPARAM_ScreenOutput, CPX_ON);
	if (inst->status) { fprintf(stderr, "[ERR][setup_CPLEX_env]: Failed to turn on screen indicator, error %d.\n", inst->status); return inst->status; }
	inst->status = CPXsetdblparam(inst->env, CPXPARAM_TimeLimit, 3600);
	if (inst->status) { fprintf(stderr, "[ERR][setup_CPLEX_env]: Failed to set time limit, error %d.\n", inst->status); return inst->status; }
	return inst->status;
}

int read_MIP_problem(instance* inst, char* filename) {
	inst->lp = CPXcreateprob(inst->env, &(inst->status), filename);
	if (inst->lp == NULL) { fprintf(stderr, "[ERR][read_MIP_problem]: Failed to create MIP.\n"); return inst->status; }
	inst->status = CPXreadcopyprob(inst->env, inst->lp, filename, NULL);
	if (inst->status) { fprintf(stderr, "[ERR][read_MIP_problem]: Failed to read and copy the problem data.\n"); return inst->status; }
	// Change problem type {CPX_MAX, CPX_MIN} as specified from command line
	// No. Let CPLEX handle it.
	// CPXchgobjsen(inst->env, inst->lp, inst->objsen);
	return inst->status;
}

void read_problem_sizes(instance* inst) {
	inst->objsen = CPXgetobjsen(inst->env, inst->lp);
	inst->cur_numrows = CPXgetnumrows(inst->env, inst->lp);
	inst->cur_numcols = CPXgetnumcols(inst->env, inst->lp);
}

int save_integer_variables(instance* inst) {
	// Get MIP variable types {CPX_CONTINUOUS, CPX_BINARY, CPX_INTEGER, CPX_SEMICONT, CPX_SEMIINT}
	inst->mip_ctype = (char*)malloc(inst->cur_numcols * sizeof(char));
	if (inst->mip_ctype == NULL) { fprintf(stderr, "[ERR][save_integer_variables]: Failed to allocate mip_ctype.\n"); return 1; }
	inst->status = CPXgetctype(inst->env, inst->lp, inst->mip_ctype, 0, inst->cur_numcols - 1);
	if (inst->status) { fprintf(stderr, "[ERR][save_integer_variables]: Failed to obtain MIP variable types.\n"); return inst->status; }
	// Remember integer variables {CPX_BINARY, CPX_INTEGER}
	inst->int_var = (int*)malloc(inst->cur_numcols * sizeof(int));
	if (inst->int_var == NULL) { fprintf(stderr, "[ERR][save_integer_variables]: Failed to allocate int_var.\n"); return 1; }
	for (int j = 0; j < inst->cur_numcols; j++) {
		inst->int_var[j] = (inst->mip_ctype[j] == CPX_BINARY) || (inst->mip_ctype[j] == CPX_INTEGER);
	}
	return inst->status;
}

int solve_continuous_relaxation(instance* inst) {
	// MIP --> continuous relaxation
	inst->status = CPXchgprobtype(inst->env, inst->lp, CPXPROB_LP);
	if (inst->status) { fprintf(stderr, "[ERR][solve_continuous_relaxation]: Failed to change problem type.\n"); return inst->status; }
	// Optimize LP
	inst->status = CPXlpopt(inst->env, inst->lp);
	if (inst->status) { fprintf(stderr, "[ERR][solve_continuous_relaxation]: Failed to optimize LP.\n"); return inst->status; }
	return inst->status;
}

int read_solution(instance* inst) {
	// Get solution
	inst->solnstat = CPXgetstat(inst->env, inst->lp);
	if (inst->solnstat == CPX_STAT_UNBOUNDED)       { fprintf(stdout, "[INFO][read_solution]: Model is unbounded.\n"); return inst->status; }
	else if (inst->solnstat == CPX_STAT_INFEASIBLE) { fprintf(stdout, "[INFO][read_solution]: Model is infeasible.\n"); return inst->status; }
	else if (inst->solnstat == CPX_STAT_INForUNBD)  { fprintf(stdout, "[INFO][read_solution]: Model is infeasible or unbounded.\n"); return inst->status; }
	inst->status = CPXsolninfo(inst->env, inst->lp, &(inst->solnmethod), &(inst->solntype), NULL, NULL);
	if (inst->status) { fprintf(stderr, "[ERR][read_solution]: Failed to obtain solution info.\n"); return inst->status; }
	if (inst->solntype == CPX_NO_SOLN) { fprintf(stderr, "[ERR][read_solution]: Solution not available.\n"); return inst->status; }
	// fprintf(stdout, "Solution status %d, solution method %d.\n", inst->solnstat, inst->solnmethod);
	inst->x = (double*)malloc(inst->cur_numcols * sizeof(double));
	if (inst->x == NULL) { fprintf(stderr, "[ERR][read_solution]: Failed to allocate solution.\n"); return 1; }
	inst->status = CPXgetx(inst->env, inst->lp, inst->x, 0, inst->cur_numcols - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_solution]: Failed to obtain primal solution.\n"); return inst->status; }
	return inst->status;
}

int read_variable_bounds(instance* inst) {
	// Get variable bounds (upper and lower)
	inst->ub = (double*)malloc(inst->cur_numcols * sizeof(double));
	inst->lb = (double*)malloc(inst->cur_numcols * sizeof(double));
	if (inst->ub == NULL || inst->lb == NULL) { fprintf(stderr, "[ERR][read_variable_bounds]: Failed to allocate variable bounds.\n"); return 1; }
	inst->status = CPXgetub(inst->env, inst->lp, inst->ub, 0, inst->cur_numcols - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_variable_bounds]: Failed to obtain upper bounds.\n"); return inst->status; }
	inst->status = CPXgetlb(inst->env, inst->lp, inst->lb, 0, inst->cur_numcols - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_variable_bounds]: Failed to obtain lower bounds.\n"); return inst->status; }
	return inst->status;
}

int read_objective_value(instance* inst) {
	// Get objective value
	inst->status = CPXgetobjval(inst->env, inst->lp, &(inst->objval));
	if (inst->status) { fprintf(stderr, "[ERR][read_objective_value]: Failed to obtain objective value.\n"); return inst->status; }
	return inst->status;
}

int read_objective_coefficients(instance* inst) {
	// Get objective coefficients
	inst->obj = (double*)malloc(inst->cur_numcols * sizeof(double));
	if (inst->obj == NULL) { fprintf(stderr, "[ERR][read_objective_coefficients]: Failed to allocate objective coefficients.\n"); return 1; }
	inst->status = CPXgetobj(inst->env, inst->lp, inst->obj, 0, inst->cur_numcols - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_objective_coefficients]: Failed to obtain objective coefficients.\n"); return inst->status; }
	return inst->status;
}

int read_constraints_coefficients(instance* inst) {
	// Get constraint matrix, both by rows and by columns
	// First, get the number of non zero coefficients of the matrix (nzcnt)
	int unused = 0;
	inst->nzcnt = CPXgetnumnz(inst->env, inst->lp);
	// Get rows
	inst->rmatbeg = (int*)malloc(inst->cur_numrows * sizeof(int));
	inst->rmatind = (int*)malloc(inst->nzcnt * sizeof(int));
	inst->rmatval = (double*)malloc(inst->nzcnt * sizeof(double));
	if (inst->rmatbeg == NULL || inst->rmatind == NULL || inst->rmatval == NULL) { fprintf(stderr, "[ERR][read_constraints_coefficients]: Failed to allocate one of rmatbeg, rmatind, rmatval.\n"); return 1; }
	inst->status = CPXgetrows(inst->env, inst->lp, &unused, inst->rmatbeg, inst->rmatind, inst->rmatval, inst->nzcnt, &unused, 0, inst->cur_numrows - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_constraints_coefficients]: Failed to obtain rows info.\n"); return inst->status; }
	// Get columns
	inst->cmatbeg = (int*)malloc(inst->cur_numcols * sizeof(int));
	inst->cmatind = (int*)malloc(inst->nzcnt * sizeof(int));
	inst->cmatval = (double*)malloc(inst->nzcnt * sizeof(double));
	if (inst->cmatbeg == NULL || inst->cmatind == NULL || inst->cmatval == NULL) { fprintf(stderr, "[ERR][read_constraints_coefficients]: Failed to allocate one of cmatbeg, cmatind, cmatval.\n"); return 1; }
	inst->status = CPXgetcols(inst->env, inst->lp, &unused, inst->cmatbeg, inst->cmatind, inst->cmatval, inst->nzcnt, &unused, 0, inst->cur_numcols - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_constraints_coefficients]: Failed to obtain columns info.\n"); return inst->status; }
	return inst->status;
}

int read_constraints_senses(instance* inst) {
	// Get constraint senses {'L','E','G'}
	inst->sense = (char*)malloc(inst->cur_numrows * sizeof(char));
	if (inst->sense == NULL) { fprintf(stderr, "[ERR][read_constraints_senses]: Failed to allocate constraint senses.\n"); return 1; }
	inst->status = CPXgetsense(inst->env, inst->lp, inst->sense, 0, inst->cur_numrows - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_constraints_senses]: Failed to obtain constraints senses.\n"); return inst->status; }
	return inst->status;
}

int read_constraints_right_hand_sides(instance* inst) {
	// Get right hand sides
	inst->rhs = (double*)malloc(inst->cur_numrows * sizeof(double));
	if (inst->rhs == NULL) { fprintf(stderr, "[ERR][read_constraints_right_hand_sides]: Failed to allocate right hand sides.\n"); return 1; }
	inst->status = CPXgetrhs(inst->env, inst->lp, inst->rhs, 0, inst->cur_numrows - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_constraints_right_hand_sides]: Failed to obtain rhs.\n"); return inst->status; }
	return inst->status;
}

int read_row_slacks(instance* inst) {
	// Get row slacks
	inst->slack = (double*)malloc(inst->cur_numrows * sizeof(double));
	if (inst->slack == NULL) { fprintf(stderr, "[ERR][read_row_slacks]: Failed to allocate slacks.\n"); return 1; }
	inst->status = CPXgetslack(inst->env, inst->lp, inst->slack, 0, inst->cur_numrows - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_row_slacks]: Failed to obtain slacks.\n"); return inst->status; }
	return inst->status;
}

int read_problem_data(instance* inst) {
	int status = inst->status;
	status = read_solution(inst);                     if (status) { fprintf(stderr, "[ERR][read_problem_data]: Error inside read_solution.\n"); return status; }
	status = read_variable_bounds(inst);              if (status) { fprintf(stderr, "[ERR][read_problem_data]: Error inside read_variable_bounds.\n"); return status; }
	status = read_objective_value(inst);              if (status) { fprintf(stderr, "[ERR][read_problem_data]: Error inside read_objective_value.\n"); return status; }
	status = read_objective_coefficients(inst);       if (status) { fprintf(stderr, "[ERR][read_problem_data]: Error inside read_objective_coefficients.\n"); return status; }
	status = read_constraints_coefficients(inst);     if (status) { fprintf(stderr, "[ERR][read_problem_data]: Error inside read_constraints_coefficients.\n"); return status; }
	status = read_constraints_senses(inst);           if (status) { fprintf(stderr, "[ERR][read_problem_data]: Error inside read_constraints_senses.\n"); return status; }
	status = read_constraints_right_hand_sides(inst); if (status) { fprintf(stderr, "[ERR][read_problem_data]: Error inside read_constraints_right_hand_sides.\n"); return status; }
	status = read_row_slacks(inst);                   if (status) { fprintf(stderr, "[ERR][read_problem_data]: Error inside read_row_slacks.\n"); return status; }
	return status;
}

//**************************************************************************************************************************************************************
//*************************************** ZI-ROUND *************************************************************************************************************
//**************************************************************************************************************************************************************
int zi_round(instance* inst) {

	int status = inst->status;
	// inst->x_prev = (double*)malloc(inst->cur_numcols * sizeof(double));
	// inst->x_updated = (double*)malloc(inst->cur_numcols * sizeof(double));
	inst->UB = (double*)malloc(inst->cur_numcols * sizeof(double));
	inst->LB = (double*)malloc(inst->cur_numcols * sizeof(double));
	if (/*inst->x_prev == NULL || inst->x_updated == NULL ||*/ inst->UB == NULL || inst->LB == NULL) { 
		fprintf(stderr, "[ERR][zi_round]: Failed to allocate x_prev, x_updated or row_infeas.\n"); return 1; 
	}

	// Outer loop (repeat until no more updates found)
	do {

		// Initialize update flag and save current solution
		inst->updated = 0;
		// clone_array(inst->x, inst->x_prev, inst->cur_numcols);

		// Inner loop (for each variable xj that was integer/binary in the original MIP)
		for (int j = 0; j < inst->cur_numcols; j++) {
			// Skip xj if it's continuous in the original MIP
			if (!(inst->int_var)) continue;

			// xj non-fractional
			if (!is_fractional(inst->x[j])) {
				if (VERBOSE > 10) fprintf(stdout, "[DEBUG][zi_round]: -x- x_%d is non-fractional\n", j + 1);
				
				// Calculate UBj and LBj (with epsilon = 1.0)
				calculate_UBjLBj(inst, j, 1.0);
				// Skip xj if both UBj and LBj are equal to zero (--> no shift necessary)
				if (fabs(inst->UB[j]) < TOLERANCE && fabs(inst->LB[j]) < TOLERANCE) continue;

				// Condition(s) for rounding of xj
				// (>= to include the case of a zero obj coefficient)
				if ((inst->obj[j] >= 0 && fabs(inst->LB[j] - 1.0) < TOLERANCE) ||
					(inst->obj[j] <= 0 && fabs(inst->UB[j] - 1.0) < TOLERANCE)) {
					// Update xj to improve objective and update slacks
					update_xj_to_improve_objective(inst, j, 0); // flag xj non-fractional (0)
				}
			} // end if xj non-fractional

			// xj fractional
			else if (is_fractional(inst->x[j])) {
				if (VERBOSE > 10) fprintf(stdout, "[DEBUG][zi_round]: -x- x_%d is fractional\n", j + 1);
				
				// Calculate UBj and LBj
				calculate_UBjLBj(inst, j, EPSILON);
				// Skip xj if both UBj and LBj are equal to zero (--> no shift necessary)
				if (fabs(inst->UB[j]) < TOLERANCE && fabs(inst->LB[j]) < TOLERANCE) continue;
				
				// ZI-Round (version 1) core (3 cases). First, calculate the fractionalities needed.
				inst->ZI = fractionality(inst->x[j]);
				inst->ZIplus = fractionality(inst->x[j] + inst->UB[j]);
				inst->ZIminus = fractionality(inst->x[j] - inst->LB[j]);

				// First case
				if (fabs(inst->ZIplus - inst->ZIminus) < TOLERANCE && 
					(inst->ZIplus - inst->ZI) < -(TOLERANCE)) { // was inst->ZIplus < inst->ZI
					// Update xj to improve objective and update slacks
					update_xj_to_improve_objective(inst, j, 1); // flag xj fractional (1)
				} // end first case

				// Second case
				else if ((inst->ZIplus - inst->ZIminus) < -(TOLERANCE) && // was inst->ZIplus < inst->ZIminus
					     (inst->ZIplus - inst->ZI) < -(TOLERANCE)) {      // was inst->ZIplus < inst->ZI
					if (VERBOSE > 10) fprintf(stdout, "[DEBUG][zi_round]: >>> Set x_%d = x_%d + UB_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], inst->UB[j], inst->x[j] + inst->UB[j]);
					inst->x[j] += inst->UB[j];
					inst->updated = 1;
					update_slacks(inst, j, inst->UB[j]);
					update_objective_value(inst, j, inst->UB[j]);
				} // end second case

				// Third case
				else if ((inst->ZIminus - inst->ZIplus) < -(TOLERANCE) && // was inst->ZIminus < inst->ZIplus
					     (inst->ZIminus - inst->ZI) < -(TOLERANCE)) {     // was inst->ZIminus < inst->ZI
					if (VERBOSE > 10) fprintf(stdout, "[DEBUG][zi_round]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], inst->LB[j], inst->x[j] - inst->LB[j]);
					inst->x[j] -= inst->LB[j];
					inst->updated = 1;
					update_slacks(inst, j, -(inst->LB[j]));
					update_objective_value(inst, j, -(inst->LB[j]));
				} // end third case

			} // end if xj fractional

		} // end inner loop

		if (VERBOSE > 10) {
			if (inst->updated) fprintf(stdout, "[DEBUG][zi_round] ...Update found, scan variables again...\n");
			else fprintf(stdout, "[DEBUG][zi_round] ...No updates found, exit outer loop...\n");
		}
		// [BRUTE FORCE] [DEBUG ONLY] Check variable bounds and constraints
		if (VERBOSE > 50) {
			status = check_bounds(inst, inst->x);      if (status) { fprintf(stderr, "[ERR][zi_round]: Error inside check_bounds.\n"); return status; }
			status = check_constraints(inst, inst->x); if (status) { fprintf(stderr, "[ERR][zi_round]: Error inside check_constraints.\n"); return status; }
		}

	} while (inst->updated); // end outer loop

	return inst->status;
}
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************

void update_xj_to_improve_objective(instance* inst, int j, int is_fractional) {

	// First, calculate obj value for both updates. (IMPROVEMENT: inst->objval is just a common offset!)
	inst->obj_plusUBj = inst->objval + (inst->obj[j] * inst->UB[j]);
	inst->obj_minusLBj = inst->objval - (inst->obj[j] * inst->LB[j]);
	/* was
	// compare delta * coef (confronto solo tra i due delta * coef_j (+ objval è solo un offset) (tenendo conto della solita float tolerance)
	clone_array(inst->x, inst->x_updated, inst->cur_numcols);
	inst->x_updated[j] = inst->x[j] + inst->UB[j];
	..inst->obj_plusUBj = dot_product(inst->obj, inst->x_updated, inst->cur_numcols); // brute force!
	inst->x_updated[j] = inst->x[j] - inst->LB[j];
	..inst->obj_minusLBj = dot_product(inst->obj, inst->x_updated, inst->cur_numcols); // brute force!*/

	// Check obj sense, then update xj, update slacks and update objective value
	switch (inst->objsen) {
		case CPX_MIN:
			if ((inst->obj_plusUBj - inst->objval) < -(TOLERANCE) &&       // was inst->obj_plusUBj < inst->objval 
				(inst->obj_plusUBj - inst->obj_minusLBj) < -(TOLERANCE)) { // was inst->obj_plusUBj < inst->obj_minusLBj
				// xj = xj + UBj (if xj is not fractional then UBj must be 1.0)
				if (VERBOSE > 50) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: _>_ UBj = %f (must be 1.0)\n", inst->UB[j]);
				if (VERBOSE > 10) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: >>> Set x_%d = x_%d + UB_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], inst->UB[j], inst->x[j] + inst->UB[j]);
				if (!is_fractional && VERBOSE > 10 && fabs(inst->UB[j] - 1.0) > TOLERANCE) {
					fprintf(stderr, "[ERR][update_xj_to_improve_objective]: UB_%d = %f (should be 1.0).\n", j + 1, inst->UB[j]);
					exit(EXIT_FAILURE);
				}
				inst->x[j] += inst->UB[j];
				inst->updated = 1;
				update_slacks(inst, j, inst->UB[j]);
				update_objective_value(inst, j, inst->UB[j]);
			}
			else if ((inst->obj_minusLBj - inst->objval) < -(TOLERANCE) &&      // was inst->obj_minusLBj < inst->objval
				     (inst->obj_minusLBj - inst->obj_plusUBj) < -(TOLERANCE)) { // was inst->obj_minusLBj < inst->obj_plusUBj
				// xj = xj - LBj (if xj is not fractional then LBj must be 1.0)
				if (VERBOSE > 50) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: _>_ LBj = %f (must be 1.0)\n", inst->LB[j]);
				if (VERBOSE > 10) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], inst->LB[j], inst->x[j] - inst->LB[j]);
				if (!is_fractional && VERBOSE > 10 && fabs(inst->LB[j] - 1.0) > TOLERANCE) {
					fprintf(stderr, "[ERR][update_xj_to_improve_objective]: LB_%d = %f (should be 1.0).\n", j + 1, inst->LB[j]);
					exit(EXIT_FAILURE);
				}
				inst->x[j] -= inst->LB[j];
				inst->updated = 1;
				update_slacks(inst, j, -(inst->LB[j]));
				update_objective_value(inst, j, -(inst->LB[j]));
			}
			break;
		// (problems from mps files should always be MIN)
		case CPX_MAX:
			if ((inst->obj_plusUBj - inst->objval) > TOLERANCE &&       // was inst->obj_plusUBj > inst->objval
				(inst->obj_plusUBj - inst->obj_minusLBj) > TOLERANCE) { // was inst->obj_plusUBj > inst->obj_minusLBj
				// xj = xj + UBj (if xj is not fractional then UBj must be 1.0)
				if (VERBOSE > 50) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: _>_ UBj = %f (must be 1.0)\n", inst->UB[j]);
				if (VERBOSE > 10) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: >>> Set x_%d = x_%d + UB_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], inst->UB[j], inst->x[j] + inst->UB[j]);
				if (!is_fractional && VERBOSE > 10 && fabs(inst->UB[j] - 1.0) > TOLERANCE) {
					fprintf(stderr, "[ERR][update_xj_to_improve_objective]: UB_%d = %f (should be 1.0).\n", j + 1, inst->UB[j]);
					exit(EXIT_FAILURE);
				}
				inst->x[j] += inst->UB[j];
				inst->updated = 1;
				update_slacks(inst, j, inst->UB[j]);
				update_objective_value(inst, j, inst->UB[j]);
			}
			else if ((inst->obj_minusLBj - inst->objval) > TOLERANCE &&      // was inst->obj_minusLBj > inst->objval
				     (inst->obj_minusLBj - inst->obj_plusUBj) > TOLERANCE) { // was inst->obj_minusLBj > inst->obj_plusUBj
				// xj = xj - LBj (if xj is not fractional then LBj must be 1.0)
				if (VERBOSE > 50) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: _>_ LBj = %f (must be 1.0)\n", inst->LB[j]);
				if (VERBOSE > 10) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], inst->LB[j], inst->x[j] - inst->LB[j]);
				if (!is_fractional && VERBOSE > 10 && fabs(inst->LB[j] - 1.0) > TOLERANCE) {
					fprintf(stderr, "[ERR][update_xj_to_improve_objective]: LB_%d = %f (should be 1.0).\n", j + 1, inst->LB[j]);
					exit(EXIT_FAILURE);
				}
				inst->x[j] -= inst->LB[j];
				inst->updated = 1;
				update_slacks(inst, j, -(inst->LB[j]));
				update_objective_value(inst, j, -(inst->LB[j]));
			}
			break;
		default:
			fprintf(stderr, "[ERR][update_xj_to_improve_objective]: Entered default case in switch.\n");
			exit(EXIT_FAILURE);
			break;
	} // end switch
}

double row_activity(instance* inst, int i, double* x) {
	double sum = 0; double aij; int col_index; 
	int end_row = (i < inst->cur_numrows - 1) ? inst->rmatbeg[i + 1] : inst->nzcnt;
	// Scan non-zero coefficients of row i
	/*
		For constraint i:
			rmatbeg[i] is the first index of rmatind and rmatval for row i
			--> Row i in range [ rmatbeg[i];rmatbeg[i+1] ) except for the last one (see end_row)
			rmatval contains coefficient values
			rmatind contains the column indices of coefficient values
	*/
	for (int k = inst->rmatbeg[i]; k < end_row; k++) {
		// Get current non-zero aij and its column index
		aij = inst->rmatval[k];
		col_index = inst->rmatind[k];
		// Update row activity
		sum += (aij * x[col_index]);
	}
	return sum;
}

// Incremental update of row slacks (only for the rows where the variable xj is involved)
// Whenever it is called, only one variable xj has been updated
void update_slacks(instance* inst, int j, double signed_delta) { 
	double aij; int row_index; char sense;
	int end_col = (j < inst->cur_numcols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;
	// Scan non-zero coefficients of column j
	/*
		For variable xj:
			cmatbeg[j] is the first index of cmatind and cmatval for column j
			--> Column j in range [ cmatbeg[j];cmatbeg[j+1] ) except for the last one (see end_col)
			cmatval contains coefficient values
			cmatind contains the row indices of coefficient values
	*/
	for (int k = inst->cmatbeg[j]; k < end_col; k++) {
		aij = inst->cmatval[k];
		row_index = inst->cmatind[k];
		sense = inst->sense[row_index];
		if (sense == 'E') {
			fprintf(stderr, "[ERR][update_slacks]: Tried to update slack for a variable involved in an equality constraint!\n");
			exit(EXIT_FAILURE);
		}
		/*
			slack = rhs - row_activity
			row_activity --> row_activity + (aij * signed_delta)
			--> slack = rhs - (row_activity + (aij * signed_delta))
			          = rhs - row_activity - (aij * signed_delta)
					  = (rhs - row_activity) - (aij * signed_delta)
			slack --> slack - (aij * signed_delta)
		*/
		inst->slack[row_index] -= (aij * signed_delta);
	}
}

// Incremental update of the objective value
// Whenever it is called, only one variable xj has been updated
void update_objective_value(instance* inst, int j, double signed_delta) {
	inst->objval += (inst->obj[j] * signed_delta);
}

void calculate_UBjLBj(instance* inst, int j, const double epsilon) {
	/*
		For 'L' (<=) constraints: (si non-negative)
			firstUBj_L = min_i{si/aij : aij > 0}
			firstLBj_L = min_i{-si/aij : aij < 0}
		For 'G' (>=) constraints: (si non-positive)
			firstUBj_G = min_i{si/aij : aij < 0}
			firstLBj_G = min_i{-si/aij : aij > 0}

		--> firstUBj = min{ firstUBj_L , firstUBj_G }
		--> firstLBj = min{ firstLBj_L , firstLBj_G }
	*/
	double firstUBj = LONG_MAX, firstLBj = LONG_MAX;
	double aij; int row_index; char sense;
	double ratioUBj, ratioLBj;
	double newUBj, newLBj;
	int end_col = (j < inst->cur_numcols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;

	double secondUBj = inst->ub[j] - inst->x[j];
	double secondLBj = inst->x[j] - inst->lb[j];
	if (VERBOSE > 10) fprintf(stdout, "[DEBUG][calculate_UBjLBj]: secondUB_%d = ub_%d - x_%d = %f - %f = %f ; secondLB_%d = x_%d - lb_%d = %f - %f = %f\n", 
						j + 1, j + 1, j + 1, inst->ub[j], inst->x[j], secondUBj, j + 1, j + 1, j + 1, inst->x[j], inst->lb[j], secondLBj);

	// Scan non-zero coefficients of column j
	/*
		For variable xj:
			cmatbeg[j] is the first index of cmatind and cmatval for column j
			--> Column j in range [ cmatbeg[j];cmatbeg[j+1] ) except for the last one (see end_col)
			cmatval contains coefficient values
			cmatind contains the row indices of coefficient values
	*/
	for (int k = inst->cmatbeg[j]; k < end_col; k++) {
		// Get current non-zero aij, its row index and constraint sense
		aij = inst->cmatval[k]; 
		row_index = inst->cmatind[k];
		sense = inst->sense[row_index];
		// Check sense, then check sign of aij, and update firstUBj, firstLBj
		switch (sense) {
			case 'L': // (slack non-negative)
				// Clip slack to zero if negative
				if (inst->slack[row_index] < -(TOLERANCE)) inst->slack[row_index] = 0.0;
				if (aij > 0) { // Update firstUBj
					if (VERBOSE > 20) fprintf(stdout, "[DEBUG][calculate_UBjLBj]: sense = L ; slack[%d] = %f ; a_%d_%d = %f ; ratioUB_%d = %f\n", 
										row_index + 1, inst->slack[row_index], row_index + 1, j + 1, aij, j + 1, inst->slack[row_index] / aij);
					ratioUBj = inst->slack[row_index] / aij;
					firstUBj = min(ratioUBj, firstUBj);
				}
				if (aij < 0) { // Update firstLBj
					if (VERBOSE > 20) fprintf(stdout, "[DEBUG][calculate_UBjLBj]: sense = L ; slack[%d] = %f ; a_%d_%d = %f ; ratioLB_%d = %f\n",
										row_index + 1, inst->slack[row_index], row_index + 1, j + 1, aij, j + 1, -(inst->slack[row_index]) / aij);
					ratioLBj = -(inst->slack[row_index]) / aij;
					firstLBj = min(ratioLBj, firstLBj);
				}
				break;
			case 'G': // (slack non-positive)
				// Clip slack to zero if positive
				if (inst->slack[row_index] > TOLERANCE) inst->slack[row_index] = 0.0;
				if (aij < 0) { // Update firstUBj
					if (VERBOSE > 20) fprintf(stdout, "[DEBUG][calculate_UBjLBj]: sense = G ; slack[%d] = %f ; a_%d_%d = %f ; ratioUB_%d = %f\n",
										row_index + 1, inst->slack[row_index], row_index + 1, j + 1, aij, j + 1, inst->slack[row_index] / aij);
					ratioUBj = inst->slack[row_index] / aij;
					firstUBj = min(ratioUBj, firstUBj);
				}
				if (aij > 0) { // Update firstLBj
					if (VERBOSE > 20) fprintf(stdout, "[DEBUG][calculate_UBjLBj]: sense = G ; slack[%d] = %f ; a_%d_%d = %f ; ratioLB_%d = %f\n",
										row_index + 1, inst->slack[row_index], row_index + 1, j + 1, aij, j + 1, -(inst->slack[row_index]) / aij);
					ratioLBj = -(inst->slack[row_index]) / aij;
					firstLBj = min(ratioLBj, firstLBj);
				}
				break;
			case 'E': // (slack zero): variable xj involved in equality constraint, thus cannot be moved
				// Set firstUBj and firstLBj to zero --> newUBj and newLBj will get value zero
				if (VERBOSE > 20) fprintf(stdout, "[DEBUG][calculate_UBjLBj]: sense = E ; Variable x_%d involved in equality constraint %d --> cannot be moved!\n", j, row_index);
				firstUBj = 0.0;
				firstLBj = 0.0;
				break;
			default:
				fprintf(stderr, "[ERR][calculate_UBjLBj]: Constraint sense not included in {'L','G','E'}.\n"); exit(EXIT_FAILURE);
				break;
		} // end switch
	} // end for

	// Results
	if (VERBOSE > 10) fprintf(stdout, "[DEBUG][calculate_UBjLBj][Results]: firstUB_%d = %f ; firstLB_%d = %f\n", j + 1, firstUBj, j + 1, firstLBj);
	newUBj = min(firstUBj, secondUBj);
	newLBj = min(firstLBj, secondLBj);
	if (VERBOSE > 10) fprintf(stdout, 
		"[DEBUG][calculate_UBjLBj][Results]: (New) UB_%d = min{firstUB_%d, secondUB_%d} = min{%f, %f} = %f ; LB_%d = min{firstLB_%d, secondLB_%d} = min{%f, %f} = %f\n", 
			j + 1, j + 1, j + 1, firstUBj, secondUBj, newUBj, j + 1, j + 1, j + 1, firstLBj, secondLBj, newLBj);
	
	// If newUBj < epsilon clip it to 0
	if (newUBj < epsilon) newUBj = 0.0;
	// If newLBj < epsilon clip it to 0
	if (newLBj < epsilon) newLBj = 0.0;
	// Update both no matter what
	inst->UB[j] = newUBj;
	inst->LB[j] = newLBj;
	/* was
	// Update UBj and LBj iff [ (they both do not fall below epsilon) || (they should both be set to zero due to an equality constraint) ]
	if (!(newUBj < epsilon && newLBj < epsilon) || (newUBj == 0.0 && newLBj == 0.0)) {
		inst->UB[j] = newUBj;
		inst->LB[j] = newLBj;
	} */
}

int check_bounds(instance* inst, double* x) {
	// Check variable bounds
	for (int j = 0; j < inst->cur_numcols; j++)
		// was x[j] < inst->lb[j] || x[j] > inst->ub[j]
		if ((x[j] - inst->lb[j]) < -(TOLERANCE) || (x[j] - inst->ub[j]) > TOLERANCE) { fprintf(stderr, "[ERR][check_bounds][!!!]: Bound %d violated!\n", j + 1); return 1; }
	return 0;
}

int check_constraints(instance* inst, double* x) {
	// Check constraints
	double* row_infeas = (double*)malloc(inst->cur_numrows * sizeof(double));
	if (row_infeas == NULL) { fprintf(stderr, "[ERR][check_constraints]: Failed to allocate row_infeas.\n"); return 1; }
	inst->status = CPXgetrowinfeas(inst->env, inst->lp, x, row_infeas, 0, inst->cur_numrows - 1);
	if (inst->status) { fprintf(stderr, "[ERR][check_constraints]: Failed to obtain row infeasibilities.\n"); free_and_null((char**)&row_infeas); return inst->status; }
	for (int i = 0; i < inst->cur_numrows; i++)
		if (fabs(row_infeas[i]) > TOLERANCE) { fprintf(stdout, "[ERR][check_constraints][!!!]: Constraint %d violated!\n", i + 1); free_and_null((char**)&row_infeas); return 1; }
	free_and_null((char**)&row_infeas);
	return inst->status;
}

int free_instance(instance* inst) {
	free_and_null((char**)&(inst->x));       free_and_null((char**)&(inst->obj));     free_and_null((char**)&(inst->lb));
	free_and_null((char**)&(inst->ub));      free_and_null((char**)&(inst->slack));   free_and_null((char**)&(inst->mip_ctype));
	free_and_null((char**)&(inst->int_var)); free_and_null((char**)&(inst->rmatbeg)); free_and_null((char**)&(inst->rmatind));
	free_and_null((char**)&(inst->rmatval)); free_and_null((char**)&(inst->cmatbeg)); free_and_null((char**)&(inst->cmatind));
	free_and_null((char**)&(inst->cmatval)); free_and_null((char**)&(inst->sense));   free_and_null((char**)&(inst->rhs));
	free_and_null((char**)&(inst->UB));      free_and_null((char**)&(inst->LB));      // free_and_null((char**)&(inst->x_prev)); 
	// free_and_null((char**)&(inst->x_updated));
	if (inst->lp != NULL) {
		inst->status = CPXfreeprob(inst->env, &(inst->lp)); 
		if (inst->status) fprintf(stderr, "[ERR][free_instance]: CPXfreeprob failed, error code %d.\n", inst->status);
	}
	if (inst->env != NULL) {
		inst->status = CPXcloseCPLEX(&(inst->env));
		if (inst->status) {
			char errmsg[CPXMESSAGEBUFSIZE];
			fprintf(stderr, "[ERR][free_instance]: Could not close CPLEX environment.\n");
			CPXgeterrorstring(inst->env, inst->status, errmsg);
			fprintf(stderr, "[ERR][free_instance]: %s", errmsg);
		}
	}
	return inst->status;
}

void print_problem_info(instance* inst, int sol_available, int to_file) {
	FILE* out_stream = (to_file) ? fopen("problem_info.txt", "w") : stdout;
	// cur_numrows, cur_numcols
	int nrows = CPXgetnumrows(inst->env, inst->lp);
	int ncols = CPXgetnumcols(inst->env, inst->lp);
	// variable bounds
	double* ub = (double*)malloc(ncols * sizeof(double));
	double* lb = (double*)malloc(ncols * sizeof(double));
	CPXgetub(inst->env, inst->lp, ub, 0, ncols - 1);
	CPXgetlb(inst->env, inst->lp, lb, 0, ncols - 1);
	// obj sense and coefficients
	int objsen = CPXgetobjsen(inst->env, inst->lp);
	double* obj = (double*)malloc(ncols * sizeof(double));
	CPXgetobj(inst->env, inst->lp, obj, 0, ncols - 1);
	// constraints coefficients, rhs and senses
	int nnz = CPXgetnumnz(inst->env, inst->lp); int unused = 0;
	int* rmatbeg = (int*)malloc(nrows * sizeof(int));
	int* rmatind = (int*)malloc(nnz * sizeof(int));
	double* rmatval = (double*)malloc(nnz * sizeof(double));
	CPXgetrows(inst->env, inst->lp, &unused, rmatbeg, rmatind, rmatval, nnz, &unused, 0, nrows - 1);
	int* cmatbeg = (int*)malloc(ncols * sizeof(int));
	int* cmatind = (int*)malloc(nnz * sizeof(int));
	double* cmatval = (double*)malloc(nnz * sizeof(double)); 
	CPXgetcols(inst->env, inst->lp, &unused, cmatbeg, cmatind, cmatval, nnz, &unused, 0, ncols - 1);
	double* rhs = (double*)malloc(nrows * sizeof(double));
	CPXgetrhs(inst->env, inst->lp, rhs, 0, nrows - 1);
	char* sense = (char*)malloc(nrows * sizeof(char));
	CPXgetsense(inst->env, inst->lp, sense, 0, nrows - 1);
	// solution and obj value
	double* x; double objval;
	if (sol_available) {
		x = (double*)malloc(ncols * sizeof(double));
		clone_array(inst->x, x, ncols);
		objval = inst->objval;
	}
	fprintf(out_stream, "\n********************************************************\n");
	fprintf(out_stream, "********************  PROBLEM INFO  ********************\n");
	fprintf(out_stream, "**** Dimensions (rows, cols): (%d, %d)\n", nrows, ncols);
	fprintf(out_stream, "**** Obj sense: %s\n", ((objsen > 0)?"MIN":"MAX"));
	fprintf(out_stream, "**** Obj coefficients:\n");
	for (int j = 0; j < ncols; j++) fprintf(out_stream, "****     x_%d: %f\n", j + 1, obj[j]);
	fprintf(out_stream, "**** Constraints: (nnz = %d)\n", nnz);
	for (int i = 0; i < nrows - 1; i++) {
		fprintf(out_stream, "****     ");
		for (int k = rmatbeg[i]; k < rmatbeg[i+1]; k++) fprintf(out_stream, "(%f)x_%d ", rmatval[k], rmatind[k]+1);
		fprintf(out_stream, " %s %f\n", ((sense[i] == 'L')?"<=":((sense[i] == 'G')?">=":"=")), rhs[i]);
	}
	fprintf(out_stream, "****     "); // last row
	for (int k = rmatbeg[nrows - 1]; k < nnz; k++) fprintf(out_stream, "(%f)x_%d ", rmatval[k], rmatind[k]+1);
	fprintf(out_stream, " %s %f\n", ((sense[nrows - 1] == 'L') ? "<=" : ((sense[nrows - 1] == 'G') ? ">=" : "=")), rhs[nrows - 1]);
	fprintf(out_stream, "**** Variable bounds:\n");
	for (int j = 0; j < ncols; j++) {
		fprintf(out_stream, "****     %f <= x_%d <= %f\n", lb[j], j + 1, ub[j]);
	}
	if (sol_available) {
		// solution and obj value
		double* x = (double*)malloc(ncols * sizeof(double));
		clone_array(inst->x, x, ncols);
		double objval = inst->objval;
		fprintf(out_stream, "**** Solution: \n");
		for (int j = 0; j < ncols; j++) fprintf(out_stream, "****     x_%d = %f\n", j + 1, x[j]);
		fprintf(out_stream, "**** Obj value: %f\n", objval);
	}
	fprintf(out_stream, "**** Integer or binary variables:\n");
	for (int j = 0; j < ncols; j++) fprintf(out_stream, "****     x_%d: %s\n", j + 1, ((inst->int_var[j])?"BIN/INT":"CONT"));
	fprintf(out_stream, "********************************************************\n");
	fprintf(out_stream, "********************************************************\n\n");
	if (to_file) fclose(out_stream);
}

void free_and_null(char** ptr) { if (*ptr != NULL) { free(*ptr); *ptr = NULL; } }

double fractionality(double xj) {
	double minusfloor = xj - floor(xj); 
	double minusceil = ceil(xj) - xj;
	return min(minusfloor, minusceil);
}

int is_fractional(double num) { return 1 - (fabs(num - round(num)) < TOLERANCE); }

double dot_product(double* coef, double* var_value, int len) {
	double dotprod = 0;
	for (int j = 0; j < len; j++) dotprod += (coef[j]) * (var_value[j]);
	return dotprod;
}

int different_arr(double* prev, double* new, int len) {
	int dif = 0;
	for (int i = 0; i < len; i++) {
		if (fabs(prev[i] - new[i]) > TOLERANCE) {
			dif = 1;
			break;
		}
	}
	return dif;
}

void clone_array(double* arr, double* clo, int len) { for (int i = 0; i < len; i++) clo[i] = arr[i]; }