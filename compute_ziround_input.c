/**
 * @file compute_ziround_input.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void setup_CPLEX_env(instance* inst) {

	// External variables
	CPXENVptr env; /**< CPLEX environment pointer. */
	// Local variables
	char* errmsg;  /**< CPLEX error message. */
	int status;	   /**< Support status flag. */

	// Initialize
	env = NULL;
	errmsg = (char*)malloc(CPXMESSAGEBUFSIZE * sizeof(char));
	status = 0;

	env = CPXopenCPLEX(&status);
	inst->env = env;
	if (env == NULL) {
		print_warning("[setup_CPLEX_env]: Could not open CPLEX environment.\n");
		CPXgeterrorstring(env, status, errmsg);
		print_error("[setup_CPLEX_env]: %s", errmsg);
	}

	// Set CPLEX parameters
	status = CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
	if (status) print_error("[setup_CPLEX_env]: Failed to turn on screen indicator, error %d.\n", status);
	status = CPXsetdblparam(env, CPXPARAM_TimeLimit, 3600);      
	if (status) print_error("[setup_CPLEX_env]: Failed to set time limit, error %d.\n", status);

	// Free
	free(errmsg);
}

void read_MIP_problem(instance* inst, char* filename) {

	// External variables
	CPXENVptr env; /**< CPLEX environment pointer. */
	CPXLPptr lp;   /**< CPLEX lp pointer. */
	// Local variables
	int status;	   /**< Support status flag. */

	// Initialize
	env = inst->env;
	lp = NULL;
	status = 0;

	// Create MIP from input file (.mps)
	lp = CPXcreateprob(env, &status, filename);
	inst->lp = lp;
	if (lp == NULL) print_error("[read_MIP_problem]: Failed to create MIP.\n");
	status = CPXreadcopyprob(env, lp, filename, NULL); 
	if (status) print_error("[read_MIP_problem]: Failed to read and copy the problem data.\n");
}

void save_integer_variables(instance* inst) {
	
	// External variables
	CPXENVptr env; /**< CPLEX environment pointer. */
	CPXLPptr lp;   /**< CPLEX lp pointer. */
	int ncols;     /**< Number of variables (columns) of the MIP problem. */
	char* vartype; /**< Variable types array. */
	int* int_var;  /**< Array of flags for integer/binary variables of the original MIP problem. */
	// Local variables
	int status;    /**< Support status flag. */

	// Allocate / Initialize
	env = inst->env;
	lp = inst->lp;
	ncols = inst->ncols;
	vartype = NULL;
	int_var = NULL;
	status = 0;

	// Allocate on the instance, assign to locals
	inst->mip_ctype = (char*)malloc(ncols * sizeof(char));
	vartype = inst->mip_ctype;
	if (vartype == NULL) print_error("[save_integer_variables]: Failed to allocate mip_ctype.\n");
	inst->int_var = (int*)malloc(ncols * sizeof(int));
	int_var = inst->int_var;
	if (int_var == NULL) print_error("[save_integer_variables]: Failed to allocate int_var.\n");

	// Get MIP variable types {CPX_CONTINUOUS, CPX_BINARY, CPX_INTEGER, CPX_SEMICONT, CPX_SEMIINT}
	status = CPXgetctype(env, lp, vartype, 0, ncols - 1);
	if (status) print_error("[save_integer_variables]: Failed to obtain MIP variable types.\n");

	// Remember integer variables {CPX_BINARY, CPX_INTEGER}
	for (int j = 0; j < ncols; j++) {
		int_var[j] = (vartype[j] == CPX_BINARY) || (vartype[j] == CPX_INTEGER);
	}
}

void solve_continuous_relaxation(instance* inst) {

	// External variables
	CPXENVptr env; /**< CPLEX environment pointer. */
	CPXLPptr lp;   /**< CPLEX lp pointer. */
	// Local variables
	int status;    /**< Support status flag. */

	// Initialize
	env = inst->env;
	lp = inst->lp;
	status = 0;

	// MIP --> continuous relaxation
	status = CPXchgprobtype(env, lp, CPXPROB_LP);
	if (status) print_error("[solve_continuous_relaxation]: Failed to change problem type.\n");

	// Optimize LP
	status = CPXlpopt(env, lp);
	if (status) print_error("[solve_continuous_relaxation]: Failed to optimize LP.\n");
}