/**
 * @file compute_ziround_input.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void setup_CPLEX_env(instance* inst) {

	char* errmsg = (char*)malloc(CPXMESSAGEBUFSIZE * sizeof(char)); if (errmsg == NULL) print_error("[setup_CPLEX_env]: Failed to allocate errmsg.\n");
	int status = 0;

	inst->env = CPXopenCPLEX(&status);
	if (inst->env == NULL) {
		print_warning("[setup_CPLEX_env]: Could not open CPLEX environment.\n");
		CPXgeterrorstring(inst->env, status, errmsg);
		print_error("[setup_CPLEX_env]: %s", errmsg);
	}

	// Set CPLEX parameters
	if (CPXsetintparam(inst->env, CPXPARAM_ScreenOutput, CPX_ON)) print_error("[setup_CPLEX_env]: Failed to turn on screen indicator.\n");     
	if (CPXsetdblparam(inst->env, CPXPARAM_TimeLimit,    3600))   print_error("[setup_CPLEX_env]: Failed to set time limit.\n");

	// Free
	free(errmsg);
}

void read_MIP_problem(instance* inst, char* filename) {
	
	int status = 0;

	// Create MIP from input file (.mps)
	inst->lp = CPXcreateprob(inst->env, &status, filename);   if (inst->lp == NULL) print_error("[read_MIP_problem]: Failed to create MIP.\n");
	if (CPXreadcopyprob(inst->env, inst->lp, filename, NULL)) print_error("[read_MIP_problem]: Failed to read and copy the problem data.\n");
}

void save_integer_variables(instance* inst) {
	
	int ncols = CPXgetnumcols(inst->env, inst->lp); assert(ncols > 0);

	// Allocate variable types
	inst->vartype = (char*)malloc(ncols * sizeof(char)); if (inst->vartype == NULL) print_error("[save_integer_variables]: Failed to allocate vartype.\n");
	inst->int_var = (int*)calloc(ncols, sizeof(int));    if (inst->int_var == NULL) print_error("[save_integer_variables]: Failed to allocate int_var.\n");

	// Get MIP variable types {CPX_CONTINUOUS, CPX_BINARY, CPX_INTEGER, CPX_SEMICONT, CPX_SEMIINT}
	if (CPXgetctype(inst->env, inst->lp, inst->vartype, 0, ncols - 1)) print_error("[save_integer_variables]: Failed to obtain MIP variable types.\n");
	assert(valid_var_types(inst->vartype, ncols));

	// Remember integer variables {CPX_BINARY, CPX_INTEGER}
	for (int j = 0; j < ncols; j++) {

		switch (inst->vartype[j]) {
			case CPX_INTEGER:
			case CPX_BINARY:
				inst->int_var[j] = 1;
				break;
			case CPX_CONTINUOUS:
				inst->int_var[j] = 0;
				break;
			default:
				print_error(100, "[save_integer_variables][x_%d]: Type '%c' not supported.\n", j + 1, inst->vartype[j]);
		}
	}
}

void solve_continuous_relaxation(instance* inst) {

	if (CPXchgprobtype(inst->env, inst->lp, CPXPROB_LP)) print_error("[solve_continuous_relaxation]: Failed to change problem type.\n");
	if (CPXlpopt(inst->env, inst->lp)) print_error("[solve_continuous_relaxation]: Failed to optimize LP.\n");
}