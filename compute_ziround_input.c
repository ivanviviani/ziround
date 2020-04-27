/**
 * @file compute_ziround_input.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void setup_CPLEX_env(instance* inst) {

	// Local variables
	char* errmsg = (char*)malloc(CPXMESSAGEBUFSIZE * sizeof(char));
	if (errmsg == NULL) print_error("[setup_CPLEX_env]: Failed to allocate errmsg.\n");
	int status = 0;

	inst->env = CPXopenCPLEX(&status);
	if (inst->env == NULL) {
		print_warning("[setup_CPLEX_env]: Could not open CPLEX environment.\n");
		CPXgeterrorstring(inst->env, status, errmsg);
		print_error("[setup_CPLEX_env]: %s", errmsg);
	}

	// Set CPLEX parameters
	status = CPXsetintparam(inst->env, CPXPARAM_ScreenOutput, CPX_ON);
	if (status) print_error("[setup_CPLEX_env]: Failed to turn on screen indicator, error %d.\n", status);
	status = CPXsetdblparam(inst->env, CPXPARAM_TimeLimit, 3600);      
	if (status) print_error("[setup_CPLEX_env]: Failed to set time limit, error %d.\n", status);

	// Free
	free(errmsg);
}

void read_MIP_problem(instance* inst, char* filename) {

	// Local variables
	int status = 0;

	// Create MIP from input file (.mps)
	inst->lp = CPXcreateprob(inst->env, &status, filename);
	if (inst->lp == NULL) print_error("[read_MIP_problem]: Failed to create MIP.\n");
	status = CPXreadcopyprob(inst->env, inst->lp, filename, NULL); 
	if (status) print_error("[read_MIP_problem]: Failed to read and copy the problem data.\n");
}

void save_integer_variables(instance* inst) {
	
	// Local variables
	int ncols = CPXgetnumcols(inst->env, inst->lp);
	int status = 0;

	// Allocate variable types
	inst->mip_ctype = (char*)malloc(ncols * sizeof(char));
	if (inst->mip_ctype == NULL) print_error("[save_integer_variables]: Failed to allocate mip_ctype.\n");
	inst->int_var = (int*)malloc(ncols * sizeof(int));
	if (inst->int_var == NULL) print_error("[save_integer_variables]: Failed to allocate int_var.\n");

	// Get MIP variable types {CPX_CONTINUOUS, CPX_BINARY, CPX_INTEGER, CPX_SEMICONT, CPX_SEMIINT}
	status = CPXgetctype(inst->env, inst->lp, inst->mip_ctype, 0, ncols - 1);
	if (status) print_error("[save_integer_variables]: Failed to obtain MIP variable types.\n");

	// Remember integer variables {CPX_BINARY, CPX_INTEGER}
	for (int j = 0; j < ncols; j++) {
		if (inst->mip_ctype[j] == CPX_INTEGER || inst->mip_ctype[j] == CPX_BINARY) {
			inst->int_var[j] = 1;
		}
		else {
			inst->int_var[j] = 0;
			if (inst->mip_ctype[j] != CPX_CONTINUOUS) print_verbose(100, "[INFO][save_integer_variables]: Variable x_%d of type '%c' found.\n", j + 1, inst->mip_ctype[j]);
		}
	}
}

void solve_continuous_relaxation(instance* inst) {

	// Local variables
	int status = 0;

	// MIP --> continuous relaxation
	status = CPXchgprobtype(inst->env, inst->lp, CPXPROB_LP);
	if (status) print_error("[solve_continuous_relaxation]: Failed to change problem type.\n");

	// Optimize LP
	status = CPXlpopt(inst->env, inst->lp);
	if (status) print_error("[solve_continuous_relaxation]: Failed to optimize LP.\n");
}