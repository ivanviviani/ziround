#include "ziround.h"

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