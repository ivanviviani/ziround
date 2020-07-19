/**
 * @file instance.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void init_inst(instance* inst) {
	
	null_all(27, inst->x, inst->obj, inst->lb, inst->ub, inst->slack, inst->vartype,
		inst->int_var, inst->row_singletons, inst->num_singletons, inst->rs_beg,
		inst->rs_coef, inst->ss_val, inst->ss_ub, inst->ss_lb, inst->rmatbeg,
		inst->rmatind, inst->rmatval, inst->cmatbeg, inst->cmatind, inst->cmatval,
		inst->sense, inst->rhs, inst->tracker_sol_frac, inst->tracker_sol_cost,
		inst->tracker_toround, inst->env, inst->lp);
	inst->objsen           = CPX_MIN;
	inst->extension        = 0;
	inst->timelimit        = 3600;
	inst->rseed            = -1;
	strcpy(inst->input_file, "NULL");
	strcpy(inst->input_folder, "NULL");
}

void free_inst(instance* inst) {

	char* errmsg = (char*)malloc(CPXMESSAGEBUFSIZE * sizeof(char)); if (errmsg == NULL) print_error("[setup_CPLEX_env]: Failed to allocate errmsg.\n");
	int status = 0;
	
	free_all(25, inst->x, inst->obj, inst->lb, inst->ub, inst->slack, inst->vartype, 
		inst->int_var, inst->row_singletons, inst->num_singletons, inst->rs_beg, 
		inst->rs_coef, inst->ss_val, inst->ss_ub, inst->ss_lb, inst->rmatbeg, inst->rmatind,
		inst->rmatval, inst->cmatbeg, inst->cmatind, inst->cmatval, inst->sense,
		inst->rhs, inst->tracker_sol_frac, inst->tracker_sol_cost, inst->tracker_toround);
	
	if (inst->lp != NULL) {
		if (CPXfreeprob(inst->env, &(inst->lp))) print_error("[free_inst]: CPXfreeprob failed, error code %d.\n", status);
	}
	if (inst->env != NULL) {
		status = CPXcloseCPLEX(&(inst->env));
		if (status) {
			print_warning("[free_inst]: Could not close CPLEX environment.\n");
			CPXgeterrorstring(inst->env, status, errmsg);
			print_error("[free_inst]: %s\n", errmsg);
		}
	}

	// Free
	free(errmsg);
}