/**
 * @file instance.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void init_inst(instance* inst) {

	inst->x = NULL; 
	inst->obj = NULL; 
	inst->lb = NULL; 
	inst->ub = NULL;
	inst->slack = NULL; 
	inst->objsen = CPX_MIN; 
	inst->mip_ctype = NULL;
	inst->int_var = NULL; 
	inst->rmatbeg = NULL; 
	inst->rmatind = NULL;
	inst->rmatval = NULL; 
	inst->cmatbeg = NULL; 
	inst->cmatind = NULL;
	inst->cmatval = NULL; 
	inst->sense = NULL; 
	inst->rhs = NULL;
	inst->env = NULL; 
	inst->lp = NULL; 
	inst->input_file = NULL;
}

void free_inst(instance* inst) {

	// Local variables
	char* errmsg = (char*)malloc(CPXMESSAGEBUFSIZE * sizeof(char));
	if (errmsg == NULL) print_error("[setup_CPLEX_env]: Failed to allocate errmsg.\n");
	int status = 0;

	free(inst->x);       
	free(inst->obj);     
	free(inst->lb);
	free(inst->ub);      
	free(inst->slack);   
	free(inst->mip_ctype);
	free(inst->int_var); 
	free(inst->rmatbeg); 
	free(inst->rmatind);
	free(inst->rmatval); 
	free(inst->cmatbeg); 
	free(inst->cmatind);
	free(inst->cmatval); 
	free(inst->sense);   
	free(inst->rhs);
	if (inst->lp != NULL) {
		status = CPXfreeprob(inst->env, &(inst->lp));
		if (status) print_error("[free_inst]: CPXfreeprob failed, error code %d.\n", status);
	}
	if (inst->env != NULL) {
		status = CPXcloseCPLEX(&(inst->env));
		if (status) {
			print_warning("[free_inst]: Could not close CPLEX environment.\n");
			CPXgeterrorstring(inst->env, status, errmsg);
			print_error(stderr, "[free_inst]: %s", errmsg);
		}
	}

	// Free
	free(errmsg);
}