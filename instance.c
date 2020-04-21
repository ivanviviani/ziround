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
	inst->status = 0; 
	inst->env = NULL; 
	inst->lp = NULL; 
	inst->input_file = NULL;
}

int free_instance(instance* inst) {

	// External variables
	CPXENVptr env; /**< CPLEX environment pointer. */
	CPXLPptr lp;   /**< CPLEX lp pointer. */
	// Local variables
	char* errmsg;
	int status;    /**< Support status flag. */

	// Allocate / Initialize
	env = inst->env;
	lp = inst->lp;
	errmsg = (char*)malloc(CPXMESSAGEBUFSIZE * sizeof(char));
	status = 0;

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
		status = CPXfreeprob(env, &lp);
		if (status) print_error("[free_instance]: CPXfreeprob failed, error code %d.\n", status);
	}
	if (env != NULL) {
		status = CPXcloseCPLEX(&env);
		if (status) {
			print_warning("[free_instance]: Could not close CPLEX environment.\n");
			CPXgeterrorstring(env, status, errmsg);
			print_error(stderr, "[free_instance]: %s", errmsg);
		}
	}
	free(errmsg);

	return status;
}