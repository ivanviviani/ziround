#include "ziround.h"

void initialize_instance(instance* inst) {
	inst->x = NULL; inst->obj = NULL; inst->lb = NULL; inst->ub = NULL;
	inst->slack = NULL; inst->objsen = CPX_MIN; inst->mip_ctype = NULL;
	inst->int_var = NULL; inst->rmatbeg = NULL; inst->rmatind = NULL;
	inst->rmatval = NULL; inst->cmatbeg = NULL; inst->cmatind = NULL;
	inst->cmatval = NULL; inst->sense = NULL; inst->rhs = NULL; inst->UB = NULL;
	inst->LB = NULL; inst->updated = 0; // inst->x_prev = NULL; inst->x_updated = NULL; 
	inst->status = 0; inst->env = NULL; inst->lp = NULL; inst->input_file = NULL;
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