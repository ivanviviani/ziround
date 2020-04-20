#include "ziround.h"

void read_problem_sizes(instance* inst) {
	inst->objsen = CPXgetobjsen(inst->env, inst->lp);
	inst->cur_numrows = CPXgetnumrows(inst->env, inst->lp);
	inst->cur_numcols = CPXgetnumcols(inst->env, inst->lp);
}

int read_solution(instance* inst) {
	// Get solution
	inst->solnstat = CPXgetstat(inst->env, inst->lp);
	if (inst->solnstat == CPX_STAT_UNBOUNDED) { fprintf(stdout, "[INFO][read_solution]: Model is unbounded.\n"); return inst->status; }
	else if (inst->solnstat == CPX_STAT_INFEASIBLE) { fprintf(stdout, "[INFO][read_solution]: Model is infeasible.\n"); return inst->status; }
	else if (inst->solnstat == CPX_STAT_INForUNBD) { fprintf(stdout, "[INFO][read_solution]: Model is infeasible or unbounded.\n"); return inst->status; }
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