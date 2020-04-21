#include "ziround.h"

void read_problem_sizes(instance* inst) {

	// External variables
	CPXENVptr env; /**< CPLEX environment pointer. */
	CPXLPptr lp;   /**< CPLEX lp pointer. */

	// Initialize
	env = inst->env;
	lp = inst->lp;

	inst->objsen = CPXgetobjsen(env, lp);
	inst->nrows = CPXgetnumrows(env, lp);
	inst->ncols = CPXgetnumcols(env, lp);
}

int read_solution(instance* inst) {

	// External variables
	CPXENVptr env; /**< CPLEX environment pointer. */
	CPXLPptr lp;   /**< CPLEX lp pointer. */
	int ncols;     /**< Number of variables (columns) of the MIP problem. */
	double* x;     /**< Initial continuous relaxation solution. */
	// Local variables
	int solstat;   /**< Solution status according to CPLEX. */
	int solmethod; /**< Solution method according to CPLEX. */
	int soltype;   /**< Solution type according to CPLEX. */
	int status;    /**< Support status flag. */

	// Initialize
	env = inst->env;
	lp = inst->lp;
	ncols = inst->ncols;
	x = NULL;
	status = 0;

	// Allocate on the instance, assign to locals
	inst->x = (double*)malloc(ncols * sizeof(double));
	x = inst->x;
	if (x == NULL) print_error("[read_solution]: Failed to allocate solution.\n");

	// Get solution status
	solstat = CPXgetstat(env, lp);
	switch (solstat) {
		case CPX_STAT_UNBOUNDED:
			print_warning("[read_solution]: Model is unbounded.\n"); 
			return 1;
		case CPX_STAT_INFEASIBLE:
			print_warning("[read_solution]: Model is infeasible.\n"); 
			return 1;
		case CPX_STAT_INForUNBD:
			print_warning("[read_solution]: Model is infeasible or unbounded.\n"); 
			return 1;
		default:
			break;
	}

	// Get solution info
	status = CPXsolninfo(env, lp, &(solmethod), &(soltype), NULL, NULL);
	if (status)                 print_error("[read_solution]: Failed to obtain solution info.\n");
	if (soltype == CPX_NO_SOLN) print_error("[read_solution]: Solution not available.\n");
	print_verbose(150, "Solution status %d, solution method %d.\n", solstat, solmethod);

	// Get solution
	status = CPXgetx(env, lp, x, 0, ncols - 1);
	if (status) print_error("[read_solution]: Failed to obtain primal solution.\n");

	return status;
}

int read_variable_bounds(instance* inst) {
	// Get variable bounds (upper and lower)
	inst->ub = (double*)malloc(inst->ncols * sizeof(double));
	inst->lb = (double*)malloc(inst->ncols * sizeof(double));
	if (inst->ub == NULL || inst->lb == NULL) { fprintf(stderr, "[ERR][read_variable_bounds]: Failed to allocate variable bounds.\n"); return 1; }
	inst->status = CPXgetub(inst->env, inst->lp, inst->ub, 0, inst->ncols - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_variable_bounds]: Failed to obtain upper bounds.\n"); return inst->status; }
	inst->status = CPXgetlb(inst->env, inst->lp, inst->lb, 0, inst->ncols - 1);
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
	inst->obj = (double*)malloc(inst->ncols * sizeof(double));
	if (inst->obj == NULL) { fprintf(stderr, "[ERR][read_objective_coefficients]: Failed to allocate objective coefficients.\n"); return 1; }
	inst->status = CPXgetobj(inst->env, inst->lp, inst->obj, 0, inst->ncols - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_objective_coefficients]: Failed to obtain objective coefficients.\n"); return inst->status; }
	return inst->status;
}

int read_constraints_coefficients(instance* inst) {
	// Get constraint matrix, both by rows and by columns
	// First, get the number of non zero coefficients of the matrix (nzcnt)
	int unused = 0;
	inst->nzcnt = CPXgetnumnz(inst->env, inst->lp);
	// Get rows
	inst->rmatbeg = (int*)malloc(inst->nrows * sizeof(int));
	inst->rmatind = (int*)malloc(inst->nzcnt * sizeof(int));
	inst->rmatval = (double*)malloc(inst->nzcnt * sizeof(double));
	if (inst->rmatbeg == NULL || inst->rmatind == NULL || inst->rmatval == NULL) { fprintf(stderr, "[ERR][read_constraints_coefficients]: Failed to allocate one of rmatbeg, rmatind, rmatval.\n"); return 1; }
	inst->status = CPXgetrows(inst->env, inst->lp, &unused, inst->rmatbeg, inst->rmatind, inst->rmatval, inst->nzcnt, &unused, 0, inst->nrows - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_constraints_coefficients]: Failed to obtain rows info.\n"); return inst->status; }
	// Get columns
	inst->cmatbeg = (int*)malloc(inst->ncols * sizeof(int));
	inst->cmatind = (int*)malloc(inst->nzcnt * sizeof(int));
	inst->cmatval = (double*)malloc(inst->nzcnt * sizeof(double));
	if (inst->cmatbeg == NULL || inst->cmatind == NULL || inst->cmatval == NULL) { fprintf(stderr, "[ERR][read_constraints_coefficients]: Failed to allocate one of cmatbeg, cmatind, cmatval.\n"); return 1; }
	inst->status = CPXgetcols(inst->env, inst->lp, &unused, inst->cmatbeg, inst->cmatind, inst->cmatval, inst->nzcnt, &unused, 0, inst->ncols - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_constraints_coefficients]: Failed to obtain columns info.\n"); return inst->status; }
	return inst->status;
}

int read_constraints_senses(instance* inst) {
	// Get constraint senses {'L','E','G'}
	inst->sense = (char*)malloc(inst->nrows * sizeof(char));
	if (inst->sense == NULL) { fprintf(stderr, "[ERR][read_constraints_senses]: Failed to allocate constraint senses.\n"); return 1; }
	inst->status = CPXgetsense(inst->env, inst->lp, inst->sense, 0, inst->nrows - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_constraints_senses]: Failed to obtain constraints senses.\n"); return inst->status; }
	return inst->status;
}

int read_constraints_right_hand_sides(instance* inst) {
	// Get right hand sides
	inst->rhs = (double*)malloc(inst->nrows * sizeof(double));
	if (inst->rhs == NULL) { fprintf(stderr, "[ERR][read_constraints_right_hand_sides]: Failed to allocate right hand sides.\n"); return 1; }
	inst->status = CPXgetrhs(inst->env, inst->lp, inst->rhs, 0, inst->nrows - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_constraints_right_hand_sides]: Failed to obtain rhs.\n"); return inst->status; }
	return inst->status;
}

int read_row_slacks(instance* inst) {
	// Get row slacks
	inst->slack = (double*)malloc(inst->nrows * sizeof(double));
	if (inst->slack == NULL) { fprintf(stderr, "[ERR][read_row_slacks]: Failed to allocate slacks.\n"); return 1; }
	inst->status = CPXgetslack(inst->env, inst->lp, inst->slack, 0, inst->nrows - 1);
	if (inst->status) { fprintf(stderr, "[ERR][read_row_slacks]: Failed to obtain slacks.\n"); return inst->status; }
	return inst->status;
}

int read_problem_data(instance* inst) {

	// Local variables
	int status; /**< Support status flag. */

	// Initialize
	status = 0;
														// TOLTE SE PRINT ERROR USATA INSIDE...
	status = read_solution(inst);                     if (status) print_error("[read_problem_data]: Error inside read_solution.\n"); return status; }
	status = read_variable_bounds(inst);              if (status) print_error("[read_problem_data]: Error inside read_variable_bounds.\n"); return status; }
	status = read_objective_value(inst);              if (status) print_error("[read_problem_data]: Error inside read_objective_value.\n"); return status; }
	status = read_objective_coefficients(inst);       if (status) print_error("[read_problem_data]: Error inside read_objective_coefficients.\n"); return status; }
	status = read_constraints_coefficients(inst);     if (status) print_error("[read_problem_data]: Error inside read_constraints_coefficients.\n"); return status; }
	status = read_constraints_senses(inst);           if (status) print_error("[read_problem_data]: Error inside read_constraints_senses.\n"); return status; }
	status = read_constraints_right_hand_sides(inst); if (status) print_error("[read_problem_data]: Error inside read_constraints_right_hand_sides.\n"); return status; }
	status = read_row_slacks(inst);                   if (status) print_error("[read_problem_data]: Error inside read_row_slacks.\n"); return status; }
	return status;
}