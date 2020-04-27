/**
 * @file read_ziround_input.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void populate_inst(instance* inst) {

	// Read problem sizes
	inst->nrows = CPXgetnumrows(inst->env, inst->lp);
	inst->ncols = CPXgetnumcols(inst->env, inst->lp);

	read_solution(inst);

	read_variable_bounds(inst);

	// Read objective info
	inst->objsen = CPXgetobjsen(inst->env, inst->lp);
	read_objective_value(inst);
	read_objective_coefficients(inst);

	// Read constraints info
	read_constraints_coefficients(inst);
	read_constraints_senses(inst);
	read_constraints_right_hand_sides(inst);
	read_row_slacks(inst);
}

void read_solution(instance* inst) {

	// Local variables
	int solstat;   /**< Solution status according to CPLEX. */
	int solmethod; /**< Solution method according to CPLEX. */
	int soltype;   /**< Solution type according to CPLEX. */
	int status = 0;

	// Allocate solution
	inst->x = (double*)malloc(inst->ncols * sizeof(double));
	if (inst->x == NULL) print_error("[read_solution]: Failed to allocate solution.\n");

	// Get solution status
	solstat = CPXgetstat(inst->env, inst->lp);

	switch (solstat) {

		case CPX_STAT_UNBOUNDED:
			print_warning("[read_solution]: Model is unbounded.\n"); 
			exit(EXIT_FAILURE);

		case CPX_STAT_INFEASIBLE:
			print_warning("[read_solution]: Model is infeasible.\n"); 
			exit(EXIT_FAILURE);

		case CPX_STAT_INForUNBD:
			print_warning("[read_solution]: Model is infeasible or unbounded.\n"); 
			exit(EXIT_FAILURE);

		default:
			break;
	}

	// Get solution info
	status = CPXsolninfo(inst->env, inst->lp, &(solmethod), &(soltype), NULL, NULL);
	if (status) print_error("[read_solution]: Failed to obtain solution info.\n");
	if (soltype == CPX_NO_SOLN) print_error("[read_solution]: Solution not available.\n");
	print_verbose(150, "Solution status %d, solution method %d.\n", solstat, solmethod);

	// Get solution
	status = CPXgetx(inst->env, inst->lp, inst->x, 0, inst->ncols - 1);
	if (status) print_error("[read_solution]: Failed to obtain primal solution.\n");
}

void read_variable_bounds(instance* inst) {

	// Local variables
	int status = 0;

	// Allocate variable bounds
	inst->ub = (double*)malloc(inst->ncols * sizeof(double));
	inst->lb = (double*)malloc(inst->ncols * sizeof(double));
	if (inst->ub == NULL || inst->lb == NULL) print_error("[read_variable_bounds]: Failed to allocate variable bounds.\n");

	// Get variable bounds (upper and lower)
	status = CPXgetub(inst->env, inst->lp, inst->ub, 0, inst->ncols - 1);
	if (status) print_error("[read_variable_bounds]: Failed to obtain upper bounds.\n");
	status = CPXgetlb(inst->env, inst->lp, inst->lb, 0, inst->ncols - 1);
	if (status) print_error("[read_variable_bounds]: Failed to obtain lower bounds.\n");
}

void read_objective_value(instance* inst) {

	// Local variables
	int status = 0;

	// Get objective value
	status = CPXgetobjval(inst->env, inst->lp, &(inst->objval));
	if (status) print_error("[read_objective_value]: Failed to obtain objective value.\n");
}

void read_objective_coefficients(instance* inst) {

	// Local variables
	int status = 0;

	// Allocate objective coefficients
	inst->obj = (double*)malloc(inst->ncols * sizeof(double));
	if (inst->obj == NULL) print_error("[read_objective_coefficients]: Failed to allocate objective coefficients.\n");

	// Get objective coefficients
	status = CPXgetobj(inst->env, inst->lp, inst->obj, 0, inst->ncols - 1);
	if (status) print_error("[read_objective_coefficients]: Failed to obtain objective coefficients.\n");
}

void read_constraints_coefficients(instance* inst) {

	// Local variables
	int unused = 0;
	int status = 0;

	// First, get the number of non zero coefficients of the matrix (nzcnt)
	inst->nzcnt = CPXgetnumnz(inst->env, inst->lp);

	// Allocate constraints info
	inst->rmatbeg = (int*)malloc(inst->nrows * sizeof(int));
	inst->rmatind = (int*)malloc(inst->nzcnt * sizeof(int));
	inst->rmatval = (double*)malloc(inst->nzcnt * sizeof(double));
	inst->cmatbeg = (int*)malloc(inst->ncols * sizeof(int));
	inst->cmatind = (int*)malloc(inst->nzcnt * sizeof(int));
	inst->cmatval = (double*)malloc(inst->nzcnt * sizeof(double));
	if (inst->rmatbeg == NULL || inst->rmatind == NULL || inst->rmatval == NULL ||
		inst->cmatbeg == NULL || inst->cmatind == NULL || inst->cmatval == NULL) {
		print_error("[read_constraints_coefficients]: Failed to allocate one of rmatbeg, rmatind, rmatval, cmatbeg, cmatind, cmatval.\n");
	}

	// Get constraint matrix, both by rows and by columns
	status = CPXgetrows(inst->env, inst->lp, &unused, inst->rmatbeg, inst->rmatind, inst->rmatval, inst->nzcnt, &unused, 0, inst->nrows - 1);
	if (status) print_error("[read_constraints_coefficients]: Failed to obtain rows info.\n");
	status = CPXgetcols(inst->env, inst->lp, &unused, inst->cmatbeg, inst->cmatind, inst->cmatval, inst->nzcnt, &unused, 0, inst->ncols - 1);
	if (status) print_error("[read_constraints_coefficients]: Failed to obtain columns info.\n");
}

void read_constraints_senses(instance* inst) {

	// Local variables
	int status = 0;

	// Allocate constraint senses
	inst->sense = (char*)malloc(inst->nrows * sizeof(char));
	if (inst->sense == NULL) print_error("[read_constraints_senses]: Failed to allocate constraint senses.\n");

	// Get constraint senses {'L','E','G'}
	status = CPXgetsense(inst->env, inst->lp, inst->sense, 0, inst->nrows - 1);
	if (status) print_error("[read_constraints_senses]: Failed to obtain constraints senses.\n");
}

void read_constraints_right_hand_sides(instance* inst) {

	// Local variables
	int status = 0;

	// Allocate constraint right hand sides
	inst->rhs = (double*)malloc(inst->nrows * sizeof(double));
	if (inst->rhs == NULL) print_error("[read_constraints_right_hand_sides]: Failed to allocate right hand sides.\n");

	// Get right hand sides
	status = CPXgetrhs(inst->env, inst->lp, inst->rhs, 0, inst->nrows - 1);
	if (status) print_error("[read_constraints_right_hand_sides]: Failed to obtain rhs.\n");
}

void read_row_slacks(instance* inst) {

	// Local variables
	int status = 0;

	// Allocate row slacks
	inst->slack = (double*)malloc(inst->nrows * sizeof(double));
	if (inst->slack == NULL) print_error("[read_row_slacks]: Failed to allocate slacks.\n");

	// Get row slacks
	status = CPXgetslack(inst->env, inst->lp, inst->slack, 0, inst->nrows - 1);
	if (status) print_error("[read_row_slacks]: Failed to obtain slacks.\n");
}