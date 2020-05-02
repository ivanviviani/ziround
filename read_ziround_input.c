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

	// Extension (if enabled)
	if (inst->extension) {
		extend_row_slacks(inst);
		extend_constraints_senses(inst);
	}
}

void read_solution(instance* inst) {

	int solstat;   /**< Solution status according to CPLEX. */
	int solmethod; /**< Solution method according to CPLEX. */
	int soltype;   /**< Solution type according to CPLEX. */
	int status = 0;

	// Allocate solution
	inst->x = (double*)malloc(inst->ncols * sizeof(double)); if (inst->x == NULL) print_error("[read_solution]: Failed to allocate solution.\n");

	// Get solution status
	solstat = CPXgetstat(inst->env, inst->lp);

	switch (solstat) {
		case CPX_STAT_UNBOUNDED:
			print_warning("[read_solution]: Model is unbounded.\n");               exit(EXIT_FAILURE);
		case CPX_STAT_INFEASIBLE:
			print_warning("[read_solution]: Model is infeasible.\n");              exit(EXIT_FAILURE);
		case CPX_STAT_INForUNBD:
			print_warning("[read_solution]: Model is infeasible or unbounded.\n"); exit(EXIT_FAILURE);
		default:
			break;
	}

	// Get solution info
	if (CPXsolninfo(inst->env, inst->lp, &(solmethod), &(soltype), NULL, NULL)) print_error("[read_solution]: Failed to obtain solution info.\n");
	if (soltype == CPX_NO_SOLN)                                                 print_error("[read_solution]: Solution not available.\n");
	print_verbose(150, "Solution status %d, solution method %d.\n", solstat, solmethod);

	// Get solution
	if (CPXgetx(inst->env, inst->lp, inst->x, 0, inst->ncols - 1))              print_error("[read_solution]: Failed to obtain primal solution.\n");
}

void read_variable_bounds(instance* inst) {

	// Allocate variable bounds
	inst->ub = (double*)malloc(inst->ncols * sizeof(double)); if (inst->ub == NULL) print_error("[read_variable_bounds]: Failed to allocate variables upper bounds.\n");
	inst->lb = (double*)malloc(inst->ncols * sizeof(double)); if (inst->lb == NULL) print_error("[read_variable_bounds]: Failed to allocate variables lower bounds.\n");

	// Get variable bounds (upper and lower)
	if (CPXgetub(inst->env, inst->lp, inst->ub, 0, inst->ncols - 1)) print_error("[read_variable_bounds]: Failed to obtain upper bounds.\n");
	if (CPXgetlb(inst->env, inst->lp, inst->lb, 0, inst->ncols - 1)) print_error("[read_variable_bounds]: Failed to obtain lower bounds.\n");
}

void read_objective_value(instance* inst) {

	// Get objective value
	if (CPXgetobjval(inst->env, inst->lp, &(inst->objval))) print_error("[read_objective_value]: Failed to obtain objective value.\n");
}

void read_objective_coefficients(instance* inst) {

	// Allocate objective coefficients
	inst->obj = (double*)malloc(inst->ncols * sizeof(double)); if (inst->obj == NULL) print_error("[read_objective_coefficients]: Failed to allocate objective coefficients.\n");

	// Get objective coefficients
	if (CPXgetobj(inst->env, inst->lp, inst->obj, 0, inst->ncols - 1)) print_error("[read_objective_coefficients]: Failed to obtain objective coefficients.\n");

	// [DEBUG ONLY] Print objective coefficients
	if (VERBOSE >= 200) {
		printf("\n\n");
		for (int j = 0; j < inst->ncols; j++) {
			printf("%f ", inst->obj[j]);
		}
		printf("\n\n");
	}
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
	if (CPXgetrows(inst->env, inst->lp, &unused, inst->rmatbeg, inst->rmatind, inst->rmatval, inst->nzcnt, &unused, 0, inst->nrows - 1)) print_error("[read_constraints_coefficients]: Failed to obtain rows info.\n");
	if (CPXgetcols(inst->env, inst->lp, &unused, inst->cmatbeg, inst->cmatind, inst->cmatval, inst->nzcnt, &unused, 0, inst->ncols - 1)) print_error("[read_constraints_coefficients]: Failed to obtain columns info.\n");
}

void read_constraints_senses(instance* inst) {

	// Local variables
	int status = 0;

	// Allocate constraint senses
	inst->sense = (char*)malloc(inst->nrows * sizeof(char)); if (inst->sense == NULL) print_error("[read_constraints_senses]: Failed to allocate constraint senses.\n");

	// Get constraint senses {'L','E','G'}
	if (CPXgetsense(inst->env, inst->lp, inst->sense, 0, inst->nrows - 1)) print_error("[read_constraints_senses]: Failed to obtain constraints senses.\n");

	// [DEBUG ONLY] Print constraints senses
	if (VERBOSE >= 120) {
		printf("\n\n");
		for (int i = 0; i < inst->nrows; i++) {
			printf("%c ", inst->sense[i]);
		}
		printf("\n\n");
	}
}

void read_constraints_right_hand_sides(instance* inst) {

	// Allocate constraint right hand sides
	inst->rhs = (double*)malloc(inst->nrows * sizeof(double)); if (inst->rhs == NULL) print_error("[read_constraints_right_hand_sides]: Failed to allocate right hand sides.\n");

	// Get right hand sides
	if (CPXgetrhs(inst->env, inst->lp, inst->rhs, 0, inst->nrows - 1)) print_error("[read_constraints_right_hand_sides]: Failed to obtain rhs.\n");
}

void read_row_slacks(instance* inst) {

	// Allocate row slacks
	inst->slack = (double*)malloc(inst->nrows * sizeof(double)); if (inst->slack == NULL) print_error("[read_row_slacks]: Failed to allocate slacks.\n");

	// Get row slacks
	if (CPXgetslack(inst->env, inst->lp, inst->slack, 0, inst->nrows - 1)) print_error("[read_row_slacks]: Failed to obtain slacks.\n");
}

void extend_row_slacks(instance* inst) {

	double aij;        /**< Current coefficient. */
	int rowind;        /**< Current row index. */
	int colend;        /**< Index of the last non-zero coefficient of the current column. */
	int eq_index = -1; /**< Index of the only equality constraint in which a continuous variable is involved. */
	double contrib;    /**< Contribution of the variable in that constraint. */

	// Allocate extended equality constraints flags (initialized to false)
	inst->eq_ext = (int*)calloc(inst->nrows, sizeof(int)); if (inst->eq_ext == NULL) print_error("[extend_row_slacks]: Failed to allocate unique equality constraints flags.\n");

	// Scan continuous variables
	for (int j = 0; j < inst->ncols; j++) {
		if (inst->vartype[j] != CPX_CONTINUOUS) continue; // Skip non-continuous variables

		colend = (j < inst->ncols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;

		// Scan equality constraints in which the j-th continuous variable is involved
		for (int k = inst->cmatbeg[j]; k < colend; k++) {

			rowind = inst->cmatind[k];
			if (inst->sense[rowind] != 'E') continue; // Skip non-equality constraints

			// Save the first row index. Skip variable if involved in more than one equality constraint.
			if (eq_index == -1) eq_index = rowind;
			else {
				eq_index = -2;
				break;
			}
		}

		// If the j-th continuous variable is involved in only one equality constraint
		if (eq_index >= 0) {
			
			// Set equality constraint flag
			inst->eq_ext[eq_index] = 1;

			// Get its contribution in that constraint
			if (CPXgetcoef(inst->env, inst->lp, eq_index, j, &aij)) print_error("[extend_row_slacks]: Failed to obtain coefficient.\n");
			contrib = aij * inst->x[j];

			// Consider its contribution as contribution to the row slack
			inst->slack[eq_index] += contrib;
		}
	}
}

void extend_constraints_senses(instance* inst) {

	// Scan unique equality constraints flags
	for (int i = 0; i < inst->nrows; i++) {
		
		if (!(inst->eq_ext[i])) continue; // Skip non-unique equality constraints

		// Set extended constraint sense
		switch (inst->sense[i]) {

			case 'E':
				if (inst->slack[i] > TOLERANCE)         inst->sense[i] = 'L';
				else if (inst->slack[i] < -(TOLERANCE)) inst->sense[i] = 'G';
				break;
			case 'L':
			case 'G':
				print_error("[extend_constraints_senses]: Tried to extend a non-equality constraint!\n");
			default:
				print_error("[extend_constraints_senses]: Constraint sense %c not supported!\n", inst->sense[i]);
		}
	}

	// [DEBUG ONLY] Print constraints senses
	if (VERBOSE >= 120) {
		printf("\n\n___ Constraints senses after extension ___\n");
		for (int i = 0; i < inst->nrows; i++) {
			printf("%c ", inst->sense[i]);
		}
		printf("\n\n");
	}
}