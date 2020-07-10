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
	assert(
		positive_integer(inst->nrows) & 
		positive_integer(inst->ncols)
	);

	read_solution(inst);

	read_variable_bounds(inst);

	// Read objective info
	inst->objsen = CPXgetobjsen(inst->env, inst->lp);
	assert(valid_obj_sense(inst->objsen));
	read_objective_value(inst);
	read_objective_coefficients(inst);

	// Read constraints info
	read_constraints_coefficients(inst);
	read_constraints_senses(inst);
	read_constraints_right_hand_sides(inst);
	read_row_slacks(inst);

	// Extension (if enabled)
	if (inst->extension) {
		find_singletons(inst);
		compute_singletons_slacks_bounds(inst);
	}
}

void read_solution(instance* inst) {

	int solstat;   /**< Solution status according to CPLEX. */
	int solmethod; /**< Solution method according to CPLEX. */
	int soltype;   /**< Solution type according to CPLEX. */

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
	if (soltype == CPX_NO_SOLN) print_error("[read_solution]: Solution not available.\n");
	print_verbose(150, "Solution status %d, solution method %d.\n", solstat, solmethod);

	// Get solution
	if (CPXgetx(inst->env, inst->lp, inst->x, 0, inst->ncols - 1)) print_error("[read_solution]: Failed to obtain primal solution.\n");
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
	if (VERBOSE >= 201) {
		printf("\n\n");
		for (int j = 0; j < inst->ncols; j++) {
			printf("%f ", inst->obj[j]);
		}
		printf("\n\n");
	}
}

void read_constraints_coefficients(instance* inst) {

	int unused = 0;

	// First, get the number of non zero coefficients of the matrix (nzcnt)
	inst->nzcnt = CPXgetnumnz(inst->env, inst->lp);
	assert(positive_integer(inst->nzcnt));

	// Allocate constraints info
	inst->rmatbeg = (int*)malloc(inst->nrows * sizeof(int));
	inst->rmatind = (int*)malloc(inst->nzcnt * sizeof(int));
	inst->rmatval = (double*)malloc(inst->nzcnt * sizeof(double));
	inst->cmatbeg = (int*)malloc(inst->ncols * sizeof(int));
	inst->cmatind = (int*)malloc(inst->nzcnt * sizeof(int));
	inst->cmatval = (double*)malloc(inst->nzcnt * sizeof(double));
	if (inst->rmatbeg == NULL || inst->rmatind == NULL || inst->rmatval == NULL || inst->cmatbeg == NULL || inst->cmatind == NULL || inst->cmatval == NULL) {
		print_error("[read_constraints_coefficients]: Failed to allocate one of rmatbeg, rmatind, rmatval, cmatbeg, cmatind, cmatval.\n");
	}

	// Get constraint matrix, both by rows and by columns
	if (CPXgetrows(inst->env, inst->lp, &unused, inst->rmatbeg, inst->rmatind, inst->rmatval, inst->nzcnt, &unused, 0, inst->nrows - 1)) print_error("[read_constraints_coefficients]: Failed to obtain rows info.\n");
	if (CPXgetcols(inst->env, inst->lp, &unused, inst->cmatbeg, inst->cmatind, inst->cmatval, inst->nzcnt, &unused, 0, inst->ncols - 1)) print_error("[read_constraints_coefficients]: Failed to obtain columns info.\n");
}

void read_constraints_senses(instance* inst) {

	// Allocate constraint senses
	inst->sense = (char*)malloc(inst->nrows * sizeof(char)); if (inst->sense == NULL) print_error("[read_constraints_senses]: Failed to allocate constraint senses.\n");

	// Get constraint senses {'L','E','G'}
	if (CPXgetsense(inst->env, inst->lp, inst->sense, 0, inst->nrows - 1)) print_error("[read_constraints_senses]: Failed to obtain constraints senses.\n");
	assert(no_ranged_constraints(inst->sense, inst->nrows));

	// [DEBUG ONLY] Print constraints senses
	if (VERBOSE >= 201) {
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
	inst->slack = (double*)malloc(inst->nrows * sizeof(double)); if (inst->slack == NULL) print_error("[read_row_slacks]: Failed to allocate row slacks.\n");

	// Get row slacks
	if (CPXgetslack(inst->env, inst->lp, inst->slack, 0, inst->nrows - 1)) print_error("[read_row_slacks]: Failed to obtain slacks.\n");
	assert(valid_row_slacks(inst->slack, inst->sense, inst->nrows));
}

// [EXTENSION]
void find_singletons(instance* inst) {

	int colend;    /**< Index of the last constraint containing variable x_j. */
	int rowind;	   /**< Index of the current constraint. */
	int* count;    /**< Support structure for counting singletons for each row. */
	int first = 0; /**< Support flag. */
	int prev = -1; /**< Index of previous row with beg index. */
	int index;     /**< Index of the current singleton. */
	int beg;       /**< Index of the first singleton of a given row. */
	int offset;    /**< Offset for a singleton of a given row. */

	// Allocate
	inst->num_singletons = (int*)calloc((size_t)inst->nrows, sizeof(int)); 
	count = (int*)calloc((size_t)inst->nrows, sizeof(int)); if (inst->num_singletons == NULL || count == NULL) print_error("[find_singletons][extension]: Failed to allocate num_singletons or count.\n");

	// Count number of singletons for each row (scan continuous variables)
	inst->rs_size = 0; // Total number of singletons
	for (int j = 0; j < inst->ncols; j++) {
		
		// Skip non-continuous variables and FIXED variables (lb = ub)
		if ((inst->vartype[j] != CPX_CONTINUOUS) || (fabs(inst->ub[j] - inst->lb[j]) < TOLERANCE)) continue;
		assert(var_type_continuous(inst->vartype[j]));

		// Row index of the last constraint in which x_j appears
		colend = (j < inst->ncols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;

		// If the variable appears in only one constraint
		if (inst->cmatbeg[j] == colend - 1) {

			rowind = inst->cmatind[inst->cmatbeg[j]];
			assert(index_in_bounds(rowind, inst->nrows));
			print_verbose(200, "[find_singletons][extension]: x_%d = %f in constraint %d ('%c')\n", j + 1, inst->x[j], rowind, inst->sense[rowind]);
			inst->num_singletons[rowind]++;
			count[rowind]++;
			inst->rs_size++;
		}
	}
	assert(non_negative_integer(inst->rs_size));
	print_verbose(120, "[find_singletons][extension]: Total number of singletons = %d\n", inst->rs_size);

	// Allocate / Initialize
	inst->row_singletons = (int*)malloc((size_t)inst->rs_size * sizeof(int));
	inst->rs_beg = (int*)malloc((size_t)inst->nrows * sizeof(int));
	inst->rs_coef = (double*)calloc((size_t)inst->rs_size, sizeof(double)); if (inst->row_singletons == NULL || inst->rs_beg == NULL || inst->rs_coef == NULL) print_error("[find_singletons]: Failed to allocate row_singletons or rs_beg or rs_coef\n");
	for (int k = 0; k < inst->rs_size; k++) inst->row_singletons[k] = -1;
	for (int i = 0; i < inst->nrows; i++) inst->rs_beg[i] = -1;

	// Populate row singletons begin indices
	for (int i = 0; i < inst->nrows; i++) {

		// Skip rows that have no singletons
		if (inst->num_singletons[i] == 0) continue;
		assert(positive_integer(inst->num_singletons[i]));

		if (first == 0) {
			inst->rs_beg[i] = 0;
			first = 1;
			prev = i;
		}
		else {
			index = inst->rs_beg[prev] + inst->num_singletons[prev]; 
			assert(index_in_bounds(index, inst->rs_size));
			inst->rs_beg[i] = index;
			prev = i;
		}
		
		// [DEBUG ONLY] Print row singletons begin indices
		print_verbose(200, "[DEBUG][find_singletons][extension]: Row %d | %d singletons | rs_beg = %d\n", i, inst->num_singletons[i], inst->rs_beg[i]);
	}
	// [DEBUG ONLY] Print size of row singletons array
	print_verbose(200, "[DEBUG][find_singletons][extension]: rs_size = %d\n", inst->rs_size);

	// Populate singleton indices and coefficients for each row
	for (int j = 0; j < inst->ncols; j++) {

		// Skip non-continuous variables and FIXED variables (lb = ub)
		if ((inst->vartype[j] != CPX_CONTINUOUS) || (fabs(inst->ub[j] - inst->lb[j]) < TOLERANCE)) continue;
		assert(var_type_continuous(inst->vartype[j]));

		colend = (j < inst->ncols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;

		// If the variable appears in only one constraint
		if (inst->cmatbeg[j] == colend - 1) {

			rowind = inst->cmatind[colend - 1];                    
			assert(index_in_bounds(rowind, inst->nrows));

			offset = inst->num_singletons[rowind] - count[rowind]; 
			assert(non_negative_integer(offset));

			beg = inst->rs_beg[rowind];                            
			assert(index_in_bounds(beg + offset, inst->rs_size));

			inst->row_singletons[beg + offset] = j;
			inst->rs_coef[beg + offset] = inst->cmatval[colend - 1];
			count[rowind]--;
		}
	}
	assert(array_of_zeros(count, inst->nrows));
	free(count);

	// [DEBUG ONLY] Print row singletons (indices)
	if (VERBOSE >= 201) {
		fprintf(stdout, "\n[DEBUG][find_singletons][extension]: Row singletons (index | coef):\n");
		for (int i = 0; i < inst->nrows; i++) {
			fprintf(stdout, "[DEBUG][find_singletons][extension]: Row %d: ", i);
			if (inst->num_singletons[i] == 0) fprintf(stdout, "-");
			beg = inst->rs_beg[i];
			for (int k = 0; k < inst->num_singletons[i]; k++) {
				fprintf(stdout, "(%d | %f) ", inst->row_singletons[beg + k], inst->rs_coef[beg + k]);
			}
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "\n");
	}
}

// [EXTENSION]
void compute_singletons_slacks_bounds(instance* inst) {

	int beg;             /**< Index of the first singleton of a given row. */
	int singleton_index; /**< Absolute index of the current singleton. */
	double coef;         /**< Coefficient of the current singleton (in its row). */

	// Allocate
	inst->ss_ub = (double*)calloc((size_t)inst->nrows, sizeof(double));
	inst->ss_lb = (double*)calloc((size_t)inst->nrows, sizeof(double)); if (inst->ss_ub == NULL || inst->ss_lb == NULL) print_error("[compute_singletons_slacks_bounds][extension]: Failed to allocate ss_ub or ss_lb.");

	// Scan constraints that have singletons
	for (int i = 0; i < inst->nrows; i++) {

		// Skip rows with no singletons
		if (inst->num_singletons[i] == 0) continue;
		assert(positive_integer(inst->num_singletons[i]));

		beg = inst->rs_beg[i]; 
		assert(index_in_bounds(beg, inst->rs_size));

		// Compute singletons slack upper and lower bound (row i)
		for (int k = 0; k < inst->num_singletons[i]; k++) {

			assert(index_in_bounds(beg + k, inst->rs_size));
			singleton_index = inst->row_singletons[beg + k]; 
			assert(index_in_bounds(singleton_index, inst->ncols));
			coef = inst->rs_coef[beg + k];

			if (coef > 0.0) {
				inst->ss_ub[i] += (coef * inst->ub[singleton_index]);
				inst->ss_lb[i] += (coef * inst->lb[singleton_index]);
			}
			if (coef < 0.0) {
				inst->ss_ub[i] += (coef * inst->lb[singleton_index]);
				inst->ss_lb[i] += (coef * inst->ub[singleton_index]);
			}
		}

		// [DEBUG ONLY] Print singletons slacks bounds
		print_verbose(200, "[DEBUG][compute_singletons_slacks_bounds][extension][row %d]: ss_lb = %f | ss_ub = %f\n", i + 1, inst->ss_lb[i], inst->ss_ub[i]);
	}
	assert(valid_bounds(inst->ss_lb, inst->ss_ub, inst->nrows));
}