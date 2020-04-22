/**
 * @file util.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void check_bounds(instance* inst, double* x) {

	// External variables
	int ncols;  /**< Number of variables (columns) of the MIP problem. */
	double* lb; /**< Lower bounds array. */
	double* ub; /**< Upper bounds array. */

	// Initialize
	ncols = inst->ncols;
	lb = inst->lb;
	ub = inst->ub;

	// Scan xj variables
	for (int j = 0; j < ncols; j++) {
		if ((x[j] - lb[j]) < -(TOLERANCE) || (x[j] - ub[j]) > TOLERANCE) {
			print_error("[check_bounds][!!!]: Bound %d violated!\n", j + 1);
		}
	}
}

void check_constraints(instance* inst, double* x) {

	// External variables
	int nrows;			/**< Number of constraints (rows) of the MIP problem.*/
	CPXLPptr lp;		/**< CPLEX lp pointer. */
	CPXENVptr env;		/**< CPLEX environment pointer. */
	// Local variables
	double* row_infeas; /**< Row infeasibilities array. */
	int status;			/**< Support status flag. */

	// Allocate / Initialize
	nrows = inst->nrows;
	lp = inst->lp;
	env = inst->env;
	row_infeas = (double*)malloc(nrows * sizeof(double));
	if (row_infeas == NULL) print_error("[check_constraints]: Failed to allocate row_infeas.\n");
	status = 0;

	// Compute row infeasibilities
	status = CPXgetrowinfeas(env, lp, x, row_infeas, 0, nrows - 1);
	if (status) print_error("[check_constraints]: Failed to obtain row infeasibilities.\n");

	// Scan rows
	for (int i = 0; i < nrows; i++) {
		if (fabs(row_infeas[i]) > TOLERANCE) {
			print_error("[check_constraints][!!!]: Constraint %d violated!\n", i + 1);
		}
	}
	
	// Free
	free(row_infeas);
}

double fractionality(double xj) {
	double minusfloor = xj - floor(xj);
	double minusceil = ceil(xj) - xj;
	return min(minusfloor, minusceil);
}

int is_fractional(double num) { return 1 - (fabs(num - round(num)) < TOLERANCE); }

double dot_product(double* coef, double* var_value, int len) {
	double dotprod = 0;
	for (int j = 0; j < len; j++) dotprod += (coef[j]) * (var_value[j]);
	return dotprod;
}

void clone_array(double* arr, double* clo, int len) { for (int i = 0; i < len; i++) clo[i] = arr[i]; }