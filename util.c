/**
 * @file util.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void check_bounds(instance* inst, double* x) {

	// Scan xj variables
	for (int j = 0; j < inst->ncols; j++) {
		if ((x[j] - inst->lb[j]) < -(TOLERANCE) || (x[j] - inst->ub[j]) > TOLERANCE) {
			print_error("[check_bounds][!!!]: Bound %d violated!\n", j + 1);
		}
	}
}

void check_constraints(instance* inst, double* x) {

	// Local variables
	double* row_infeas = (double*)malloc(inst->nrows * sizeof(double));
	if (row_infeas == NULL) print_error("[check_constraints]: Failed to allocate row_infeas.\n");
	int status = 0;

	// Compute row infeasibilities
	status = CPXgetrowinfeas(inst->env, inst->lp, x, row_infeas, 0, inst->nrows - 1);
	if (status) print_error("[check_constraints]: Failed to obtain row infeasibilities.\n");

	// Scan rows
	for (int i = 0; i < inst->nrows; i++) {
		if (fabs(row_infeas[i]) > TOLERANCE) {
			print_error("[check_constraints]: Constraint %d violated!\n", i + 1);
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