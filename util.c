/**
 * @file util.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void check_bounds(instance* inst, double* x) {

	// Scan xj variables
	for (int j = 0; j < inst->ncols; j++) {
		if ((x[j] < inst->lb[j] - TOLERANCE) || (x[j] > inst->ub[j] + TOLERANCE)) {
			print_error("[check_bounds]: Bound %d violated: %f <= x_%d %f <= %f\n", j + 1, inst->lb[j], j + 1, inst->x[j], inst->ub[j]);
		}
	}
}

void check_constraints(instance* inst, double* x) {

	// Scan constraints
	for (int i = 0; i < inst->nrows; i++) {

		int rowend = (i < inst->nrows - 1) ? inst->rmatbeg[i + 1] : inst->nzcnt;
		double rowact = 0.0;

		// Scan non-zero coefficients of the constraint and compute row activity
		for (int k = inst->rmatbeg[i]; k < rowend; k++) {

			int varind = inst->rmatind[k];
			assert(index_in_bounds(varind, inst->ncols));
			rowact += (inst->rmatval[k] * x[varind]);
		}

		// Check compliance with the constraint sense
		switch (inst->sense[i]) {
			case 'L':
				if (rowact > inst->rhs[i] + TOLERANCE) print_error("[check_constraints]: Constraint %d violated: rowact %f <= rhs %f\n", i + 1, rowact, inst->rhs[i]);
				break;
			case 'G':
				if (rowact < inst->rhs[i] - TOLERANCE) print_error("[check_constraints]: Constraint %d violated: rowact %f >= rhs %f\n", i + 1, rowact, inst->rhs[i]);
				break;
			case 'E':
				if (fabs(rowact - inst->rhs[i]) > TOLERANCE) print_error("[check_constraints]: Constraint %d violated: rowact %f = rhs %f\n", i + 1, rowact, inst->rhs[i]);
				break;
			default:
				print_error("[check_constraints]: Constraint sense '%c' not supported.\n", inst->sense[i]);
		}
	}
	print_verbose(100, "[check_constraints][OK]: Constraints satisfied.\n");
}

void check_rounding(instance* inst) {

	int rounded = 1;

	// Scan integer/binary variables
	for (int j = 0; j < inst->ncols; j++) {
		if (!(inst->int_var[j])) continue;
		assert(var_type_integer_or_binary(inst->vartype[j]));

		if (is_fractional(inst->x[j])) {
			print_warning("[check_rounding]: Variable (type '%c') x_%d = %f has not been rounded!\n", inst->vartype[j], j + 1, inst->x[j]);
			rounded = 0;
		}
	}

	if (!rounded) print_error("[check_rounding]: ... Failed to round all integer/binary variables of the MIP ...\n");
}

double fractionality(double xj) {
	double minusfloor = xj - floor(xj);
	double minusceil = ceil(xj) - xj;
	return min(minusfloor, minusceil);
}

double sol_fractionality(double* x, int* int_var, int len) {
	double solfrac = 0.0;
	for (int j = 0; j < len; j++) {
		if (!(int_var[j])) continue;
		solfrac += fractionality(x[j]);
	}
	return solfrac;
}

int is_fractional(double num) { return 1 - (fabs(num - round(num)) < TOLERANCE); }

double dot_product(double* coef, double* var_value, int len) {
	double dotprod = 0;
	for (int j = 0; j < len; j++) dotprod += (coef[j]) * (var_value[j]);
	return dotprod;
}

void clone_array(double* arr, double* clo, int len) { for (int i = 0; i < len; i++) clo[i] = arr[i]; }