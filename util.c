/**
 * @file util.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void check_bounds(double* x, double* lb, double* ub, int ncols) {

	// Scan xj variables
	for (int j = 0; j < ncols; j++) {
		if ((x[j] < lb[j] - TOLERANCE) || (x[j] > ub[j] + TOLERANCE)) {
			print_error("[check_bounds]: Bound %d violated: %f <= x_%d %f <= %f\n", j + 1, lb[j], j + 1, x[j], ub[j]);
		}
	}
}

void check_constraints(double* x, int ncols, int nrows, int nzcnt, int* rmatbeg, int* rmatind, double* rmatval, char* sense, double* rhs) {

	int rowend;    /**< Last variable index of the current row. */
	double rowact; /**< Current row activity. */
	int varind;    /**< Current variable index. */
	int violated;  /**< Violated constraints flag. */

	violated = 0;

	// Scan constraints
	for (int i = 0; i < nrows; i++) {

		rowend = (i < nrows - 1) ? rmatbeg[i + 1] : nzcnt;
		rowact = 0.0;

		// Scan non-zero coefficients of the constraint and compute row activity
		for (int k = rmatbeg[i]; k < rowend; k++) {

			varind = rmatind[k];
			assert(index_in_bounds(varind, ncols));
			rowact += (rmatval[k] * x[varind]);
		}

		// Check compliance with the constraint sense
		switch (sense[i]) {
			case 'L':
				if (rowact > rhs[i] + TOLERANCE) violated = 1; // print_error("[check_constraints]: Constraint %d violated: rowact %f <= rhs %f\n", i + 1, rowact, rhs[i]);
				break;
			case 'G':
				if (rowact < rhs[i] - TOLERANCE) violated = 1; // print_error("[check_constraints]: Constraint %d violated: rowact %f >= rhs %f\n", i + 1, rowact, rhs[i]);
				break;
			case 'E':
				if (fabs(rowact - rhs[i]) > TOLERANCE) violated = 1; // print_error("[check_constraints]: Constraint %d violated: rowact %f = rhs %f\n", i + 1, rowact, rhs[i]);
				break;
			default:
				print_error("[check_constraints]: Constraint sense '%c' not supported.\n", sense[i]);
		}

		// Terminate program at the first violated constraint
		if (violated) print_error("[check_constraints]: Some constraints are not satisfied!\n");
	}
	print_verbose(100, "[check_constraints][OK]: Constraints satisfied.\n");
}

int check_rounding(double* x, int ncols, int* int_var, char* vartype) {

	// Scan integer/binary variables
	for (int j = 0; j < ncols; j++) {
		if (!(int_var[j])) continue;
		assert(var_type_integer_or_binary(vartype[j]));

		if (is_fractional(x[j])) return 0;
	}

	return 1;
}

int count_rounded(double* x, int ncols, int* int_var, char* vartype) {

	int count = 0;

	// Scan integer/binary variables
	for (int j = 0; j < ncols; j++) {
		if (!(int_var[j])) continue;
		assert(var_type_integer_or_binary(vartype[j]));

		count += (!is_fractional(x[j]));
	}

	return count;
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

void free_all(int count, ...) {
	va_list args;
	va_start(args, count);

	for (int i = 0; i < count; i++) {
		void* elem = va_arg(args, void*);
		free(elem);
	}

	va_end(args);
	fflush(NULL);
}