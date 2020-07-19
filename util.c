/**
 * @file util.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void check_bounds(double* x, double* lb, double* ub, int ncols) {

	// Scan xj variables
	for (int j = 0; j < ncols; j++) {

		// Terminate program at the first violated bounds
		if (!var_in_bounds(x[j], lb[j], ub[j])) print_error("[check_bounds]: Some variable bounds are violated!\n");
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
				if (rowact > rhs[i] + TOLERANCE) violated = 1;
				break;
			case 'G':
				if (rowact < rhs[i] - TOLERANCE) violated = 1;
				break;
			case 'E':
				if (fabs(rowact - rhs[i]) > TOLERANCE) violated = 1;
				break;
			default:
				print_error("[check_constraints]: Constraint sense '%c' not supported.\n", sense[i]);
		}

		// Terminate program at the first violated constraint
		if (violated) print_error("[check_constraints]: Some constraints are violated!\n");
	}

	print_verbose(100, "[check_constraints][OK]: Constraints satisfied.\n");
}

int check_rounding(double* x, int ncols, int* int_var, char* vartype) {

	// Scan integer/binary variables
	for (int j = 0; j < ncols; j++) {
		if (!(int_var[j])) continue;
		assert(var_type_integer_or_binary(vartype[j]));

		// Stop at the first integer variable not rounded
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

	return min(xj - floor(xj), ceil(xj) - xj);
}

double sol_fractionality(double* x, int* int_var, int len) {

	double solfrac = 0.0;

	for (int j = 0; j < len; j++) {
		if (!(int_var[j])) continue;

		solfrac += fractionality(x[j]);
	}

	return solfrac;
}

int is_fractional(double num) { 

	return 1 - (fabs(num - round(num)) < TOLERANCE); 
}

double dot_product(double* coef, double* var_value, int len) {

	double dotprod = 0;

	for (int j = 0; j < len; j++) dotprod += (coef[j]) * (var_value[j]);

	return dotprod;
}

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