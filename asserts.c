/**
 * @file asserts.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

int positive_integer(int num) {

	return (num > 0);
}

int non_negative_integer(int num) {

	return (num >= 0);
}

int non_negative_double(double num) {

	return (num >= 0.0 - TOLERANCE);
}

int non_positive_double(double num) {

	return (num <= 0.0 + TOLERANCE);
}

int zero_double(double num) {

	return (fabs(num) < TOLERANCE);
}

int equals_double(double x, double y) {

	return zero_double(x - y);
}

int index_in_bounds(int ind, int len) {

	return ((ind >= 0) && (ind < len));
}

int valid_obj_sense(int objsen) {

	return ((objsen == CPX_MIN) || (objsen == CPX_MAX));
}

int no_ranged_constraints(char* sense, int nrows) {

	for (int i = 0; i < nrows; i++) {
		if (sense[i] == 'R') return (0);
	}
	return (1);
}

int valid_row_slacks(double* slack, char* sense, int nrows) {

	for (int i = 0; i < nrows; i++) {
		switch (sense[i]) {
			case 'L': // row slack must be non-negative
				if (slack[i] < -(TOLERANCE)) return (0);
				break;
			case 'G': // row slack must be non-positive
				if (slack[i] > TOLERANCE) return (0);
				break;
			case 'E': // row slack must be zero
				if (fabs(slack[i]) > TOLERANCE) return (0);
				break;
			default:
				return (0);
		}
	}
	return (1);
}

int valid_var_types(char* vartype, int ncols) {

	for (int j = 0; j < ncols; j++) {
		if (vartype[j] == CPX_SEMICONT || vartype[j] == CPX_SEMIINT) return (0);
	}
	return (1);
}

int var_type_integer_or_binary(char vartype) {

	return ((vartype == CPX_INTEGER) || (vartype == CPX_BINARY));
}

int var_type_continuous(char vartype) {

	return (vartype == CPX_CONTINUOUS);
}

int array_of_zeros(int* arr, int len) {

	for (int i = 0; i < len; i++) {
		if (arr[i] != 0) return (0);
	}
	return (1);
}

int valid_bounds(double* lower, double* upper, int nvars) {

	for (int i = 0; i < nvars; i++) {
		if (lower[i] > upper[i]) return (0);
	}
	return (1);
}

int var_in_bounds(double var, double lb, double ub) {

	return ((var > lb - TOLERANCE) && (var < ub + TOLERANCE));
}