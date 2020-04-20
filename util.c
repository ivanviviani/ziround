#include "ziround.h"

int check_bounds(instance* inst, double* x) {
	// Check variable bounds
	for (int j = 0; j < inst->cur_numcols; j++)
		// was x[j] < inst->lb[j] || x[j] > inst->ub[j]
		if ((x[j] - inst->lb[j]) < -(TOLERANCE) || (x[j] - inst->ub[j]) > TOLERANCE) { fprintf(stderr, "[ERR][check_bounds][!!!]: Bound %d violated!\n", j + 1); return 1; }
	return 0;
}

int check_constraints(instance* inst, double* x) {
	// Check constraints
	double* row_infeas = (double*)malloc(inst->cur_numrows * sizeof(double));
	if (row_infeas == NULL) { fprintf(stderr, "[ERR][check_constraints]: Failed to allocate row_infeas.\n"); return 1; }
	inst->status = CPXgetrowinfeas(inst->env, inst->lp, x, row_infeas, 0, inst->cur_numrows - 1);
	if (inst->status) { fprintf(stderr, "[ERR][check_constraints]: Failed to obtain row infeasibilities.\n"); free_and_null((char**)&row_infeas); return inst->status; }
	for (int i = 0; i < inst->cur_numrows; i++)
		if (fabs(row_infeas[i]) > TOLERANCE) { fprintf(stdout, "[ERR][check_constraints][!!!]: Constraint %d violated!\n", i + 1); free_and_null((char**)&row_infeas); return 1; }
	free_and_null((char**)&row_infeas);
	return inst->status;
}

double row_activity(instance* inst, int i, double* x) {
	double sum = 0; double aij; int col_index;
	int end_row = (i < inst->cur_numrows - 1) ? inst->rmatbeg[i + 1] : inst->nzcnt;
	// Scan non-zero coefficients of row i
	/*
		For constraint i:
			rmatbeg[i] is the first index of rmatind and rmatval for row i
			--> Row i in range [ rmatbeg[i];rmatbeg[i+1] ) except for the last one (see end_row)
			rmatval contains coefficient values
			rmatind contains the column indices of coefficient values
	*/
	for (int k = inst->rmatbeg[i]; k < end_row; k++) {
		// Get current non-zero aij and its column index
		aij = inst->rmatval[k];
		col_index = inst->rmatind[k];
		// Update row activity
		sum += (aij * x[col_index]);
	}
	return sum;
}

void free_and_null(char** ptr) { if (*ptr != NULL) { free(*ptr); *ptr = NULL; } }

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

int different_arr(double* prev, double* new, int len) {
	int dif = 0;
	for (int i = 0; i < len; i++) {
		if (fabs(prev[i] - new[i]) > TOLERANCE) {
			dif = 1;
			break;
		}
	}
	return dif;
}

void clone_array(double* arr, double* clo, int len) { for (int i = 0; i < len; i++) clo[i] = arr[i]; }