/**
 * @file print.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void print_warning(const char* warn, ...) {
	printf("\n\nWARNING: ");
	va_list args;
	va_start(args, warn);
	vprintf(warn, args);
	va_end(args);
	fflush(NULL);
}

void print_error(const char* err, ...) {
	printf("\n\nERROR: ");
	va_list args;
	va_start(args, err);
	vprintf(err, args);
	va_end(args);
	fflush(NULL);

	exit(EXIT_FAILURE);
}

void print_verbose(int msg_verb, const char* format, ...) {
	if (VERBOSE >= msg_verb) {
		va_list args;
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
		fflush(NULL);
	}
}

void print_problem_info(instance* inst, int sol_available, int to_file) {
	FILE* out_stream = (to_file) ? fopen("problem_info.txt", "w") : stdout;
	// nrows, ncols
	int nrows = CPXgetnumrows(inst->env, inst->lp);
	int ncols = CPXgetnumcols(inst->env, inst->lp);
	// variable bounds
	double* ub = (double*)malloc(ncols * sizeof(double));
	double* lb = (double*)malloc(ncols * sizeof(double));
	CPXgetub(inst->env, inst->lp, ub, 0, ncols - 1);
	CPXgetlb(inst->env, inst->lp, lb, 0, ncols - 1);
	// obj sense and coefficients
	int objsen = CPXgetobjsen(inst->env, inst->lp);
	double* obj = (double*)malloc(ncols * sizeof(double));
	CPXgetobj(inst->env, inst->lp, obj, 0, ncols - 1);
	// constraints coefficients, rhs and senses
	int nnz = CPXgetnumnz(inst->env, inst->lp); int unused = 0;
	int* rmatbeg = (int*)malloc(nrows * sizeof(int));
	int* rmatind = (int*)malloc(nnz * sizeof(int));
	double* rmatval = (double*)malloc(nnz * sizeof(double));
	CPXgetrows(inst->env, inst->lp, &unused, rmatbeg, rmatind, rmatval, nnz, &unused, 0, nrows - 1);
	int* cmatbeg = (int*)malloc(ncols * sizeof(int));
	int* cmatind = (int*)malloc(nnz * sizeof(int));
	double* cmatval = (double*)malloc(nnz * sizeof(double));
	CPXgetcols(inst->env, inst->lp, &unused, cmatbeg, cmatind, cmatval, nnz, &unused, 0, ncols - 1);
	double* rhs = (double*)malloc(nrows * sizeof(double));
	CPXgetrhs(inst->env, inst->lp, rhs, 0, nrows - 1);
	char* sense = (char*)malloc(nrows * sizeof(char));
	CPXgetsense(inst->env, inst->lp, sense, 0, nrows - 1);
	// solution and obj value
	double* x; double objval;
	if (sol_available) {
		x = (double*)malloc(ncols * sizeof(double));
		clone_array(inst->x, x, ncols);
		objval = inst->objval;
	}
	fprintf(out_stream, "\n********************************************************\n");
	fprintf(out_stream, "********************  PROBLEM INFO  ********************\n");
	fprintf(out_stream, "**** Dimensions (rows, cols): (%d, %d)\n", nrows, ncols);
	fprintf(out_stream, "**** Obj sense: %s\n", ((objsen > 0) ? "MIN" : "MAX"));
	fprintf(out_stream, "**** Obj coefficients:\n");
	for (int j = 0; j < ncols; j++) fprintf(out_stream, "****     x_%d: %f\n", j + 1, obj[j]);
	fprintf(out_stream, "**** Constraints: (nnz = %d)\n", nnz);
	for (int i = 0; i < nrows - 1; i++) {
		fprintf(out_stream, "****     ");
		for (int k = rmatbeg[i]; k < rmatbeg[i + 1]; k++) fprintf(out_stream, "(%f)x_%d ", rmatval[k], rmatind[k] + 1);
		fprintf(out_stream, " %s %f\n", ((sense[i] == 'L') ? "<=" : ((sense[i] == 'G') ? ">=" : "=")), rhs[i]);
	}
	fprintf(out_stream, "****     "); // last row
	for (int k = rmatbeg[nrows - 1]; k < nnz; k++) fprintf(out_stream, "(%f)x_%d ", rmatval[k], rmatind[k] + 1);
	fprintf(out_stream, " %s %f\n", ((sense[nrows - 1] == 'L') ? "<=" : ((sense[nrows - 1] == 'G') ? ">=" : "=")), rhs[nrows - 1]);
	fprintf(out_stream, "**** Variable bounds:\n");
	for (int j = 0; j < ncols; j++) {
		fprintf(out_stream, "****     %f <= x_%d <= %f\n", lb[j], j + 1, ub[j]);
	}
	if (sol_available) {
		// solution and obj value
		double* x = (double*)malloc(ncols * sizeof(double));
		clone_array(inst->x, x, ncols);
		double objval = inst->objval;
		fprintf(out_stream, "**** Solution: \n");
		for (int j = 0; j < ncols; j++) fprintf(out_stream, "****     x_%d = %f\n", j + 1, x[j]);
		fprintf(out_stream, "**** Obj value: %f\n", objval);
	}
	fprintf(out_stream, "**** Integer or binary variables:\n");
	for (int j = 0; j < ncols; j++) {
		fprintf(out_stream, "****     x_%d: ", j + 1);
		if (inst->int_var[j]) fprintf(out_stream, "BIN/INT\n");
		else fprintf(out_stream, "CONT\n");
	}
	fprintf(out_stream, "********************************************************\n");
	fprintf(out_stream, "********************************************************\n\n");
	if (to_file) fclose(out_stream);
}