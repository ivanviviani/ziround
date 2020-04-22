/**
 * @file read_ziround_input.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void read_problem_sizes(instance* inst) {

	// External variables
	CPXENVptr env; /**< CPLEX environment pointer. */
	CPXLPptr lp;   /**< CPLEX lp pointer. */

	// Initialize
	env = inst->env;
	lp = inst->lp;

	inst->objsen = CPXgetobjsen(env, lp);
	inst->nrows = CPXgetnumrows(env, lp);
	inst->ncols = CPXgetnumcols(env, lp);
}

void populate_inst(instance* inst) {

	read_solution(inst);
	read_variable_bounds(inst);
	read_objective_value(inst);
	read_objective_coefficients(inst);
	read_constraints_coefficients(inst);
	read_constraints_senses(inst);
	read_constraints_right_hand_sides(inst);
	read_row_slacks(inst);
}

void read_solution(instance* inst) {

	// External variables
	CPXENVptr env; /**< CPLEX environment pointer. */
	CPXLPptr lp;   /**< CPLEX lp pointer. */
	int ncols;     /**< Number of variables (columns) of the MIP problem. */
	double* x;     /**< Initial continuous relaxation solution. */
	// Local variables
	int solstat;   /**< Solution status according to CPLEX. */
	int solmethod; /**< Solution method according to CPLEX. */
	int soltype;   /**< Solution type according to CPLEX. */
	int status;    /**< Support status flag. */

	// Initialize
	env = inst->env;
	lp = inst->lp;
	ncols = inst->ncols;
	x = NULL;
	status = 0;

	// Allocate on the instance, assign to locals
	inst->x = (double*)malloc(ncols * sizeof(double));
	x = inst->x;
	if (x == NULL) print_error("[read_solution]: Failed to allocate solution.\n");

	// Get solution status
	solstat = CPXgetstat(env, lp);
	switch (solstat) {
		case CPX_STAT_UNBOUNDED:
			print_warning("[read_solution]: Model is unbounded.\n"); 
			return 1;
		case CPX_STAT_INFEASIBLE:
			print_warning("[read_solution]: Model is infeasible.\n"); 
			return 1;
		case CPX_STAT_INForUNBD:
			print_warning("[read_solution]: Model is infeasible or unbounded.\n"); 
			return 1;
		default:
			break;
	}

	// Get solution info
	status = CPXsolninfo(env, lp, &(solmethod), &(soltype), NULL, NULL);
	if (status) print_error("[read_solution]: Failed to obtain solution info.\n");
	if (soltype == CPX_NO_SOLN) print_error("[read_solution]: Solution not available.\n");
	print_verbose(150, "Solution status %d, solution method %d.\n", solstat, solmethod);

	// Get solution
	status = CPXgetx(env, lp, x, 0, ncols - 1);
	if (status) print_error("[read_solution]: Failed to obtain primal solution.\n");
}

void read_variable_bounds(instance* inst) {

	// External variables
	CPXENVptr env; /**< CPLEX environment pointer. */
	CPXLPptr lp;   /**< CPLEX lp pointer. */
	int ncols;     /**< Number of variables (columns) of the MIP problem. */
	double* ub;	   /**< Variable upper bounds array. */
	double* lb;	   /**< Variable lower bounds array. */
	// Local variables
	int status;    /**< Support status flag. */

	// Initialize
	env = inst->env;
	lp = inst->lp;
	ncols = inst->ncols;
	ub = NULL;
	lb = NULL;
	status = 0;

	// Allocate on the instance, assign to locals
	inst->ub = (double*)malloc(inst->ncols * sizeof(double));
	ub = inst->ub;
	inst->lb = (double*)malloc(inst->ncols * sizeof(double));
	lb = inst->lb;
	if (ub == NULL || lb == NULL) print_error("[read_variable_bounds]: Failed to allocate variable bounds.\n");

	// Get variable bounds (upper and lower)
	status = CPXgetub(env, lp, ub, 0, ncols - 1);
	if (status) print_error("[read_variable_bounds]: Failed to obtain upper bounds.\n");
	status = CPXgetlb(env, lp, lb, 0, ncols - 1);
	if (status) print_error("[read_variable_bounds]: Failed to obtain lower bounds.\n");
}

void read_objective_value(instance* inst) {

	// External variables
	CPXENVptr env;      /**< CPLEX environment pointer. */
	CPXLPptr lp;        /**< CPLEX lp pointer. */
	int ncols;          /**< Number of variables (columns) of the MIP problem. */
	double* objval_ptr; /**< Pointer to the initial objective value of the solution of the continuous relaxation. */
	// Local variables
	int status;         /**< Support status flag. */

	// Initialize
	env = inst->env;
	lp = inst->lp;
	ncols = inst->ncols;
	objval_ptr = &(inst->objval);
	status = 0;

	// Get objective value
	status = CPXgetobjval(env, lp, objval_ptr);
	if (status) print_error("[read_objective_value]: Failed to obtain objective value.\n");
}

void read_objective_coefficients(instance* inst) {

	// External variables
	CPXENVptr env;      /**< CPLEX environment pointer. */
	CPXLPptr lp;        /**< CPLEX lp pointer. */
	int ncols;          /**< Number of variables (columns) of the MIP problem. */
	double* obj;        /**< Objective coefficients. */
	// Local variables
	int status;         /**< Support status flag. */

	// Initialize
	env = inst->env;
	lp = inst->lp;
	ncols = inst->ncols;
	obj = NULL;
	status = 0;

	// Allocate on the instance, assign to locals
	inst->obj = (double*)malloc(ncols * sizeof(double));
	obj = inst->obj;
	if (obj == NULL) print_error("[read_objective_coefficients]: Failed to allocate objective coefficients.\n");

	// Get objective coefficients
	status = CPXgetobj(env, lp, obj, 0, ncols - 1);
	if (status) print_error("[read_objective_coefficients]: Failed to obtain objective coefficients.\n");
}

void read_constraints_coefficients(instance* inst) {

	// External variables
	CPXENVptr env;  /**< CPLEX environment pointer. */
	CPXLPptr lp;    /**< CPLEX lp pointer. */
	int ncols;      /**< Number of variables (columns) of the MIP problem. */
	int nrows;      /**< Number of constraints (rows) of the MIP problem. */
	int nnz; 		/**< Number of non-zero coefficients. */
	int* rowbeg; 	/**< Begin row indices of non-zero coefficients for rmatind and rmatval. */
	int* rowind; 	/**< Column indices of non-zero coefficients. */
	double* rowval; /**< Non-zero coefficients (row major). */
	int* colbeg; 	/**< Begin column indices of non-zero coefficients for cmatind and cmatval. */
	int* colind; 	/**< Row indices of non-zero coefficients. */
	double* colval; /**< Non-zero coefficients (column major). */
	// Local variables
	int unused;     /**< Unused parameter. */
	int status;		/**< Support status flag. */

	// Initialize
	env = inst->env;
	lp = inst->lp;
	ncols = inst->ncols;
	nrows = inst->nrows;
	nnz = 0;
	rowbeg = NULL;
	rowind = NULL;
	rowval = NULL;
	colbeg = NULL;
	colind = NULL;
	colval = NULL;
	unused = 0;
	status = 0;

	// First, get the number of non zero coefficients of the matrix (nzcnt)
	inst->nzcnt = CPXgetnumnz(env, lp);
	nnz = inst->nzcnt;

	// Allocate on the instance, assign to locals
	inst->rmatbeg = (int*)malloc(nrows * sizeof(int));
	rowbeg = inst->rmatbeg;
	inst->rmatind = (int*)malloc(nnz * sizeof(int));
	rowind = inst->rmatind;
	inst->rmatval = (double*)malloc(nnz * sizeof(double));
	rowval = inst->rmatval;
	inst->cmatbeg = (int*)malloc(ncols * sizeof(int));
	colbeg = inst->cmatbeg;
	inst->cmatind = (int*)malloc(nnz * sizeof(int));
	colind = inst->cmatind;
	inst->cmatval = (double*)malloc(nnz * sizeof(double));
	colval = inst->cmatval;
	if (rowbeg == NULL || rowind == NULL || rowval == NULL || 
		colbeg == NULL || colind == NULL || colval == NULL) {
		print_error("[read_constraints_coefficients]: Failed to allocate one of rmatbeg, rmatind, rmatval, cmatbeg, cmatind, cmatval.\n");
	}

	// Get constraint matrix, both by rows and by columns
	status = CPXgetrows(env, lp, &unused, rowbeg, rowind, rowval, nnz, &unused, 0, nrows - 1);
	if (status) print_error("[read_constraints_coefficients]: Failed to obtain rows info.\n");
	status = CPXgetcols(env, lp, &unused, colbeg, colind, colval, nnz, &unused, 0, ncols - 1);
	if (status) print_error("[read_constraints_coefficients]: Failed to obtain columns info.\n");
}

void read_constraints_senses(instance* inst) {

	// External variables
	CPXENVptr env; /**< CPLEX environment pointer. */
	CPXLPptr lp;   /**< CPLEX lp pointer. */
	int nrows;     /**< Number of constraints (rows) of the MIP problem. */
	char* sense;   /**< Constraint senses array. */
	// Local variables
	int status;    /**< Support status flag. */

	// Initialize
	env = inst->env;
	lp = inst->lp;
	nrows = inst->nrows;
	sense = NULL;
	status = 0;

	// Allocate on the instance, assign to locals
	inst->sense = (char*)malloc(nrows * sizeof(char));
	sense = inst->sense;
	if (sense == NULL) print_error("[read_constraints_senses]: Failed to allocate constraint senses.\n");

	// Get constraint senses {'L','E','G'}
	status = CPXgetsense(env, lp, sense, 0, nrows - 1);
	if (status) print_error("[read_constraints_senses]: Failed to obtain constraints senses.\n");
}

void read_constraints_right_hand_sides(instance* inst) {

	// External variables
	CPXENVptr env; /**< CPLEX environment pointer. */
	CPXLPptr lp;   /**< CPLEX lp pointer. */
	int nrows;     /**< Number of constraints (rows) of the MIP problem. */
	double* rhs;   /**< Constraint right hand sides array. */
	// Local variables
	int status;    /**< Support status flag. */

	// Initialize
	env = inst->env;
	lp = inst->lp;
	nrows = inst->nrows;
	rhs = NULL;
	status = 0;

	// Allocate on the instance, assign to locals
	inst->rhs = (double*)malloc(nrows * sizeof(double));
	rhs = inst->rhs;
	if (rhs == NULL) print_error("[read_constraints_right_hand_sides]: Failed to allocate right hand sides.\n");

	// Get right hand sides
	status = CPXgetrhs(env, lp, rhs, 0, nrows - 1);
	if (status) print_error("[read_constraints_right_hand_sides]: Failed to obtain rhs.\n");
}

void read_row_slacks(instance* inst) {

	// External variables
	CPXENVptr env; /**< CPLEX environment pointer. */
	CPXLPptr lp;   /**< CPLEX lp pointer. */
	int nrows;     /**< Number of constraints (rows) of the MIP problem. */
	double* slack; /**< Row slacks array. */
	// Local variables
	int status;    /**< Support status flag. */

	// Initialize
	env = inst->env;
	lp = inst->lp;
	nrows = inst->nrows;
	slack = NULL;
	status = 0;

	// Allocate on the instance, assign to locals
	inst->slack = (double*)malloc(nrows * sizeof(double));
	slack = inst->slack;
	if (slack == NULL) print_error("[read_row_slacks]: Failed to allocate slacks.\n");

	// Get row slacks
	status = CPXgetslack(env, lp, slack, 0, nrows - 1);
	if (status) print_error("[read_row_slacks]: Failed to obtain slacks.\n");
}