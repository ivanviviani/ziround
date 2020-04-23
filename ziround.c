/**
 * @file ziround.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

//**************************************************************************************************************************************************************
//*************************************** ZI-ROUND *************************************************************************************************************
//**************************************************************************************************************************************************************
void zi_round(instance* inst) {

	// External variables
	int ncols;			 /**< Number of variables (columns) of the MIP problem. */
	double objval;		 /**< Objective value of the initial continuous relaxation solution. */
	double* obj;		 /**< Objective coefficients. */
	double* x;           /**< Initial continuous relaxation solution. */
	int* int_var;		 /**< Array of flags for integer/binary variables of the original MIP problem. */
	// Local variables
	double* delta_up; 	 /**< Maximum variable up-shifts. */
	double* delta_down;  /**< Maximum variable down-shifts. */
	double ZI; 			 /**< Fractionality of a variable (used in function zi_round). */
	double ZIplus; 		 /**< Fractionality of a shifted up variable (used in function zi_round). */
	double ZIminus; 	 /**< Fractionality of a shifted down variable (used in fucntion zi_round). */
	double obj_plusUBj;  /**< Objective value for a shifted up variable (used in function zi_round). */
	double obj_minusLBj; /**< Objective value for a shifted down variable (used in function zi_round). */
	int updated;		 /**< Flag set to 1 when at least one variable shift has been made. */

	// Allocate / Initialize
	ncols = inst->ncols;
	objval = inst->objval;
	obj = inst->obj;
	x = inst->x;
	int_var = inst->int_var;
	delta_up = (double*)malloc(ncols * sizeof(double));
	delta_down = (double*)malloc(ncols * sizeof(double));
	if (delta_up == NULL || delta_down == NULL) print_error("[zi_round]: Failed to allocate delta arrays.\n");
	updated = 0;

	// Outer loop (repeat until no more updates found)
	do {
		updated = 0;

		// Inner loop (for each variable xj that was integer/binary in the original MIP)
		for (int j = 0; j < ncols; j++) {
			if (!(int_var[j])) continue;

			switch (is_fractional(x[j])) {

				// xj non fractional
				case 0:

					print_verbose(150, "[zi_round]: -x- x_%d is non-fractional\n", j + 1);

					// Calculate deltas (with epsilon = 1.0)
					delta_updown(inst, j, delta_up, delta_down, 1.0);

					// Skip xj if both deltas are equal to zero (no shift necessary)
					if (fabs(delta_up[j]) < TOLERANCE && fabs(delta_down[j]) < TOLERANCE) continue;

					// Condition(s) for rounding of xj (>= to include the case of a zero obj coefficient)
					if ((obj[j] >= 0 && fabs(delta_down[j] - 1.0) < TOLERANCE) ||
						(obj[j] <= 0 && fabs(delta_up[j] - 1.0) < TOLERANCE)) {

						// Round xj to improve objective and update slacks
						updated = round_xj_bestobj(inst, &objval, j, delta_up, delta_down, 0); // flag xj non-fractional (0)
					}

					break;

				// xj fractional
				case 1:

					print_verbose(150, "[zi_round]: -x- x_%d is fractional\n", j + 1);

					// Calculate deltas
					delta_updown(inst, j, delta_up, delta_down, EPSILON);

					// Skip xj if both UBj and LBj are equal to zero (no shift necessary)
					if (fabs(delta_up[j]) < TOLERANCE && fabs(delta_down[j]) < TOLERANCE) continue;

					ZI = fractionality(x[j]);
					ZIplus = fractionality(x[j] + delta_up[j]);
					ZIminus = fractionality(x[j] - delta_down[j]);

					// First case
					if (fabs(ZIplus - ZIminus) < TOLERANCE && ZIplus - ZI < -(TOLERANCE)) {

						// Round xj to improve objective and update slacks
						updated = round_xj_bestobj(inst, &objval, j, delta_up, delta_down, 1); // flag xj fractional (1)
					}

					// Second case
					else if (ZIplus - ZIminus < -(TOLERANCE) && ZIplus - ZI < -(TOLERANCE)) {

						print_verbose(100, "[zi_round]: >>> Set x_%d = x_%d + delta_up_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, x[j], delta_up[j], x[j] + delta_up[j]);
						
						// Round UP
						x[j] += delta_up[j];

						updated = 1;
						update_slacks(inst, j, delta_up[j]);
						update_objval(inst, &objval, j, delta_up[j]);
					}

					// Third case
					else if (ZIminus - ZIplus < -(TOLERANCE) && ZIminus - ZI < -(TOLERANCE)) {
						
						print_verbose(100, "[zi_round]: >>> Set x_%d = x_%d - delta_down_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, x[j], delta_down[j], x[j] - delta_down[j]);
						
						// Round DOWN
						x[j] -= delta_down[j];

						updated = 1;
						update_slacks(inst, j, -(delta_down[j]));
						update_objval(inst, &objval, j, -(delta_down[j]));
					}

					break;

				default:
					print_error(" in function is_fractional.\n");
			}
		} // end inner loop

		
		if (updated) print_verbose(100, "[zi_round] ...Update found, scan variables again...\n");
		else print_verbose(100, "[zi_round] ...No updates found, exit outer loop...\n");
		
		// [BRUTE FORCE] [DEBUG ONLY] Check variable bounds and constraints
		if (VERBOSE >= 150) {
			check_bounds(inst, x);
			check_constraints(inst, x);
		}

	} while (updated); // end outer loop

	// Free
	free(delta_up);
	free(delta_down);
}
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************

int round_xj_bestobj(instance* inst, double* objval, int j, double* delta_up, double* delta_down, int xj_fractional) {

	// External variables
	int objsen;			   /**< Objective sense. */
	double* obj;           /**< Objective coefficients. */
	double* x;			   /**< Current solution (rounding in progress). */
	// Local variables
	double obj_deltaplus;  /**< Delta obj if xj is shifted up. */
	double obj_deltaminus; /**< Delta obj if xj is shifted down. */
	int updated;           /**< Flag set to 1 when at least one variable shift has been made. */

	// Initialize
	objsen = inst->objsen;
	obj = inst->obj;
	x = inst->x;
	obj_deltaplus = 0.0 + obj[j] * delta_up[j];
	obj_deltaminus = 0.0 - obj[j] * delta_down[j];
	updated = 0;

	// Check obj sense, then update xj, update slacks and update objective value
	switch (objsen) {

		case CPX_MIN:

			if (obj_deltaplus < -(TOLERANCE) && obj_deltaplus - obj_deltaminus < -(TOLERANCE)) {

				print_verbose(100, "[round_xj_bestobj]: >>> Set x_%d = x_%d + delta_up_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, x[j], delta_up[j], x[j] + delta_up[j]);
				if (!xj_fractional && VERBOSE >= 150 && fabs(delta_up[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: delta_up_%d = %f (should be 1.0).\n", j + 1, delta_up[j]);
				}

				// xj = xj + delta_up[j] (if xj is not fractional then delta_up[j] must be 1.0)
				x[j] += delta_up[j];
				updated = 1;
				update_slacks(inst, j, delta_up[j]);
				update_objval(inst, objval, j, delta_up[j]);
			}

			else if (obj_deltaminus < -(TOLERANCE) && obj_deltaminus - obj_deltaplus < -(TOLERANCE)) {

				print_verbose(100, "[round_xj_bestobj]: >>> Set x_%d = x_%d - delta_down_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, x[j], delta_down[j], x[j] - delta_down[j]);
				if (!xj_fractional && VERBOSE >= 150 && fabs(delta_down[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: delta_down_%d = %f (should be 1.0).\n", j + 1, delta_down[j]);
				}

				// xj = xj - delta_down[j] (if xj is not fractional then delta_down[j] must be 1.0)
				x[j] -= delta_down[j];
				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				update_objval(inst, objval, j, -(delta_down[j]));
			}

			break;

		// (problems from mps files should always be MIN)
		case CPX_MAX:

			if (obj_deltaplus > TOLERANCE && obj_deltaplus - obj_deltaminus > TOLERANCE) {
				
				print_verbose(100, "[round_xj_bestobj]: >>> Set x_%d = x_%d + delta_up_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, x[j], delta_up[j], x[j] + delta_up[j]);
				if (!xj_fractional && VERBOSE >= 150 && fabs(delta_up[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: delta_up_%d = %f (should be 1.0).\n", j + 1, delta_up[j]);
				}

				// xj = xj + delta_up[j] (if xj is not fractional then delta_up[j] must be 1.0)
				x[j] += delta_up[j];
				updated = 1;
				update_slacks(inst, j, delta_up[j]);
				update_objval(inst, objval, j, delta_up[j]);
			}

			else if (obj_deltaminus > TOLERANCE && obj_deltaminus - obj_deltaplus > TOLERANCE) {
				
				print_verbose(100, "[round_xj_bestobj]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, x[j], delta_down[j], x[j] - delta_down[j]);
				if (!xj_fractional && VERBOSE >= 150 && fabs(delta_down[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: delta_down_%d = %f (should be 1.0).\n", j + 1, delta_down[j]);
				}

				// xj = xj - delta_down[j] (if xj is not fractional then delta_down[j] must be 1.0)
				x[j] -= delta_down[j];
				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				update_objval(inst, objval, j, -(delta_down[j]));
			}

			break;

		default:
			print_error("[round_xj_bestobj]: Objective sense not supported.\n");
	}

	return updated;
}

// Incremental update of row slacks (only for the rows where the variable xj is involved)
// Whenever it is called, only one variable xj has been updated
void update_slacks(instance* inst, int j, double signed_delta) {

	// External variables
	int ncols;      /**< Number of variables (columns) of the MIP problem. */
	int nnz;        /**< Number of non-zero coefficients of the problem matrix. */
	int* colbeg;    /**< CPLEX cmatbeg. */
	int* colind;    /**< CPLEX cmatind. */
	double* colval; /**< CPLEX cmatval. */
	double* slack;  /**< Row slacks array. */
	// Local Variables
	double aij;     /**< Current coefficient. */
	int rowind;     /**< Current row index. */
	char* sense;    /**< Constraint senses array. */
	int colend;     /**< Index of the last non-zero coefficient of the current column. */

	// Initialize
	nnz = inst->nzcnt;
	ncols = inst->ncols;
	colbeg = inst->cmatbeg;
	colind = inst->cmatind;
	colval = inst->cmatval;
	slack = inst->slack;
	sense = inst->sense;
	colend = (j < ncols - 1) ? colbeg[j + 1] : nnz;

	// Scan non-zero coefficients of column j
	for (int k = colbeg[j]; k < colend; k++) {

		aij = colval[k];
		rowind = colind[k];

		switch (sense[rowind]) {

			case 'L':
			case 'G':
				slack[rowind] -= aij * signed_delta;
				break;
			case 'E':
				print_error("[update_slacks]: Tried to update slack for a variable involved in an equality constraint!\n");
			default:
				print_error("[update_slacks]: Constraint sense %c not supported!\n", sense[rowind]);
		}
	}
}

// Incremental update of the objective value
// Whenever it is called, only one variable xj has been updated
void update_objval(instance* inst, double* objval, int j, double signed_delta) {

	// External variables
	double* obj; /**< Objective coefficients. */

	// Initialize
	obj = inst->obj;

	*objval += (obj[j] * signed_delta);
}

/*
	For 'L' (<=) constraints: (si non-negative)
		delta_up1_L = min_i{si/aij : aij > 0}
		delta_down1_L = min_i{-si/aij : aij < 0}
	For 'G' (>=) constraints: (si non-positive)
		delta_up1_G = min_i{si/aij : aij < 0}
		delta_down1_G = min_i{-si/aij : aij > 0}

	--> delta_up1 = min{ delta_up1_L , delta_up1_G }
	--> delta_down1 = min{ delta_down1_L , delta_down1_G }
*/
void delta_updown(instance* inst, int j, double* delta_up, double* delta_down, const double epsilon) {

	// External variables
	int ncols;				/**< Number of variables (columns) of the MIP problem. */
	double* ub;				/**< Variable upper bounds array. */
	double* lb;				/**< Variable lower bounds array. */
	double* x;				/**< Current solution (rounding in progress). */
	int nnz;				/**< Number of non-zero coefficients of the problem matrix. */
	int* colbeg;			/**< CPLEX cmatbeg. */
	int* colind;			/**< CPLEX cmatind. */
	double* colval;			/**< CPLEX cmatval. */
	char* sense;			/**< Constraint senses array. */
	double* slack;			/**< Row slacks array. */
	// Local variables
	double delta_up1;		/**< First delta_up[j] major candidate. */
	double delta_down1;		/**< First delta_down[j] major candidate. */
	double delta_up2;		/**< Second delta_up[j] major candidate. */
	double delta_down2;		/**< Second delta_down[j] major candidate. */
	double candidate_up1;	/**< Current delta_up[j] minor candidate. */
	double candidate_down1; /**< Current delta_down[j] minor candidate. */
	double new_delta_up;    /**< Final delta_up[j] winner. */
	double new_delta_down;  /**< Final delta_down[j] winner. */
	double aij;				/**< Current coefficient. */
	int rowind;				/**< Current row index. */
	int colend;				/**< Index of the last non-zero coefficient of the current column. */
	

	// Initialize
	ncols = inst->ncols;
	ub = inst->ub;
	lb = inst->lb;
	x = inst->x;
	nnz = inst->nzcnt;
	colbeg = inst->cmatbeg;
	colind = inst->cmatind;
	colval = inst->cmatval;
	sense = inst->sense;
	slack = inst->slack;
	delta_up1 = LONG_MAX;
	delta_down1 = LONG_MAX;
	delta_up2 = ub[j] - x[j];
	delta_down2 = x[j] - lb[j];
	colend = (j < ncols - 1) ? colbeg[j + 1] : nnz;
	
	print_verbose(200, "[delta_updown]: delta_up2_%d = ub_%d - x_%d = %f - %f = %f ; delta_down2_%d = x_%d - lb_%d = %f - %f = %f\n",
		j + 1, j + 1, j + 1, ub[j], x[j], delta_up2, j + 1, j + 1, j + 1, x[j], lb[j], delta_down2);

	// Scan non-zero coefficients of column j
	for (int k = colbeg[j]; k < colend; k++) {

		aij = colval[k];
		rowind = colind[k];
		
		// Check sense, then check sign of aij, and update delta_up1, delta_down1
		switch (sense[rowind]) {

			case 'L': // (slack non-negative)

				// Clip slack to zero if negative
				if (slack[rowind] < -(TOLERANCE)) slack[rowind] = 0.0;

				if (aij > 0) { 
					
					print_verbose(200, "[delta_updown]: sense = L ; slack[%d] = %f ; a_%d_%d = %f ; ratioUB_%d = %f\n",
						rowind + 1, slack[rowind], rowind + 1, j + 1, aij, j + 1, slack[rowind] / aij);
					
					// Update delta_up1
					candidate_up1 = slack[rowind] / aij;
					delta_up1 = min(candidate_up1, delta_up1);
				}
				if (aij < 0) {

					print_verbose(200, "[delta_updown]: sense = L ; slack[%d] = %f ; a_%d_%d = %f ; ratioLB_%d = %f\n",
						rowind + 1, slack[rowind], rowind + 1, j + 1, aij, j + 1, -(slack[rowind]) / aij);

					// Update delta_down1
					candidate_down1 = -(slack[rowind]) / aij;
					delta_down1 = min(candidate_down1, delta_down1);
				}

				break;

			case 'G': // (slack non-positive)

				// Clip slack to zero if positive
				if (slack[rowind] > TOLERANCE) slack[rowind] = 0.0;

				if (aij < 0) {

					print_verbose(200, "[delta_updown]: sense = G ; slack[%d] = %f ; a_%d_%d = %f ; ratioUB_%d = %f\n",
						rowind + 1, slack[rowind], rowind + 1, j + 1, aij, j + 1, slack[rowind] / aij);

					// Update delta_up1
					candidate_up1 = slack[rowind] / aij;
					delta_up1 = min(candidate_up1, delta_up1);
				}
				if (aij > 0) {

					print_verbose(200, "[delta_updown]: sense = G ; slack[%d] = %f ; a_%d_%d = %f ; ratioLB_%d = %f\n",
						rowind + 1, slack[rowind], rowind + 1, j + 1, aij, j + 1, -(slack[rowind]) / aij);

					// Update delta_down1
					candidate_down1 = -(slack[rowind]) / aij;
					delta_down1 = min(candidate_down1, delta_down1);
				}

				break;

			case 'E': // (slack zero): variable xj involved in equality constraint, thus cannot be moved

				print_verbose(200, "[delta_updown]: sense = E ; Variable x_%d involved in equality constraint %d --> cannot be moved!\n", j + 1, rowind);
				
				// Set delta_up1 and delta_down1 to zero --> new_delta_up and new_delta_down will get value zero
				delta_up1 = 0.0;
				delta_down1 = 0.0;

				break;

			default:
				print_error("[delta_updown]: Constraint sense not included in {'L','G','E'}.\n");
		}
	} // end for

	print_verbose(200, "[delta_updown][candidates]: delta_up1_%d = %f ; delta_down1_%d = %f\n", j + 1, delta_up1, j + 1, delta_down1);

	// Results
	new_delta_up = min(delta_up1, delta_up2);
	new_delta_down = min(delta_down1, delta_down2);
	print_verbose(200, "[delta_updown][results]: (NEW) delta_up_%d = min{%f, %f} = %f ; delta_down_%d = min{%f, %f} = %f\n",
		j + 1, delta_up1, delta_up2, new_delta_up, j + 1, delta_down1, delta_down2, new_delta_down);

	// Update deltas
	if (new_delta_up < epsilon) new_delta_up = 0.0;
	if (new_delta_down < epsilon) new_delta_down = 0.0;
	delta_up[j] = new_delta_up;
	delta_down[j] = new_delta_down;
}