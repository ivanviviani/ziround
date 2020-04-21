#include "ziround.h"

//**************************************************************************************************************************************************************
//*************************************** ZI-ROUND *************************************************************************************************************
//**************************************************************************************************************************************************************
int zi_round(instance* inst) {

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
	int status;			 /**< Support status flag. */

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
	status = 0;

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

						print_verbose(100, "[zi_round]: >>> Set x_%d = x_%d + UB_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, x[j], delta_up[j], x[j] + delta_up[j]);
						
						// Round UP
						x[j] += delta_up[j];

						updated = 1;
						update_slacks(inst, j, delta_up[j]);
						update_objval(inst, &objval, j, delta_up[j]);
					}

					// Third case
					else if (ZIminus - ZIplus < -(TOLERANCE) && ZIminus - ZI < -(TOLERANCE)) {
						
						print_verbose(100, "[zi_round]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, x[j], delta_down[j], x[j] - delta_down[j]);
						
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
		if (VERBOSE > 150) {
			check_bounds(inst, x);
			check_constraints(inst, x);
		}

	} while (updated); // end outer loop

	// Free
	free(delta_up);
	free(delta_down);

	return status;
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

				print_verbose(150, "[round_xj_bestobj]: _>_ UBj = %f (must be 1.0)\n", delta_up[j]);
				print_verbose(100, "[round_xj_bestobj]: >>> Set x_%d = x_%d + UB_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, x[j], delta_up[j], x[j] + delta_up[j]);
				if (!xj_fractional && VERBOSE > 10 && fabs(delta_up[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: UB_%d = %f (should be 1.0).\n", j + 1, delta_up[j]);
				}

				// xj = xj + delta_up[j] (if xj is not fractional then delta_up[j] must be 1.0)
				x[j] += delta_up[j];
				updated = 1;
				update_slacks(inst, j, delta_up[j]);
				update_objval(inst, objval, j, delta_up[j]);
			}

			else if (obj_deltaminus < -(TOLERANCE) && obj_deltaminus - obj_deltaplus < -(TOLERANCE)) {

				print_verbose(150, "[round_xj_bestobj]: _>_ LBj = %f (must be 1.0)\n", delta_down[j]);
				print_verbose(100, "[round_xj_bestobj]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, x[j], delta_down[j], x[j] - delta_down[j]);
				if (!xj_fractional && VERBOSE > 10 && fabs(delta_down[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: LB_%d = %f (should be 1.0).\n", j + 1, delta_down[j]);
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
				
				print_verbose(150, "[round_xj_bestobj]: _>_ UBj = %f (must be 1.0)\n", delta_up[j]);
				print_verbose(100, "[round_xj_bestobj]: >>> Set x_%d = x_%d + UB_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, x[j], delta_up[j], x[j] + delta_up[j]);
				if (!xj_fractional && VERBOSE > 10 && fabs(delta_up[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: UB_%d = %f (should be 1.0).\n", j + 1, delta_up[j]);
				}

				// xj = xj + delta_up[j] (if xj is not fractional then delta_up[j] must be 1.0)
				x[j] += delta_up[j];
				updated = 1;
				update_slacks(inst, j, delta_up[j]);
				update_objval(inst, objval, j, delta_up[j]);
			}

			else if (obj_deltaminus > TOLERANCE && obj_deltaminus - obj_deltaplus > TOLERANCE) {
				
				print_verbose(150, "[round_xj_bestobj]: _>_ LBj = %f (must be 1.0)\n", delta_down[j]);
				print_verbose(100, "[round_xj_bestobj]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, x[j], delta_down[j], x[j] - delta_down[j]);
				if (!xj_fractional && VERBOSE > 10 && fabs(delta_down[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: LB_%d = %f (should be 1.0).\n", j + 1, delta_down[j]);
				}

				// xj = xj - delta_down[j] (if xj is not fractional then delta_down[j] must be 1.0)
				x[j] -= delta_down[j];
				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				update_objval(inst, objval, j, -(delta_down[j]));
			}

			break;

		default:
			print_error("[round_xj_bestobj]: Entered default case in switch.\n");
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
	char* sense;    /**< Current row (constraint) sense. */
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
	*objval += (inst->obj[j] * signed_delta);
}

void delta_updown(instance* inst, int j, double* delta_up, double* delta_down, const double epsilon) {
	/*
		For 'L' (<=) constraints: (si non-negative)
			firstUBj_L = min_i{si/aij : aij > 0}
			firstLBj_L = min_i{-si/aij : aij < 0}
		For 'G' (>=) constraints: (si non-positive)
			firstUBj_G = min_i{si/aij : aij < 0}
			firstLBj_G = min_i{-si/aij : aij > 0}

		--> firstUBj = min{ firstUBj_L , firstUBj_G }
		--> firstLBj = min{ firstLBj_L , firstLBj_G }
	*/
	double firstUBj = LONG_MAX, firstLBj = LONG_MAX;
	double aij; int row_index; char sense;
	double ratioUBj, ratioLBj;
	double newUBj, newLBj;
	int end_col = (j < inst->ncols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;

	double secondUBj = inst->ub[j] - inst->x[j];
	double secondLBj = inst->x[j] - inst->lb[j];
	if (VERBOSE > 10) fprintf(stdout, "[DEBUG][delta_updown]: secondUB_%d = ub_%d - x_%d = %f - %f = %f ; secondLB_%d = x_%d - lb_%d = %f - %f = %f\n",
		j + 1, j + 1, j + 1, inst->ub[j], inst->x[j], secondUBj, j + 1, j + 1, j + 1, inst->x[j], inst->lb[j], secondLBj);

	// Scan non-zero coefficients of column j
	/*
		For variable xj:
			cmatbeg[j] is the first index of cmatind and cmatval for column j
			--> Column j in range [ cmatbeg[j];cmatbeg[j+1] ) except for the last one (see end_col)
			cmatval contains coefficient values
			cmatind contains the row indices of coefficient values
	*/
	for (int k = inst->cmatbeg[j]; k < end_col; k++) {
		// Get current non-zero aij, its row index and constraint sense
		aij = inst->cmatval[k];
		row_index = inst->cmatind[k];
		sense = inst->sense[row_index];
		// Check sense, then check sign of aij, and update firstUBj, firstLBj
		switch (sense) {
		case 'L': // (slack non-negative)
			// Clip slack to zero if negative
			if (inst->slack[row_index] < -(TOLERANCE)) inst->slack[row_index] = 0.0;
			if (aij > 0) { // Update firstUBj
				if (VERBOSE > 20) fprintf(stdout, "[DEBUG][delta_updown]: sense = L ; slack[%d] = %f ; a_%d_%d = %f ; ratioUB_%d = %f\n",
					row_index + 1, inst->slack[row_index], row_index + 1, j + 1, aij, j + 1, inst->slack[row_index] / aij);
				ratioUBj = inst->slack[row_index] / aij;
				firstUBj = min(ratioUBj, firstUBj);
			}
			if (aij < 0) { // Update firstLBj
				if (VERBOSE > 20) fprintf(stdout, "[DEBUG][delta_updown]: sense = L ; slack[%d] = %f ; a_%d_%d = %f ; ratioLB_%d = %f\n",
					row_index + 1, inst->slack[row_index], row_index + 1, j + 1, aij, j + 1, -(inst->slack[row_index]) / aij);
				ratioLBj = -(inst->slack[row_index]) / aij;
				firstLBj = min(ratioLBj, firstLBj);
			}
			break;
		case 'G': // (slack non-positive)
			// Clip slack to zero if positive
			if (inst->slack[row_index] > TOLERANCE) inst->slack[row_index] = 0.0;
			if (aij < 0) { // Update firstUBj
				if (VERBOSE > 20) fprintf(stdout, "[DEBUG][delta_updown]: sense = G ; slack[%d] = %f ; a_%d_%d = %f ; ratioUB_%d = %f\n",
					row_index + 1, inst->slack[row_index], row_index + 1, j + 1, aij, j + 1, inst->slack[row_index] / aij);
				ratioUBj = inst->slack[row_index] / aij;
				firstUBj = min(ratioUBj, firstUBj);
			}
			if (aij > 0) { // Update firstLBj
				if (VERBOSE > 20) fprintf(stdout, "[DEBUG][delta_updown]: sense = G ; slack[%d] = %f ; a_%d_%d = %f ; ratioLB_%d = %f\n",
					row_index + 1, inst->slack[row_index], row_index + 1, j + 1, aij, j + 1, -(inst->slack[row_index]) / aij);
				ratioLBj = -(inst->slack[row_index]) / aij;
				firstLBj = min(ratioLBj, firstLBj);
			}
			break;
		case 'E': // (slack zero): variable xj involved in equality constraint, thus cannot be moved
			// Set firstUBj and firstLBj to zero --> newUBj and newLBj will get value zero
			if (VERBOSE > 20) fprintf(stdout, "[DEBUG][delta_updown]: sense = E ; Variable x_%d involved in equality constraint %d --> cannot be moved!\n", j, row_index);
			firstUBj = 0.0;
			firstLBj = 0.0;
			break;
		default:
			fprintf(stderr, "[ERR][delta_updown]: Constraint sense not included in {'L','G','E'}.\n"); exit(EXIT_FAILURE);
			break;
		} // end switch
	} // end for

	// Results
	if (VERBOSE > 10) fprintf(stdout, "[DEBUG][delta_updown][Results]: firstUB_%d = %f ; firstLB_%d = %f\n", j + 1, firstUBj, j + 1, firstLBj);
	newUBj = min(firstUBj, secondUBj);
	newLBj = min(firstLBj, secondLBj);
	if (VERBOSE > 10) fprintf(stdout,
		"[DEBUG][delta_updown][Results]: (New) UB_%d = min{firstUB_%d, secondUB_%d} = min{%f, %f} = %f ; LB_%d = min{firstLB_%d, secondLB_%d} = min{%f, %f} = %f\n",
		j + 1, j + 1, j + 1, firstUBj, secondUBj, newUBj, j + 1, j + 1, j + 1, firstLBj, secondLBj, newLBj);

	// If newUBj < epsilon clip it to 0
	if (newUBj < epsilon) newUBj = 0.0;
	// If newLBj < epsilon clip it to 0
	if (newLBj < epsilon) newLBj = 0.0;
	// Update both no matter what
	delta_up[j] = newUBj;
	delta_down[j] = newLBj;
}