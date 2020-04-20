#include "ziround.h"

//**************************************************************************************************************************************************************
//*************************************** ZI-ROUND *************************************************************************************************************
//**************************************************************************************************************************************************************
int zi_round(instance* inst) {
	int status = 0;

	// Local variables
	double* UB; 		 /**< Maximum variable up-shifts. */
	double* LB; 		 /**< Maximum variable down-shifts. */
	double ZI; 			 /**< Fractionality of a variable (used in function zi_round). */
	double ZIplus; 		 /**< Fractionality of a shifted up variable (used in function zi_round). */
	double ZIminus; 	 /**< Fractionality of a shifted down variable (used in fucntion zi_round). */
	double objval;		 /**< Objective value of the initial continuous relaxation solution. */
	double obj_plusUBj;  /**< Objective value for a shifted up variable (used in function zi_round). */
	double obj_minusLBj; /**< Objective value for a shifted down variable (used in function zi_round). */
	int updated;		 /**< Flag set to 1 when at least one variable shift has been made. */
	
	// Allocate / Initialize
	UB = (double*)malloc(inst->ncols * sizeof(double));
	LB = (double*)malloc(inst->ncols * sizeof(double));
	objval = inst->objval;
	updated = 0;

	// inst->x_prev = (double*)malloc(inst->ncols * sizeof(double));
	// inst->x_updated = (double*)malloc(inst->ncols * sizeof(double));
	
	if (/*inst->x_prev == NULL || inst->x_updated == NULL ||*/ UB == NULL || LB == NULL) {
		fprintf(stderr, "[ERR][zi_round]: Failed to allocate x_prev, x_updated or row_infeas.\n"); return 1;
	}

	// Outer loop (repeat until no more updates found)
	do {

		// Initialize update flag and save current solution
		updated = 0;
		// clone_array(inst->x, inst->x_prev, inst->ncols);

		// Inner loop (for each variable xj that was integer/binary in the original MIP)
		for (int j = 0; j < inst->ncols; j++) {
			// Skip xj if it's continuous in the original MIP
			if (!(inst->int_var)) continue;

			// xj non-fractional
			if (!is_fractional(inst->x[j])) {
				if (VERBOSE > 10) fprintf(stdout, "[DEBUG][zi_round]: -x- x_%d is non-fractional\n", j + 1);

				// Calculate UBj and LBj (with epsilon = 1.0)
				calculate_UBjLBj(inst, j, UB, LB, 1.0);
				// Skip xj if both UBj and LBj are equal to zero (--> no shift necessary)
				if (fabs(UB[j]) < TOLERANCE && fabs(LB[j]) < TOLERANCE) continue;

				// Condition(s) for rounding of xj
				// (>= to include the case of a zero obj coefficient)
				if ((inst->obj[j] >= 0 && fabs(LB[j] - 1.0) < TOLERANCE) ||
					(inst->obj[j] <= 0 && fabs(UB[j] - 1.0) < TOLERANCE)) {
					// Update xj to improve objective and update slacks
					update_xj_to_improve_objective(inst, &objval, j, UB, LB, &updated, 0); // flag xj non-fractional (0)
				}
			} // end if xj non-fractional

			// xj fractional
			else if (is_fractional(inst->x[j])) {
				if (VERBOSE > 10) fprintf(stdout, "[DEBUG][zi_round]: -x- x_%d is fractional\n", j + 1);

				// Calculate UBj and LBj
				calculate_UBjLBj(inst, j, UB, LB, EPSILON);
				// Skip xj if both UBj and LBj are equal to zero (--> no shift necessary)
				if (fabs(UB[j]) < TOLERANCE && fabs(LB[j]) < TOLERANCE) continue;

				// ZI-Round (version 1) core (3 cases). First, calculate the fractionalities needed.
				ZI = fractionality(inst->x[j]);
				ZIplus = fractionality(inst->x[j] + UB[j]);
				ZIminus = fractionality(inst->x[j] - LB[j]);

				// First case
				if (fabs(ZIplus - ZIminus) < TOLERANCE &&
					(ZIplus - ZI) < -(TOLERANCE)) { // was inst->ZIplus < inst->ZI
					// Update xj to improve objective and update slacks
					update_xj_to_improve_objective(inst, &objval, j, UB, LB, &updated, 1); // flag xj fractional (1)
				} // end first case

				// Second case
				else if ((ZIplus - ZIminus) < -(TOLERANCE) && // was inst->ZIplus < inst->ZIminus
					(ZIplus - ZI) < -(TOLERANCE)) {      // was inst->ZIplus < inst->ZI
					if (VERBOSE > 10) fprintf(stdout, "[DEBUG][zi_round]: >>> Set x_%d = x_%d + UB_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], UB[j], inst->x[j] + UB[j]);
					inst->x[j] += UB[j];
					updated = 1;
					update_slacks(inst, j, UB[j]);
					update_objective_value(inst, &objval, j, UB[j]);
				} // end second case

				// Third case
				else if ((ZIminus - ZIplus) < -(TOLERANCE) && // was inst->ZIminus < inst->ZIplus
					(ZIminus - ZI) < -(TOLERANCE)) {     // was inst->ZIminus < inst->ZI
					if (VERBOSE > 10) fprintf(stdout, "[DEBUG][zi_round]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], LB[j], inst->x[j] - LB[j]);
					inst->x[j] -= LB[j];
					updated = 1;
					update_slacks(inst, j, -(LB[j]));
					update_objective_value(inst, &objval, j, -(LB[j]));
				} // end third case

			} // end if xj fractional

		} // end inner loop

		if (VERBOSE > 10) {
			if (updated) fprintf(stdout, "[DEBUG][zi_round] ...Update found, scan variables again...\n");
			else fprintf(stdout, "[DEBUG][zi_round] ...No updates found, exit outer loop...\n");
		}
		// [BRUTE FORCE] [DEBUG ONLY] Check variable bounds and constraints
		if (VERBOSE > 50) {
			status = check_bounds(inst, inst->x);      if (status) { fprintf(stderr, "[ERR][zi_round]: Error inside check_bounds.\n"); return status; }
			status = check_constraints(inst, inst->x); if (status) { fprintf(stderr, "[ERR][zi_round]: Error inside check_constraints.\n"); return status; }
		}

	} while (updated); // end outer loop

	// Free
	free(UB);
	free(LB);

	return inst->status;
}
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************

void update_xj_to_improve_objective(instance* inst, double* objval, int j, double* delta_up, double* delta_down, int* updated, int is_fractional) {

	// First, calculate obj value for both updates. (IMPROVEMENT: inst->objval is just a common offset!)
	double obj_plusUBj = *objval + (inst->obj[j] * delta_up[j]);
	double obj_minusLBj = *objval - (inst->obj[j] * delta_down[j]);
	/* was
	// compare delta * coef (confronto solo tra i due delta * coef_j (+ objval è solo un offset) (tenendo conto della solita float tolerance)
	clone_array(inst->x, inst->x_updated, inst->ncols);
	inst->x_updated[j] = inst->x[j] + inst->UB[j];
	..inst->obj_plusUBj = dot_product(inst->obj, inst->x_updated, inst->ncols); // brute force!
	inst->x_updated[j] = inst->x[j] - inst->LB[j];
	..inst->obj_minusLBj = dot_product(inst->obj, inst->x_updated, inst->ncols); // brute force!*/

	// Check obj sense, then update xj, update slacks and update objective value
	switch (inst->objsen) {
	case CPX_MIN:
		if ((obj_plusUBj - *objval) < -(TOLERANCE) &&       // was inst->obj_plusUBj < inst->objval 
			(obj_plusUBj - obj_minusLBj) < -(TOLERANCE)) { // was inst->obj_plusUBj < inst->obj_minusLBj
			// xj = xj + UBj (if xj is not fractional then UBj must be 1.0)
			if (VERBOSE > 50) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: _>_ UBj = %f (must be 1.0)\n", delta_up[j]);
			if (VERBOSE > 10) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: >>> Set x_%d = x_%d + UB_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
			if (!is_fractional && VERBOSE > 10 && fabs(delta_up[j] - 1.0) > TOLERANCE) {
				fprintf(stderr, "[ERR][update_xj_to_improve_objective]: UB_%d = %f (should be 1.0).\n", j + 1, delta_up[j]);
				exit(EXIT_FAILURE);
			}
			inst->x[j] += delta_up[j];
			*updated = 1;
			update_slacks(inst, j, delta_up[j]);
			update_objective_value(inst, objval, j, delta_up[j]);
		}
		else if ((obj_minusLBj - *objval) < -(TOLERANCE) &&      // was inst->obj_minusLBj < inst->objval
			(obj_minusLBj - obj_plusUBj) < -(TOLERANCE)) { // was inst->obj_minusLBj < inst->obj_plusUBj
	   // xj = xj - LBj (if xj is not fractional then LBj must be 1.0)
			if (VERBOSE > 50) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: _>_ LBj = %f (must be 1.0)\n", delta_down[j]);
			if (VERBOSE > 10) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);
			if (!is_fractional && VERBOSE > 10 && fabs(delta_down[j] - 1.0) > TOLERANCE) {
				fprintf(stderr, "[ERR][update_xj_to_improve_objective]: LB_%d = %f (should be 1.0).\n", j + 1, delta_down[j]);
				exit(EXIT_FAILURE);
			}
			inst->x[j] -= delta_down[j];
			*updated = 1;
			update_slacks(inst, j, -(delta_down[j]));
			update_objective_value(inst, objval, j, -(delta_down[j]));
		}
		break;
		// (problems from mps files should always be MIN)
	case CPX_MAX:
		if ((obj_plusUBj - *objval) > TOLERANCE &&       // was inst->obj_plusUBj > inst->objval
			(obj_plusUBj - obj_minusLBj) > TOLERANCE) { // was inst->obj_plusUBj > inst->obj_minusLBj
			// xj = xj + UBj (if xj is not fractional then UBj must be 1.0)
			if (VERBOSE > 50) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: _>_ UBj = %f (must be 1.0)\n", delta_up[j]);
			if (VERBOSE > 10) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: >>> Set x_%d = x_%d + UB_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
			if (!is_fractional && VERBOSE > 10 && fabs(delta_up[j] - 1.0) > TOLERANCE) {
				fprintf(stderr, "[ERR][update_xj_to_improve_objective]: UB_%d = %f (should be 1.0).\n", j + 1, delta_up[j]);
				exit(EXIT_FAILURE);
			}
			inst->x[j] += delta_up[j];
			*updated = 1;
			update_slacks(inst, j, delta_up[j]);
			update_objective_value(inst, objval, j, delta_up[j]);
		}
		else if ((obj_minusLBj - *objval) > TOLERANCE &&      // was inst->obj_minusLBj > inst->objval
			(obj_minusLBj - obj_plusUBj) > TOLERANCE) { // was inst->obj_minusLBj > inst->obj_plusUBj
	   // xj = xj - LBj (if xj is not fractional then LBj must be 1.0)
			if (VERBOSE > 50) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: _>_ LBj = %f (must be 1.0)\n", delta_down[j]);
			if (VERBOSE > 10) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);
			if (!is_fractional && VERBOSE > 10 && fabs(delta_down[j] - 1.0) > TOLERANCE) {
				fprintf(stderr, "[ERR][update_xj_to_improve_objective]: LB_%d = %f (should be 1.0).\n", j + 1, delta_down[j]);
				exit(EXIT_FAILURE);
			}
			inst->x[j] -= delta_down[j];
			*updated = 1;
			update_slacks(inst, j, -(delta_down[j]));
			update_objective_value(inst, objval, j, -(delta_down[j]));
		}
		break;
	default:
		fprintf(stderr, "[ERR][update_xj_to_improve_objective]: Entered default case in switch.\n");
		exit(EXIT_FAILURE);
		break;
	} // end switch
}

// Incremental update of row slacks (only for the rows where the variable xj is involved)
// Whenever it is called, only one variable xj has been updated
void update_slacks(instance* inst, int j, double signed_delta) {
	double aij; int row_index; char sense;
	int end_col = (j < inst->ncols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;
	// Scan non-zero coefficients of column j
	/*
		For variable xj:
			cmatbeg[j] is the first index of cmatind and cmatval for column j
			--> Column j in range [ cmatbeg[j];cmatbeg[j+1] ) except for the last one (see end_col)
			cmatval contains coefficient values
			cmatind contains the row indices of coefficient values
	*/
	for (int k = inst->cmatbeg[j]; k < end_col; k++) {
		aij = inst->cmatval[k];
		row_index = inst->cmatind[k];
		sense = inst->sense[row_index];
		if (sense == 'E') {
			fprintf(stderr, "[ERR][update_slacks]: Tried to update slack for a variable involved in an equality constraint!\n");
			exit(EXIT_FAILURE);
		}
		/*
			slack = rhs - row_activity
			row_activity --> row_activity + (aij * signed_delta)
			--> slack = rhs - (row_activity + (aij * signed_delta))
					  = rhs - row_activity - (aij * signed_delta)
					  = (rhs - row_activity) - (aij * signed_delta)
			slack --> slack - (aij * signed_delta)
		*/
		inst->slack[row_index] -= (aij * signed_delta);
	}
}

// Incremental update of the objective value
// Whenever it is called, only one variable xj has been updated
void update_objective_value(instance* inst, double* objval, int j, double signed_delta) {
	*objval += (inst->obj[j] * signed_delta);
}

void calculate_UBjLBj(instance* inst, int j, double* delta_up, double* delta_down, const double epsilon) {
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
	if (VERBOSE > 10) fprintf(stdout, "[DEBUG][calculate_UBjLBj]: secondUB_%d = ub_%d - x_%d = %f - %f = %f ; secondLB_%d = x_%d - lb_%d = %f - %f = %f\n",
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
				if (VERBOSE > 20) fprintf(stdout, "[DEBUG][calculate_UBjLBj]: sense = L ; slack[%d] = %f ; a_%d_%d = %f ; ratioUB_%d = %f\n",
					row_index + 1, inst->slack[row_index], row_index + 1, j + 1, aij, j + 1, inst->slack[row_index] / aij);
				ratioUBj = inst->slack[row_index] / aij;
				firstUBj = min(ratioUBj, firstUBj);
			}
			if (aij < 0) { // Update firstLBj
				if (VERBOSE > 20) fprintf(stdout, "[DEBUG][calculate_UBjLBj]: sense = L ; slack[%d] = %f ; a_%d_%d = %f ; ratioLB_%d = %f\n",
					row_index + 1, inst->slack[row_index], row_index + 1, j + 1, aij, j + 1, -(inst->slack[row_index]) / aij);
				ratioLBj = -(inst->slack[row_index]) / aij;
				firstLBj = min(ratioLBj, firstLBj);
			}
			break;
		case 'G': // (slack non-positive)
			// Clip slack to zero if positive
			if (inst->slack[row_index] > TOLERANCE) inst->slack[row_index] = 0.0;
			if (aij < 0) { // Update firstUBj
				if (VERBOSE > 20) fprintf(stdout, "[DEBUG][calculate_UBjLBj]: sense = G ; slack[%d] = %f ; a_%d_%d = %f ; ratioUB_%d = %f\n",
					row_index + 1, inst->slack[row_index], row_index + 1, j + 1, aij, j + 1, inst->slack[row_index] / aij);
				ratioUBj = inst->slack[row_index] / aij;
				firstUBj = min(ratioUBj, firstUBj);
			}
			if (aij > 0) { // Update firstLBj
				if (VERBOSE > 20) fprintf(stdout, "[DEBUG][calculate_UBjLBj]: sense = G ; slack[%d] = %f ; a_%d_%d = %f ; ratioLB_%d = %f\n",
					row_index + 1, inst->slack[row_index], row_index + 1, j + 1, aij, j + 1, -(inst->slack[row_index]) / aij);
				ratioLBj = -(inst->slack[row_index]) / aij;
				firstLBj = min(ratioLBj, firstLBj);
			}
			break;
		case 'E': // (slack zero): variable xj involved in equality constraint, thus cannot be moved
			// Set firstUBj and firstLBj to zero --> newUBj and newLBj will get value zero
			if (VERBOSE > 20) fprintf(stdout, "[DEBUG][calculate_UBjLBj]: sense = E ; Variable x_%d involved in equality constraint %d --> cannot be moved!\n", j, row_index);
			firstUBj = 0.0;
			firstLBj = 0.0;
			break;
		default:
			fprintf(stderr, "[ERR][calculate_UBjLBj]: Constraint sense not included in {'L','G','E'}.\n"); exit(EXIT_FAILURE);
			break;
		} // end switch
	} // end for

	// Results
	if (VERBOSE > 10) fprintf(stdout, "[DEBUG][calculate_UBjLBj][Results]: firstUB_%d = %f ; firstLB_%d = %f\n", j + 1, firstUBj, j + 1, firstLBj);
	newUBj = min(firstUBj, secondUBj);
	newLBj = min(firstLBj, secondLBj);
	if (VERBOSE > 10) fprintf(stdout,
		"[DEBUG][calculate_UBjLBj][Results]: (New) UB_%d = min{firstUB_%d, secondUB_%d} = min{%f, %f} = %f ; LB_%d = min{firstLB_%d, secondLB_%d} = min{%f, %f} = %f\n",
		j + 1, j + 1, j + 1, firstUBj, secondUBj, newUBj, j + 1, j + 1, j + 1, firstLBj, secondLBj, newLBj);

	// If newUBj < epsilon clip it to 0
	if (newUBj < epsilon) newUBj = 0.0;
	// If newLBj < epsilon clip it to 0
	if (newLBj < epsilon) newLBj = 0.0;
	// Update both no matter what
	delta_up[j] = newUBj;
	delta_down[j] = newLBj;
	/* was
	// Update UBj and LBj iff [ (they both do not fall below epsilon) || (they should both be set to zero due to an equality constraint) ]
	if (!(newUBj < epsilon && newLBj < epsilon) || (newUBj == 0.0 && newLBj == 0.0)) {
		inst->UB[j] = newUBj;
		inst->LB[j] = newLBj;
	} */
}