#include "ziround.h"

//**************************************************************************************************************************************************************
//*************************************** ZI-ROUND *************************************************************************************************************
//**************************************************************************************************************************************************************
int zi_round(instance* inst) {

	int status = inst->status;
	// inst->x_prev = (double*)malloc(inst->cur_numcols * sizeof(double));
	// inst->x_updated = (double*)malloc(inst->cur_numcols * sizeof(double));
	inst->UB = (double*)malloc(inst->cur_numcols * sizeof(double));
	inst->LB = (double*)malloc(inst->cur_numcols * sizeof(double));
	if (/*inst->x_prev == NULL || inst->x_updated == NULL ||*/ inst->UB == NULL || inst->LB == NULL) {
		fprintf(stderr, "[ERR][zi_round]: Failed to allocate x_prev, x_updated or row_infeas.\n"); return 1;
	}

	// Outer loop (repeat until no more updates found)
	do {

		// Initialize update flag and save current solution
		inst->updated = 0;
		// clone_array(inst->x, inst->x_prev, inst->cur_numcols);

		// Inner loop (for each variable xj that was integer/binary in the original MIP)
		for (int j = 0; j < inst->cur_numcols; j++) {
			// Skip xj if it's continuous in the original MIP
			if (!(inst->int_var)) continue;

			// xj non-fractional
			if (!is_fractional(inst->x[j])) {
				if (VERBOSE > 10) fprintf(stdout, "[DEBUG][zi_round]: -x- x_%d is non-fractional\n", j + 1);

				// Calculate UBj and LBj (with epsilon = 1.0)
				calculate_UBjLBj(inst, j, 1.0);
				// Skip xj if both UBj and LBj are equal to zero (--> no shift necessary)
				if (fabs(inst->UB[j]) < TOLERANCE && fabs(inst->LB[j]) < TOLERANCE) continue;

				// Condition(s) for rounding of xj
				// (>= to include the case of a zero obj coefficient)
				if ((inst->obj[j] >= 0 && fabs(inst->LB[j] - 1.0) < TOLERANCE) ||
					(inst->obj[j] <= 0 && fabs(inst->UB[j] - 1.0) < TOLERANCE)) {
					// Update xj to improve objective and update slacks
					update_xj_to_improve_objective(inst, j, 0); // flag xj non-fractional (0)
				}
			} // end if xj non-fractional

			// xj fractional
			else if (is_fractional(inst->x[j])) {
				if (VERBOSE > 10) fprintf(stdout, "[DEBUG][zi_round]: -x- x_%d is fractional\n", j + 1);

				// Calculate UBj and LBj
				calculate_UBjLBj(inst, j, EPSILON);
				// Skip xj if both UBj and LBj are equal to zero (--> no shift necessary)
				if (fabs(inst->UB[j]) < TOLERANCE && fabs(inst->LB[j]) < TOLERANCE) continue;

				// ZI-Round (version 1) core (3 cases). First, calculate the fractionalities needed.
				inst->ZI = fractionality(inst->x[j]);
				inst->ZIplus = fractionality(inst->x[j] + inst->UB[j]);
				inst->ZIminus = fractionality(inst->x[j] - inst->LB[j]);

				// First case
				if (fabs(inst->ZIplus - inst->ZIminus) < TOLERANCE &&
					(inst->ZIplus - inst->ZI) < -(TOLERANCE)) { // was inst->ZIplus < inst->ZI
					// Update xj to improve objective and update slacks
					update_xj_to_improve_objective(inst, j, 1); // flag xj fractional (1)
				} // end first case

				// Second case
				else if ((inst->ZIplus - inst->ZIminus) < -(TOLERANCE) && // was inst->ZIplus < inst->ZIminus
					(inst->ZIplus - inst->ZI) < -(TOLERANCE)) {      // was inst->ZIplus < inst->ZI
					if (VERBOSE > 10) fprintf(stdout, "[DEBUG][zi_round]: >>> Set x_%d = x_%d + UB_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], inst->UB[j], inst->x[j] + inst->UB[j]);
					inst->x[j] += inst->UB[j];
					inst->updated = 1;
					update_slacks(inst, j, inst->UB[j]);
					update_objective_value(inst, j, inst->UB[j]);
				} // end second case

				// Third case
				else if ((inst->ZIminus - inst->ZIplus) < -(TOLERANCE) && // was inst->ZIminus < inst->ZIplus
					(inst->ZIminus - inst->ZI) < -(TOLERANCE)) {     // was inst->ZIminus < inst->ZI
					if (VERBOSE > 10) fprintf(stdout, "[DEBUG][zi_round]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], inst->LB[j], inst->x[j] - inst->LB[j]);
					inst->x[j] -= inst->LB[j];
					inst->updated = 1;
					update_slacks(inst, j, -(inst->LB[j]));
					update_objective_value(inst, j, -(inst->LB[j]));
				} // end third case

			} // end if xj fractional

		} // end inner loop

		if (VERBOSE > 10) {
			if (inst->updated) fprintf(stdout, "[DEBUG][zi_round] ...Update found, scan variables again...\n");
			else fprintf(stdout, "[DEBUG][zi_round] ...No updates found, exit outer loop...\n");
		}
		// [BRUTE FORCE] [DEBUG ONLY] Check variable bounds and constraints
		if (VERBOSE > 50) {
			status = check_bounds(inst, inst->x);      if (status) { fprintf(stderr, "[ERR][zi_round]: Error inside check_bounds.\n"); return status; }
			status = check_constraints(inst, inst->x); if (status) { fprintf(stderr, "[ERR][zi_round]: Error inside check_constraints.\n"); return status; }
		}

	} while (inst->updated); // end outer loop

	return inst->status;
}
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************

void update_xj_to_improve_objective(instance* inst, int j, int is_fractional) {

	// First, calculate obj value for both updates. (IMPROVEMENT: inst->objval is just a common offset!)
	inst->obj_plusUBj = inst->objval + (inst->obj[j] * inst->UB[j]);
	inst->obj_minusLBj = inst->objval - (inst->obj[j] * inst->LB[j]);
	/* was
	// compare delta * coef (confronto solo tra i due delta * coef_j (+ objval è solo un offset) (tenendo conto della solita float tolerance)
	clone_array(inst->x, inst->x_updated, inst->cur_numcols);
	inst->x_updated[j] = inst->x[j] + inst->UB[j];
	..inst->obj_plusUBj = dot_product(inst->obj, inst->x_updated, inst->cur_numcols); // brute force!
	inst->x_updated[j] = inst->x[j] - inst->LB[j];
	..inst->obj_minusLBj = dot_product(inst->obj, inst->x_updated, inst->cur_numcols); // brute force!*/

	// Check obj sense, then update xj, update slacks and update objective value
	switch (inst->objsen) {
	case CPX_MIN:
		if ((inst->obj_plusUBj - inst->objval) < -(TOLERANCE) &&       // was inst->obj_plusUBj < inst->objval 
			(inst->obj_plusUBj - inst->obj_minusLBj) < -(TOLERANCE)) { // was inst->obj_plusUBj < inst->obj_minusLBj
			// xj = xj + UBj (if xj is not fractional then UBj must be 1.0)
			if (VERBOSE > 50) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: _>_ UBj = %f (must be 1.0)\n", inst->UB[j]);
			if (VERBOSE > 10) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: >>> Set x_%d = x_%d + UB_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], inst->UB[j], inst->x[j] + inst->UB[j]);
			if (!is_fractional && VERBOSE > 10 && fabs(inst->UB[j] - 1.0) > TOLERANCE) {
				fprintf(stderr, "[ERR][update_xj_to_improve_objective]: UB_%d = %f (should be 1.0).\n", j + 1, inst->UB[j]);
				exit(EXIT_FAILURE);
			}
			inst->x[j] += inst->UB[j];
			inst->updated = 1;
			update_slacks(inst, j, inst->UB[j]);
			update_objective_value(inst, j, inst->UB[j]);
		}
		else if ((inst->obj_minusLBj - inst->objval) < -(TOLERANCE) &&      // was inst->obj_minusLBj < inst->objval
			(inst->obj_minusLBj - inst->obj_plusUBj) < -(TOLERANCE)) { // was inst->obj_minusLBj < inst->obj_plusUBj
	   // xj = xj - LBj (if xj is not fractional then LBj must be 1.0)
			if (VERBOSE > 50) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: _>_ LBj = %f (must be 1.0)\n", inst->LB[j]);
			if (VERBOSE > 10) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], inst->LB[j], inst->x[j] - inst->LB[j]);
			if (!is_fractional && VERBOSE > 10 && fabs(inst->LB[j] - 1.0) > TOLERANCE) {
				fprintf(stderr, "[ERR][update_xj_to_improve_objective]: LB_%d = %f (should be 1.0).\n", j + 1, inst->LB[j]);
				exit(EXIT_FAILURE);
			}
			inst->x[j] -= inst->LB[j];
			inst->updated = 1;
			update_slacks(inst, j, -(inst->LB[j]));
			update_objective_value(inst, j, -(inst->LB[j]));
		}
		break;
		// (problems from mps files should always be MIN)
	case CPX_MAX:
		if ((inst->obj_plusUBj - inst->objval) > TOLERANCE &&       // was inst->obj_plusUBj > inst->objval
			(inst->obj_plusUBj - inst->obj_minusLBj) > TOLERANCE) { // was inst->obj_plusUBj > inst->obj_minusLBj
			// xj = xj + UBj (if xj is not fractional then UBj must be 1.0)
			if (VERBOSE > 50) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: _>_ UBj = %f (must be 1.0)\n", inst->UB[j]);
			if (VERBOSE > 10) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: >>> Set x_%d = x_%d + UB_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], inst->UB[j], inst->x[j] + inst->UB[j]);
			if (!is_fractional && VERBOSE > 10 && fabs(inst->UB[j] - 1.0) > TOLERANCE) {
				fprintf(stderr, "[ERR][update_xj_to_improve_objective]: UB_%d = %f (should be 1.0).\n", j + 1, inst->UB[j]);
				exit(EXIT_FAILURE);
			}
			inst->x[j] += inst->UB[j];
			inst->updated = 1;
			update_slacks(inst, j, inst->UB[j]);
			update_objective_value(inst, j, inst->UB[j]);
		}
		else if ((inst->obj_minusLBj - inst->objval) > TOLERANCE &&      // was inst->obj_minusLBj > inst->objval
			(inst->obj_minusLBj - inst->obj_plusUBj) > TOLERANCE) { // was inst->obj_minusLBj > inst->obj_plusUBj
	   // xj = xj - LBj (if xj is not fractional then LBj must be 1.0)
			if (VERBOSE > 50) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: _>_ LBj = %f (must be 1.0)\n", inst->LB[j]);
			if (VERBOSE > 10) fprintf(stdout, "[DEBUG][update_xj_to_improve_objective]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], inst->LB[j], inst->x[j] - inst->LB[j]);
			if (!is_fractional && VERBOSE > 10 && fabs(inst->LB[j] - 1.0) > TOLERANCE) {
				fprintf(stderr, "[ERR][update_xj_to_improve_objective]: LB_%d = %f (should be 1.0).\n", j + 1, inst->LB[j]);
				exit(EXIT_FAILURE);
			}
			inst->x[j] -= inst->LB[j];
			inst->updated = 1;
			update_slacks(inst, j, -(inst->LB[j]));
			update_objective_value(inst, j, -(inst->LB[j]));
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
	int end_col = (j < inst->cur_numcols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;
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
void update_objective_value(instance* inst, int j, double signed_delta) {
	inst->objval += (inst->obj[j] * signed_delta);
}

void calculate_UBjLBj(instance* inst, int j, const double epsilon) {
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
	int end_col = (j < inst->cur_numcols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;

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
	inst->UB[j] = newUBj;
	inst->LB[j] = newLBj;
	/* was
	// Update UBj and LBj iff [ (they both do not fall below epsilon) || (they should both be set to zero due to an equality constraint) ]
	if (!(newUBj < epsilon && newLBj < epsilon) || (newUBj == 0.0 && newLBj == 0.0)) {
		inst->UB[j] = newUBj;
		inst->LB[j] = newLBj;
	} */
}