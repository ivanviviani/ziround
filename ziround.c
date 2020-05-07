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

	// Local variables
	double* delta_up; 	 /**< Maximum variable up-shifts. */
	double* delta_down;  /**< Maximum variable down-shifts. */
	double ZI; 			 /**< Fractionality of a variable (used in function zi_round). */
	double ZIplus; 		 /**< Fractionality of a shifted up variable (used in function zi_round). */
	double ZIminus; 	 /**< Fractionality of a shifted down variable (used in fucntion zi_round). */
	double obj_plusUBj;  /**< Objective value for a shifted up variable (used in function zi_round). */
	double obj_minusLBj; /**< Objective value for a shifted down variable (used in function zi_round). */
	int updated;	     /**< Flag set to 1 when at least one variable shift has been made. */

	// Allocate
	delta_up   = (double*)malloc(inst->ncols * sizeof(double));
	delta_down = (double*)malloc(inst->ncols * sizeof(double));
	if (delta_up == NULL || delta_down == NULL) print_error("[zi_round]: Failed to allocate delta arrays.\n");

	// Outer loop (repeat until no more updates found)
	do {
		updated = 0; 

		// Inner loop (for each variable xj that was integer/binary in the original MIP)
		for (int j = 0; j < inst->ncols; j++) {

			if (!(inst->int_var[j])) continue; // Skip non-integer variables

			switch (is_fractional(inst->x[j])) {

				// xj non fractional
				case 0:

					print_verbose(150, "[zi_round]: -x- x_%d is non-fractional\n", j + 1);

					// Calculate deltas (with epsilon = 1.0)
					delta_updown(inst, j, delta_up, delta_down, 1.0);

					// Skip xj if both deltas are equal to zero (no shift necessary)
					if (fabs(delta_up[j] - 0.0) < TOLERANCE && fabs(delta_down[j] - 0.0) < TOLERANCE) continue;

					// Condition(s) for rounding of xj (>= to include the case of a zero obj coefficient)
					if ((inst->obj[j] > 0 && fabs(delta_down[j] - 1.0) < TOLERANCE) ||
						(inst->obj[j] < 0 && fabs(delta_up[j] - 1.0) < TOLERANCE)) {

						// Round xj to improve objective and update slacks
						updated = round_xj_bestobj(inst, j, delta_up, delta_down, 0); // flag xj non-fractional (0)
					}

					break;

				// xj fractional
				case 1:

					print_verbose(150, "[zi_round]: -x- x_%d is fractional\n", j + 1);

					// Calculate deltas
					delta_updown(inst, j, delta_up, delta_down, EPSILON);

					// Skip xj if both UBj and LBj are equal to zero (no shift necessary)
					if (fabs(delta_up[j] - 0.0) < TOLERANCE && fabs(delta_down[j] - 0.0) < TOLERANCE) continue;

					ZI      = fractionality(inst->x[j]);
					ZIplus  = fractionality(inst->x[j] + delta_up[j]);
					ZIminus = fractionality(inst->x[j] - delta_down[j]);

					// First case
					if (fabs(ZIplus - ZIminus) < TOLERANCE && ZIplus < ZI - TOLERANCE) {

						// Round xj to improve objective and update slacks
						updated = round_xj_bestobj(inst, j, delta_up, delta_down, 1); // flag xj fractional (1)
					}

					// Second case
					else if (ZIplus < ZIminus - TOLERANCE && ZIplus < ZI - TOLERANCE) {

						print_verbose(100, "[zi_round][FRA]        : >>> Set x_%d = x_%d + delta_up_%d = %f + %f = %f\n", 
							j + 1, j + 1, j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
						
						// Round UP
						inst->x[j] = inst->x[j] + delta_up[j];

						updated = 1;
						update_slacks(inst, j, delta_up[j]);
						inst->objval = inst->objval + (inst->obj[j] * delta_up[j]);
					}

					// Third case
					else if (ZIminus < ZIplus - TOLERANCE && ZIminus < ZI - TOLERANCE) {
						
						print_verbose(100, "[zi_round][FRA]        : >>> Set x_%d = x_%d - delta_down_%d = %f - %f = %f\n", 
							j + 1, j + 1, j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);
						
						// Round DOWN
						inst->x[j] = inst->x[j] - delta_down[j];

						updated = 1;
						update_slacks(inst, j, -(delta_down[j]));
						inst->objval = inst->objval - (inst->obj[j] * delta_down[j]);
					}

					break;

				default:
					print_error(" in function is_fractional.\n");
			}
		} // end inner loop

		
		if (updated) print_verbose(100, "[zi_round] ...Update found, scan variables again...\n");
		else print_verbose(100, "[zi_round] ...No updates found, exit outer loop...\n");
		
		// [DEBUG ONLY] [BRUTE FORCE]  Check variable bounds and constraints
		if (VERBOSE >= 150) {
			check_bounds(inst, inst->x);
			check_constraints(inst, inst->x);
		}

	} while (updated); // end outer loop

	// Free
	free(delta_up);
	free(delta_down);
}
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************

int round_xj_bestobj(instance* inst, int j, double* delta_up, double* delta_down, int xj_fractional) {

	// Local variables
	double obj_deltaplus;  /**< Delta obj if xj is shifted up. */
	double obj_deltaminus; /**< Delta obj if xj is shifted down. */
	int updated = 0;       /**< Flag set to 1 when at least one variable shift has been made. */
	char vtype[3];		   /**< Support string. */

	// Initialize
	obj_deltaplus  = 0.0 + (inst->obj[j] * delta_up[j]);
	obj_deltaminus = 0.0 - (inst->obj[j] * delta_down[j]);
	sprintf(vtype, ((xj_fractional)?"FRA":"INT"));

	// Check obj sense, then update xj, update slacks and update objective value
	switch (inst->objsen) {

		case CPX_MIN:

			if (obj_deltaplus < 0.0 - TOLERANCE && obj_deltaplus < obj_deltaminus - TOLERANCE) {

				print_verbose(100, "[round_xj_bestobj][%s]: >>> Set x_%d = x_%d + delta_up_%d = %f + %f = %f\n", vtype, j + 1, j + 1, j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
				if (!xj_fractional && VERBOSE >= 150 && fabs(delta_up[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: delta_up_%d = %f (should be 1.0).\n", j + 1, delta_up[j]);
				}

				// Round UP (if xj is not fractional then delta_up[j] must be 1.0)
				inst->x[j] = inst->x[j] + delta_up[j];

				updated = 1;
				update_slacks(inst, j, delta_up[j]);
				inst->objval = inst->objval + (inst->obj[j] * delta_up[j]);
			}

			else if (obj_deltaminus < 0.0 - TOLERANCE && obj_deltaminus < obj_deltaplus - TOLERANCE) {

				print_verbose(100, "[round_xj_bestobj][%s]: >>> Set x_%d = x_%d - delta_down_%d = %f - %f = %f\n", vtype, j + 1, j + 1, j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);
				if (!xj_fractional && VERBOSE >= 150 && fabs(delta_down[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: delta_down_%d = %f (should be 1.0).\n", j + 1, delta_down[j]);
				}

				// Round DOWN (if xj is not fractional then delta_down[j] must be 1.0)
				inst->x[j] = inst->x[j] - delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval = inst->objval - (inst->obj[j] * delta_down[j]);
			}

			else print_verbose(120, "[INFO][round_xj_bestobj]: obj_deltaplus == obj_deltaminus.\n");

			break;

		// (problems from mps files should always be MIN)
		case CPX_MAX:

			if (obj_deltaplus > 0.0 + TOLERANCE && obj_deltaplus > obj_deltaminus + TOLERANCE) {
				
				print_verbose(100, "[round_xj_bestobj][%s]: >>> Set x_%d = x_%d + delta_up_%d = %f + %f = %f\n", 
					vtype, j + 1, j + 1, j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
				if (!xj_fractional && VERBOSE >= 150 && fabs(delta_up[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: delta_up_%d = %f (should be 1.0).\n", j + 1, delta_up[j]);
				}

				// Round UP (if xj is not fractional then delta_up[j] must be 1.0)
				inst->x[j] = inst->x[j] + delta_up[j];

				updated = 1;
				update_slacks(inst, j, delta_up[j]);
				inst->objval = inst->objval + (inst->obj[j] * delta_up[j]);
			}

			else if (obj_deltaminus > 0.0 + TOLERANCE && obj_deltaminus > obj_deltaplus + TOLERANCE) {
				
				print_verbose(100, "[round_xj_bestobj][%s]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", 
					vtype, j + 1, j + 1, j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);
				if (!xj_fractional && VERBOSE >= 150 && fabs(delta_down[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: delta_down_%d = %f (should be 1.0).\n", j + 1, delta_down[j]);
				}

				// Round UP (if xj is not fractional then delta_down[j] must be 1.0)
				inst->x[j] = inst->x[j] - delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval = inst->objval - (inst->obj[j] * delta_down[j]);
			}

			else print_verbose(120, "[INFO][round_xj_bestobj]: obj_deltaplus == obj_deltaminus.\n");

			break;

		default:
			print_error("[round_xj_bestobj]: Objective sense '%d' not supported.\n", inst->objsen);
	}

	return updated;
}

// Incremental update of row slacks (only for the rows where the variable xj is involved)
// Whenever it is called, only one variable xj has been updated
void update_slacks(instance* inst, int j, double signed_delta) {

	// Local Variables
	//double aij; /**< Current coefficient. */
	//int rowind; /**< Current row index. */
	int colend; /**< Index of the last non-zero coefficient of the current column. */

	// Initialize
	colend = (j < inst->ncols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;

	// Scan constraints of variable j
	for (int k = inst->cmatbeg[j]; k < colend; k++) {

		double aij = inst->cmatval[k];
		int rowind = inst->cmatind[k];

		switch (inst->sense[rowind]) {

			case 'L':
			case 'G':
				inst->slack[rowind] = inst->slack[rowind] - aij * signed_delta;
				break;
			case 'E':
				// Distinguish "special" equality constraints (if extension enabled)
				if (inst->extension && inst->row_slack_var[rowind] != -1) {

					print_verbose(120, "[update_slacks]: ENTERED SPECIAL EQUALITY CASE.\n");

					int slack_var = inst->row_slack_var[rowind];
					double slack_coef = 0.0;

					// Update (continuous) slack variable AND objective value
					int rowend = (rowind < inst->nrows - 1) ? inst->rmatbeg[rowind] : inst->nzcnt;
					double sum = 0.0;
					for (int h = inst->rmatbeg[rowind]; h < rowend; h++) {

						int varind = inst->rmatind[h];

						if (varind == slack_var) {
							slack_coef = inst->rmatval[h];
							continue;
						}
						
						sum += inst->rmatval[h] * inst->x[varind];
					}
					inst->objval = inst->objval - (inst->obj[slack_var] * inst->x[slack_var]);
					inst->x[slack_var] = (inst->rhs[rowind] - sum) / slack_coef;
					inst->objval = inst->objval + (inst->obj[slack_var] * inst->x[slack_var]);
				}
				else {
					print_error("[update_slacks]: Tried to update slack of an equality constraint (with no slack variable)!\n");
				}
				break;
			default:
				print_error("[update_slacks]: Constraint sense %c not supported!\n", inst->sense[rowind]);
		}
	}
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
	delta_up1   = LONG_MAX;
	delta_down1 = LONG_MAX;
	delta_up2   = inst->ub[j] - inst->x[j];
	delta_down2 = inst->x[j] - inst->lb[j];
	colend      = (j < inst->ncols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;
	
	print_verbose(200, "[delta_updown]: delta_up2_%d = ub_%d - x_%d = %f - %f = %f ; delta_down2_%d = x_%d - lb_%d = %f - %f = %f\n",
		j + 1, j + 1, j + 1, inst->ub[j], inst->x[j], delta_up2, j + 1, j + 1, j + 1, inst->x[j], inst->lb[j], delta_down2);

	// Scan non-zero coefficients of column j
	for (int k = inst->cmatbeg[j]; k < colend; k++) {

		aij = inst->cmatval[k];
		rowind = inst->cmatind[k];
		
		// Check sense, then check sign of aij, and update delta_up1, delta_down1
		switch (inst->sense[rowind]) {

			case 'L': // (slack non-negative)

				// Clip slack to zero if negative
				if (inst->slack[rowind] < 0.0 - TOLERANCE) inst->slack[rowind] = 0.0;

				if (aij > 0.0) { 
					
					print_verbose(200, "[delta_updown]: sense = L ; slack[%d] = %f ; a_%d_%d = %f ; candidate_up1_%d = %f\n",
						rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1, inst->slack[rowind] / aij);
					
					// Update delta_up1
					candidate_up1 = inst->slack[rowind] / aij;
					delta_up1 = min(candidate_up1, delta_up1);
				}
				if (aij < 0.0) {

					print_verbose(200, "[delta_updown]: sense = L ; slack[%d] = %f ; a_%d_%d = %f ; candidate_down1_%d = %f\n",
						rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1, -(inst->slack[rowind]) / aij);

					// Update delta_down1
					candidate_down1 = -(inst->slack[rowind]) / aij;
					delta_down1 = min(candidate_down1, delta_down1);
				}

				break;

			case 'G': // (slack non-positive)

				// Clip slack to zero if positive
				if (inst->slack[rowind] > 0.0 + TOLERANCE) inst->slack[rowind] = 0.0;

				if (aij < 0.0) {

					print_verbose(200, "[delta_updown]: sense = G ; slack[%d] = %f ; a_%d_%d = %f ; candidate_up1_%d = %f\n",
						rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1, inst->slack[rowind] / aij);

					// Update delta_up1
					candidate_up1 = inst->slack[rowind] / aij;
					delta_up1 = min(candidate_up1, delta_up1);
				}
				if (aij > 0.0) {

					print_verbose(200, "[delta_updown]: sense = G ; slack[%d] = %f ; a_%d_%d = %f ; candidate_down1_%d = %f\n",
						rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1, -(inst->slack[rowind]) / aij);

					// Update delta_down1
					candidate_down1 = -(inst->slack[rowind]) / aij;
					delta_down1 = min(candidate_down1, delta_down1);
				}

				break;

			case 'E': // (slack zero)

				// Distinguish "special" equality constraints (if extension enabled)
				if (inst->extension && inst->row_slack_var[rowind] != -1) {

					int slack_var = inst->row_slack_var[rowind];
					double slack_coef;
					CPXgetcoef(inst->env, inst->lp, rowind, slack_var, &slack_coef);
					double slack = slack_coef * inst->x[slack_var];

					// DEBUG
					int rowend = (rowind < inst->nrows - 1) ? inst->rmatbeg[rowind + 1] : inst->nzcnt;
					double rowact = 0.0;
					for (int i = inst->rmatbeg[rowind]; i < rowend; i++) {
						int varind = inst->rmatind[i];
						rowact += inst->rmatval[i] * inst->x[varind];
					}
					print_verbose(120, "[delta_updown][DEBUG]: rowact %f = rhs %f | slack_var x_%d = %f | slack %f\n", rowact, inst->rhs[rowind], slack_var, inst->x[slack_var], slack);
					// DEBUG

					// Check slack sign
					if (slack > 0.0 + TOLERANCE) {

						// Without the slack, the constraint would be 'L' (slack non-negative)
						if (aij > 0.0) {

							print_verbose(200, "[delta_updown]: sense = L ; slack[%d] = %f ; a_%d_%d = %f ; candidate_up1_%d = %f\n",
								rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1, inst->slack[rowind] / aij);

							// Update delta_up1
							candidate_up1 = slack / aij;
							delta_up1 = min(candidate_up1, delta_up1);
						}
						if (aij < 0.0) {

							print_verbose(200, "[delta_updown]: sense = L ; slack[%d] = %f ; a_%d_%d = %f ; candidate_down1_%d = %f\n",
								rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1, -(inst->slack[rowind]) / aij);

							// Update delta_down1
							candidate_down1 = -(slack) / aij;
							delta_down1 = min(candidate_down1, delta_down1);
						}
					}
					else if (slack < 0.0 - TOLERANCE) {

						// Without the slack, the constraint would be 'G' (slack non-positive)
						if (aij < 0.0) {

							print_verbose(200, "[delta_updown]: sense = G ; slack[%d] = %f ; a_%d_%d = %f ; candidate_up1_%d = %f\n",
								rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1, inst->slack[rowind] / aij);

							// Update delta_up1
							candidate_up1 = slack / aij;
							delta_up1 = min(candidate_up1, delta_up1);
						}
						if (aij > 0.0) {

							print_verbose(200, "[delta_updown]: sense = G ; slack[%d] = %f ; a_%d_%d = %f ; candidate_down1_%d = %f\n",
								rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1, -(inst->slack[rowind]) / aij);

							// Update delta_down1
							candidate_down1 = -(slack) / aij;
							delta_down1 = min(candidate_down1, delta_down1);
						}
					}
					else {
						print_verbose(200, "[delta_updown][extension]: sense = E ; Variable x_%d involved in equality constraint %d with slack zero --> cannot be moved!\n", j + 1, rowind);

						// Set delta_up1 and delta_down1 to zero --> new_delta_up and new_delta_down will get value zero
						delta_up1 = 0.0;
						delta_down1 = 0.0;
					}
				}
				else {
					print_verbose(200, "[delta_updown]: sense = E ; Variable x_%d involved in equality constraint %d with slack zero --> cannot be moved!\n", j + 1, rowind);

					// Set delta_up1 and delta_down1 to zero --> new_delta_up and new_delta_down will get value zero
					delta_up1 = 0.0;
					delta_down1 = 0.0;
				}

				break;

			default:
				print_error("[delta_updown]: Constraint sense '%c' not included in {'L','G','E'}.\n", inst->sense[rowind]);
		}
	} // end for

	print_verbose(200, "[delta_updown][candidates]: delta_up1_%d = %f ; delta_down1_%d = %f\n", j + 1, delta_up1, j + 1, delta_down1);

	// Results
	new_delta_up = min(delta_up1, delta_up2);
	new_delta_down = min(delta_down1, delta_down2);
	print_verbose(200, "[delta_updown][results]: (NEW) delta_up_%d = min{%f, %f} = %f ; delta_down_%d = min{%f, %f} = %f\n",
		j + 1, delta_up1, delta_up2, new_delta_up, j + 1, delta_down1, delta_down2, new_delta_down);

	// Update deltas
	if (new_delta_up < epsilon - TOLERANCE && new_delta_down < epsilon - TOLERANCE) {
		new_delta_up = 0.0;
		new_delta_down = 0.0;
	}
	delta_up[j] = new_delta_up;
	delta_down[j] = new_delta_down;
}