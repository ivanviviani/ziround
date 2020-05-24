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

					// Skip xj if both deltas are equal to zero (no shift necessary)
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

						// Check whether all affected constraints have enough slack for a ROUND UP of xj
						check_slacks(inst, j, delta_up, delta_down, 'U');
						
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

						// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
						check_slacks(inst, j, delta_up, delta_down, 'D');
						
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
		
		if (updated) { print_verbose(100, "[zi_round] ...Update found, scan variables again...\n"); }
		else { print_verbose(100, "[zi_round] ...No updates found, exit outer loop...\n"); }
		
		// [DEBUG ONLY] [BRUTE FORCE]  Check variable bounds and constraints
		if (VERBOSE >= 120) {
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

void check_slacks(instance* inst, int j, double* delta_up, double* delta_down, const char round_updown) {

	if (round_updown != 'U' && round_updown != 'D') print_error("[check_slacks]: Rounding sense '%c' undefined.\n", round_updown);

	// Check whether all affected constraints have enough slack for a ROUND UP/DOWN of xj
	int colend = (j < inst->ncols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;
	for (int k = inst->cmatbeg[j]; k < colend; k++) {

		int rowind = inst->cmatind[k];
		double aij = inst->cmatval[k];
		double curr_slack = 0.0;
		double singletons_slack = 0.0;
		double delta_slack = 0.0;
		double new_slack = 0.0;
		char msg[17];
		int enough_slack = 0;

		switch (inst->sense[rowind]) {

			case 'L': // (slack non-negative)

				curr_slack = inst->slack[rowind];
				delta_slack = (round_updown == 'U') ? (aij * delta_up[j]) : (aij * (-delta_down[j]));

				// Row slack after rounding (negative iff also singletons slack should be used)
				new_slack = curr_slack - delta_slack;
				print_verbose(20, "[check_slacks][x_%d aij %f][row %d '%c']: new_slack = %f ", j + 1, rowind + 1, inst->sense[rowind], new_slack);
				if (new_slack < -(TOLERANCE)) print_verbose(20, "(need singletons slack!)\n");
				else print_verbose(20, "\n");

				// [EXTENSION] Distinguish inequality constraints with singletons (if extension enabled)
				if (inst->extension && inst->num_singletons[rowind] > 0) {

					// If the new row slack is negative then the remaining amount must be covered by the singletons slack
					if (new_slack < -(TOLERANCE)) {

						// Compute singletons slack of constraint rowind and get bounds
						double ss_lb = inst->ss_lb[rowind];
						double ss_ub = inst->ss_ub[rowind];
						singletons_slack = compute_singletons_slack(inst, rowind);

						// Singletons slack after rounding
						double delta_ss = new_slack;
						double new_ss = singletons_slack + delta_ss; // + because signed delta

						// New singletons slack must stay within its bounds
						sprintf(msg, (!(new_ss < ss_lb - TOLERANCE || new_ss > ss_ub + TOLERANCE)) ? "Enough slack [OK]" : "NOT enough slack!");
					}
				}
				else {
					// Extension disabled
					if (!(inst->extension)) sprintf(msg, (new_slack > -(TOLERANCE)) ? "Enough slack [OK]" : "NOT enough slack!");

					// If new_slack was negative then the extension must be enabled, otherwise an unknown error occured
					if (new_slack < -(TOLERANCE)) print_error("[check_slacks][x_%d aij %f][row %d '%c']: new_slack negative but no singletons to cover it.", j + 1, rowind + 1, inst->sense[rowind]);
				}

				print_verbose(120, "[check_slacks][x_%d aij %f][row %d '%c']: curr_slack %f (-) delta_slack = %f * %f = %f | %s\n", j + 1, aij, rowind + 1, inst->sense[rowind], curr_slack, aij, (round_updown == 'U') ? delta_up[j] : (-delta_down[j]), delta_slack, msg);

				if (msg[0] == 'N') print_error("[check_slacks][x_%d][row %d '%c']: After rounding, invalid slack. Found %f < 0.\n", j + 1, rowind + 1, inst->sense[rowind], new_slack);

				break;

			case 'G': // (slack non-positive)
				
				curr_slack = inst->slack[rowind];
				delta_slack = (round_updown == 'U') ? (aij * delta_up[j]) : (aij * (-delta_down[j]));

				// Row slack after rounding (positive iff also singletons slack should be used)
				new_slack = curr_slack - delta_slack;
				print_verbose(20, "[check_slacks][x_%d aij %f][row %d '%c']: new_slack = %f ", j + 1, rowind + 1, inst->sense[rowind], new_slack);
				if (new_slack > TOLERANCE) print_verbose(20, "(need singletons slack!)\n");
				else print_verbose(20, "\n");

				// [EXTENSION] Distinguish inequality constraints with singletons (if extension enabled)
				if (inst->extension && inst->num_singletons[rowind] > 0) {

					// If the new row slack is positive then the remaining amount must be covered by the singletons slack
					if (new_slack > TOLERANCE) {

						// Compute singletons slack of constraint rowind and get bounds
						double ss_lb = inst->ss_lb[rowind];
						double ss_ub = inst->ss_ub[rowind];
						singletons_slack = compute_singletons_slack(inst, rowind);

						// Singletons slack after rounding
						double delta_ss = new_slack;
						double new_ss = singletons_slack + delta_ss; // + because signed delta

						// New singletons slack must stay within its bounds
						sprintf(msg, (!(new_ss < ss_lb - TOLERANCE || new_ss > ss_ub + TOLERANCE)) ? "Enough slack [OK]" : "NOT enough slack!");
					}
				}
				else {
					// Extension disabled
					if (!(inst->extension)) sprintf(msg, (new_slack < TOLERANCE) ? "Enough slack [OK]" : "NOT enough slack!");

					// If new_slack was positive then the extension must be enabled, otherwise an unknown error occured
					if (new_slack > TOLERANCE) print_error("[check_slacks][x_%d aij %f][row %d '%c']: new_slack positive but no singletons to cover it.", j + 1, rowind + 1, inst->sense[rowind]);
				}

				print_verbose(120, "[check_slacks][x_%d aij %f][row %d '%c']: curr_slack %f (-) delta_slack = %f * %f = %f | %s\n", j + 1, aij, rowind + 1, inst->sense[rowind], curr_slack, aij, (round_updown == 'U') ? delta_up[j] : (-delta_down[j]), delta_slack, msg);

				if (msg[0] == 'N') print_error("[check_slacks][x_%d][row %d '%c']: After rounding, invalid slack. Found %f > 0.\n", j + 1, rowind + 1, inst->sense[rowind], new_slack);

				break;

			case 'E': // (slack zero if extension disabled)

				// [EXTENSION] Distinguish equality constraints with singletons (if extension enabled)
				if (inst->extension && inst->num_singletons[rowind] > 0) {

					// Compute singletons slack of constraint rowind (with bounds)
					double ss_lb = inst->ss_lb[rowind];
					double ss_ub = inst->ss_ub[rowind];
					curr_slack = compute_singletons_slack(inst, rowind);
					delta_slack = (round_updown == 'U') ? (aij * delta_up[j]) : (aij * (-delta_down[j]));

					// Singletons slack after rounding
					double new_ss = curr_slack - delta_slack;
					sprintf(msg, (!(new_ss < ss_lb - TOLERANCE || new_ss > ss_ub + TOLERANCE)) ? "Enough slack [OK]" : "NOT enough slack!");
						
					print_verbose(120, "[check_slacks][extension][x_%d aij %f][row %d '%c']: curr_slack %f (-) delta_slack = %f * %f = %f | %s\n", j + 1, aij, rowind + 1, inst->sense[rowind], curr_slack, aij, (round_updown == 'U') ? delta_up[j] : (-delta_down[j]), delta_slack, msg);

					if (msg[0] == 'N') print_error("[check_slacks][extension][x_%d][row %d '%c']: After rounding, singletons slack out of bounds. Found %f <= %f <= %f.\n", j + 1, rowind + 1, inst->sense[rowind], ss_lb, new_ss, ss_ub);
				}
				else {
					// Extension disabled
					print_error("[check_slacks][x_%d][row %d '%c']: Constraint has no singletons --> slack ZERO --> x_%d cannot be rounded.\n", j + 1, rowind + 1, inst->sense[rowind], j + 1);
				}

				break;

			default:
				print_error("[check_slacks]: Constraint sense '%c' not included in {'L','G','E'}.\n", inst->sense[rowind]);
		}

	}
}

int round_xj_bestobj(instance* inst, int j, double* delta_up, double* delta_down, int xj_fractional) {

	// Local variables
	double obj_deltaplus;  /**< Delta obj if xj is shifted up. */
	double obj_deltaminus; /**< Delta obj if xj is shifted down. */
	int updated = 0;       /**< Flag set to 1 when at least one variable shift has been made. */
	char vtype[3];		   /**< Support string. */

	// Initialize
	obj_deltaplus  = (inst->obj[j] * delta_up[j]);
	obj_deltaminus = -(inst->obj[j] * delta_down[j]);
	sprintf(vtype, ((xj_fractional)?"FRA":"INT"));

	print_verbose(120, "[round_xj_bestobj]: %f <= x_%d = %f <= %f | delta_up = %f delta down = %f objcoef = %f | obj_deltaplus = %f obj_deltaminus = %f\n", 
		inst->lb[j], j + 1, inst->x[j], inst->ub[j], delta_up[j], delta_down[j], inst->obj[j], obj_deltaplus, obj_deltaminus);

	// Check obj sense, then update xj, update slacks and update objective value
	switch (inst->objsen) {

		case CPX_MIN:

			if (obj_deltaplus < 0.0 - TOLERANCE && obj_deltaplus < obj_deltaminus - TOLERANCE) {

				print_verbose(100, "[round_xj_bestobj][%s]: >>> Set x_%d = x_%d + delta_up_%d = %f + %f = %f\n", vtype, j + 1, j + 1, j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
				if (!xj_fractional && VERBOSE >= 150 && fabs(delta_up[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: delta_up_%d = %f (should be 1.0).\n", j + 1, delta_up[j]);
				}

				// Check whether all affected constraints have enough slack for a ROUND UP of xj
				check_slacks(inst, j, delta_up, delta_down, 'U');

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

				// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
				check_slacks(inst, j, delta_up, delta_down, 'D');

				// Round DOWN (if xj is not fractional then delta_down[j] must be 1.0)
				inst->x[j] = inst->x[j] - delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval = inst->objval - (inst->obj[j] * delta_down[j]);
			}

			else if (obj_deltaminus < 0.0 - TOLERANCE && obj_deltaplus < 0.0 - TOLERANCE && fabs(obj_deltaminus - obj_deltaplus) < TOLERANCE) {

				print_verbose(120, "[round_xj_bestobj]: obj_deltaplus = %f = %f = obj_deltaminus. >>> Round x_%d arbitrarily (DOWN).\n", obj_deltaplus, obj_deltaminus, j + 1);

				// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
				check_slacks(inst, j, delta_up, delta_down, 'D');

				// Round arbitrarily (DOWN)
				inst->x[j] = inst->x[j] - delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval = inst->objval - (inst->obj[j] * delta_down[j]);
			}

			break;

		// (problems from mps files should always be MIN)
		case CPX_MAX:

			if (obj_deltaplus > 0.0 + TOLERANCE && obj_deltaplus > obj_deltaminus + TOLERANCE) {
				
				print_verbose(100, "[round_xj_bestobj][%s]: >>> Set x_%d = x_%d + delta_up_%d = %f + %f = %f\n", 
					vtype, j + 1, j + 1, j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
				if (!xj_fractional && VERBOSE >= 150 && fabs(delta_up[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: delta_up_%d = %f (should be 1.0).\n", j + 1, delta_up[j]);
				}

				// Check whether all affected constraints have enough slack for a ROUND UP of xj
				check_slacks(inst, j, delta_up, delta_down, 'U');

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

				// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
				check_slacks(inst, j, delta_up, delta_down, 'D');

				// Round DOWN (if xj is not fractional then delta_down[j] must be 1.0)
				inst->x[j] = inst->x[j] - delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval = inst->objval - (inst->obj[j] * delta_down[j]);
			}

			else if (obj_deltaminus > 0.0 + TOLERANCE && obj_deltaplus > 0.0 + TOLERANCE && fabs(obj_deltaminus - obj_deltaplus) < TOLERANCE) {
				
				print_verbose(120, "[round_xj_bestobj]: obj_deltaplus = %f = %f = obj_deltaminus. >>> Round x_%d arbitrarily (DOWN).\n", obj_deltaplus, obj_deltaminus, j + 1);

				// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
				check_slacks(inst, j, delta_up, delta_down, 'D');

				// Round arbitrarily (DOWN)
				inst->x[j] = inst->x[j] - delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval = inst->objval - (inst->obj[j] * delta_down[j]);
			}

			break;

		default:
			print_error("[round_xj_bestobj]: Objective sense '%d' not supported.\n", inst->objsen);
	}

	return updated;
}

// Incremental update of row slacks (only for the rows where the variable xj is involved)
// Whenever it is called, only one variable xj has been updated
void update_slacks(instance* inst, int j, double signed_delta) {

	int colend = (j < inst->ncols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;

	// Scan constraints of variable j
	for (int k = inst->cmatbeg[j]; k < colend; k++) {

		double aij = inst->cmatval[k];
		int rowind = inst->cmatind[k];
		double delta_slack = aij * signed_delta;
		double curr_slack = inst->slack[rowind];
		double temp_slack = 0.0;

		switch (inst->sense[rowind]) {

			case 'L': // (slack non-negative)

				// [EXTENSION] Distinguish inequality constraints with singletons (if extension enabled)
				if (inst->extension && inst->num_singletons[rowind] > 0) {

					// First, use at most all the row slack available to cover delta_slack
					temp_slack = curr_slack - delta_slack;
					// Update row slack
					inst->slack[rowind] = max(0.0, temp_slack);

					// If not enough row slack (temp_slack negative), resort to singletons slack
					if (temp_slack < -(TOLERANCE)) {

						// Delta singletons slack to distribute among the singletons [new_ss = ss + delta_ss (+ because signed delta)]
						double delta_ss = temp_slack;

						// Distribute delta among the singletons, stop when done (delta_ss negative --> singletons slack must decrease)
						int beg = inst->rs_beg[rowind];
						for (int k = 0; k < inst->num_singletons[rowind]; k++) {

							// Stop updating the singletons when delta singletons slack has been covered
							if (delta_ss > -(TOLERANCE)) break;
							print_verbose(120, "[update_slacks][extension][row %d '%c']: Remaining delta singletons slack to distribute: %f.\n", rowind + 1, inst->sense[rowind], delta_ss);

							// Singleton info
							int singleton_index = inst->row_singletons[beg + k];
							double coef = inst->rs_coef[beg + k];
							double s_lb = inst->lb[singleton_index];
							double s_ub = inst->ub[singleton_index];
							double s_val = inst->x[singleton_index];
							double covered_delta_ss = 0.0;
							double max_s_delta = 0.0;
							double s_delta = 0.0;

							// Compute covered delta of the singleton
							if (coef > 0.0) {
								max_s_delta = s_val - s_lb;
								covered_delta_ss = max(delta_ss, -coef * max_s_delta);
							}
							if (coef < 0.0) {
								max_s_delta = s_ub - s_val;
								covered_delta_ss = max(delta_ss, coef * max_s_delta);
							}
							// Update remaining delta to be covered by the next singletons
							delta_ss = delta_ss - covered_delta_ss;

							// Compute singleton delta
							s_delta = covered_delta_ss / coef;
							// Update singleton
							inst->x[singleton_index] = s_val + s_delta;

							// Update objective value
							inst->objval = inst->objval + (inst->obj[singleton_index] * s_delta);
						} // end for
					}
					else {
						// Enough row slack, already updated above
					}
				}
				else {
					// Extension disabled
					if (!(inst->extension)) {
						// Just update row slack
						print_verbose(120, "[update_slacks][x_%d][row %d '%c']: slack = %f - (%f * %f) = %f.\n", j + 1, rowind + 1, inst->sense[rowind], inst->slack[rowind], aij, signed_delta, inst->slack[rowind] - aij * signed_delta);
						inst->slack[rowind] = inst->slack[rowind] - aij * signed_delta;
					}

					// If temp_slack was negative then the extension must be enabled, otherwise an unknown error occured
					if (temp_slack < -(TOLERANCE)) print_error("[update_slacks][x_%d aij %f][row %d '%c']: temp_slack negative but no singletons to cover it.", j + 1, rowind + 1, inst->sense[rowind]);
				}

				break;

			case 'G': // (slack non-positive)

				// [EXTENSION] Distinguish inequality constraints with singletons (if extension enabled)
				if (inst->extension && inst->num_singletons[rowind] > 0) {

					// First, use at most all the row slack available to cover delta_slack
					temp_slack = curr_slack - delta_slack;
					// Update row slack
					inst->slack[rowind] = min(0.0, temp_slack);

					// If not enough row slack (temp_slack positive), resort to singletons slack
					if (temp_slack > TOLERANCE) {

						// Delta singletons slack to distribute among the singletons [new_ss = ss + delta_ss (+ because signed delta)]
						double delta_ss = temp_slack;

						// Distribute delta among the singletons, stop when done (delta_ss positive --> singletons slack must increase)
						int beg = inst->rs_beg[rowind];
						for (int k = 0; k < inst->num_singletons[rowind]; k++) {

							// Stop updating the singletons when delta singletons slack has been covered
							if (delta_ss < TOLERANCE) break;
							print_verbose(120, "[update_slacks][extension][row %d '%c']: Remaining delta singletons slack to distribute: %f.\n", rowind + 1, inst->sense[rowind], delta_ss);

							// Singleton info
							int singleton_index = inst->row_singletons[beg + k];
							double coef = inst->rs_coef[beg + k];
							double s_lb = inst->lb[singleton_index];
							double s_ub = inst->ub[singleton_index];
							double s_val = inst->x[singleton_index];
							double covered_delta_ss = 0.0;
							double max_s_delta = 0.0;
							double s_delta = 0.0;

							// Compute covered delta of the singleton
							if (coef > 0.0) {
								max_s_delta = s_ub - s_val;
								covered_delta_ss = min(delta_ss, coef * max_s_delta);
							}
							if (coef < 0.0) {
								max_s_delta = s_val - s_lb;
								covered_delta_ss = min(delta_ss, -coef * max_s_delta);
							}
							// Update remaining delta to be covered by the next singletons
							delta_ss = delta_ss - covered_delta_ss;

							// Compute singleton delta
							s_delta = covered_delta_ss / coef;
							// Update singleton
							inst->x[singleton_index] = s_val + s_delta;

							// Update objective value
							inst->objval = inst->objval + (inst->obj[singleton_index] * s_delta);
						} // end for
					}
					else {
						// Enough row slack, already updated above
					}
				}
				else {
					// Extension disabled
					if (!(inst->extension)) {
						// Just update row slack
						print_verbose(120, "[update_slacks][x_%d][row %d '%c']: slack = %f - (%f * %f) = %f.\n", j + 1, rowind + 1, inst->sense[rowind], inst->slack[rowind], aij, signed_delta, inst->slack[rowind] - aij * signed_delta);
						inst->slack[rowind] = inst->slack[rowind] - aij * signed_delta;
					}

					// If temp_slack was negative then the extension must be enabled, otherwise an unknown error occured
					if (temp_slack < -(TOLERANCE)) print_error("[update_slacks][x_%d aij %f][row %d '%c']: temp_slack negative but no singletons to cover it.", j + 1, rowind + 1, inst->sense[rowind]);
				}

				break;

			case 'E':

				// [EXTENSION] Distinguish equality constraints with singletons (if extension enabled)
				if (inst->extension && inst->num_singletons[rowind] > 0) {

					// Compute singletons slack of constraint rowind (with bounds)
					double ss_lb = inst->ss_lb[rowind];
					double ss_ub = inst->ss_ub[rowind];
					double singletons_slack = compute_singletons_slack(inst, rowind);
					
					// Compute new singletons slack (after rounding xj):
					// xj contribution increases (decreases) --> singletons slack decreases (increases) of the same amount
					double delta_slack = aij * signed_delta;
					if (fabs(delta_slack) < TOLERANCE) print_error("[update_slacks][extension][x_%d][row %d '%c']: Found delta_slack ZERO.\n");
					double new_ss = singletons_slack - delta_slack;
					
					// Update singletons
					if (inst->num_singletons[rowind] == 1) {

						int beg = inst->rs_beg[rowind];
						int singleton_index = inst->row_singletons[beg];
						double coef = inst->rs_coef[beg];
						double new_singleton_val = new_ss / coef;

						// Check bounds of new singleton value
						if (!(new_singleton_val > inst->lb[singleton_index] - TOLERANCE &&
							  new_singleton_val < inst->ub[singleton_index] + TOLERANCE))
							print_error("[update_slacks][extension][singleton x_%d][row %d '%c']: Updated singleton out of bounds. Found %f <= %f <= %f\n",
								singleton_index + 1, rowind + 1, inst->sense[rowind], inst->lb[singleton_index], new_singleton_val, inst->ub[singleton_index]);

						// Update the only singleton of the equality
						double signed_delta = new_singleton_val - inst->x[singleton_index];
						inst->x[singleton_index] = new_singleton_val;

						// Update objective value
						inst->objval = inst->objval + (inst->obj[singleton_index] * signed_delta);
					}
					else if (inst->num_singletons[rowind] > 1) {

						// New singletons slack to distribute among the singletons
						double remaining = new_ss;
						print_verbose(120, "[update_slacks][extension][row %d '%c']: Singletons slack to distribute %f\n", rowind + 1, inst->sense[rowind], remaining);

						int beg = inst->rs_beg[rowind];

						// Positive singletons slack to distribute
						if (new_ss > TOLERANCE) {

							// Scan singletons
							for (int k = 0; k < inst->num_singletons[rowind]; k++) {

								print_verbose(120, "[update_slacks][extension][row %d '%c']: Remaining singletons slack to distribute %f.\n",
									rowind + 1, inst->sense[rowind], remaining);

								// Singleton info
								int singleton_index = inst->row_singletons[beg + k];
								double coef = inst->rs_coef[beg + k];
								double s_lb = inst->lb[singleton_index];
								double s_ub = inst->ub[singleton_index];

								print_verbose(120, "[update_slacks][extension][row %d '%c'][singleton x_%d]: %f <= ss <= %f | ",
									rowind + 1, inst->sense[rowind], singleton_index + 1, s_lb, s_ub);

								// If singletons slack distributed already, try to set the following singletons to ZERO
								if (fabs(remaining) < TOLERANCE) {

									if (s_lb < TOLERANCE && s_ub > -(TOLERANCE)) {

										double new_singleton_val = 0.0;

										// Update singleton
										double signed_delta = new_singleton_val - inst->x[singleton_index];
										inst->x[singleton_index] = new_singleton_val;

										// Update objective value
										inst->objval = inst->objval + (inst->obj[singleton_index] * signed_delta);
									}
									else print_error("[update_slacks][extension][row %d '%c'][singleton x_%d]: Failed to set singleton to ZERO.\n",
										rowind + 1, inst->sense[rowind], singleton_index + 1);
								}
								// Remaining singletons slack to distribute
								else {

									double new_singleton_val = 0.0; // initialized

									// Distribute max possible amount of singletons slack
									if (coef > 0.0 && s_ub > -(TOLERANCE)) {

										new_singleton_val = min(s_ub, remaining / coef);
									}
									if (coef < 0.0 && s_lb < TOLERANCE) {

										new_singleton_val = max(s_lb, remaining / coef);
									}

									// Check bounds of new singleton value
									if (!(new_singleton_val > inst->lb[singleton_index] - TOLERANCE &&
										new_singleton_val < inst->ub[singleton_index] + TOLERANCE))
										print_error("[update_slacks][extension][singleton x_%d][row %d '%c']: Updated singleton out of bounds. Found %f <= %f <= %f\n",
											singleton_index + 1, rowind + 1, inst->sense[rowind], inst->lb[singleton_index], new_singleton_val, inst->ub[singleton_index]);

									// Update remaining singletons slack to distribute
									remaining = remaining - (coef * new_singleton_val);
									print_verbose(120, "[update_slacks][extension][row %d '%c']: Remaining singletons slack to distribute %f.\n",
										rowind + 1, inst->sense[rowind], remaining);

									// Update singleton
									double signed_delta = new_singleton_val - inst->x[singleton_index];
									inst->x[singleton_index] = new_singleton_val;

									// Update objective value
									inst->objval = inst->objval + (inst->obj[singleton_index] * signed_delta);
								}
							}
						}
						// Negative singletons slack to distribute
						else if (new_ss < -(TOLERANCE)) {

							// Scan singletons
							for (int k = 0; k < inst->num_singletons[rowind]; k++) {

								print_verbose(120, "[update_slacks][extension][row %d '%c']: Remaining singletons slack to distribute %f.\n",
									rowind + 1, inst->sense[rowind], remaining);

								// Singleton info
								int singleton_index = inst->row_singletons[beg + k];
								double coef = inst->rs_coef[beg + k];
								double s_lb = inst->lb[singleton_index];
								double s_ub = inst->ub[singleton_index];

								print_verbose(120, "[update_slacks][extension][row %d '%c'][singleton x_%d]: %f <= ss <= %f | ",
									rowind + 1, inst->sense[rowind], singleton_index + 1, s_lb, s_ub);

								// If singletons slack distributed already, try to set the following singletons to ZERO
								if (fabs(remaining) < TOLERANCE) {

									if (s_lb < TOLERANCE && s_ub > -(TOLERANCE)) {

										double new_singleton_val = 0.0;

										// Update singleton
										double signed_delta = new_singleton_val - inst->x[singleton_index];
										inst->x[singleton_index] = new_singleton_val;

										// Update objective value
										inst->objval = inst->objval + (inst->obj[singleton_index] * signed_delta);
									}
									else print_error("[update_slacks][extension][row %d '%c'][singleton x_%d]: Failed to set singleton to ZERO.\n", 
											rowind + 1, inst->sense[rowind], singleton_index + 1);
								}
								// Remaining singletons slack to distribute
								else {

									double new_singleton_val = 0.0; // initialized

									// Distribute max possible amount of singletons slack
									if (coef < 0.0 && s_ub > -(TOLERANCE)) {

										new_singleton_val = min(s_ub, remaining / coef);
									}
									if (coef > 0.0 && s_lb < TOLERANCE) {

										new_singleton_val = max(s_lb, remaining / coef);
									}

									// Check bounds of new singleton value
									if (!(new_singleton_val > inst->lb[singleton_index] - TOLERANCE &&
										new_singleton_val < inst->ub[singleton_index] + TOLERANCE))
										print_error("[update_slacks][extension][singleton x_%d][row %d '%c']: Updated singleton out of bounds. Found %f <= %f <= %f\n",
											singleton_index + 1, rowind + 1, inst->sense[rowind], inst->lb[singleton_index], new_singleton_val, inst->ub[singleton_index]);

									// Update remaining singletons slack to distribute
									remaining = remaining - (coef * new_singleton_val);
									print_verbose(120, "[update_slacks][extension][row %d '%c']: Remaining singletons slack to distribute %f.\n",
										rowind + 1, inst->sense[rowind], remaining);

									// Update singleton
									double signed_delta = new_singleton_val - inst->x[singleton_index];
									inst->x[singleton_index] = new_singleton_val;

									// Update objective value
									inst->objval = inst->objval + (inst->obj[singleton_index] * signed_delta);
								}
							}
						}
						// Zero singletons slack to distribute
						else if (fabs(new_ss) < TOLERANCE) {

							// Try to set all the singletons to ZERO
							// Scan singletons
							for (int k = 0; k < inst->num_singletons[rowind]; k++) {

								// Singleton info
								int singleton_index = inst->row_singletons[beg + k];
								double s_lb = inst->lb[singleton_index];
								double s_ub = inst->ub[singleton_index];

								if (s_lb < TOLERANCE && s_ub > -(TOLERANCE)) {

									double new_singleton_val = 0.0;

									// Update singleton
									double signed_delta = new_singleton_val - inst->x[singleton_index];
									inst->x[singleton_index] = new_singleton_val;

									// Update objective value
									inst->objval = inst->objval + (inst->obj[singleton_index] * signed_delta);
								}
								else print_error("[update_slacks][extension][row %d '%c'][singleton x_%d]: Failed to set singleton to ZERO.\n",
									rowind + 1, inst->sense[rowind], singleton_index + 1);
							}
						}
					}
				}
				else {
					// Extension disabled
					print_error("[update_slacks]: Tried to update slack of an equality constraint with no singletons!\n");
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
	int colend;				/**< Index of the last non-zero coefficient of column j. */

	// Initialize
	delta_up1   = LONG_MAX;
	delta_down1 = LONG_MAX;
	delta_up2   = inst->ub[j] - inst->x[j];
	delta_down2 = inst->x[j] - inst->lb[j];
	colend      = (j < inst->ncols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;
	
	print_verbose(200, "[delta_updown]: delta_up2_%d = ub_%d - x_%d = %f - %f = %f ; delta_down2_%d = x_%d - lb_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, inst->ub[j], inst->x[j], delta_up2, j + 1, j + 1, j + 1, inst->x[j], inst->lb[j], delta_down2);

	// Scan constraints of variable xj
	for (int k = inst->cmatbeg[j]; k < colend; k++) {

		double aij = inst->cmatval[k];      // Coefficient of xj in the constraint
		int rowind = inst->cmatind[k];      // Constraint index
		double slack = inst->slack[rowind]; // Constraint row slack (no singleton slack included)

		// [EXTENSION] Get singletons slack info (if any)
		double ss_lb = LONG_MAX;
		double ss_ub = LONG_MIN;
		double singletons_slack = 0.0;
		double ss_delta_up = 0.0;
		double ss_delta_down = 0.0;
		if (inst->extension && inst->num_singletons[rowind] > 0) {

			// Compute singletons slack of constraint rowind and get bounds
			ss_lb = inst->ss_lb[rowind];
			ss_ub = inst->ss_ub[rowind];
			singletons_slack = compute_singletons_slack(inst, rowind);

			// Conpute singletons slack deltas (clip to zero if non-positive)
			ss_delta_up = ss_ub - singletons_slack;
			ss_delta_down = singletons_slack - ss_lb;
			if (ss_delta_up < 0.0) ss_delta_up = 0.0;
			if (ss_delta_down < 0.0) ss_delta_down = 0.0;
		}
		
		// Check sense, then check sign of aij, and update delta_up1, delta_down1
		switch (inst->sense[rowind]) {

			case 'L': // (slack non-negative)

				// [TEMPORARY]
				if (inst->extension && inst->num_singletons[rowind] > 0) print_error("[delta_updown][extension][x_%d][row %d '%c']: Inequality constraint has %d singletons.\n", j + 1, rowind + 1, inst->sense[rowind], inst->num_singletons[rowind]);

				// Clip slack to zero if negative
				if (inst->slack[rowind] < 0.0 - TOLERANCE) { 
					inst->slack[rowind] = 0.0; 
					slack = inst->slack[rowind]; 
				}

				// [EXTENSION] Update available slack: 'L' constraint --> singletons slack (if any) should decrease
				if (inst->extension && inst->num_singletons[rowind] > 0) slack = slack + ss_delta_down; // overall slack increases

				if (aij > 0.0) { 
					
					print_verbose(200, "[delta_updown]: sense = L ; slack[%d] = %f ; a_%d_%d = %f (--> candidate_up1)\n", rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1);

					// Update delta_up1
					candidate_up1 = slack / aij;
					delta_up1 = min(candidate_up1, delta_up1);
				}
				if (aij < 0.0) {

					print_verbose(200, "[delta_updown]: sense = L ; slack[%d] = %f ; a_%d_%d = %f (--> candidate_down1)\n", rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1);

					// Update delta_down1
					candidate_down1 = -(slack) / aij;
					delta_down1 = min(candidate_down1, delta_down1);
				}

				break;

			case 'G': // (slack non-positive)

				// [TEMPORARY]
				if (inst->extension && inst->num_singletons[rowind] > 0) print_error("[delta_updown][extension][x_%d][row %d '%c']: Inequality constraint has %d singletons.\n", j + 1, rowind + 1, inst->sense[rowind], inst->num_singletons[rowind]);

				// Clip slack to zero if positive
				if (inst->slack[rowind] > 0.0 + TOLERANCE) { 
					inst->slack[rowind] = 0.0;
					slack = inst->slack[rowind];
				}

				// [EXTENSION] Update available slack: 'G' constraint --> singletons slack (if any) should increase
				if (inst->extension && inst->num_singletons[rowind] > 0) slack = slack - ss_delta_up; // overall slack decreases (increases in absolute value)

				if (aij < 0.0) {

					print_verbose(200, "[delta_updown]: sense = G ; slack[%d] = %f ; a_%d_%d = %f (--> candidate_up1)\n", rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1);

					// Update delta_up1
					candidate_up1 = slack / aij;
					delta_up1 = min(candidate_up1, delta_up1);
				}
				if (aij > 0.0) {

					print_verbose(200, "[delta_updown]: sense = G ; slack[%d] = %f ; a_%d_%d = %f (--> candidate_down1)\n", rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1);

					// Update delta_down1
					candidate_down1 = -(slack) / aij;
					delta_down1 = min(candidate_down1, delta_down1);
				}

				break;

			case 'E': // (slack zero if extension disabled)

				// [EXTENSION] Distinguish equality constraints with singletons (if extension enabled)
				if (inst->extension && inst->num_singletons[rowind] > 0) {

					print_verbose(120, "[delta_updown][extension][x_%d][row %d '%c']: %d singletons. ", j + 1, rowind + 1, inst->sense[rowind], inst->num_singletons[rowind]);

					// Compute singletons slack of constraint rowind and get bounds (done above)
					print_verbose(120, "[delta_updown][extension][x_%d][row %d '%c']: Singletons slack = %f. Bounds %f <= ss <= %f\n", j + 1, rowind + 1, inst->sense[rowind], singletons_slack, ss_lb, ss_ub);
					// Compute singletons slack deltas (done above)

					// Update candidate deltas
					if (aij > 0.0) {

						candidate_down1 = ss_delta_up / aij;
						candidate_up1 = ss_delta_down / aij;
						// Clip candidates to zero if inaccurate
						if (fabs(candidate_down1) < TOLERANCE) candidate_down1 = 0.0;
						if (fabs(candidate_up1) < TOLERANCE) candidate_up1 = 0.0;

						// Check candidates signs
						if (candidate_down1 < -(TOLERANCE) || candidate_up1 < -(TOLERANCE))
							print_error("[delta_updown][extension][x_%d][row %d '%c']: Negative candidate deltas. Found %f and %f.", j + 1, rowind + 1, inst->sense[rowind], candidate_down1, candidate_up1);

						delta_down1 = min(candidate_down1, delta_down1);
						delta_up1 = min(candidate_up1, delta_up1);
					}
					if (aij < 0.0) {

						candidate_up1 = -(ss_delta_up) / aij;
						candidate_down1 = -(ss_delta_down) / aij;
						// Clip candidates to zero if inaccurate
						if (fabs(candidate_down1) < TOLERANCE) candidate_down1 = 0.0;
						if (fabs(candidate_up1) < TOLERANCE) candidate_up1 = 0.0;

						// Check candidates signs
						if (candidate_down1 < -(TOLERANCE) || candidate_up1 < -(TOLERANCE))
							print_error("[delta_updown][extension][x_%d][row %d '%c']: Negative candidate deltas. Found %f and %f.", j + 1, rowind + 1, inst->sense[rowind], candidate_down1, candidate_up1);

						delta_down1 = min(candidate_down1, delta_down1);
						delta_up1 = min(candidate_up1, delta_up1);
					}
				}
				else {
					// Extension disabled
					print_verbose(200, "[delta_updown][x_%d][row %d '%c']: Slack ZERO (no singletons) --> x_%d cannot be moved!\n", j + 1, rowind + 1, inst->sense[rowind], j + 1);

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
	if (new_delta_up < (epsilon - TOLERANCE) && new_delta_down < (epsilon - TOLERANCE)) {
		new_delta_up = 0.0;
		new_delta_down = 0.0;
	}
	delta_up[j] = new_delta_up;
	delta_down[j] = new_delta_down;
}

// [EXTENSION]
double compute_singletons_slack(instance* inst, int rowind) {

	if (inst->num_singletons[rowind] <= 0) print_error("[compute_singletons_slack][extension]: Tried to compute singletons slack of row %d with no singletons.\n", rowind + 1);

	// Compute singletons slack
	double singletons_slack = 0.0;
	int beg = inst->rs_beg[rowind];
	for (int k = 0; k < inst->num_singletons[rowind]; k++) {

		int singleton_index = inst->row_singletons[beg + k];
		double coef = inst->rs_coef[beg + k];
		singletons_slack += (coef * inst->x[singleton_index]);
	}

	return singletons_slack;
}