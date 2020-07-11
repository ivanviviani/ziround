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

	double* delta_up; 	 /**< Maximum variable up-shifts. */
	double* delta_down;  /**< Maximum variable down-shifts. */
	double ZI; 			 /**< Fractionality of a variable (used in function zi_round). */
	double ZIplus; 		 /**< Fractionality of a shifted up variable (used in function zi_round). */
	double ZIminus; 	 /**< Fractionality of a shifted down variable (used in fucntion zi_round). */
	double obj_plusUBj;  /**< Objective value for a shifted up variable (used in function zi_round). */
	double obj_minusLBj; /**< Objective value for a shifted down variable (used in function zi_round). */
	int updated;	     /**< Flag set to 1 when at least one variable shift has been made. */
	double curr_frac;	 /**< Current solution fractionality. */
	double frac[2];      /**< Circular buffer for current solution fractionality. */
	int bufind = 0;      /**< Current index in the circular buffer. */

	// Allocate / Initialize
	delta_up   = (double*)malloc(inst->ncols * sizeof(double));
	delta_down = (double*)malloc(inst->ncols * sizeof(double)); if (delta_up == NULL || delta_down == NULL) print_error("[zi_round]: Failed to allocate delta arrays.\n");
	frac[0] = 0.0;
	frac[1] = 0.0;

	// Plotting variables ---------------------------------------------------------------------------------------------------
	int size_frac;            /**< Actual current size of the solution fractionality tracker array. */
	int size_cost;            /**< Actual current size of the solution cost tracker array. */
	int len_frac;             /**< Maximum length of the solution fractionality tracker array (resizable). */
	int len_cost;             /**< Maximum length of the solution cost tracker array (resizable). */
	double* tracker_sol_frac; /**< Tracker of solution fractionality. */
	double* tracker_sol_cost; /**< Tracker of solution cost. */
	if (VERBOSE >= 10) {
		size_frac = 0;
		size_cost = 0;
		len_frac = 10;
		len_cost = 10;
		tracker_sol_frac = (double*)calloc(len_frac, sizeof(double));
		tracker_sol_cost = (double*)calloc(len_cost, sizeof(double)); 
		if (tracker_sol_frac == NULL || tracker_sol_cost == NULL) print_error("[ziround]: Failed to allocate trackers.\n");
	}
	// ----------------------------------------------------------------------------------------------------------------------

	// [DEBUG ONLY]: Print solution fractionality and cost ------------------------------------------------------------------
	print_verbose(10, "* Sol.fract. | Objval *\n");
	curr_frac = sol_fractionality(inst->x, inst->int_var, inst->ncols);
	frac[bufind] = curr_frac;
	if (fabs(frac[bufind] - frac[!bufind]) > TOLERANCE) print_verbose(10, "* %f | %f *\n", frac[bufind], inst->objval);
	if (VERBOSE >= 10) {
		add_point_single_tracker(curr_frac, &tracker_sol_frac, &len_frac, &size_frac);
		add_point_single_tracker(inst->objval, &tracker_sol_cost, &len_cost, &size_cost);
	}
	bufind = !bufind;
	// ----------------------------------------------------------------------------------------------------------------------

	// Outer loop (repeat until no more updates found)
	do {
		updated = 0;

		// Inner loop (for each variable xj that was integer/binary in the original MIP)
		for (int j = 0; j < inst->ncols; j++) {

			// Skip non-integer variables
			if (!(inst->int_var[j])) continue;
			assert(var_type_integer_or_binary(inst->vartype[j]));

			switch (is_fractional(inst->x[j])) {

				// xj non fractional
				case 0:

					print_verbose(150, "[zi_round]: -x- x_%d is non-fractional\n", j + 1);

					// Calculate deltas (with epsilon = 1.0)
					delta_updown(inst, j, delta_up, delta_down, 1.0);
					assert(
						var_in_bounds(inst->x[j] + delta_up[j], inst->lb[j], inst->ub[j]) & 
						var_in_bounds(inst->x[j] - delta_down[j], inst->lb[j], inst->ub[j])
					);

					// Skip xj if both deltas are equal to zero (no shift necessary)
					if ((fabs(delta_up[j]) < TOLERANCE) && (fabs(delta_down[j]) < TOLERANCE)) continue;

					// Condition(s) for rounding of xj (>= to include the case of a zero obj coefficient)
					if ((inst->obj[j] >= 0 && fabs(delta_down[j] - 1.0) < TOLERANCE) ||
						(inst->obj[j] <= 0 && fabs(delta_up[j] - 1.0) < TOLERANCE)) {

						// Round xj to improve objective and update slacks
						updated = updated | round_xj_bestobj(inst, j, delta_up, delta_down, 0); // flag xj non-fractional (0)
					}

					break;

				// xj fractional
				case 1:

					print_verbose(150, "[zi_round]: -x- x_%d is fractional\n", j + 1);

					// Calculate deltas
					delta_updown(inst, j, delta_up, delta_down, EPSILON);
					assert(
						var_in_bounds(inst->x[j] + delta_up[j], inst->lb[j], inst->ub[j]) & 
						var_in_bounds(inst->x[j] - delta_down[j], inst->lb[j], inst->ub[j])
					);

					// Skip xj if both deltas are equal to zero (no shift necessary)
					if ((fabs(delta_up[j]) < TOLERANCE) && (fabs(delta_down[j]) < TOLERANCE)) continue;

					ZI      = fractionality(inst->x[j]);
					ZIplus  = fractionality(inst->x[j] + delta_up[j]);
					ZIminus = fractionality(inst->x[j] - delta_down[j]);

					// First case: ZIplus = ZIminus && both < ZI --> Round to improve objective
					if ((fabs(ZIplus - ZIminus) < TOLERANCE) && (ZIplus < ZI - TOLERANCE)) {

						// Round xj to improve objective and update slacks
						updated = updated | round_xj_bestobj(inst, j, delta_up, delta_down, 1); // flag xj fractional (1)
					}

					// Second case: ZIplus < ZIminus && ZIplus < ZI --> Round UP
					else if ((ZIplus < ZIminus - TOLERANCE) && (ZIplus < ZI - TOLERANCE)) {

						// Skip variable if delta_up = 0
						if (fabs(delta_up[j]) < TOLERANCE) continue;

						print_verbose(100, "[zi_round][FRA]        : Round UP >>> Set x_%d = x_%d + delta_up_%d = %f + %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);

						// Check whether all affected constraints have enough slack for a ROUND UP of xj
						check_slacks(inst, j, delta_up, delta_down, 'U');
						
						// Round UP
						inst->x[j] = inst->x[j] + delta_up[j];

						updated = 1;
						update_slacks(inst, j, delta_up[j]);
						inst->objval += (inst->obj[j] * delta_up[j]);
					}

					// Third case: ZIminus < ZIplus && ZIminus < ZI --> Round DOWN
					else if ((ZIminus < ZIplus - TOLERANCE) && (ZIminus < ZI - TOLERANCE)) {

						// Skip variable if delta_down = 0
						if (fabs(delta_down[j]) < TOLERANCE) continue;
						
						print_verbose(100, "[zi_round][FRA]        : Round DOWN >>> Set x_%d = x_%d - delta_down_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);

						// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
						check_slacks(inst, j, delta_up, delta_down, 'D');
						
						// Round DOWN
						inst->x[j] = inst->x[j] - delta_down[j];

						updated = 1;
						update_slacks(inst, j, -(delta_down[j]));
						inst->objval -= (inst->obj[j] * delta_down[j]);
					}

					break;

				default:
					print_error(" in function is_fractional.\n");
			}

			// [DEBUG ONLY]: Print solution fractionality and cost ---------------------------------------------------------------
			curr_frac = sol_fractionality(inst->x, inst->int_var, inst->ncols);
			frac[bufind] = curr_frac;
			if (fabs(frac[bufind] - frac[!bufind]) > TOLERANCE) print_verbose(10, "* %f | %f *\n", frac[bufind], inst->objval);
			if (VERBOSE >= 10) {
				add_point_single_tracker(curr_frac, &tracker_sol_frac, &len_frac, &size_frac);
				add_point_single_tracker(inst->objval, &tracker_sol_cost, &len_cost, &size_cost);
			}
			bufind = !bufind;
			// -------------------------------------------------------------------------------------------------------------------

		} // end inner loop

		if (updated) { print_verbose(10, "[zi_round]: ... Update found, scan variables again ...\n"); }
		else { print_verbose(10, "[zi_round]: ... No updates found, exit outer loop ...\n"); }

		// [DEBUG ONLY] Pause after each inner loop execution
		system("pause");
		
		// [DEBUG ONLY] (BRUTE FORCE)  Check variable bounds and constraints
		if (VERBOSE >= 1) {
			check_bounds(inst, inst->x);
			check_constraints(inst, inst->x);
		}

	} while (updated); // end outer loop

	// [DEBUG ONLY]: Plot solution fractionality and cost trackers -------------------------------------------------
	if (VERBOSE >= 10) {
		char** labels_sol_frac = (char**)calloc(2, sizeof(char*));
		labels_sol_frac[0] = (char*)calloc(20, sizeof(char));
		labels_sol_frac[1] = (char*)calloc(20, sizeof(char));
		char** labels_sol_cost = (char**)calloc(2, sizeof(char*));
		labels_sol_cost[0] = (char*)calloc(20, sizeof(char));
		labels_sol_cost[1] = (char*)calloc(20, sizeof(char));
		char* name = (char*)calloc(20, sizeof(char));
		sprintf(labels_sol_frac[0], "Round");
		sprintf(labels_sol_frac[1], "Fractionality");
		sprintf(labels_sol_cost[0], "Round");
		sprintf(labels_sol_cost[1], "Cost");
		sprintf(name, "Solution");
		plot_tracker(tracker_sol_frac, name, labels_sol_frac, size_frac, NULL);
		plot_tracker(tracker_sol_cost, name, labels_sol_cost, size_cost, NULL);
		free_all(9, tracker_sol_frac, tracker_sol_cost, labels_sol_frac[0], labels_sol_frac[1],
				labels_sol_cost[0], labels_sol_cost[1], labels_sol_frac, labels_sol_cost, name);
	}
	// -------------------------------------------------------------------------------------------------------------

	// Free
	free(delta_up);
	free(delta_down);
}
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************

void check_slacks(instance* inst, int j, double* delta_up, double* delta_down, const char round_updown) {

	if (round_updown != 'U' && round_updown != 'D') print_error("[check_slacks]: Rounding sense '%c' undefined.\n", round_updown);

	int colend;              /**< Index of the last constraint containing variable \p j. */
	int rowind;              /**< Current row index. */
	double aij;              /**< Current constraint coefficient of variable \p j. */
	double curr_slack;       /**< Slack of the current constraint. */
	double singletons_slack; /**< Singletons slack of the current constraint. */
	double delta_slack;      /**< Delta slack of the current constraint (to be distributed). */
	double new_slack;        /**< Current row slack after rounding. */
	int enough_slack;        /**< Support flag. */
	double ss_lb;            /**< Lower bound of current singletons slack (in its row). */
	double ss_ub;            /**< Upper bound of current singletons slack (in its row). */
	double delta_ss;         /**< Delta singletons slack after rounding. */
	double new_ss;           /**< Singletons slack after rounding. */

	colend = (j < inst->ncols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;

	// Check whether all affected constraints have enough slack for a ROUND UP/DOWN of xj
	for (int k = inst->cmatbeg[j]; k < colend; k++) {

		rowind = inst->cmatind[k];
		assert(index_in_bounds(rowind, inst->nrows));
		aij = inst->cmatval[k];
		curr_slack = 0.0;
		singletons_slack = 0.0;
		delta_slack = 0.0;
		new_slack = 0.0;
		enough_slack = 0;

		switch (inst->sense[rowind]) {

			case 'L': // (slack non-negative)

				curr_slack = inst->slack[rowind];
				assert(non_negative_double(curr_slack));
				delta_slack = (round_updown == 'U') ? (aij * delta_up[j]) : (aij * (-delta_down[j]));

				// Row slack after rounding (negative iff also singletons slack should be used)
				new_slack = curr_slack - delta_slack;
				print_verbose(200, "[check_slacks][x_%d aij %f][row %d '%c']: new_slack = %f ", j + 1, aij, rowind + 1, inst->sense[rowind], new_slack);
				if (new_slack < -(TOLERANCE)) print_verbose(200, "(need singletons slack!)\n");
				else print_verbose(200, "\n");

				// [EXTENSION] Distinguish inequality constraints with singletons (if extension enabled)
				if ((inst->extension) && (inst->num_singletons[rowind] > 0)) {

					// If the new row slack is negative then the remaining amount must be covered by the singletons slack
					if (new_slack < -(TOLERANCE)) {

						// Compute singletons slack of constraint rowind and get bounds
						ss_lb = inst->ss_lb[rowind];
						ss_ub = inst->ss_ub[rowind];
						singletons_slack = compute_singletons_slack(inst, rowind);
						assert(var_in_bounds(singletons_slack, ss_lb, ss_ub));

						// Singletons slack after rounding
						delta_ss = new_slack;
						new_ss = singletons_slack + delta_ss; // + because signed delta

						// New singletons slack must stay within its bounds
						enough_slack = !((new_ss < ss_lb - TOLERANCE) || (new_ss > ss_ub + TOLERANCE));
					}
					else {
						// Extension enabled, row has singletons, but new_slack is non-negative --> no need to use singletons
						enough_slack = 1;
					}
				}
				else {
					// Extension disabled OR enabled but no singletons
					enough_slack = (new_slack > -(TOLERANCE));					
				}

				if (!enough_slack) print_error("[check_slacks][x_%d][row %d '%c']: After rounding, invalid slack.\n", j + 1, rowind + 1, inst->sense[rowind]);

				break;

			case 'G': // (slack non-positive)
				
				curr_slack = inst->slack[rowind];
				assert(non_positive_double(curr_slack));
				delta_slack = (round_updown == 'U') ? (aij * delta_up[j]) : (aij * (-delta_down[j]));

				// Row slack after rounding (positive iff also singletons slack should be used)
				new_slack = curr_slack - delta_slack;
				print_verbose(200, "[check_slacks][x_%d aij %f][row %d '%c']: new_slack = %f ", j + 1, aij, rowind + 1, inst->sense[rowind], new_slack);
				if (new_slack > TOLERANCE) print_verbose(200, "(need singletons slack!)\n");
				else print_verbose(200, "\n");

				// [EXTENSION] Distinguish inequality constraints with singletons (if extension enabled)
				if ((inst->extension) && (inst->num_singletons[rowind] > 0)) {

					// If the new row slack is positive then the remaining amount must be covered by the singletons slack
					if (new_slack > TOLERANCE) {

						// Compute singletons slack of constraint rowind and get bounds
						ss_lb = inst->ss_lb[rowind];
						ss_ub = inst->ss_ub[rowind];
						singletons_slack = compute_singletons_slack(inst, rowind);
						assert(var_in_bounds(singletons_slack, ss_lb, ss_ub));

						// Singletons slack after rounding
						delta_ss = new_slack;
						new_ss = singletons_slack + delta_ss; // + because signed delta

						// New singletons slack must stay within its bounds
						enough_slack = !((new_ss < ss_lb - TOLERANCE) || (new_ss > ss_ub + TOLERANCE));
					}
					else {
						// Extension enabled, row has singletons, but new_slack is non-positive --> no need to use singletons
						enough_slack = 1;
					}
				}
				else {
					// Extension disabled OR enabled but no singletons
					enough_slack = (new_slack < TOLERANCE);
				}

				if (!enough_slack) print_error("[check_slacks][x_%d][row %d '%c']: After rounding, invalid slack.\n", j + 1, rowind + 1, inst->sense[rowind]);

				break;

			case 'E': // (slack zero if extension disabled)

				// [EXTENSION] Distinguish equality constraints with singletons (if extension enabled)
				if (inst->extension && inst->num_singletons[rowind] > 0) {

					// Compute singletons slack of constraint rowind (with bounds)
					ss_lb = inst->ss_lb[rowind];
					ss_ub = inst->ss_ub[rowind];
					curr_slack = compute_singletons_slack(inst, rowind);
					assert(var_in_bounds(curr_slack, ss_lb, ss_ub));
					delta_slack = (round_updown == 'U') ? (aij * delta_up[j]) : (aij * (-delta_down[j]));

					// Singletons slack after rounding
					new_ss = curr_slack - delta_slack;

					// New singletons slack must stay within its bounds
					enough_slack = !((new_ss < ss_lb - TOLERANCE) || (new_ss > ss_ub + TOLERANCE));
					
					if (!enough_slack) print_error("[check_slacks][extension][x_%d][row %d '%c']: After rounding, singletons slack out of bounds. Found %f <= %f <= %f.\n", j + 1, rowind + 1, inst->sense[rowind], ss_lb, new_ss, ss_ub);
				}
				else {
					// Extension disabled OR enabled but no singletons
					print_error("[check_slacks][x_%d][row %d '%c']: Extension disabled OR constraint has no singletons --> slack ZERO --> x_%d cannot be rounded.\n", j + 1, rowind + 1, inst->sense[rowind], j + 1);
				}

				break;

			default:
				print_error("[check_slacks]: Constraint sense '%c' not included in {'L','G','E'}.\n", inst->sense[rowind]);
		} // end switch
	} // end for
}

int round_xj_bestobj(instance* inst, int j, double* delta_up, double* delta_down, int xj_fractional) {

	double obj_deltaplus = 0.0;  /**< Delta obj if xj is shifted up. */
	double obj_deltaminus = 0.0; /**< Delta obj if xj is shifted down. */
	int updated = 0;             /**< Flag set to 1 when at least one variable shift has been made. */

	if (fabs(inst->obj[j]) >= TOLERANCE) {
		obj_deltaplus = (inst->obj[j] * delta_up[j]);
		obj_deltaminus = -(inst->obj[j] * delta_down[j]);
	}

	print_verbose(201, "[round_xj_bestobj]: %f <= x_%d = %f <= %f | delta_up = %f delta down = %f objcoef = %f | obj_deltaplus = %f obj_deltaminus = %f\n", inst->lb[j], j + 1, inst->x[j], inst->ub[j], delta_up[j], delta_down[j], inst->obj[j], obj_deltaplus, obj_deltaminus);

	// Check obj sense, then update xj, update slacks and update objective value
	switch (inst->objsen) {

		case CPX_MIN:

			// [] Adding delta_up to x_j improves objval
			if ((obj_deltaplus < -(TOLERANCE)) && (obj_deltaplus < obj_deltaminus - TOLERANCE)) {

				// Skip variable if delta_up = 0
				if (fabs(delta_up[j]) < TOLERANCE) return 0;

				// print_verbose(100, "[round_xj_bestobj][frac?%d]: >>> Set x_%d = x_%d + delta_up_%d = %f + %f = %f\n", xj_fractional, j + 1, j + 1, j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
				print_verbose(10, "[round_xj_bestobj]: >>> Round x_%d = %f + %f = %f\n", j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
				if ((!(xj_fractional)) && (VERBOSE >= 150) && (fabs(delta_up[j] - 1.0) > TOLERANCE)) {
					print_error("[round_xj_bestobj]: delta_up_%d = %f (should be 1.0).\n", j + 1, delta_up[j]);
				}

				// Check whether all affected constraints have enough slack for a ROUND UP of xj
				check_slacks(inst, j, delta_up, delta_down, 'U');

				// Round UP (if xj is not fractional then delta_up[j] must be 1.0)
				inst->x[j] += delta_up[j];

				updated = 1;
				update_slacks(inst, j, delta_up[j]);
				inst->objval += obj_deltaplus;
			}
			// [] Adding -delta_down to x_j improves objval
			else if ((obj_deltaminus < -(TOLERANCE)) && (obj_deltaminus < obj_deltaplus - TOLERANCE)) {

				// Skip variable if delta_down = 0
				if (fabs(delta_down[j]) < TOLERANCE) return 0;

				// print_verbose(100, "[round_xj_bestobj][frac?%d]: >>> Set x_%d = x_%d - delta_down_%d = %f - %f = %f\n", xj_fractional, j + 1, j + 1, j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);
				print_verbose(10, "[round_xj_bestobj]: >>> Round x_%d = %f - %f = %f\n", j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);
				if ((!(xj_fractional)) && (VERBOSE >= 150) && (fabs(delta_down[j] - 1.0) > TOLERANCE)) {
					print_error("[round_xj_bestobj]: delta_down_%d = %f (should be 1.0).\n", j + 1, delta_down[j]);
				}

				// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
				check_slacks(inst, j, delta_up, delta_down, 'D');

				// Round DOWN (if xj is not fractional then delta_down[j] must be 1.0)
				inst->x[j] -= delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval += obj_deltaminus;
			}
			// [] Both deltas improve objval of the same amount < 0 --> Round arbitrarily (DOWN)
			else if ((obj_deltaminus < -(TOLERANCE)) && (obj_deltaplus < -(TOLERANCE)) && (fabs(obj_deltaminus - obj_deltaplus) < TOLERANCE)) {

				// Skip variable if delta_down = 0
				if (fabs(delta_down[j]) < TOLERANCE) return 0;

				// print_verbose(120, "[round_xj_bestobj][frac?%d]: obj_deltaplus = %f = %f = obj_deltaminus. >>> Round x_%d arbitrarily.\n", xj_fractional, obj_deltaplus, obj_deltaminus, j + 1);
				print_verbose(10, "[round_xj_bestobj]: >>> Round x_%d = %f - %f = %f\n", j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);

				// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
				check_slacks(inst, j, delta_up, delta_down, 'D');

				// Round arbitrarily (DOWN)
				inst->x[j] -= delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval += obj_deltaminus;
			}
			// [] Both deltas do not change objval (both = 0) --> Round arbitrarily (UP)
			else if ((fabs(obj_deltaplus) < TOLERANCE) && (fabs(obj_deltaplus - obj_deltaminus) < TOLERANCE)) {

				// Skip variable if delta_up = 0
				if (fabs(delta_up[j]) < TOLERANCE) return 0;

				// print_verbose(120, "[round_xj_bestobj][frac?%d]: obj_deltaplus = %f = %f = obj_deltaminus. >>> Round x_%d arbitrarily.\n", xj_fractional, obj_deltaplus, obj_deltaminus, j + 1);
				print_verbose(10, "[round_xj_bestobj]: >>> Round x_%d = %f + %f = %f\n", j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
				
				// Check whether all affected constraints have enough slack for a ROUND UP of xj
				check_slacks(inst, j, delta_up, delta_down, 'U');

				// Round arbitrarily (UP)
				inst->x[j] += delta_up[j];
				updated = 1;
				update_slacks(inst, j, delta_up[j]);
				// no objval update
			}

			break;

		// (problems from mps files should always be MIN)
		case CPX_MAX:

			// [] Adding delta_up to x_j improves objval
			if ((obj_deltaplus > TOLERANCE) && (obj_deltaplus > obj_deltaminus + TOLERANCE)) {
				
				// print_verbose(100, "[round_xj_bestobj][frac?%d]: >>> Set x_%d = x_%d + delta_up_%d = %f + %f = %f\n", xj_fractional, j + 1, j + 1, j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
				print_verbose(10, "[round_xj_bestobj]: >>> Round x_%d = %f + %f = %f\n", j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
				if ((!(xj_fractional)) && (VERBOSE >= 150) && (fabs(delta_up[j] - 1.0) > TOLERANCE)) {
					print_error("[round_xj_bestobj]: delta_up_%d = %f (should be 1.0).\n", j + 1, delta_up[j]);
				}

				// Check whether all affected constraints have enough slack for a ROUND UP of xj
				check_slacks(inst, j, delta_up, delta_down, 'U');

				// Round UP (if xj is not fractional then delta_up[j] must be 1.0)
				inst->x[j] += delta_up[j];

				updated = 1;
				update_slacks(inst, j, delta_up[j]);
				inst->objval += obj_deltaplus;
			}
			// [] Adding delta_down to x_j improves objval
			else if ((obj_deltaminus > TOLERANCE) && (obj_deltaminus > obj_deltaplus + TOLERANCE)) {
				
				// print_verbose(100, "[round_xj_bestobj][frac?%d]: >>> Set x_%d = x_%d - LB_%d = %f - %f = %f\n", xj_fractional, j + 1, j + 1, j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);
				print_verbose(10, "[round_xj_bestobj]: >>> Round x_%d = %f - %f = %f\n", j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);
				if (!xj_fractional && VERBOSE >= 150 && fabs(delta_down[j] - 1.0) > TOLERANCE) {
					print_error("[round_xj_bestobj]: delta_down_%d = %f (should be 1.0).\n", j + 1, delta_down[j]);
				}

				// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
				check_slacks(inst, j, delta_up, delta_down, 'D');

				// Round DOWN (if xj is not fractional then delta_down[j] must be 1.0)
				inst->x[j] -= delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval += obj_deltaminus;
			}
			// [] Both deltas improve objval of the same amount > 0 --> Round arbitrarily (DOWN)
			else if ((obj_deltaminus > TOLERANCE) && (obj_deltaplus > TOLERANCE) && (fabs(obj_deltaminus - obj_deltaplus) < TOLERANCE)) {
				
				print_verbose(120, "[round_xj_bestobj]: obj_deltaplus = %f = %f = obj_deltaminus. >>> Round x_%d arbitrarily.\n", obj_deltaplus, obj_deltaminus, j + 1);
				print_verbose(10, "[round_xj_bestobj]: >>> Round x_%d = %f - %f = %f\n", j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);

				// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
				check_slacks(inst, j, delta_up, delta_down, 'D');

				// Round arbitrarily (DOWN)
				inst->x[j] -= delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval += obj_deltaminus;
			}
			// [] Both deltas do not change objval (both = 0) --> Round arbitrarily (UP)
			else if ((fabs(obj_deltaplus) < TOLERANCE) && (fabs(obj_deltaplus - obj_deltaminus) < TOLERANCE)) {

				// Skip variable if delta_up = 0
				if (fabs(delta_up[j]) < TOLERANCE) return 0;

				print_verbose(120, "[round_xj_bestobj][frac?%d]: obj_deltaplus = %f = %f = obj_deltaminus. >>> Round x_%d arbitrarily.\n", xj_fractional, obj_deltaplus, obj_deltaminus, j + 1);
				print_verbose(10, "[round_xj_bestobj]: >>> Round x_%d = %f + %f = %f\n", j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
				
				// Check whether all affected constraints have enough slack for a ROUND UP of xj
				check_slacks(inst, j, delta_up, delta_down, 'U');

				// Round arbitrarily (UP)
				inst->x[j] += delta_up[j];
				updated = 1;
				update_slacks(inst, j, delta_up[j]);
				// no objval update
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

	int colend;         /**< Index of the last constraint containing variable \p j. */
	double aij;         /**< Current constraint coefficient of variable \p j. */
	int rowind;         /**< Current row index. */
	double delta_slack; /**< Delta slack of the current constraint (to be distributed). */
	double curr_slack;  /**< Slack of the current constraint. */
	double temp_slack;  /**< Support variable for \p delta_slack distribution. */
	double delta_ss;    /**< Delta singletons slack of the current constraint (to be distributed). */

	colend = (j < inst->ncols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;

	// Scan constraints of variable j
	for (int k = inst->cmatbeg[j]; k < colend; k++) {

		aij = inst->cmatval[k];
		rowind = inst->cmatind[k];
		delta_slack = aij * signed_delta;
		curr_slack = inst->slack[rowind];
		temp_slack = 0.0;

		switch (inst->sense[rowind]) {

			case 'L': // (slack non-negative)
			case 'G': // (slack non-positive)

				// [EXTENSION] Distinguish inequality constraints with singletons (if extension enabled)
				if (inst->extension && inst->num_singletons[rowind] > 0) {

					// First, use at most all the row slack available to cover delta_slack
					temp_slack = curr_slack - delta_slack;
					// Update row slack
					inst->slack[rowind] = (inst->sense[rowind] == 'L') ? max(0.0, temp_slack) : min(0.0, temp_slack);

					// If not enough row slack (temp_slack negative for 'L', positive for 'G' constraints), resort to singletons slack
					if ((inst->sense[rowind] == 'L' && temp_slack < -(TOLERANCE)) || 
						(inst->sense[rowind] == 'R' && temp_slack > TOLERANCE)) {

						// Delta singletons slack to distribute among the singletons [new_ss = ss + delta_ss (+ because signed delta)]
						delta_ss = temp_slack; // negative for 'L', positive for 'G' constraints

						// Distribute delta among the singletons, stop when done (delta_ss negative --> singletons slack must decrease)
						update_singletons(inst, rowind, delta_ss);
					}
					else {
						// Enough row slack, update it
						assert(
							(inst->sense[rowind] == 'L' && non_negative_double(inst->slack[rowind] - delta_slack)) ||
							(inst->sense[rowind] == 'G' && non_positive_double(inst->slack[rowind] - delta_slack))
						);
						inst->slack[rowind] -= delta_slack;
					}
				}
				else {
					// Extension disabled OR enabled but zero singletons
					// Just update row slack
					print_verbose(201, "[update_slacks][x_%d][row %d '%c']: slack = %f - (%f * %f) = %f\n", j + 1, rowind + 1, inst->sense[rowind], inst->slack[rowind], aij, signed_delta, inst->slack[rowind] - delta_slack);
					assert( 
						(inst->sense[rowind] == 'L' && non_negative_double(inst->slack[rowind] - delta_slack)) || 
						(inst->sense[rowind] == 'G' && non_positive_double(inst->slack[rowind] - delta_slack))
					);
					inst->slack[rowind] -= delta_slack;
				}

				break;

			case 'E':

				// [EXTENSION] Distinguish equality constraints with singletons (if extension enabled)
				if (inst->extension && inst->num_singletons[rowind] > 0) {

					// Equality constraint --> row slack is always zero
					
					// Delta singletons slack to distribute among the singletons (could be positive or negative)
					delta_ss = -(delta_slack);
					
					// Distribute delta among the singletons, stop when done
					update_singletons(inst, rowind, delta_ss);
				}
				else {
					// Extension disabled OR enabled but zero singletons
					print_error("[update_slacks]: Tried to update slack of an equality constraint with extension disabled or zero singletons!\n");
				}

				break;

			default:
				print_error("[update_slacks]: Constraint sense %c not supported!\n", inst->sense[rowind]);
		} // end switch
	} // end for
}

// [EXTENSION]
void update_singletons(instance* inst, int rowind, double delta_ss) {

	int beg = inst->rs_beg[rowind]; /**< Begin index of the singletons for constraint rowind. */
	int singleton_index;            /**< Current singleton index. */
	double coef;                    /**< Current singleton coefficient. */
	double s_lb;                    /**< Singleton lower bound. */
	double s_ub;                    /**< Singleton upper bound. */
	double s_val;                   /**< Current singleton value. */
	double covered_delta_ss;        /**< Delta covered by the current singleton. */
	double max_s_delta;             /**< Maximum delta coverable by the current singleton. */
	double s_delta;                 /**< Delta of the current singleton (to be updated). */
	int s_slack_increase;           /**< Flag set to 1 iff singletons slack should increase, 0 otherwise. */

	s_slack_increase = (delta_ss >= 0.0) ? 1 : 0;

	// Distribute delta among the singletons, stop when done (delta_ss positive(negative) --> singletons slack must increase(decrease))
	for (int k = 0; k < inst->num_singletons[rowind]; k++) {

		// Stop updating the singletons when delta singletons slack has been covered
		if ((s_slack_increase && delta_ss < TOLERANCE) || 
			(!s_slack_increase && delta_ss > -(TOLERANCE))) {
			print_verbose(200, "[update_singletons][extension][row %d '%c']: delta_ss covered, found %f\n", rowind + 1, inst->sense[rowind], delta_ss);
			break;
		}
		assert(
			(s_slack_increase && non_negative_double(delta_ss)) || 
			(!s_slack_increase && non_positive_double(delta_ss))
		);
		print_verbose(120, "[update_slacks][extension][row %d '%c']: Remaining delta singletons slack to distribute: %f.\n", rowind + 1, inst->sense[rowind], delta_ss);

		// Singleton info
		assert(index_in_bounds(beg + k, inst->rs_size));
		singleton_index = inst->row_singletons[beg + k];
		assert(index_in_bounds(singleton_index, inst->ncols));
		coef = inst->rs_coef[beg + k];
		s_lb = inst->lb[singleton_index];
		s_ub = inst->ub[singleton_index];
		s_val = inst->x[singleton_index];
		assert(var_in_bounds(s_val, s_lb, s_ub));
		covered_delta_ss = 0.0;
		max_s_delta = 0.0;
		s_delta = 0.0;

		// Compute covered delta of the singleton
		if (coef > 0.0) {
			if (s_slack_increase) {
				// Singletons slack increase
				max_s_delta = s_ub - s_val;
				covered_delta_ss = min(delta_ss, coef * max_s_delta);
			}
			else {
				// Singletons slack decrease
				max_s_delta = s_val - s_lb;
				covered_delta_ss = max(delta_ss, -coef * max_s_delta);
			}
		}
		if (coef < 0.0) {
			if (s_slack_increase) {
				// Singletons slack increase
				max_s_delta = s_val - s_lb;
				covered_delta_ss = min(delta_ss, -coef * max_s_delta);
			}
			else {
				// Singletons slack decrease
				max_s_delta = s_ub - s_val;
				covered_delta_ss = max(delta_ss, coef * max_s_delta);
			}
		}
		// Update remaining delta to be covered by the next singletons
		delta_ss -= covered_delta_ss;

		// Compute singleton delta
		s_delta = covered_delta_ss / coef;
		// Update singleton
		assert(var_in_bounds(s_val + s_delta, s_lb, s_ub));
		inst->x[singleton_index] = s_val + s_delta;

		// Update objective value
		inst->objval += (inst->obj[singleton_index] * s_delta);
	} // end for

	// Delta slack must have been distributed among the singletons
	assert(zero_double(delta_ss));
	print_verbose(120, "[update_singletons][extension][row %d '%c']: delta_ss distributed, remaining %f\n", rowind + 1, inst->sense[rowind], delta_ss);
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

	double delta_up1;		 /**< First delta_up[j] major candidate. */
	double delta_down1;		 /**< First delta_down[j] major candidate. */
	double delta_up2;		 /**< Second delta_up[j] major candidate. */
	double delta_down2;		 /**< Second delta_down[j] major candidate. */
	double candidate_up1;	 /**< Current delta_up[j] minor candidate. */
	double candidate_down1;  /**< Current delta_down[j] minor candidate. */
	double new_delta_up;     /**< Final delta_up[j] winner. */
	double new_delta_down;   /**< Final delta_down[j] winner. */
	int colend;				 /**< Index of the last constraint containing variable x_j. */
	double aij;              /**< Coefficient of xj in the constraint. */
	int rowind;              /**< Constraint index. */
	double slack;            /**< Row slack (no singleton slack included). */
	double ss_lb;            /**< Lower bound of singletons slack (in its row). */
	double ss_ub;            /**< Upper bound of singletons slack (in its row). */
	double singletons_slack; /**< Singletons slack value (in its row). */
	double ss_delta_up;      /**< Maximum delta up for current singletons slack. */
	double ss_delta_down;    /**< Maximum delta down for current singletons slack. */

	delta_up1   = LONG_MAX;
	delta_down1 = LONG_MAX;
	delta_up2   = inst->ub[j] - inst->x[j];
	delta_down2 = inst->x[j] - inst->lb[j];
	assert(
		non_negative_double(delta_up2) & 
		non_negative_double(delta_down2)
	);
	colend = (j < inst->ncols - 1) ? inst->cmatbeg[j + 1] : inst->nzcnt;
	
	print_verbose(201, "[delta_updown]: delta_up2_%d = ub_%d - x_%d = %f - %f = %f ; delta_down2_%d = x_%d - lb_%d = %f - %f = %f\n", j + 1, j + 1, j + 1, inst->ub[j], inst->x[j], delta_up2, j + 1, j + 1, j + 1, inst->x[j], inst->lb[j], delta_down2);

	// Scan constraints of variable xj
	for (int k = inst->cmatbeg[j]; k < colend; k++) {

		aij = inst->cmatval[k];
		rowind = inst->cmatind[k];
		assert(index_in_bounds(rowind, inst->nrows));
		slack = inst->slack[rowind]; // (no singleton slack included)

		// [EXTENSION] Get singletons slack info (if any)
		ss_lb = LONG_MAX;
		ss_ub = LONG_MIN;
		singletons_slack = 0.0;
		ss_delta_up = 0.0;
		ss_delta_down = 0.0;
		if (inst->extension && inst->num_singletons[rowind] > 0) {

			// Compute singletons slack of constraint rowind and get bounds
			ss_lb = inst->ss_lb[rowind];
			ss_ub = inst->ss_ub[rowind];
			singletons_slack = compute_singletons_slack(inst, rowind);
			assert(var_in_bounds(singletons_slack, ss_lb, ss_ub));

			// Conpute singletons slack deltas (clip to zero if slightly non-positive)
			ss_delta_up = ss_ub - singletons_slack;
			ss_delta_down = singletons_slack - ss_lb;
			if ((ss_delta_up < 0.0) && (ss_delta_up > -(TOLERANCE))) ss_delta_up = 0.0;
			if ((ss_delta_down < 0.0) && (ss_delta_down > -(TOLERANCE))) ss_delta_down = 0.0;
			assert(
				non_negative_double(ss_delta_up) & 
				non_negative_double(ss_delta_down)
			);
		}
		
		// Check sense, then check sign of aij, and update delta_up1, delta_down1
		switch (inst->sense[rowind]) {

			case 'L': // (slack non-negative)

				if (inst->slack[rowind] < -(TOLERANCE)) print_error("[delta_updown][row %d 'L']: Found negative row slack = %f\n", rowind + 1, inst->slack[rowind]);

				// Clip slack to zero if slightly negative
				if ((inst->slack[rowind] < 0.0) && (inst->slack[rowind] > -(TOLERANCE))) { 
					inst->slack[rowind] = 0.0; 
					slack = inst->slack[rowind]; 
				}
				assert(equals_double(slack, inst->slack[rowind]));

				// [EXTENSION] Update available slack: 'L' constraint --> singletons slack (if any) should decrease
				if (inst->extension && inst->num_singletons[rowind] > 0) slack = slack + ss_delta_down; // overall slack increases

				if (aij > 0.0) { 
					
					print_verbose(201, "[delta_updown]: sense = L ; slack[%d] = %f ; a_%d_%d = %f (--> candidate_up1)\n", rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1);

					// Update delta_up1
					candidate_up1 = slack / aij;
					assert(non_negative_double(candidate_up1));
					delta_up1 = min(candidate_up1, delta_up1);
					assert(non_negative_double(delta_up1));
				}
				if (aij < 0.0) {

					print_verbose(201, "[delta_updown]: sense = L ; slack[%d] = %f ; a_%d_%d = %f (--> candidate_down1)\n", rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1);

					// Update delta_down1
					candidate_down1 = -(slack) / aij;
					assert(non_negative_double(candidate_down1));
					delta_down1 = min(candidate_down1, delta_down1);
					assert(non_negative_double(delta_down1));
				}

				break;

			case 'G': // (slack non-positive)

				if (inst->slack[rowind] > TOLERANCE) print_error("[delta_updown][row %d 'G']: Found positive row slack = %f\n", rowind + 1, inst->slack[rowind]);

				// Clip slack to zero if slightly positive
				if ((inst->slack[rowind] > 0.0) && (inst->slack[rowind] < TOLERANCE)) { 
					inst->slack[rowind] = 0.0;
					slack = inst->slack[rowind];
				}
				assert(equals_double(slack, inst->slack[rowind]));

				// [EXTENSION] Update available slack: 'G' constraint --> singletons slack (if any) should increase
				if (inst->extension && inst->num_singletons[rowind] > 0) slack = slack - ss_delta_up; // overall slack decreases (increases in absolute value)

				if (aij < 0.0) {

					print_verbose(201, "[delta_updown]: sense = G ; slack[%d] = %f ; a_%d_%d = %f (--> candidate_up1)\n", rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1);

					// Update delta_up1
					candidate_up1 = slack / aij;
					assert(non_negative_double(candidate_up1));
					delta_up1 = min(candidate_up1, delta_up1);
					assert(non_negative_double(delta_up1));
				}
				if (aij > 0.0) {

					print_verbose(201, "[delta_updown]: sense = G ; slack[%d] = %f ; a_%d_%d = %f (--> candidate_down1)\n", rowind + 1, inst->slack[rowind], rowind + 1, j + 1, aij, j + 1);

					// Update delta_down1
					candidate_down1 = -(slack) / aij;
					assert(non_negative_double(candidate_down1));
					delta_down1 = min(candidate_down1, delta_down1);
					assert(non_negative_double(delta_down1));
				}

				break;

			case 'E': // (slack zero if extension disabled)

				// [EXTENSION] Distinguish equality constraints with singletons (if extension enabled)
				if (inst->extension && inst->num_singletons[rowind] > 0) {

					print_verbose(201, "[delta_updown][extension][x_%d][row %d '%c']: %d singletons. ", j + 1, rowind + 1, inst->sense[rowind], inst->num_singletons[rowind]);

					// Compute singletons slack of constraint rowind and get bounds (done above)
					print_verbose(201, "Singletons slack = %f. Bounds %f <= ss <= %f\n", singletons_slack, ss_lb, ss_ub);
					// Compute singletons slack deltas (done above)

					// Update candidate deltas
					if (aij > 0.0) {

						candidate_down1 = ss_delta_up / aij;
						candidate_up1 = ss_delta_down / aij;
						assert(
							non_negative_double(candidate_down1) & 
							non_negative_double(candidate_up1)
						);
						// Clip candidates to zero if inaccurate
						if (fabs(candidate_down1) < TOLERANCE) candidate_down1 = 0.0;
						if (fabs(candidate_up1) < TOLERANCE) candidate_up1 = 0.0;

						// Check candidates signs
						// if (candidate_down1 < -(TOLERANCE) || candidate_up1 < -(TOLERANCE)) print_error("[delta_updown][extension][x_%d][row %d '%c']: Negative candidate deltas. Found %f and %f.", j + 1, rowind + 1, inst->sense[rowind], candidate_down1, candidate_up1);

						delta_down1 = min(candidate_down1, delta_down1);
						delta_up1 = min(candidate_up1, delta_up1);
					}
					if (aij < 0.0) {

						candidate_up1 = -(ss_delta_up) / aij;
						candidate_down1 = -(ss_delta_down) / aij;
						assert(
							non_negative_double(candidate_up1) & 
							non_negative_double(candidate_down1)
						);
						// Clip candidates to zero if inaccurate
						if (fabs(candidate_down1) < TOLERANCE) candidate_down1 = 0.0;
						if (fabs(candidate_up1) < TOLERANCE) candidate_up1 = 0.0;

						// Check candidates signs
						// if (candidate_down1 < -(TOLERANCE) || candidate_up1 < -(TOLERANCE)) print_error("[delta_updown][extension][x_%d][row %d '%c']: Negative candidate deltas. Found %f and %f.", j + 1, rowind + 1, inst->sense[rowind], candidate_down1, candidate_up1);

						delta_down1 = min(candidate_down1, delta_down1);
						delta_up1 = min(candidate_up1, delta_up1);
					}
					assert(
						non_negative_double(delta_down1) & 
						non_negative_double(delta_up1)
					);
				}
				else {
					// Extension disabled OR enabled but zero singletons
					print_verbose(201, "[delta_updown][x_%d][row %d '%c']: Slack ZERO (no singletons) --> x_%d cannot be moved!\n", j + 1, rowind + 1, inst->sense[rowind], j + 1);

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
	assert(
		non_negative_double(delta_up1) & 
		non_negative_double(delta_down1)
	);
	new_delta_up = min(delta_up1, delta_up2);
	new_delta_down = min(delta_down1, delta_down2);
	assert(
		non_negative_double(new_delta_up) & 
		non_negative_double(new_delta_down)
	);
	print_verbose(201, "[delta_updown][results]: (NEW) delta_up_%d = min{%f, %f} = %f ; delta_down_%d = min{%f, %f} = %f\n", j + 1, delta_up1, delta_up2, new_delta_up, j + 1, delta_down1, delta_down2, new_delta_down);

	// Update deltas
	if (new_delta_up < (epsilon - TOLERANCE) && new_delta_down < (epsilon - TOLERANCE)) {
		new_delta_up = 0.0;
		new_delta_down = 0.0;
	}
	delta_up[j] = new_delta_up;
	delta_down[j] = new_delta_down;
	assert(
		equals_double(delta_up[j], new_delta_up) & 
		equals_double(delta_down[j], new_delta_down)
	);
}

// [EXTENSION]
double compute_singletons_slack(instance* inst, int rowind) {

	assert(index_in_bounds(rowind, inst->nrows));
	if (inst->num_singletons[rowind] <= 0) print_error("[compute_singletons_slack][extension]: Tried to compute singletons slack of row %d with no singletons.\n", rowind + 1);

	double singletons_slack; /**< Current singletons slack value. */
	int beg;                 /**< Begin index of singleton indices for row \p rowind. */
	int singleton_index;     /**< Current singleton index. */
	double coef;             /**< Current singleton coefficient in row \p rowind. */

	// Compute singletons slack
	singletons_slack = 0.0;
	beg = inst->rs_beg[rowind];
	assert(index_in_bounds(beg, inst->rs_size));
	for (int k = 0; k < inst->num_singletons[rowind]; k++) {

		assert(index_in_bounds(beg + k, inst->rs_size));
		singleton_index = inst->row_singletons[beg + k];
		assert(index_in_bounds(singleton_index, inst->ncols));
		coef = inst->rs_coef[beg + k];
		singletons_slack += (coef * inst->x[singleton_index]);
	}
	assert(var_in_bounds(singletons_slack, inst->ss_lb[rowind], inst->ss_ub[rowind]));

	return singletons_slack;
}