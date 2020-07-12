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
	double frac[2];      /**< Circular buffer for current solution fractionality. */
	double objval[2];    /**< Circular buffer for current objective value. */
	int bufind = 0;      /**< Current index in the circular buffer. */

	// Allocate / Initialize
	delta_up   = (double*)malloc(inst->ncols * sizeof(double));
	delta_down = (double*)malloc(inst->ncols * sizeof(double)); if (delta_up == NULL || delta_down == NULL) print_error("[zi_round]: Failed to allocate delta arrays.\n");
	frac[0] = 0.0;
	frac[1] = 0.0;

	// Plotting variables ---------------------------------------------------------------------------------------------------
	int size_frac;            /**< Actual current size of the solution fractionality tracker array. */
	int size_cost;            /**< Actual current size of the solution cost tracker array. */
	int size_rounded;         /**< Actual current size of the number of rounded variables array. */
	int len_frac;             /**< Maximum length of the solution fractionality tracker array (resizable). */
	int len_cost;             /**< Maximum length of the solution cost tracker array (resizable). */
	int len_rounded;          /**< Maximum length of the number of rounded variables array (resizable). */
	double* tracker_sol_frac; /**< Tracker of solution fractionality. */
	double* tracker_sol_cost; /**< Tracker of solution cost. */
	double* tracker_rounded;  /**< Tracker of number of rounded variables. */
	if (VERBOSE >= 10) {
		if (PLOT_SOL_FRAC) {
			size_frac = 0;
			len_frac = 10;
			tracker_sol_frac = (double*)calloc(len_frac, sizeof(double)); if (tracker_sol_frac == NULL) print_error("[ziround]: Failed to allocate solution fractionality tracker.\n");
		}
		if (PLOT_SOL_COST) {
			size_cost = 0;
			len_cost = 10;
			tracker_sol_cost = (double*)calloc(len_cost, sizeof(double)); if (tracker_sol_cost == NULL) print_error("[ziround]: Failed to allocate solution cost tracker.\n");
		}
		if (PLOT_NUM_ROUNDED_VARS) {
			size_rounded = 0;
			len_rounded = 10;
			tracker_rounded = (double*)calloc(len_rounded, sizeof(double)); if (tracker_rounded == NULL) print_error("[ziround]: Failed to allocate number of rounded variables tracker.\n");
		}
	}
	// ----------------------------------------------------------------------------------------------------------------------

	// [DEBUG ONLY]: Print solution fractionality and cost ------------------------------------------------------------------
	print_verbose(10, "* Sol.fract. | Objval *\n");
	frac[bufind] = sol_fractionality(inst->x, inst->int_var, inst->ncols);
	objval[bufind] = inst->objval;
	if (not_equals(frac[bufind], frac[!bufind]) || not_equals(objval[bufind], objval[!bufind])) print_verbose(10, "* %f | %f *\n", frac[bufind], objval[bufind]);
	if (VERBOSE >= 10) {
		if (PLOT_SOL_FRAC) add_point_single_tracker(frac[bufind], &tracker_sol_frac, &len_frac, &size_frac);
		if (PLOT_SOL_COST) add_point_single_tracker(objval[bufind], &tracker_sol_cost, &len_cost, &size_cost);
		if (PLOT_NUM_ROUNDED_VARS) add_point_single_tracker(count_rounded(inst->x, inst->ncols, inst->int_var, inst->vartype), &tracker_rounded, &len_rounded, &size_rounded);
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

					// Calculate deltas (with epsilon = 1.0)
					delta_updown(inst, j, delta_up, delta_down, 1.0);
					assert(
						var_in_bounds(inst->x[j] + delta_up[j], inst->lb[j], inst->ub[j]) & 
						var_in_bounds(inst->x[j] - delta_down[j], inst->lb[j], inst->ub[j])
					);

					// Skip xj if both deltas are equal to zero (no shift necessary)
					if (zero(delta_up[j]) && zero(delta_down[j])) continue;

					// Condition(s) for rounding of xj (>= to include the case of a zero obj coefficient)
					if ((inst->obj[j] >= 0 && equals(delta_down[j], 1.0)) ||
						(inst->obj[j] <= 0 && equals(delta_up[j], 1.0))) {

						// Round xj to improve objective and update slacks
						updated = updated | round_xj_bestobj(inst, j, delta_up, delta_down, 0); // flag xj non-fractional (0)
					}

					break;

				// xj fractional
				case 1:

					// Calculate deltas
					delta_updown(inst, j, delta_up, delta_down, EPSILON);
					assert(
						var_in_bounds(inst->x[j] + delta_up[j], inst->lb[j], inst->ub[j]) & 
						var_in_bounds(inst->x[j] - delta_down[j], inst->lb[j], inst->ub[j])
					);

					// Skip xj if both deltas are equal to zero (no shift necessary)
					if (zero(delta_up[j]) && zero(delta_down[j])) continue;

					ZI      = fractionality(inst->x[j]);
					ZIplus  = fractionality(inst->x[j] + delta_up[j]);
					ZIminus = fractionality(inst->x[j] - delta_down[j]);

					// First case: ZIplus = ZIminus && both < ZI --> Round to improve objective
					if (equals(ZIplus, ZIminus) && less_than(ZIplus, ZI)) {

						// Round xj to improve objective and update slacks
						updated = updated | round_xj_bestobj(inst, j, delta_up, delta_down, 1); // flag xj fractional (1)
					}

					// Second case: ZIplus < ZIminus && ZIplus < ZI --> Round UP
					else if (less_than(ZIplus, ZIminus) && less_than(ZIplus, ZI)) {

						// Skip variable if delta_up = 0
						if (zero(delta_up[j])) continue;

						print_verbose(20, "[ziround]: >>> Round x_%d = %f + %f = %f\n", j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);

						// Check whether all affected constraints have enough slack for a ROUND UP of xj
						check_slacks(inst, j, delta_up[j], delta_down[j], 'U');
						
						// Round UP
						inst->x[j] += delta_up[j];

						updated = 1;
						update_slacks(inst, j, delta_up[j]);
						inst->objval += (inst->obj[j] * delta_up[j]);
					}

					// Third case: ZIminus < ZIplus && ZIminus < ZI --> Round DOWN
					else if (less_than(ZIminus, ZIplus) && less_than(ZIminus, ZI)) {

						// Skip variable if delta_down = 0
						if (zero(delta_down[j])) continue;
						
						print_verbose(20, "[ziround]: >>> Round x_%d = %f - %f = %f\n", j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);

						// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
						check_slacks(inst, j, delta_up[j], delta_down[j], 'D');
						
						// Round DOWN
						inst->x[j] -= delta_down[j];

						updated = 1;
						update_slacks(inst, j, -(delta_down[j]));
						inst->objval -= (inst->obj[j] * delta_down[j]);
					}

					break;

				default:
					print_error(" in function is_fractional.\n");
			}

			// [DEBUG ONLY]: Print solution fractionality and cost ---------------------------------------------------------------
			frac[bufind] = sol_fractionality(inst->x, inst->int_var, inst->ncols);
			objval[bufind] = inst->objval;
			if (not_equals(frac[bufind], frac[!bufind]) || not_equals(objval[bufind], objval[!bufind])) print_verbose(10, "* %f | %f *\n", frac[bufind], objval[bufind]);
			if (VERBOSE >= 10) {
				if (PLOT_SOL_FRAC) add_point_single_tracker(frac[bufind], &tracker_sol_frac, &len_frac, &size_frac);
				if (PLOT_SOL_COST) add_point_single_tracker(objval[bufind], &tracker_sol_cost, &len_cost, &size_cost);
				if (PLOT_NUM_ROUNDED_VARS) add_point_single_tracker(count_rounded(inst->x, inst->ncols, inst->int_var, inst->vartype), &tracker_rounded, &len_rounded, &size_rounded);
			}
			bufind = !bufind;
			// -------------------------------------------------------------------------------------------------------------------

		} // end inner loop

		if (updated) { print_verbose(10, "[zi_round]: ... Some roundings occured, scan variables again ...\n"); }
		else { print_verbose(10, "[zi_round]: ... No roundings, exit outer loop ...\n"); }

		// [DEBUG ONLY] Pause after each inner loop execution
		if (VERBOSE >= 100) system("pause");
		
		// [DEBUG ONLY] (BRUTE FORCE)  Check variable bounds and constraints
		if (VERBOSE >= 201) {
			check_bounds(inst->x, inst->lb, inst->ub, inst->ncols);
			check_constraints(inst->x, inst->ncols, inst->nrows, inst->nzcnt, inst->rmatbeg, inst->rmatind, inst->rmatval, inst->sense, inst->rhs);
		}

	} while (updated); // end outer loop

	// [DEBUG ONLY]: Plot solution fractionality and cost trackers -------------------------------------------------
	if (VERBOSE >= 10) {
		if (PLOT_SOL_FRAC) {
			char** labels = (char**)calloc(2, sizeof(char*));
			labels[0] = (char*)calloc(20, sizeof(char));
			labels[1] = (char*)calloc(20, sizeof(char));
			char* name = (char*)calloc(20, sizeof(char));
			sprintf(labels[0], "Round");
			sprintf(labels[1], "Fractionality");
			sprintf(name, "Solution");
			plot_tracker(tracker_sol_frac, name, labels, size_frac, NULL);
			free_all(5, tracker_sol_frac, labels[0], labels[1], labels, name);
		}
		if (PLOT_SOL_COST) {
			char** labels = (char**)calloc(2, sizeof(char*));
			labels[0] = (char*)calloc(20, sizeof(char));
			labels[1] = (char*)calloc(20, sizeof(char));
			char* name = (char*)calloc(20, sizeof(char));
			sprintf(labels[0], "Round");
			sprintf(labels[1], "Cost");
			sprintf(name, "Solution");
			plot_tracker(tracker_sol_cost, name, labels, size_cost, NULL);
			free_all(5, tracker_sol_cost, labels[0], labels[1], labels, name);
		}
		if (PLOT_NUM_ROUNDED_VARS) {
			char** labels = (char**)calloc(2, sizeof(char*));
			labels[0] = (char*)calloc(20, sizeof(char));
			labels[1] = (char*)calloc(20, sizeof(char));
			char* name = (char*)calloc(20, sizeof(char));
			sprintf(labels[0], "Round");
			sprintf(labels[1], "Rounded Variables");
			sprintf(name, "Solution");
			plot_tracker(tracker_rounded, name, labels, size_rounded, NULL);
			free_all(5, tracker_rounded, labels[0], labels[1], labels, name);
		}
	}
	// -------------------------------------------------------------------------------------------------------------

	// Free
	free(delta_up);
	free(delta_down);
}
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************
//**************************************************************************************************************************************************************

void check_slacks(instance* inst, int j, double delta_up, double delta_down, const char round_updown) {

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
			case 'G': // (slack non-positive)

				curr_slack = inst->slack[rowind];
				(inst->sense[rowind] == 'L') ? assert(non_negative(curr_slack)) : assert(non_positive(curr_slack));
				delta_slack = (round_updown == 'U') ? (aij * delta_up) : (aij * (-delta_down));

				// Row slack after rounding (negative for 'L', positive for 'G' constraints iff also singletons slack should be used)
				new_slack = curr_slack - delta_slack;
				print_verbose(200, "[check_slacks][x_%d aij %f][row %d '%c']: new_slack = %f\n", j + 1, aij, rowind + 1, inst->sense[rowind], new_slack);

				// [EXTENSION] Distinguish inequality constraints with singletons (if extension enabled)
				if ((inst->extension) && (inst->num_singletons[rowind] > 0)) {

					// If the new row slack is negative for 'L', positive for 'G' constraints, then the remaining amount must be covered by the singletons slack
					if (negative(new_slack) || positive(new_slack)) {

						(inst->sense[rowind] == 'L') ? assert(negative(new_slack)) : assert(positive(new_slack));

						// Compute singletons slack of constraint rowind and get bounds
						ss_lb = inst->ss_lb[rowind];
						ss_ub = inst->ss_ub[rowind];
						singletons_slack = compute_singletons_slack(inst, rowind);
						assert(var_in_bounds(singletons_slack, ss_lb, ss_ub));

						// Singletons slack after rounding
						delta_ss = new_slack;
						new_ss = singletons_slack + delta_ss; // + because signed delta
						assert(var_in_bounds(new_ss, ss_lb, ss_ub));

						// New singletons slack must stay within its bounds
						enough_slack = var_in_bounds(new_ss, ss_lb, ss_ub); // !((new_ss < ss_lb - TOLERANCE) || (new_ss > ss_ub + TOLERANCE));
					}
					else {
						// Extension enabled, row has singletons, but new_slack is non-negative --> no need to use singletons
						enough_slack = 1;
					}
				}
				else {
					// Extension disabled OR enabled but no singletons
					enough_slack = (inst->sense[rowind] == 'L') ? non_negative(new_slack) : non_positive(new_slack);
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
					delta_slack = (round_updown == 'U') ? (aij * delta_up) : (aij * (-delta_down));

					// Singletons slack after rounding
					new_ss = curr_slack - delta_slack;
					assert(var_in_bounds(new_ss, ss_lb, ss_ub));

					// New singletons slack must stay within its bounds
					enough_slack = var_in_bounds(new_ss, ss_lb, ss_ub); // !((new_ss < ss_lb - TOLERANCE) || (new_ss > ss_ub + TOLERANCE));
					
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

	if (!zero(inst->obj[j])) {
		obj_deltaplus = (inst->obj[j] * delta_up[j]);
		obj_deltaminus = -(inst->obj[j] * delta_down[j]);
	}

	// Check obj sense, then update xj, update slacks and update objective value
	switch (inst->objsen) {

		case CPX_MIN:

			// [] Adding delta_up to x_j improves objval
			if (negative(obj_deltaplus) && less_than(obj_deltaplus, obj_deltaminus)) {

				// Skip variable if delta_up = 0
				if (zero(delta_up[j])) return 0;

				print_verbose(20, "[round_xj_bestobj]: >>> Round x_%d = %f + %f = %f\n", j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);

				// If xj is integer (not fractional) then delta_up[j] should be 1.0
				assert(xj_fractional || equals(delta_up[j], 1.0));

				// Check whether all affected constraints have enough slack for a ROUND UP of xj
				check_slacks(inst, j, delta_up[j], delta_down[j], 'U');

				// Round UP (if xj is not fractional then delta_up[j] must be 1.0)
				inst->x[j] += delta_up[j];

				updated = 1;
				update_slacks(inst, j, delta_up[j]);
				inst->objval += obj_deltaplus;
			}
			// [] Adding -delta_down to x_j improves objval
			else if (negative(obj_deltaminus) && less_than(obj_deltaminus, obj_deltaplus)) {

				// Skip variable if delta_down = 0
				if (zero(delta_down[j])) return 0;

				print_verbose(20, "[round_xj_bestobj]: >>> Round x_%d = %f - %f = %f\n", j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);
				
				// If xj is integer (not fractional) then delta_down[j] should be 1.0
				assert(xj_fractional || equals(delta_down[j], 1.0));

				// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
				check_slacks(inst, j, delta_up[j], delta_down[j], 'D');

				// Round DOWN (if xj is not fractional then delta_down[j] must be 1.0)
				inst->x[j] -= delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval += obj_deltaminus;
			}
			// [] Both deltas improve objval of the same amount < 0 --> Round arbitrarily (DOWN)
			else if (equals(obj_deltaminus, obj_deltaplus) && negative(obj_deltaminus)) {

				// Skip variable if delta_down = 0
				if (zero(delta_down[j])) return 0;

				print_verbose(20, "[round_xj_bestobj]: >>> Round x_%d = %f - %f = %f\n", j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);

				// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
				check_slacks(inst, j, delta_up[j], delta_down[j], 'D');

				// Round arbitrarily (DOWN)
				inst->x[j] -= delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval += obj_deltaminus;
			}
			// [] Both deltas do not change objval (both = 0) --> Round arbitrarily (UP)
			else if (zero(obj_deltaplus) && zero(obj_deltaminus)) {

				// Skip variable if delta_up = 0
				if (zero(delta_up[j])) return 0;

				print_verbose(20, "[round_xj_bestobj]: >>> Round x_%d = %f + %f = %f\n", j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
				
				// Check whether all affected constraints have enough slack for a ROUND UP of xj
				check_slacks(inst, j, delta_up[j], delta_down[j], 'U');

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
			if (positive(obj_deltaplus) && greater_than(obj_deltaplus, obj_deltaminus)) {

				// Skip variable if delta_up = 0
				if (zero(delta_up[j])) return 0;
				
				print_verbose(20, "[round_xj_bestobj]: >>> Round x_%d = %f + %f = %f\n", j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
				
				// If xj is integer (not fractional) then delta_up[j] should be 1.0
				assert(xj_fractional || equals(delta_up[j], 1.0));

				// Check whether all affected constraints have enough slack for a ROUND UP of xj
				check_slacks(inst, j, delta_up[j], delta_down[j], 'U');

				// Round UP (if xj is not fractional then delta_up[j] must be 1.0)
				inst->x[j] += delta_up[j];

				updated = 1;
				update_slacks(inst, j, delta_up[j]);
				inst->objval += obj_deltaplus;
			}
			// [] Adding delta_down to x_j improves objval
			else if (positive(obj_deltaminus) && greater_than(obj_deltaminus, obj_deltaplus)) {

				// Skip variable if delta_down = 0
				if (zero(delta_down[j])) return 0;
				
				print_verbose(20, "[round_xj_bestobj]: >>> Round x_%d = %f - %f = %f\n", j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);
				
				// If xj is integer (not fractional) then delta_down[j] should be 1.0
				assert(xj_fractional || equals(delta_down[j], 1.0));

				// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
				check_slacks(inst, j, delta_up[j], delta_down[j], 'D');

				// Round DOWN (if xj is not fractional then delta_down[j] must be 1.0)
				inst->x[j] -= delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval += obj_deltaminus;
			}
			// [] Both deltas improve objval of the same amount > 0 --> Round arbitrarily (DOWN)
			else if (equals(obj_deltaminus, obj_deltaplus) && positive(obj_deltaminus)) {

				// Skip variable if delta_down = 0
				if (zero(delta_down[j])) return 0;
				
				print_verbose(20, "[round_xj_bestobj]: >>> Round x_%d = %f - %f = %f\n", j + 1, inst->x[j], delta_down[j], inst->x[j] - delta_down[j]);

				// Check whether all affected constraints have enough slack for a ROUND DOWN of xj
				check_slacks(inst, j, delta_up[j], delta_down[j], 'D');

				// Round arbitrarily (DOWN)
				inst->x[j] -= delta_down[j];

				updated = 1;
				update_slacks(inst, j, -(delta_down[j]));
				inst->objval += obj_deltaminus;
			}
			// [] Both deltas do not change objval (both = 0) --> Round arbitrarily (UP)
			else if (zero(obj_deltaplus) && zero(obj_deltaminus)) {

				// Skip variable if delta_up = 0
				if (zero(delta_up[j])) return 0;

				print_verbose(20, "[round_xj_bestobj]: >>> Round x_%d = %f + %f = %f\n", j + 1, inst->x[j], delta_up[j], inst->x[j] + delta_up[j]);
				
				// Check whether all affected constraints have enough slack for a ROUND UP of xj
				check_slacks(inst, j, delta_up[j], delta_down[j], 'U');

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
					if (negative(temp_slack) || positive(temp_slack)) {

						(inst->sense[rowind] == 'L') ? assert(negative(temp_slack)) : assert(positive(temp_slack));

						// Delta singletons slack to distribute among the singletons [new_ss = ss + delta_ss (+ because signed delta)]
						delta_ss = temp_slack; // negative for 'L', positive for 'G' constraints

						// Distribute delta among the singletons, stop when done (delta_ss negative --> singletons slack must decrease)
						update_singletons(inst, rowind, delta_ss);
					}
					else {
						// Enough row slack, update it
						(inst->sense[rowind] == 'L') ? assert(non_negative(inst->slack[rowind] - delta_slack)) : assert(non_positive(inst->slack[rowind] - delta_slack));
						inst->slack[rowind] -= delta_slack;
					}
				}
				else {
					// Extension disabled OR enabled but zero singletons
					// Just update row slack
					print_verbose(201, "[update_slacks][x_%d][row %d '%c']: slack = %f - (%f * %f) = %f\n", j + 1, rowind + 1, inst->sense[rowind], inst->slack[rowind], aij, signed_delta, inst->slack[rowind] - delta_slack);
					(inst->sense[rowind] == 'L') ? assert(non_negative(inst->slack[rowind] - delta_slack)) : assert(non_positive(inst->slack[rowind] - delta_slack));
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

	s_slack_increase = (delta_ss >= 0.0);

	// Distribute delta among the singletons, stop when done (delta_ss positive(negative) --> singletons slack must increase(decrease))
	for (int k = 0; k < inst->num_singletons[rowind]; k++) {

		// Stop updating the singletons when delta singletons slack has been covered (s_slack_increase in the two conditions is necessary...)
		if ((s_slack_increase && non_positive(delta_ss)) || (!s_slack_increase &&  non_negative(delta_ss))) {
			print_verbose(200, "[update_singletons][extension][row %d '%c']: delta_ss covered, found %f\n", rowind + 1, inst->sense[rowind], delta_ss);
			break;
		}
		(s_slack_increase) ? assert(non_negative(delta_ss)) : assert(non_positive(delta_ss));
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
	assert(zero(delta_ss));
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
		non_negative(delta_up2) & 
		non_negative(delta_down2)
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
				non_negative(ss_delta_up) & 
				non_negative(ss_delta_down)
			);
		}
		
		// Check sense, then check sign of aij, and update delta_up1, delta_down1
		switch (inst->sense[rowind]) {

			case 'L': // (slack non-negative)

				if (negative(inst->slack[rowind])) print_error("[delta_updown][row %d 'L']: Found negative row slack = %f\n", rowind + 1, inst->slack[rowind]);

				// Clip slack to zero if slightly negative
				if ((inst->slack[rowind] < 0.0) && (inst->slack[rowind] > -(TOLERANCE))) { 
					inst->slack[rowind] = 0.0; 
					slack = inst->slack[rowind]; 
				}
				assert(equals(slack, inst->slack[rowind]));

				// [EXTENSION] Update available slack: 'L' constraint --> singletons slack (if any) should decrease
				if (inst->extension && inst->num_singletons[rowind] > 0) slack += ss_delta_down; // overall slack increases

				if (aij > 0.0) { 
					
					// Update delta_up1
					candidate_up1 = slack / aij;
					assert(non_negative(candidate_up1));
					delta_up1 = min(candidate_up1, delta_up1);
					assert(non_negative(delta_up1));
				}
				if (aij < 0.0) {

					// Update delta_down1
					candidate_down1 = -(slack) / aij;
					assert(non_negative(candidate_down1));
					delta_down1 = min(candidate_down1, delta_down1);
					assert(non_negative(delta_down1));
				}

				break;

			case 'G': // (slack non-positive)

				if (positive(inst->slack[rowind])) print_error("[delta_updown][row %d 'G']: Found positive row slack = %f\n", rowind + 1, inst->slack[rowind]);

				// Clip slack to zero if slightly positive
				if ((inst->slack[rowind] > 0.0) && (inst->slack[rowind] < TOLERANCE)) { 
					inst->slack[rowind] = 0.0;
					slack = inst->slack[rowind];
				}
				assert(equals(slack, inst->slack[rowind]));

				// [EXTENSION] Update available slack: 'G' constraint --> singletons slack (if any) should increase
				if (inst->extension && inst->num_singletons[rowind] > 0) slack -= ss_delta_up; // overall slack decreases (increases in absolute value)

				if (aij < 0.0) {

					// Update delta_up1
					candidate_up1 = slack / aij;
					assert(non_negative(candidate_up1));
					delta_up1 = min(candidate_up1, delta_up1);
					assert(non_negative(delta_up1));
				}
				if (aij > 0.0) {

					// Update delta_down1
					candidate_down1 = -(slack) / aij;
					assert(non_negative(candidate_down1));
					delta_down1 = min(candidate_down1, delta_down1);
					assert(non_negative(delta_down1));
				}

				break;

			case 'E': // (slack zero if extension disabled)

				// [EXTENSION] Distinguish equality constraints with singletons (if extension enabled)
				if (inst->extension && inst->num_singletons[rowind] > 0) {

					// Compute singletons slack of constraint rowind and get bounds (done above)
					print_verbose(201, "Singletons slack = %f. Bounds %f <= ss <= %f\n", singletons_slack, ss_lb, ss_ub);
					// Compute singletons slack deltas (done above)

					// Update candidate deltas
					if (aij > 0.0) {

						candidate_down1 = ss_delta_up / aij;
						candidate_up1 = ss_delta_down / aij;
						assert(
							non_negative(candidate_down1) & 
							non_negative(candidate_up1)
						);
						// Clip candidates to zero if inaccurate
						if (zero(candidate_down1)) candidate_down1 = 0.0;
						if (zero(candidate_up1)) candidate_up1 = 0.0;

						delta_down1 = min(candidate_down1, delta_down1);
						delta_up1 = min(candidate_up1, delta_up1);
					}
					if (aij < 0.0) {

						candidate_up1 = -(ss_delta_up) / aij;
						candidate_down1 = -(ss_delta_down) / aij;
						assert(
							non_negative(candidate_up1) & 
							non_negative(candidate_down1)
						);
						// Clip candidates to zero if inaccurate
						if (zero(candidate_down1)) candidate_down1 = 0.0;
						if (zero(candidate_up1)) candidate_up1 = 0.0;

						delta_down1 = min(candidate_down1, delta_down1);
						delta_up1 = min(candidate_up1, delta_up1);
					}
					assert(
						non_negative(delta_down1) & 
						non_negative(delta_up1)
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

	print_verbose(201, "[delta_updown][candidates]: delta_up1_%d = %f ; delta_down1_%d = %f\n", j + 1, delta_up1, j + 1, delta_down1);

	// Results
	assert(
		non_negative(delta_up1) & 
		non_negative(delta_down1)
	);
	new_delta_up = min(delta_up1, delta_up2);
	new_delta_down = min(delta_down1, delta_down2);
	assert(
		non_negative(new_delta_up) & 
		non_negative(new_delta_down)
	);
	print_verbose(201, "[delta_updown][results]: (NEW) delta_up_%d = min{%f, %f} = %f ; delta_down_%d = min{%f, %f} = %f\n", j + 1, delta_up1, delta_up2, new_delta_up, j + 1, delta_down1, delta_down2, new_delta_down);

	// Update deltas (clip them to zero if they are both less than epsilon)
	if (less_than(new_delta_up, epsilon) && less_than(new_delta_down, epsilon)) {
		new_delta_up = 0.0;
		new_delta_down = 0.0;
	}
	delta_up[j] = new_delta_up;
	delta_down[j] = new_delta_down;
	assert(
		equals(delta_up[j], new_delta_up) & 
		equals(delta_down[j], new_delta_down)
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