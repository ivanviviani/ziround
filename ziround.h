/**
 * @file ziround.h
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#pragma once

/*  README.md
    ZI-Round MIP Rounding Heuristic (version 2)

    Usage: ziroundv2 --input-mip [filename.mps] --obj-sense [max/min]

    Variables:

        Primal feasible point x = [x1, x2, ..., xj, ...]
        Upper bounds of all xj
        Lower bounds of all xj
        Slacks si
        Constraint matrix coefficients aij (both by rows and by columns)
        Constraints senses {'L','E','G'}
        Objective coefficients
        Problem type {CPX_MAX,CPX_MIN}
        Threshold epsilon = 0.00001

    Formulas:

        Fractionality of xj: ZI(xj) = min{xj - floor(xj), ceil(xj) - xj}
        Integer infeasibility of x[]: ZI(x) = sum{ZI(xj)}
        For 'L' (<=) constraints: (si non-negative)
            ub_j = min_i{si/aij : aij > 0}
            lb_j = min_i{-si/aij : aij < 0}
        For 'G' (>=) constraints: (si non-positive)
            ub_j = min_i{si/aij : aij < 0}
            lb_j = min_i{-si/aij : aij > 0}
        UBj = min{ub_j, ub(xj) - xj}
        LBj = min{lb_j, xj - lb(xj)}

    ZI-Round Algorithm (version 2):
        Input: Primal feasible point x

        Repeat
        Loop For each integer variable xj
            If xj non-fractional
            Then
                Calculate UBj,LBj,1
                If cj > 0 and LBj = 1 or cj < 0 and UB = 1
                Then Update xj to improve objective and Update slacks
            If xj fractional
            Then
                Calculate UBj,LBj,epsilon
                If ZI(xj + UBj) = ZI(xj - LBj) and ZI(xj + UBj) < ZI(xj)
                Then Update xj to improve objective and Update slacks
                Else If ZI(xj + UBj) < ZI(xj - LBj) and ZI(xj + UBj) < ZI(xj)
                Then Set xj = xj + UBj and Update slacks
                Else If ZI(xj - LBj) < ZI(xj + UBj) and ZI(xj - LBj) < ZI(xj)
                Then Set xj = xj - LBj and Update slacks
        Until No more updates found

        Detail: stop calculating UBj,LBj if they both fall below epsilon
*/
#include <cplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

/**
 * @brief Verbosity level.
 */
#define VERBOSE 5

/**
 * @brief Tolerance for non-integer numbers as considered by CPLEX.
 */
#define TOLERANCE 1e-6

/**
 * @brief Threshold used in the ZI-Round heuristic.
 */
#define EPSILON 1e-5

/**
 * @brief Structure holding local and global information about a
 * 	      problem instance, parameters included.
 */
typedef struct {

    // Variables
    int nrows; 		        /**< Number of rows of the coefficients matrix. */
    int ncols; 		        /**< Number of variables, also columns of the coefficients matrix. */
    double* x; 				/**< Current problem solution. Will be rounded. */
    double* obj; 			/**< Objective function coefficients. */
    double* lb; 			/**< Variable lower bounds. */
    double* ub; 			/**< Variable upper bounds. */
    double* slack; 			/**< Row (constraint) slacks, defined as right hand side minus row activity. */
    double objval; 			/**< Current objective value (for current problem solution). */
    int objsen; 			/**< Objective function sense, CPX_MIN (default) or CPX_MAX (specified from command line). */
    char* mip_ctype; 		/**< Variable types (before converting MIP to LP), integer/binary or continuous. */
    int* int_var; 			/**< Flags array that keeps track of integer/binary (value 1) and continuous (value 0) variables. */

    // Coefficient matrix
    int nzcnt; 				/**< Number of non-zero coefficients. */
    int* rmatbeg; 			/**< Begin row indices of non-zero coefficients for rmatind and rmatval. */
    int* rmatind; 			/**< Column indices of non-zero coefficients. */
    double* rmatval; 		/**< Non-zero coefficients (row major). */
    int* cmatbeg; 			/**< Begin column indices of non-zero coefficients for cmatind and cmatval. */
    int* cmatind; 			/**< Row indices of non-zero coefficients. */
    double* cmatval; 		/**< Non-zero coefficients (column major). */
    char* sense; 			/**< Constraint (row) senses, 'L' (<=) or 'G' (>=) or 'E' (=). */
    double* rhs; 			/**< Constraint right hand sides (rhs). */

    // Parameters
    CPXENVptr env; 			/**< CPLEX environment pointer. */
    CPXLPptr lp; 			/**< CPLEX lp pointer. */
    char* input_file; 		/**< Input filename (mps format, specified from command line). */

} instance;

/**
 * @brief Initialize the appropriate fields of the instance.
 *
 * @param inst Pointer to the instance.
 */
void init_inst(instance* inst);

/**
 * @brief Parse the arguments from the command line and populate the
 *        appropriate fields of the already initialized instance.
 *
 * @param argc Number of command line arguments (includes the program name).
 * @param argv Program arguments (strings).
 * @param inst Pointer to the instance.
 */
void parse_cmd(int argc, char** argv, instance* inst);

/**
 * @brief Setup the CPLEX environment for the problem represented by the instance.
 * 		  Also turn on CPLEX screen output.
 *
 * @param inst Pointer to the already initialized and populated instance.
 */
void setup_CPLEX_env(instance* inst);

/**
 * @brief Create the lp into the CPLEX env of the instance and populate the lp
 * 		  with problem data read from a file (mps format).
 *
 * @param inst Pointer to the already populated instance.
 * @param filename Name of the input file (mps format).
 */
void read_MIP_problem(instance* inst, char* filename);

/**
 * @brief Read the current number of rows and columns of the coefficients matrix
 * 	 	  and populate the corresponding fields of the instance.
 *
 * @param inst Pointer to the already populated instance.
 */
void read_problem_sizes(instance* inst);

/**
 * @brief Read and save variable types from the MIP and save the integer/binary ones in a
 * 		  flags array. Both arrays are fields of the instance.
 *
 * @param inst Pointer to the already populated instance.
 */
void save_integer_variables(instance* inst);

/**
 * @brief Change the problem type from MIP to LP and solve its continuous relaxation.
 *
 * @param inst Pointer to the already populated instance.
 */
void solve_continuous_relaxation(instance* inst);

/**
 * @brief Read the continuous relaxation solution from the CPLEX lp
 *        and populate the corresponding array field of the instance.
 *
 * @param inst Pointer to the already populated instance.
 */
void read_solution(instance* inst);

/**
 * @brief Read the variable bounds from the CPLEX lp
 * 		  and populate the corresponding array fields of the instance.
 *
 * @param inst Pointer to the already populated instance.
 */
void read_variable_bounds(instance* inst);

/**
 * @brief Read the objective value from the CPLEX lp
 * 		  and populate the corresponding field of the instance.
 *
 * @param inst Pointer to the already populated instance.
 */
void read_objective_value(instance* inst);

/**
 * @brief Read the objective function coefficients from the CPLEX lp
 * 		  and populate the corresponding array field of the instance.
 *
 * @param inst Pointer to the already populated instance.
 */
void read_objective_coefficients(instance* inst);

/**
 * @brief Read the non-zero constraint coefficients both by row and column
 * 		  from the CPLEX lp, according to the data structures used by CPLEX,
 * 	      and populate the corresponding array fields of the instance.
 *
 * @param inst Pointer to the already populated instance.
 */
void read_constraints_coefficients(instance* inst);

/**
 * @brief Read the constraint senses from the CPLEX lp,
 * 		  and populate the corresponding array field of the instance.
 *
 * @param inst Pointer to the already populated instance.
 */
void read_constraints_senses(instance* inst);

/**
 * @brief Read the constraint right hand sides from the CPLEX lp,
 * 		  and populate the corresponding array field of the instance.
 *
 * @param inst Pointer to the already populated instance.
 */
void read_constraints_right_hand_sides(instance* inst);

/**
 * @brief Read the constraint (row) slacks for the continuous relaxation
 * 		  solution from the CPLEX lp, and populate the corresponding
 * 		  array field of the instance. They will then be updated during
 * 		  the rounding phase of zi_round.
 *
 * @param inst Pointer to the already populated instance.
 */
void read_row_slacks(instance* inst);

/**
 * @brief Read the problem data from the CPLEX lp using all the read functions,
 * 		  and populate the corresponding fields of the instance.
 *
 * @param inst Pointer to the already populated instance.
 */
void populate_inst(instance* inst);

/**
 * @brief Use the ZI-Round heuristic to round the continuous relaxation solution
 * 	      until an integer feasible solution is found or the procedure terminates
 * 		  unsuccessfully.
 *
 * @param inst Pointer to the already populated instance.
 */
void zi_round(instance* inst);

/**
 * @brief Update variable xj to improve objective, according to the values
 * of UBj and LBj calculated before. Also update the row slacks.
 *
 * @param inst Pointer to the already populated instance.
 * @param objval Pointer to the objective value.
 * @param j Column (variable xj) index.
 * @param delta_up Array of possible up-shifts.
 * @param delta_down Array of possible down-shifts.
 * @param is_fractional Flag that indicates whether xj is fractional or integer.
 * @return 1 if the variable xj has been updated, 0 otherwise.
 */
int round_xj_bestobj(instance* inst, double* objval, int j, double* delta_up, double* delta_down, int is_fractional);

/**
 * @brief Update the slack array field of the instance (incrementally)
 *        after a rounding of the current solution (x) field of the instance.
 *
 * @details Whenever it is called, only one variable xj has been updated.
 *
 * @param inst Pointer to the already populated instance.
 * @param j Index of the (only) variable just updated.
 * @param signed_delta Signed delta xj.
 */
void update_slacks(instance* inst, int j, double signed_delta);

/**
 * @brief Update the objective value (incrementally) for the current updated solution x.
 *
 * @details Whenever it is called, only one variable xj has been updated.
 *
 * @param inst Pointer to the already populated instance.
 * @param objval Pointer to the objective value.
 * @param j Index of the (only) variable just updated.
 * @param signed_delta Signed delta xj.
 */
void update_objval(instance* inst, double* objval, int j, double signed_delta);

/**
 * @brief Compute the j-th entries of the arrays of possible up-shifts and down-shifts
 *        according to the ZI-Round heuristic specifications.
 *
 * @param inst Pointer to the already populated instance.
 * @param j Variable index.
 * @param delta_up Array of possible up-shifts.
 * @param delta_down Array of possible down-shifts.
 * @param epsilon Tolerance for the computed up/down shifts.
 */
void delta_updown(instance* inst, int j, double* delta_up, double* delta_down, const double epsilon);

/**
 * @brief Check that all the variable bounds are satisfied, for the
 * 		  given solution x.
 *
 * @param inst Pointer to the already populated instance.
 * @param x Solution to be used for evaluating bounds satisfiability.
 */
void check_bounds(instance* inst, double* x);

/**
 * @brief Check that all the constraints are satisfied, for the
 * 		  given solution x.
 *
 * @param inst Pointer to the already populated instance.
 * @param x Solution to be used for evaluating constraints satisfiability.
 */
void check_constraints(instance* inst, double* x);

/**
 * @brief Deallocate all the fields of the instance.
 *
 * @param inst Pointer to the already populated instance.
 */
void free_inst(instance* inst);

/**
 * @brief [DEBUG] Print problem information to standard output or file.
 *
 * @param inst Pointer to the already populated instance.
 * @param sol_available Flag set to 1 if the instance contains a solution, 0 otherwise.
 * @param to_file Flag set to 1 if the problem info is to be printed to file, 0 otherwise.
 */
void print_problem_info(instance* inst, int sol_available, int to_file);

/**
 * @brief Calculate the fractionality of the given value.
 *
 * @param xj Value.
 * @return Fractionality of the value.
 */
double fractionality(double xj);

/**
 * @brief Check the integrality of the given value.
 *
 * @param num Value.
 * @return 1 if the value is fractional, 0 otherwise.
 */
int is_fractional(double num);

/**
 * @brief Calculate the dot product of two arrays (of the same length).
 *
 * @param coef First array.
 * @param var_value Second array.
 * @param len Length of both arrays.
 * @return Dot product of the two arrays.
 */
double dot_product(double* coef, double* var_value, int len);

/**
 * @brief Clone an array into a copy (both already allocated).
 *
 * @param arr The array to be cloned.
 * @param clo Clone copy of the given array.
 * @param len Length of the array.
 */
void clone_array(double* arr, double* clo, int len);

/**
 *	@brief Print a warning message. The message is preceeded by `"\n\nWARNING: "`.
 *
 *	@param warn The warning message to print, including the format for the ellipsis parameters.
 *	@param ... The multiple parameters.
 */
void print_warning(const char* warn, ...);

/**
 *	@brief Print an error message and exit. The message is preceeded by `"\n\nERROR: "`.
 *
 *	@param err The error message to print, including the format for the ellipsis parameters.
 *	@param ... The multiple parameters.
 */
void print_error(const char* err, ...);

/**
 *	@brief Print the passed message iff the verbosity level of the message is lower
 *	than the global parameter `VERBOSE`.
 *
 *	@param msg_verb The level of verbosity of this message.
 *	@param format The message to print, including the format for the ellipsis parameters.
 *	@param ... The multiple parameters.
 */
void print_verbose(int msg_verb, const char* format, ...);