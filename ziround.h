/**
 * @file ziround.h
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#pragma once
#ifndef ZIROUND_H_
#define ZIROUND_H_

#include <cplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#include <time.h>

/**
 * @brief Verbosity level.
 */
#define VERBOSE 10

/**
 * @brief Plot solution fractionality flag (0-1).
 */
#define PLOT_SOL_FRAC 1

/**
* @brief Plot solution cost flag (0-1).
*/
#define PLOT_SOL_COST 1

/**
 * @brief Plot number of rounded variables flag (0-1).
 */
#define PLOT_NUM_VARS_TO_ROUND 1

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
    int nrows;                /**< Number of rows of the coefficients matrix. */
    int ncols;                /**< Number of variables, also columns of the coefficients matrix. */
    double* x;                /**< Current problem solution. Will be rounded. */
    double* obj;              /**< Objective function coefficients. */
    double* lb;               /**< Variable lower bounds. */
    double* ub;               /**< Variable upper bounds. */
    double* slack;            /**< Row (constraint) slacks, defined as right hand side minus row activity. */
    double objval;            /**< Current objective value (for current problem solution). */
    int objsen;               /**< Objective function sense, CPX_MIN (default) or CPX_MAX (specified from command line). */
    char* vartype;            /**< Variable types (before converting MIP to LP), integer/binary or continuous. */
    int* int_var;             /**< Flags array that keeps track of integer/binary (value 1) and continuous (value 0) variables. */
    int num_int_vars;         /**< Number of integer/binary variables to round. */
    double solfrac;           /**< Solution fractionality. */

    // Singletons
    int* row_singletons;      /**< Singleton indices. */
    int rs_size;              /**< Total number of singletons. */
    int* rs_beg;              /**< Begin index of the singleton indices of each row that contains at least one. */
    double* rs_coef;          /**< Coefficients of the singletons. */
    int* num_singletons;      /**< Number of singletons for each row. */
    double* ss_val;           /**< Singleton slack values for each row. */
    double* ss_ub;            /**< Upper bounds of the singletons slack for each row. */
    double* ss_lb;            /**< Lower bounds of the singletons slack for each row. */

    // Constraints
    int nzcnt;                /**< Number of non-zero coefficients. */
    int* rmatbeg;             /**< Begin row indices of non-zero coefficients for rmatind and rmatval. */
    int* rmatind;             /**< Column indices of non-zero coefficients. */
    double* rmatval;          /**< Non-zero coefficients (row major). */
    int* cmatbeg;             /**< Begin column indices of non-zero coefficients for cmatind and cmatval. */
    int* cmatind;             /**< Row indices of non-zero coefficients. */
    double* cmatval;          /**< Non-zero coefficients (column major). */
    char* sense;              /**< Constraint (row) senses, 'L' (<=) or 'G' (>=) or 'E' (=). */
    double* rhs;              /**< Constraint right hand sides (rhs). */

    // Plotting variables
    int size_frac;            /**< Actual current size of the solution fractionality tracker array. */
    int size_cost;            /**< Actual current size of the solution cost tracker array. */
    int size_toround;         /**< Actual current size of the number of variables to round array. */
    int len_frac;             /**< Maximum length of the solution fractionality tracker array (resizable). */
    int len_cost;             /**< Maximum length of the solution cost tracker array (resizable). */
    int len_toround;          /**< Maximum length of the number of rounded variables array (resizable). */
    double* tracker_sol_frac; /**< Tracker of solution fractionality. */
    double* tracker_sol_cost; /**< Tracker of solution cost. */
    double* tracker_toround;  /**< Tracker of number of variables to round. */

    // Parameters
    CPXENVptr env;            /**< CPLEX environment pointer. */
    CPXLPptr lp;              /**< CPLEX lp pointer. */
    char input_file[30];      /**< Input filename (mps format, specified from command line). */
    char input_folder[30];    /**< Input folder. */
    int singletons;           /**< Flag for the use of singletons in ZI-Round (default 1 = ON). */
    int timelimit;            /**< Time limit in seconds. */
    int rseed;                /**< Random seed. */

} instance;

// MAIN.C ----------------------------------------------------------------------------------------------

/**
 * @brief Test ZI-Round on a single instance.
 *
 * @param inst Pointer to the instance.
 */
void test_instance(instance* inst);

/**
 * @brief Test ZI-Round on a folder of instances.
 *
 * @param inst Pointer to the instance.
 */
void test_folder(instance* inst);

/**
 * @brief Plot the trackers implemented according to the PLOT_* macros.
 *
 * @param inst Pointer to the instance.
 */
void plot(instance* inst);
// -----------------------------------------------------------------------------------------------------

// INSTANCE.C ------------------------------------------------------------------------------------------

/**
 * @brief Initialize the appropriate fields of the instance.
 *
 * @param inst Pointer to the instance.
 */
void init_inst(instance* inst);

/**
 * @brief Deallocate all the fields of the instance.
 *
 * @param inst Pointer to the already populated instance.
 */
void free_inst(instance* inst);
// -----------------------------------------------------------------------------------------------------

// CMD_INTERFACE.C -------------------------------------------------------------------------------------

/**
 * @brief Parse the arguments from the command line and populate the
 *        appropriate fields of the already initialized instance.
 *
 * @param argc Number of command line arguments (includes the program name).
 * @param argv Program arguments (strings).
 * @param inst Pointer to the instance.
 */
void parse_cmd(int argc, char** argv, instance* inst);
// -----------------------------------------------------------------------------------------------------

// COMPUTE_ZIROUND_INPUT.C -----------------------------------------------------------------------------

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
// -----------------------------------------------------------------------------------------------------

// READ_ZIROUND_INPUT.C --------------------------------------------------------------------------------

/**
 * @brief Read the problem data from the CPLEX lp using all the read functions,
 * 		  and populate the corresponding fields of the instance.
 *
 * @param inst Pointer to the already populated instance.
 */
void populate_inst(instance* inst);

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
 * @brief Find singletons of the problem, i.e. continuous variables that
 * appear in only one constraint, and populate the corresponding data structures.
 *
 * @param inst Pointer to the instance.
 */
void find_singletons(instance* inst);

/**
 * @brief Compute values and bounds of the singletons slacks, seen as a single variable 
 * (sum of each singleton contribution in the row).
 *
 * @param inst Pointer to the instance.
 */
void compute_singletons_slacks(instance* inst);

/**
 * @brief Sort singleton indices and coefficients of a single row by lowest objective function coefficient.
 *
 * @details The row is represented by a \p start and \p end index of the arrays \p rs_ind and \p rs_coef.
 *
 * @param start Start index of the singleton indices and coefficients of the row.
 * @param end End index of the singleton indices and coefficients of the row.
 * @param rs_ind Singleton indices.
 * @param rs_coef Singleton coefficients.
 * @param obj Objective function coefficients.
 */
void sort_singletons(int start, int end, int* rs_ind, double* rs_coef, double* obj);
// -----------------------------------------------------------------------------------------------------

// ZIROUND.C -------------------------------------------------------------------------------------------

/**
 * @brief Use the ZI-Round heuristic to round the continuous relaxation solution
 * 	      until an integer feasible solution is found or the procedure terminates
 * 		  unsuccessfully.
 *
 * @param inst Pointer to the already populated instance.
 * @param numrounds Number of rounds (outer loops) of ZI-Round.
 */
void zi_round(instance* inst, int* numrounds);

/**
 * @brief Check whether all constraints affected by a round up/down of xj have enough slack for it.
 *
 * @param inst Pointer to the instance.
 * @param j Variable index.
 * @param delta_up Delta up of variable \p j.
 * @param delta_down Delta down of variable \p j.
 * @param round_updown
 */
void check_slacks(instance* inst, int j, double delta_up, double delta_down, const char round_updown);

/**
 * @brief Update variable \p j to improve objective, according to the values
 * of \p delta_up and \p delta_down calculated before. Also update the row slacks.
 *
 * @param inst Pointer to the already populated instance.
 * @param j Column (variable xj) index.
 * @param objcoef Objective coefficient of variable \p j.
 * @param delta_up Candidate up-shift of variable \p j.
 * @param delta_down Candidate down-shift of variable \p j.
 * @param is_fractional Flag that indicates whether variable \p j is fractional or integer.
 * @param solfrac Pointer to the solution fractionality.
 * @param num_toround Pointer to the number of variables to round.
 * @return 1 if variable \p j has been updated, 0 otherwise.
 */
int round_xj_bestobj(instance* inst, int j, double objcoef, double delta_up, double delta_down, int is_fractional, double* solfrac, int* num_toround);

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
 * @brief Update singletons of the constraint \p rowind by distributing \p delta_ss.
 *
 * @param inst Pointer to the instance.
 * @param rowind Index of the constraint.
 * @param delta_ss Delta singletons slack to distribute among the singletons of the constraint.
 */
void update_singletons(instance* inst, int rowind, double delta_ss);

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
 * @brief Compute singletons slack of a given constraint (row).
 *
 * @param inst Pointer to the instance.
 * @param rowind Index of the constraint.
 * @return Singletons slack of the constraint.
 */
double compute_ss_val(instance* inst, int rowind);
// -----------------------------------------------------------------------------------------------------

// UTIL.C ----------------------------------------------------------------------------------------------

/**
 * @brief Check that all the variable bounds are satisfied, for the given solution \p x.
 *
 * @param x Solution to be used for evaluating bounds satisfiability.
 * @param lb Lower bounds of the variables.
 * @param ub Upper bounds of the variables.
 * @param ncols Number of variables.
 */
void check_bounds(double* x, double* lb, double* ub, int ncols);

/**
 * @brief Check whether all the constraints are satisfied, for the given solution \p x.
 *
 * @param x Solution to be used for evaluating constraints satisfiability.
 * @param ncols Number of variables.
 * @param nrows Number of constraints.
 * @param nzcnt Number of non-zero coefficients in the constraints.
 * @param rmatbeg Constraints begin indices structure.
 * @param rmatind Constraints column indices structure.
 * @param rmatval Constraint coefficients.
 * @param sense Constraint senses.
 * @param rhs Constraint right hand sides.
 */
void check_constraints(double* x, int ncols, int nrows, int nzcnt, int* rmatbeg, int* rmatind, double* rmatval, char* sense, double* rhs);

/**
 * @brief Check whether all the integer variables of the original MIP have been rounded.
 *
 * @param x Current solution.
 * @param ncols Number of variables.
 * @param int_var Array of flags for binary/integer variables in the original MIP.
 * @param vartype Array of variable types.
 * @return 1 iff all integer/binary variables have been rounded, 0 otherwise.
 */
int check_rounding(double* x, int ncols, int* int_var, char* vartype);

/**
 * @brief Count number of variables rounded so far.
 *
 * @param x Current solution.
 * @param ncols Number of variables.
 * @param int_var Array of flags for binary/integer variables in the original MIP.
 * @param vartype Array of variable types.
 * @return Number of variables rounded so far.
 */
int count_rounded(double* x, int ncols, int* int_var, char* vartype);

/**
 * @brief Calculate the fractionality of the given value.
 *
 * @param xj Value.
 * @return Fractionality of the value.
 */
double fractionality(double xj);

/**
 * @brief Calculate the fractionality of a solution.
 *
 * @param x Solution.
 * @param int_var Flags array for integer variables.
 * @param len Length of the array.
 * @return Fractionality of the solution.
 */
double sol_fractionality(double* x, int* int_var, int len);

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
 * @brief Free memory allocated to multiple pointers.
 *
 * @details It does not check the type of the passed ellipsis parameters, so
 * be sure that actual pointers are passed. Use with care.
 *
 * @param count The number of pointers.
 * @param ... The actual pointers.
 */
void free_all(int count, ...);
// -----------------------------------------------------------------------------------------------------

// PLOT.C ----------------------------------------------------------------------------------------------

/**
 * @brief Given a tracker, add a new point to it.
 *
 * @param point The point to add to the \p tracker.
 * @param tracker The tracking array.
 * @param len The actual length of the \p tracker.
 * @param size The number of \p point elements currently in the tracker.
 */
void add_point_single_tracker(double point, double** tracker, int* len, int* size);

/**
 *  @brief Given a multivariate tracker, add a new point to it.
 *
 *  @param point   The new point to add to the multivariate tracker, represented by an array.
 *  @param tracker The tracking matrix.
 *  @param num     The number of dimensions of the multivariate tracker.
 *  @param len     The size of the tracker.
 *  @param size    The number of elements currently in each tracker.
 */
void add_point_multivariate_tracker(double* point, double** tracker, int num, int* len, int* size);

/**
 * @brief Plot the results of a single tracker.
 *
 * @param tracker The array with the tracking information.
 * @param name The name of the tracker.
 * @param label The labels of x and y axes.
 * @param size The size of the array.
 * @param filename The name of the file, where the plot is saved, this file
 *        should have the `.png` extension. If this parameter is `NULL` the plot
 *        is displayed on the screen.
 */
void plot_tracker(double* tracker, char* name, char** label, int size, char* filename);

/**
 * @brief Plot the results of a pair of trackers with two different y-axis.
 *
 * @param first_tracker The array with the first tracker information.
 * @param second_tracker The array with the second tracker information.
 * @param name The names of the two trackers.
 * @param label The labels of the two trackers (x, y1, y2).
 * @param size The size of both arrays.
 * @param filename The name of the file, where the plot is saved, this file
 *        should have the `.png` extension. If this parameter is `NULL` the plot
 *        is displayed on the screen.
 */
void plot_tracker_pair(double* first_tracker, double* second_tracker, char** name, char** label, int size, char* filename);

/**
 * @brief Plot the results of a multivariate tracker on the same plot.
 *
 * @param tracker  The matrix with tracking information.
 * @param name     The names of the tracker dimensions.
 * @param label    The labels of x and y axes.
 * @param num      The number of dimensions of the tracker.
 * @param size     The size of the tracker.
 * @param filename The name of the file, where the plot is saved, this file
 *        should have the .png extension. If this parameter is null, the plot
 *        is displayed on the screen.
 */
void plot_multivariate_tracker(double** tracker, char** name, char** label, int num, int size, char* filename);

/**
 * @brief Open a pipe, and set basic properties.
 *
 * @details The pipe is opened and the FILE* of the pipe is inserted in the
 * pointer at address \p pointer_to_g_plot_pipe. \p filename is used to
 * choose if the plot must be displayed or saved as an image: according to it,
 * gnuplot is properly set. Finally a set of `6` styles is defined: `(1,2,...6)`.
 *
 * @param pointer_to_g_plot_pipe The pointer to the pointer to the file that
 *        represents the pipe. This indirection is required to modify the
 *        value of the pointer to the file.
 * @param filename The name of the file, with the meaning described in the
 *        main plot functions.
 */
static void open_pipe(FILE** pointer_to_g_plot_pipe, char* filename);

/**
 * @brief Close the pipe, forcing a flush.
 *
 * @param g_plot_pipe The pipe to close.
 */
static void close_pipe(FILE* g_plot_pipe);
// -----------------------------------------------------------------------------------------------------

// PRINT.C ---------------------------------------------------------------------------------------------

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
// -----------------------------------------------------------------------------------------------------

// ASSERTS.C -------------------------------------------------------------------------------------------

/**
 * @brief Check whether an integer number is positive.
 *
 * @param num Integer number.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int positive_integer(int num);

/**
 * @brief Check whether an integer number is non-negative.
 *
 * @param num Integer number.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int non_negative_integer(int num);

/**
 * @brief Check whether a floating point number is non-negative.
 *
 * @param num Number.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int non_negative(double num);

/**
 * @brief Check whether a floating point number is non-positive.
 *
 * @param num Number.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int non_positive(double num);

/**
 * @brief Check whether a floating point number is negative.
 *
 * @param num Number.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int negative(double num);

/**
 * @brief Check whether a floating point number is positive.
 *
 * @param num Number.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int positive(double num);

/**
 * @brief Check whether a floating point number is equal to zero.
 *
 * @param num Number.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int zero(double num);

/**
 * @brief Check whether two floating point numbers are equal.
 *
 * @param x First number.
 * @param y Second number.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int equals(double x, double y);

/**
 * @brief Check whether two floating point numbers are not equal.
 *
 * @param x First number.
 * @param y Second number.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int not_equals(double x, double y);

/**
 * @brief Check whether for two floating point numbers it holds that the
 * first is less than the second.
 *
 * @param x First number.
 * @param y Second number.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int less_than(double x, double y);

/**
 * @brief Check whether for two floating point numbers it holds that the
 * first is greater than the second.
 *
 * @param x First number.
 * @param y Second number.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int greater_than(double x, double y);

/**
 * @brief Check whether the numeric code for the objective sense is valid
 * according to CPLEX.
 *
 * @param objsen Objective sense code.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int valid_obj_sense(int objsen);

/**
 * @brief Check whether there is a ranged constraint in an array of
 * constraint senses.
 *
 * @param sense Array of constraint senses.
 * @param nrows Size of the array.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int no_ranged_constraints(char* sense, int nrows);

/**
 * @brief Check whether row slacks are valid, i.e. have the correct sign.
 *
 * @param slack Array of row slacks.
 * @param sense Array of constraint senses.
 * @param nrows Size of the two arrays.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int valid_row_slacks(double* slack, char* sense, int nrows);

/**
 * @brief Check whether variable types are supported in the program.
 *
 * @param vartype Array of variable types.
 * @param ncols Size of the array.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int valid_var_types(char* vartype, int ncols);

/**
 * @brief Check whether a variable type is integer or binary.
 *
 * @param vartype Variable type.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int var_type_integer_or_binary(char vartype);

/**
 * @brief Check whether a variable type is continuous.
 *
 * @param vartype Variable type.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int var_type_continuous(char vartype);

/**
 * @brief Check whether an index is within its bounds 
 * (it is assumed to have a lower bound of zero).
 *
 * @param ind Index.
 * @param len Upper bound.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int index_in_bounds(int ind, int len);

/**
 * @brief Check whether an array of integers is all zeros.
 *
 * @param arr Array of integers.
 * @param len Size of the array.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int array_of_zeros(int* arr, int len);

/**
 * @brief Check whether a set of bounds is valid, i.e. lower bounds
 * less than upper bounds.
 *
 * @param lower Array of lower bounds.
 * @param upper Array of upper bounds.
 * @param nvars Size of the two arrays.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int valid_bounds(double* lower, double* upper, int nvars);

/**
 * @brief Check whether a variable is within its bounds.
 *
 * @param var Variable value.
 * @param lb Lower bound.
 * @param ub Upper bound.
 * @return 1 If the assert succeeds, 0 otherwise.
 */
int var_in_bounds(double var, double lb, double ub);
// -----------------------------------------------------------------------------------------------------

#endif /* ZIROUND_H_ */