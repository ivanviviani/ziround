/**
 * @file plot.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void plot(instance* inst) {

	// Plot trackers of solution fractionality, cost or number of variables to round

	// Try to plot solution fractionality and number of variables to round as a pair first
	if (PLOT_SOL_FRAC && PLOT_NUM_VARS_TO_ROUND) {
		char** label = (char**)calloc(3, sizeof(char*));
		label[0] = (char*)calloc(10, sizeof(char));
		label[1] = (char*)calloc(30, sizeof(char));
		label[2] = (char*)calloc(30, sizeof(char));
		char** name = (char**)calloc(2, sizeof(char*));
		name[0] = (char*)calloc(20, sizeof(char));
		name[1] = (char*)calloc(20, sizeof(char));
		char* temp = (char*)calloc(20, sizeof(char));
		sprintf(label[0], "Iteration");
		sprintf(label[1], "Solution Fractionality (SF)");
		sprintf(label[2], "#Variables to Round (#VR)");
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name[0], "%s (SF)", strtok(NULL, "/."));
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name[1], "%s (#VR)", strtok(NULL, "/."));
		assert(inst->size_frac == inst->size_toround);
		plot_tracker_pair(inst->tracker_sol_frac, inst->tracker_toround, name, label, inst->size_frac, NULL);
		free_all(8, label[0], label[1], label[2], label, name[0], name[1], name, temp);
	}
	else if (PLOT_SOL_FRAC) {
		char** label = (char**)calloc(2, sizeof(char*));
		label[0] = (char*)calloc(10, sizeof(char));
		label[1] = (char*)calloc(30, sizeof(char));
		char* name = (char*)calloc(20, sizeof(char));
		char* temp = (char*)calloc(20, sizeof(char));
		sprintf(label[0], "Iteration");
		sprintf(label[1], "Solution Fractionality");
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name, strtok(NULL, "/."));
		plot_tracker(inst->tracker_sol_frac, name, label, inst->size_frac, NULL);
		free_all(5, label[0], label[1], label, name, temp);
	}
	else if (PLOT_NUM_VARS_TO_ROUND) {
		char** label = (char**)calloc(2, sizeof(char*));
		label[0] = (char*)calloc(10, sizeof(char));
		label[1] = (char*)calloc(30, sizeof(char));
		char* name = (char*)calloc(20, sizeof(char));
		char* temp = (char*)calloc(20, sizeof(char));
		sprintf(label[0], "Iteration");
		sprintf(label[1], "#Variables to Round");
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name, strtok(NULL, "/."));
		plot_tracker(inst->tracker_toround, name, label, inst->size_toround, NULL);
		free_all(5, label[0], label[1], label, name, temp);
	}

	// Try to plot solution fractionality and solution cost as a pair first
	if (PLOT_SOL_FRAC && PLOT_SOL_COST) {
		char** label = (char**)calloc(3, sizeof(char*));
		label[0] = (char*)calloc(10, sizeof(char));
		label[1] = (char*)calloc(30, sizeof(char));
		label[2] = (char*)calloc(30, sizeof(char));
		char** name = (char**)calloc(2, sizeof(char*));
		name[0] = (char*)calloc(20, sizeof(char));
		name[1] = (char*)calloc(20, sizeof(char));
		char* temp = (char*)calloc(20, sizeof(char));
		sprintf(label[0], "Iteration");
		sprintf(label[1], "Solution Fractionality (SF)");
		sprintf(label[2], "Solution Cost (SC)");
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name[0], "%s (SF)", strtok(NULL, "/."));
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name[1], "%s (SC)", strtok(NULL, "/."));
		assert(inst->size_frac == inst->size_cost);
		plot_tracker_pair(inst->tracker_sol_frac, inst->tracker_sol_cost, name, label, inst->size_frac, NULL);
		free_all(8, label[0], label[1], label[2], label, name[0], name[1], name, temp);
	}
	else if (PLOT_SOL_FRAC) {
		char** label = (char**)calloc(2, sizeof(char*));
		label[0] = (char*)calloc(10, sizeof(char));
		label[1] = (char*)calloc(30, sizeof(char));
		char* name = (char*)calloc(20, sizeof(char));
		char* temp = (char*)calloc(20, sizeof(char));
		sprintf(label[0], "Iteration");
		sprintf(label[1], "Solution Fractionality");
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name, strtok(NULL, "/."));
		plot_tracker(inst->tracker_sol_frac, name, label, inst->size_frac, NULL);
		free_all(5, label[0], label[1], label, name, temp);
	}
	else if (PLOT_SOL_COST) {
		char** label = (char**)calloc(2, sizeof(char*));
		label[0] = (char*)calloc(10, sizeof(char));
		label[1] = (char*)calloc(30, sizeof(char));
		char* name = (char*)calloc(20, sizeof(char));
		char* temp = (char*)calloc(20, sizeof(char));
		sprintf(label[0], "Iteration");
		sprintf(label[1], "Solution Cost");
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name, strtok(NULL, "/."));
		plot_tracker(inst->tracker_sol_cost, name, label, inst->size_cost, NULL);
		free_all(5, label[0], label[1], label, name, temp);
	}

	// Try to plot solution cost and number of variables to round as a pair first
	/*
	if (PLOT_SOL_COST && PLOT_NUM_VARS_TO_ROUND) {
		char** label = (char**)calloc(3, sizeof(char*));
		label[0] = (char*)calloc(10, sizeof(char));
		label[1] = (char*)calloc(30, sizeof(char));
		label[2] = (char*)calloc(30, sizeof(char));
		char** name = (char**)calloc(2, sizeof(char*));
		name[0] = (char*)calloc(20, sizeof(char));
		name[1] = (char*)calloc(20, sizeof(char));
		char* temp = (char*)calloc(20, sizeof(char));
		sprintf(label[0], "Iteration");
		sprintf(label[1], "Solution Cost (SC)");
		sprintf(label[2], "#Variables to Round (#VR)");
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name[0], "%s (SC)", strtok(NULL, "/."));
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name[1], "%s (#VR)", strtok(NULL, "/."));
		assert(inst->size_cost == inst->size_toround);
		plot_tracker_pair(inst->tracker_sol_cost, inst->tracker_toround, name, label, inst->size_cost, NULL);
		free_all(8, label[0], label[1], label[2], label, name[0], name[1], name, temp);
	}
	else if (PLOT_SOL_COST) {
		char** label = (char**)calloc(2, sizeof(char*));
		label[0] = (char*)calloc(10, sizeof(char));
		label[1] = (char*)calloc(30, sizeof(char));
		char* name = (char*)calloc(20, sizeof(char));
		char* temp = (char*)calloc(20, sizeof(char));
		sprintf(label[0], "Iteration");
		sprintf(label[1], "Solution Cost");
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name, strtok(NULL, "/."));
		plot_tracker(inst->tracker_sol_cost, name, label, inst->size_cost, NULL);
		free_all(5, label[0], label[1], label, name, temp);
	}
	else if (PLOT_NUM_VARS_TO_ROUND) {
		char** label = (char**)calloc(2, sizeof(char*));
		label[0] = (char*)calloc(10, sizeof(char));
		label[1] = (char*)calloc(30, sizeof(char));
		char* name = (char*)calloc(20, sizeof(char));
		char* temp = (char*)calloc(20, sizeof(char));
		sprintf(label[0], "Iteration");
		sprintf(label[1], "#Variables to Round");
		sprintf(temp, inst->input_file);
		strtok(temp, "/.");
		sprintf(name, strtok(NULL, "/."));
		plot_tracker(inst->tracker_toround, name, label, inst->size_toround, NULL);
		free_all(5, label[0], label[1], label, name, temp);
	}
	*/
}

void add_point_single_tracker(double point, double** tracker, int* len, int* size) {

	size_t length = *len;

	// Resize the array if necessary
	if ((*len) == (*size)) {
		double* new_tracker = (double*)calloc(length * 2, sizeof(double));
		for (int i = 0; i < (*len); i++) new_tracker[i] = (*tracker)[i];
		free(*tracker);
		*len = (*len) * 2;
		*tracker = new_tracker;
	}

	// Add a new point
	(*tracker)[*size] = point;
	*size = *size + 1;
}

void add_point_multivariate_tracker(double* point, double** tracker, int num, int* len, int* size) {

	size_t length = *len;

	// Resize the array if necessary
	if ((*len) == (*size)) {
		for (int i = 0; i < num; i++) {
			double* new_tracker = (double*)calloc(length * 2, sizeof(double));
			for (int j = 0; j < (*len); j++) new_tracker[j] = tracker[i][j];
			free(tracker[i]);
			tracker[i] = new_tracker;
		}
		*len = (*len) * 2;
	}

	// Add a new point
	for (int i = 0; i < num; i++) {
		tracker[i][*size] = point[i];
	}
	*size = (*size) + 1;
}

void plot_tracker(double* tracker, char* name, char** label, int size, char* filename) {
	
	plot_multivariate_tracker(&tracker, &name, label, 1, size, filename);
}

void plot_tracker_pair(double* first_tracker, double* second_tracker, char** name, char** label, int size, char* filename) {
	
	FILE* g_plot_pipe;
	open_pipe(&g_plot_pipe, filename);

	fprintf(g_plot_pipe, "set xlabel \"%s\" \nset ylabel \"%s\" \nset y2label \"%s\" \nset ytics nomirror \nset y2tics nomirror \n", label[0], label[1], label[2]);
	fprintf(g_plot_pipe, "plot ");
	fprintf(g_plot_pipe, "'-' axis x1y1 with lines linestyle 1 title \"%s\"", name[0]);
	fprintf(g_plot_pipe, ", ");
	fprintf(g_plot_pipe, "'-' axis x1y2 with lines linestyle 2 title \"%s\" \n", name[1]);
	
	// Print first tracker
	for (int i = 0; i < size; i++) fprintf(g_plot_pipe, "%f %f\n", (double)i, first_tracker[i]);
	fprintf(g_plot_pipe, "e\n");

	// Print second tracker
	for (int i = 0; i < size; i++) fprintf(g_plot_pipe, "%f %f\n", (double)i, second_tracker[i]);
	fprintf(g_plot_pipe, "e\n");

	close_pipe(g_plot_pipe);
}

void plot_multivariate_tracker(double** tracker, char** name, char** label, int num, int size, char* filename) {
	FILE* g_plot_pipe;
	open_pipe(&g_plot_pipe, filename);

	fprintf(g_plot_pipe, "set xlabel \"%s\" \nset ylabel \"%s\" \n", label[0], label[1]);

	fprintf(g_plot_pipe, "plot ");
	for (int i = 0; i < num; i++) {
		fprintf(g_plot_pipe, "'-' with lines linestyle %d title \"%s\"", (i % 6) + 1, name[i]);
		if (i < num - 1) fprintf(g_plot_pipe, ", ");
	}
	fprintf(g_plot_pipe, "\n");

	// Print every tracker
	for (int i = 0; i < num; i++) {

		// Print every point of the tracker
		for (int j = 0; j < size; j++) fprintf(g_plot_pipe, "%f %f\n", (double)j, tracker[i][j]);
		fprintf(g_plot_pipe, "e\n");
	}

	close_pipe(g_plot_pipe);
}

static void open_pipe(FILE** pointer_to_g_plot_pipe, char* filename) {

	// Note: on Windows NT use _popen instead of popen
	*pointer_to_g_plot_pipe = _popen("gnuplot -persistent \n", "w");

	// Create styles for lines and points
	fprintf(*pointer_to_g_plot_pipe, "set style line 1 linecolor rgb '#0000C0' linewidth 2 pointtype 7 pointsize 1 \n");
	fprintf(*pointer_to_g_plot_pipe, "set style line 2 linecolor rgb '#00C000' linewidth 2 pointtype 7 pointsize 1 \n");
	fprintf(*pointer_to_g_plot_pipe, "set style line 3 linecolor rgb '#C00000' linewidth 2 pointtype 7 pointsize 1 \n");
	fprintf(*pointer_to_g_plot_pipe, "set style line 4 linecolor rgb '#00C0C0' linewidth 2 pointtype 7 pointsize 1 \n");
	fprintf(*pointer_to_g_plot_pipe, "set style line 5 linecolor rgb '#C000C0' linewidth 2 pointtype 7 pointsize 1 \n");
	fprintf(*pointer_to_g_plot_pipe, "set style line 6 linecolor rgb '#C0C000' linewidth 2 pointtype 7 pointsize 1 \n");
	fprintf(*pointer_to_g_plot_pipe, "set grid ytics lc rgb '#bbbbbb' lw 1 lt 0 \n");
	fprintf(*pointer_to_g_plot_pipe, "set grid xtics lc rgb '#bbbbbb' lw 1 lt 0 \n");

	// fprintf(*pointer_to_g_plot_pipe, "set key left top \n");

	// Save the image or display it on screen
	if (filename != NULL) {
		fprintf(*pointer_to_g_plot_pipe, "set terminal png size 800, 600 \n");
		char command[180];
		sprintf(command, "set output '%s'\n", filename);
		fprintf(*pointer_to_g_plot_pipe, command);
	}
}

static void close_pipe(FILE* g_plot_pipe) {

	// Write buffered data 
	fflush(g_plot_pipe);

	// Note: on Windows NT use _pclose instead of pclose
	_pclose(g_plot_pipe);
}