/**
 * @file plot.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

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

void add_point_multiple_trackers(double* point, double** tracker, int num, int* len, int* size) {

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

	fprintf(*pointer_to_g_plot_pipe, "set key left top \n");

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