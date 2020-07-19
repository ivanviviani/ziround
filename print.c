/**
 * @file print.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void print_warning(const char* warn, ...) {
	printf("\n\nWARNING: ");
	va_list args;
	va_start(args, warn);
	vprintf(warn, args);
	va_end(args);
	fflush(NULL);
}

void print_error(const char* err, ...) {
	printf("\n\nERROR: ");
	va_list args;
	va_start(args, err);
	vprintf(err, args);
	va_end(args);
	fflush(NULL);

	exit(EXIT_FAILURE);
}

void print_verbose(int msg_verb, const char* format, ...) {
	if (VERBOSE >= msg_verb) {
		va_list args;
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
		fflush(NULL);
	}
}