/**
 * @file cmd_interface.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void parse_cmd(int argc, char** argv, instance* inst) {

	int help = (argc < 2) ? 1 : 0;

	// Parse command line arguments
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--input-mip") == 0) { inst->input_file = argv[++i];      continue; }
		if (strcmp(argv[i], "--ext") == 0)       { inst->extension = atoi(argv[++i]); continue; }
		if (strcmp(argv[i], "-ext") == 0)        { inst->extension = atoi(argv[++i]); continue; }
		if (strcmp(argv[i], "--help") == 0)      { help = 1;                          continue; }
		if (strcmp(argv[i], "-help") == 0)       { help = 1;                          continue; }
		help = 1;
	}

	// Print chosen parameters
	print_verbose(10, "CHOSEN PARAMETERS -------------------------------------------------------------\n");
	print_verbose(10, "[] input-mip %s\n", inst->input_file);
	print_verbose(10, "[] ext %d\n",       inst->extension);
	print_verbose(10, "--------------------------------------------------------------------------------\n\n");

	// Help menu
	if (help) {
		print_verbose(10, "HELP MENU ----------------------------------------------------------------------\n");
		print_verbose(10, "[] --input-mip <path/filename.mps>: Input MIP problem.\n");
		print_verbose(10, "[] --ext [1|0]: Flag for activating the ZI-Round extension (default 0).\n");
		print_verbose(10, "[] --help: Show help menu.\n");
		print_verbose(10, "--------------------------------------------------------------------------------\n\n");
		exit(EXIT_FAILURE);
	}
}