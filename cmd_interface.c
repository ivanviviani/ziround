/**
 * @file cmd_interface.c
 * @author Ivan Viviani
 * @copyright Copyright (c) 2020
 */

#include "ziround.h"

void parse_cmd(int argc, char** argv, instance* inst) {

	int help = (argc < 2);

	// Parse command line arguments
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-input"))		{ strcpy(inst->input_file, argv[++i]);   continue; }
		if (!strcmp(argv[i], "-folder"))    { strcpy(inst->input_folder, argv[++i]); continue; }
		if (!strcmp(argv[i], "-ext"))       { inst->extension = atoi(argv[++i]);     continue; }
		if (!strcmp(argv[i], "-timelimit")) { inst->timelimit = atoi(argv[++i]);     continue; }
		if (!strcmp(argv[i], "-rseed"))     { inst->rseed = atoi(argv[++i]);	     continue; }
		if (!strcmp(argv[i], "--help"))     { help = 1;                              continue; }
		if (!strcmp(argv[i], "-help"))      { help = 1;                              continue; }
		if (!strcmp(argv[i], "-h"))			{ help = 1;                              continue; }
		print_warning("Invalid command: %s.\n", argv[i]);
		help = 1;
	}

	// Print chosen parameters
	print_verbose(10, "CHOSEN PARAMETERS -------------------------------------------------------------\n");
	print_verbose(10, "[] input %s\n",  inst->input_file);
	print_verbose(10, "[] folder %s\n", inst->input_folder);
	print_verbose(10, "[] ext %d\n",    inst->extension);
	print_verbose(10, "[] timelimit\n");
	print_verbose(10, "--------------------------------------------------------------------------------\n\n");

	// Help menu
	if (help) {
		print_verbose(10, "HELP MENU ----------------------------------------------------------------------\n");
		print_verbose(10, "[] -input <path/filename.mps>: Input MIP problem.\n");
		print_verbose(10, "[] -folder <foldername>:       Input folder. \n");
		print_verbose(10, "[] -ext [1|0]:                 Flag for activating the ZI-Round extension (default 0).\n");
		print_verbose(10, "[] -timelimit <seconds>:       Execution time limit (default 3600).\n");
		print_verbose(10, "[] -rseed <integer>:           Random seed (default -1). \n");
		print_verbose(10, "[] -help, --help, -h:          Show help menu.\n");
		print_verbose(10, "--------------------------------------------------------------------------------\n\n");
		exit(EXIT_FAILURE);
	}
}