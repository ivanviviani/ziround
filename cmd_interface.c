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
		if (!strcmp(argv[i], "-input"))		  { strcpy(inst->input_file, argv[++i]);       continue; }
		if (!strcmp(argv[i], "-folder"))      { strcpy(inst->input_folder, argv[++i]);     continue; }
		if (!strcmp(argv[i], "-singletons"))  { inst->singletons = atoi(argv[++i]);        continue; }
		if (!strcmp(argv[i], "-nonfracvars")) { inst->shift_nonfracvars = atoi(argv[++i]); continue; }
		if (!strcmp(argv[i], "-sortsinglet")) { inst->sort_singletons = atoi(argv[++i]);   continue; }
		if (!strcmp(argv[i], "-after0frac"))  { inst->after0frac = atoi(argv[++i]);        continue; }
		if (!strcmp(argv[i], "-timelimit"))   { inst->timelimit = atoi(argv[++i]);         continue; }
		if (!strcmp(argv[i], "-rseed"))       { inst->rseed = atoi(argv[++i]);	           continue; }
		if (!strcmp(argv[i], "--help"))       { help = 1;                                  continue; }
		if (!strcmp(argv[i], "-help"))        { help = 1;                                  continue; }
		if (!strcmp(argv[i], "-h"))			  { help = 1;                                  continue; }
		print_warning("Invalid command: %s.\n", argv[i]);
		help = 1;
	}

	// Print chosen parameters
	print_verbose(10, "CHOSEN PARAMETERS -------------------------------------------------------------\n");
	print_verbose(10, "[] input %s\n",       inst->input_file);
	print_verbose(10, "[] folder %s\n",      inst->input_folder);
	print_verbose(10, "[] singletons %d\n",  inst->singletons);
	print_verbose(10, "[] nonfracvars %d\n", inst->shift_nonfracvars);
	print_verbose(10, "[] sortsinglet %d\n", inst->sort_singletons);
	print_verbose(10, "[] after0frac %d\n", inst->after0frac);
	print_verbose(10, "[] timelimit %d\n",   inst->timelimit);
	print_verbose(10, "[] rseed %d\n",       inst->rseed);
	print_verbose(10, "--------------------------------------------------------------------------------\n\n");

	// Help menu
	if (help) {
		print_verbose(10, "HELP MENU ----------------------------------------------------------------------\n");
		print_verbose(10, "[] -input <path/filename.mps>: Input MIP problem.\n");
		print_verbose(10, "[] -folder <foldername>:       Input folder. \n");
		print_verbose(10, "[] -singletons [1|0]:          Flag for controlling the use of singletons in ZI-Round (default 1 = ON).\n");
		print_verbose(10, "[] -nonfracvars [1|0]:         Flag for controlling the shifting of also non-fractional integer variables in ZI-Round (default 1 = ON).\n");
		print_verbose(10, "[] -sortsinglet [1|0]:         Flag for controlling the sorting of the singletons in ascending order of objective coefficients (default 0 = OFF).\n");
		print_verbose(10, "[] -after0frac [1|0]:          Flag for activating the shifting of also non-fractional integer variables in ZI-Round only when fractionality reaches zero (default 0 = OFF).\n");
		print_verbose(10, "[] -timelimit <seconds>:       Execution time limit in seconds (default 300).\n");
		print_verbose(10, "[] -rseed <integer>:           Random seed (default -1). \n");
		print_verbose(10, "[] -help, --help, -h:          Show help menu.\n");
		print_verbose(10, "--------------------------------------------------------------------------------\n\n");
		exit(EXIT_FAILURE);
	}
}