#include "ziround.h"

void parse_command_line(int argc, char** argv, instance* inst) {
	// Parse command line arguments
	int help = (argc < 2) ? 1 : 0;
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--input-mip") == 0) { inst->input_file = argv[++i]; continue; }
		// if (strcmp(argv[i], "--obj-sense") == 0) { inst->objsen = (strcmp(argv[++i], "min") == 0) ? CPX_MIN : CPX_MAX; continue; }
		if (strcmp(argv[i], "--help") == 0) { help = 1; continue; }
		if (strcmp(argv[i], "-help") == 0) { help = 1; continue; }
		help = 1;
	}
	// Print chosen parameters
	if (help || VERBOSE > 10) {
		fprintf(stdout, "\n\n**** Parameters ***************************\n");
		fprintf(stdout, "**   --input-mip %s\n", inst->input_file);
		// fprintf(stdout, "**   --obj-sense %s\n", ((inst->objsen == CPX_MIN) ? "min" : "max"));
		fprintf(stdout, "*******************************************\n");
	}
	if (help) exit(EXIT_FAILURE);
}