#include "../include/utils.h"

void parse_command_line(int argc, char **argv, instance inst) {

    // Setting default values
    inst.model_type = 0;
    strcpy(inst.input_file, "NOT DEFINED YET");
    inst.time_limit = CPX_INFBOUND;

    if (argc < 2) {
        fprintf(stderr,"Usage: %s -h for help\n", argv[0]);
        exit(1);

    } else if (argc == 2 && strcmp(argv[1], "-h") == 0) {

        print_help();

    } else {

        for (int i = 1; i < argc; i++) {

            if (strcmp(argv[i], "-f") == 0) {   // Input file
                strcpy(inst.input_file, argv[++i]);
                continue;
            }

            if (strcmp(argv[i], "-t") == 0) {   // Time limit
                inst.time_limit = strtof(argv[++i], NULL);
                continue;
            }

            if (strcmp(argv[i], "-m") == 0) {   // Model type
                inst.model_type = strtol(argv[++i], NULL, 10);
                continue;
            }

            if (strcmp(argv[i], "-v") == 0) {   // Verbosity
                switch (strtol(argv[++i], NULL, 10)) {
                    case NORMAL:
                        verbose = NORMAL;
                        break;
                    case VERBOSE:
                        verbose = VERBOSE;
                        break;
                    case NERD:
                        verbose = NERD;
                        break;
                    case DEBUG:
                        verbose = DEBUG;
                        break;
                    default:
                        verbose = NORMAL;
                        break;
                }
                continue;
            }

        }

        print_command_line(inst);

    }

}

void parse_instance(instance inst) {

//    TODO write the read_input function

}

void print_command_line(instance inst) {
    printf("\nPARAMETERS ---------------------------------------------\n");
    printf("-f %s\n", inst.input_file);
    printf("-t %f\n", inst.time_limit);
    printf("-m %d\n", inst.model_type);
    printf("-v %s\n", verbose_name[verbose]);
    printf("--------------------------------------------------------\n\n");
}

void print_help() {
    printf("\nHELP ---------------------------------------------------\n");
    printf("-f <path>  : used to pass the relative instance path \n");
    printf("-t <time>  : used to pass the total running time allowed\n");
    printf("-m <model> : used to set up the model type\n");
    printf("-v <value> : used to set up the verbosity, from NORMAL (0) up to DEBUG (3)\n");
    printf("--------------------------------------------------------\n\n");
    exit(1);
}

void print_error(const char *err) {
    fprintf(stderr, "\nERROR: %s \n\n", err);
    fflush(NULL);
    exit(1);
}