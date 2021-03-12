#include "../include/utils.h"

int main(int argc, char **argv) {

    time_t start, end;
    instance inst;          // current instance of the problem

    initialize_instance(&inst);
    //print_instance(&inst); // TODO remove after check
    parse_command_line(argc, argv, &inst);

    parse_instance(&inst);
    print_instance(&inst);

    start = time(NULL);

    TSPopt(&inst) ? print_error("TSPopt() ") : print_message("All went good inside TSPopt()");

    end = time(NULL);

    if (verbose >= QUIET) printf("... TSP solved in %ld sec\n", end - start);

    // Free the memory used by the instance
    free_instance(&inst);

    return 0;
}
