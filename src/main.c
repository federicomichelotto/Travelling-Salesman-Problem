#include "../include/utils.h"
#include "../include/tsp.h"

int main(int argc, char **argv) {

    time_t start, end;
    instance inst;          // current instance of the problem

    parse_command_line(argc, argv, inst);
    parse_instance(inst);

    start = time(NULL);

    // todo solve the input instance

    end = time(NULL);

    if (verbose) printf("... TSP solved in %ld sec\n", end - start);

    // Free the memory used by the instance
    free(&inst);

    return 0;
}
