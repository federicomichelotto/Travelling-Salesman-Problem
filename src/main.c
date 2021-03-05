#include "../include/utils.h"
#include <time.h>
int main(int argc, char **argv) {

    time_t start, end;
    instance inst;          // current instance of the problem

//    if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }
//    if ( VERBOSE >= 2 ) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }

    start = time(NULL);

    parse_command_line(argc,argv, &inst);
    read_input(&inst);

    // todo solve the input instance
    for (int a = 0; a < 1000000; a++) printf("%d\n", a);

    end = time(NULL);

    if ( VERBOSE >= 1 ) printf("... VRP problem solved in %ld sec\n", end - start);

    // Free the memory used by the instance
    free(&inst);

    return 0;
}
