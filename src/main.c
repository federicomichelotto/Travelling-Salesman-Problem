#include "../include/utils.h"
#include <sys/time.h>

int main(int argc, char **argv)
{

    time_t start, end;
    instance inst; // current instance of the problem

    initialize_instance(&inst);
    parse_command_line(argc, argv, &inst);

    parse_instance(&inst);
    print_instance(&inst);

    struct timeval tval_before, tval_after, tval_result;

    gettimeofday(&tval_before, NULL);

    TSPopt(&inst) ? print_error("TSPopt() ") : print_message("All went good inside TSPopt()");

    gettimeofday(&tval_after, NULL);

    timersub(&tval_after, &tval_before, &tval_result);

    if (verbose >= QUIET)
        printf("TSP solved in %ld.%06ld\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

    if (plot_solution(&inst)) print_error("plot_solution() error");
    //plot_solution_edges(inst.n_edges, inst.nodes, inst.edges) ? print_error("plot_solution_edges() error") : printf("... gnuplot ok\n");

    generate_csv_record(inst.param.name, inst.param.seed, inst.model_type, inst.z_best, (long int)tval_result.tv_sec, (long int)tval_result.tv_usec, inst.param.run);

    // Free the memory used by the instance
    free_instance(&inst);

    return 0;
}
