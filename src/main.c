#include "../include/utils.h"

int main(int argc, char **argv)
{

    double time_elapsed;
    instance inst; // current instance of the problem

    initialize_instance(&inst);
    parse_command_line(argc, argv, &inst);

    parse_instance(&inst);
    print_instance(&inst);

    TSPopt(&inst) ? print_error("TSPopt() ") : print_message("All went good inside TSPopt()");
    time_elapsed = inst.timestamp_finish - inst.timestamp_start;
    if (verbose >= QUIET)
        printf("TSP solved in %f\n seconds", time_elapsed);
    if (plot_solution(&inst)) print_error("plot_solution() error");
    //plot_solution_edges(inst.n_edges, inst.nodes, inst.edges) ? print_error("plot_solution_edges() error") : printf("... gnuplot ok\n");

    generate_csv_record(inst.param.name, inst.param.seed, inst.model_type, inst.z_best, time_elapsed, inst.param.run);

    // Free the memory used by the instance
    free_instance(&inst);

    return 0;
}
