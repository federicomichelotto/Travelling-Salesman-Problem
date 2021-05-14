#include "../include/utils.h"

int main(int argc, char **argv)
{

    double time_elapsed;
    instance inst; // current instance of the problem

    initialize_instance(&inst);
    parse_command_line(argc, argv, &inst);

    parse_instance(&inst);
    print_instance(&inst);

    switch (inst.param.solver)
    {
    case 0:
        optimal_solver(&inst) ? print_error("optimal_solver() ") : print_message("All went good inside optimal_solver()");
        break;
    case 1:
        math_solver(&inst) ? print_error("math_solver() ") : print_message("All went good inside math_solver()");
        break;
    case 2:
        heuristic_solver(&inst) ? print_error("heuristic_solver() ") : print_message("All went good inside heuristic_solver()");
        break;
    default:
        print_error("No implemented solver selected");
    }
    save_and_plot_solution(&inst, -1);

    time_elapsed = inst.timestamp_finish - inst.timestamp_start;
    if (inst.param.verbose >= QUIET)
        inst.param.ticks ? printf("TSP solved in %f ticks\n", time_elapsed) : printf("TSP solved in %f seconds\n", time_elapsed);

    switch (inst.param.solver)
    {
    case 0:
        generate_csv_record(inst.param.name, inst.param.seed, optimal_model_name[inst.model_type], inst.z_best, time_elapsed, inst.param.run);
        break;
    case 1:
        generate_csv_record(inst.param.name, inst.param.seed, math_model_name[inst.model_type], inst.z_best, time_elapsed, inst.param.run);
        break;
    case 2:
        generate_csv_record(inst.param.name, inst.param.seed, heuristic_model_name[inst.model_type], inst.z_best, time_elapsed, inst.param.run);
        break;
    default:
        print_error("No implemented solver selected");
    }

    // Free the memory used by the instance
    free_instance(&inst);
    return 0;
}
