#include "../include/utils.h"
#include "../include/genetic.h"

int main(int argc, char **argv)
{
    double time_elapsed;
    instance inst;         // FSTSP-VSD instance of the problem
    instance_TSP inst_tsp; // TSP instance of the problem

    initialize_instance(&inst);

    struct timespec timestamp;
    if (clock_gettime(CLOCK_REALTIME, &timestamp) == -1)
        print_error("Error clock_gettime");
    inst.timestamp_start = timestamp.tv_sec + timestamp.tv_nsec * pow(10, -9);

    parse_command_line(argc, argv, &inst);
    parse_instance(&inst);
    print_instance(&inst);

    // create the folders related to this run
    createInstanceFolders(&inst);

    if (!inst.param.heur) // exact method
    {
        // generate a TSP instance from the FSTSP instance
        generate_TSP_instance(&inst, &inst_tsp);
        // compute the optimal TSP solution
        tsp_solver(&inst_tsp);
        // compute the optimal FSTPS-VSD solution
        optimal_solver(&inst, &inst_tsp) ? print_error("optimal_solver() ") : print_message("All went good inside optimal_solver()");
    }
    else // heuristic method
    {

        generate_TSP_instance(&inst, &inst_tsp);
        // compute the optimal TSP solution
        tsp_solver(&inst_tsp);
        int seq[inst.dimension];
        seq[0] = 0;
        seq[inst.dimension - 1] = inst.dimension - 1;
        int i = 1;
        for (int i = 1; i < inst.dimension; i++)
        {
            seq[i] = inst_tsp.succ[seq[i - 1]];
        }

        int truck_seq[inst.dimension];
        int drone_seq[inst.dimension];

        double tsp_score = compute_sequence_min_time(seq, &inst, truck_seq, drone_seq);
        printf("TSP derived score : %f\n", tsp_score);
        save_and_plot_solution_general(&inst, truck_seq, drone_seq, -1);

        //int seq[] = {0, 7, 8, 2, 3, 5, 6, 9, 10, 4, 1, 11};
        double obj = DBL_MAX;
        for (int i = 0; i < 10000; i++)
        {
            //generate_random_sequence(seq, &inst);
            nearest_neighbours(&inst, seq, 2);
            double score = compute_sequence_min_time(seq, &inst, truck_seq, drone_seq);

            //double score2 = compute_score(truck_seq, drone_seq, &inst);
            //double delta = reverse_drone_leg(truck_seq, drone_seq, &inst);
            // if (delta < 0)
            //     print_error("DELTA < 0\n");
            // double score3 = compute_score(truck_seq, drone_seq, &inst);

            // if (score2 - score3 - delta > 1e-4)
            // {
            //     printf("*** ACHTUNG ***\n");
            //     printf("SEQUENCE: ");
            //     for (int i = 0; i < inst.dimension; i++)
            //         printf("%d, ", seq[i]);
            //     printf("\n");
            //     printf("TRUCK SEQUENCE: ");
            //     for (int i = 0; i < inst.dimension; i++)
            //     {
            //         printf("%d, ", truck_seq[i]);
            //         if (truck_seq[i] == inst.dimension - 1)
            //             break;
            //     }
            //     printf("\n");
            //     printf("DRONE SEQUENCE: ");
            //     for (int i = 0; i < inst.dimension; i++)
            //     {
            //         printf("%d, ", drone_seq[i]);
            //         if (drone_seq[i] == inst.dimension - 1)
            //             break;
            //     }
            //     printf("\n");
            //     printf("initial score = %f\n", score);
            //     printf("check score = %f\n", score2);
            //     printf("final score = %f\n", score3);
            //     printf("delta = %f\n", delta);
            //     printf("check delta = %f\n", score2 - score3);
            // }
            double delta = 0.0;
            if (score < obj * 1.1)
            {
                delta = reverse_drone_leg(truck_seq, drone_seq, &inst);
                score -= delta;
            }
            if (score < obj)
            {
                obj = score;
                if (delta > 0)
                    printf("delta > 0\n");
                printf("new obj found : %f\n", score);
                int truck_succ[inst.dimension];
                int drone_succ[inst.dimension];
                convert_seq2succ(&inst, truck_seq, truck_succ);
                convert_seq2succ(&inst, drone_seq, drone_succ);
                save_and_plot_solution_general(&inst, truck_succ, drone_succ, -1);
            }
        }
        
        (tsp_score < obj) ? printf("\nBest score: %f\n", tsp_score) : printf("\nBest score: %f\n", obj);

        // printf("TRUCK SUCCESSOR:\n");
        // for (int i = 0; i < inst.dimension; i++)
        //     printf("%d -> %d\n", i, truck_succ[i]);
        // printf("DRONE SUCCESSOR:\n");
        // for (int i = 0; i < inst.dimension; i++)
        //     printf("%d -> %d\n", i, drone_succ[i]);
    }

    if (clock_gettime(CLOCK_REALTIME, &timestamp) == -1)
        print_error("Error clock_gettime");
    inst.timestamp_finish = timestamp.tv_sec + timestamp.tv_nsec * pow(10, -9);

    time_elapsed = inst.timestamp_finish - inst.timestamp_start;
    inst.param.ticks ? printf("Elapsed time in ticks: %f\n", time_elapsed) : printf("Elapsed time: %f seconds\n", time_elapsed);

    char score_path[] = "../output/scores.csv";
    generate_csv_record(score_path, inst.instance_name, inst.param.seed, inst.model_type, inst.param.run, inst.z_best, time_elapsed, inst.param.ticks);

    // Free the memory used by the instance
    //free_instance(&inst);
    return 0;
}
