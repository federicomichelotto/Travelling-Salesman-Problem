double optimal_solver(instance *inst, instance_TSP *inst_tsp)
{
    printf("Computing the optimal FSTSP-VDS solution...\n");
    // Open CPLEX model
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "FSTSP-VDS");
    // get timestamp
    inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_start) : getTimeStamp(&inst->timestamp_start);

    // build model
    build_model(env, lp, inst, inst_tsp->z_best);
    // path for LP and log files
    char lp_path[1000];
    sprintf(lp_path, "../output/%s/seed_%d/run_%d/model.lp", inst->instance_name, inst->param.seed, inst->param.run);
    char log_path[1000];
    sprintf(log_path, "../output/%s/seed_%d/run_%d/log.txt", inst->instance_name, inst->param.seed, inst->param.run);

    if (CPXwriteprob(env, lp, lp_path, NULL)) // write LP model
        print_error("Error in CPXwriteprob() (ftsp-vds)\n");
    if (CPXsetlogfilename(env, log_path, "w")) // set log file path
        print_error("Error in CPXsetlogfilename() (ftsp-vds)\n");
    inst->cols = CPXgetnumcols(env, lp);
    inst->best_sol = (double *)calloc(inst->cols, sizeof(double));

    post_tsp_sol(env, lp, inst_tsp, inst);

    // ****** CPLEX's parameter setting ******
    if (CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->param.seed)) // Set seed
        print_error("CPX_PARAM_RANDOMSEED error");

    if (inst->param.ticks)
    {
        if (CPXsetdblparam(env, CPX_PARAM_DETTILIM, inst->time_limit))
            print_error("CPX_PARAM_DETTILIM error");
    }
    else
    {
        if (CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit))
            print_error("CPX_PARAM_TILIM error");
    }

    if (CPXsetintparam(env, CPXPARAM_Parallel, CPX_PARALLEL_OPPORTUNISTIC)) // Set opportunistic mode
        print_error("CPXPARAM_Parallel error");

    // if (CPXsetintparam(env, CPX_PARAM_THREADS, 1)) // Set one thread
    //     print_error("CPX_PARAM_THREADS error");

    if (CPXsetintparam(env, CPX_PARAM_CLONELOG, -1)) // CPLEX does not clone log files. (off)
        print_error("CPXPARAM_Output_CloneLog error");

    // CPLEX's precision setting
    if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) // very important if big-M is present
        print_error("CPX_PARAM_EPINT error");
    if (CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9))
        print_error("CPX_PARAM_EPRHS error");
    if (CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-5)) // abort Cplex when relative gap below this value
        print_error("CPX_PARAM_EPGAP error");

    // if (inst->model_type == 11)
    // { // callback method
    //     CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION;
    //     if (CPXcallbacksetfunc(env, lp, contextid, callback_driver, inst))
    //         print_error("CPXcallbacksetfunc() error");
    // }

    if (CPXmipopt(env, lp))
        print_error("CPXmipopt() error");

    // solution status of the problem
    int lpstat = CPXgetstat(env, lp);
    printf("CPLEX status: %d\n", lpstat);
    if (lpstat == 108)
        print_error("Time limit exceeded; no integer solution found.");

    // Use the optimal solution found by CPLEX
    if (CPXgetx(env, lp, inst->best_sol, 0, inst->cols - 1))
        print_error("CPXgetx() error");
    CPXgetobjval(env, lp, &inst->z_best);      // Best objective value
    CPXgetbestobjval(env, lp, &inst->best_lb); // Best lower bound

    //printf("\nSOLUTION -----------------------------------------------\n");
    //printf("\nRUNNING : %s\n", optimal_model_full_name[inst->model_type]);

    // directed graph
    gather_solution(inst, inst->best_sol);

    printf("\nObjective value: %lf\n", inst->z_best);
    printf("Lower bound: %lf\n", inst->best_lb);

    // get timestamp
    inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_finish) : getTimeStamp(&inst->timestamp_finish);

    // Plot optimal solution
    save_and_plot_solution(inst, 0);

    printf("\nTRUCK EDGES:\n");
    printf("  [i]    [j]   [departure time]     [arrival time]     [waiting time at node j]\n");
    {
        int node = 0;
        int n_edges = 0;
        while (node != inst->dimension - 1)
        {
            int succ = -1;
            double max = 0.0;
            for (int j = 1; j < inst->dimension; j++)
            {
                if (node == j)
                    continue;
                if (inst->best_sol[xTruck_pos(node, j, inst)] > max) // succ[i] != node ??
                {
                    max = inst->best_sol[xTruck_pos(node, j, inst)];
                    succ = j;
                }
            }
            if (succ == -1)
                print_error("Number of selected truck edges too small.");
            printf("  %2d  -> %2d | %16.6f   %16.6f        %16.6f \n", node, succ, inst->best_sol[a_pos(node, inst)], inst->best_sol[a_pos(node, inst)] + inst->truck_times[node][succ], fabs(inst->best_sol[a_pos(succ, inst)] - inst->best_sol[a_pos(node, inst)] - inst->truck_times[node][succ]));
            node = succ;
            n_edges++;
            if (n_edges > inst->dimension)
                print_error("Truck path is not a tour.");
        }
    }

    // printf("\nTRUCK EDGES:\n");
    // printf("  [i]    [j]   [departure time]     [arrival time]     [waiting time at node j]\n");
    // for (int i = 0; i < inst->dimension; i++)
    // {
    //     int max_i = -1;
    //     int max_j = -1;
    //     double max = 0.0;
    //     for (int j = 0; j < inst->dimension; j++)
    //     {
    //         if (inst->best_sol[xTruck_pos(i, j, inst)] > max) // succ[i] != node ??
    //         {
    //             max = inst->best_sol[xDrone_pos(i, j, inst)];
    //             max_i = i;
    //             max_j = j;
    //         }
    //     }
    //     if (max_i == -1)
    //         continue;
    //     printf("  %2d  -> %2d | %16.6f   %16.6f        %16.6f \n", max_i, max_j, inst->best_sol[a_pos(max_i, inst)], inst->best_sol[a_pos(max_i, inst)] + inst->truck_times[max_i][max_j], fabs(inst->best_sol[a_pos(max_j, inst)] - inst->best_sol[a_pos(max_i, inst)] - inst->truck_times[max_i][max_j]));
    // }

    // {
    //     int flag = 0;
    //     printf("\nDRONE EDGES:\n");
    //     printf("  [i]    [j]   [departure time]     [arrival time] \n");
    //     int node = 0;
    //     int n_edges = 0;
    //     while (node != inst->dimension - 1)
    //     {
    //         int succ = -1;
    //         double max = 0.0;
    //         for (int j = 1; j < inst->dimension; j++)
    //         {
    //             if (node == j)
    //                 continue;
    //             if (inst->best_sol[xDrone_pos(node, j, inst)] > max) // succ[i] != node ??
    //             {
    //                 max = inst->best_sol[xDrone_pos(node, j, inst)];
    //                 succ = j;
    //             }
    //         }
    //         if (succ == -1)
    //             print_error("Number of selected truck edges too small.");
    //         // if (inst->best_sol[yD_pos(succ, inst)] > 0.5)
    //         //     printf("  %2d  -> %2d | %16.6f            - \n", node, succ, inst->best_sol[a_pos(node, inst)]);
    //         // else if (inst->best_sol[yD_pos(node, inst)] > 0.5)
    //         //     printf("  %2d  -> %2d |          -         %16.6f  \n", node, succ, inst->best_sol[a_pos(succ, inst)]);
    //         // else
    //         printf("  %2d  -> %2d | %16.6f  %16.6f \n", node, succ, inst->best_sol[a_pos(node, inst)], inst->best_sol[a_pos(succ, inst)]);
    //         node = succ;
    //         n_edges++;
    //         if (n_edges > inst->dimension)
    //             print_error("Drone path is not a tour.");
    //     }
    // }

    {
        printf("\nDRONE EDGES:\n");
        printf("  [i]    [j]   [departure time]     [arrival time] \n");

        for (int i = 0; i < inst->dimension; i++)
        {

            for (int j = 0; j < inst->dimension; j++)
            {
                if (inst->best_sol[xDrone_pos(i, j, inst)] > 0.5)
                    printf("  %2d  -> %2d | %16.6f  %16.6f \n", i, j, inst->best_sol[a_pos(i, inst)], inst->best_sol[a_pos(j, inst)]);
            }
        }
    }

    {
        int flag = 0;
        printf("\nyT:");
        for (int i = 0; i < inst->dimension; i++)
        {
            if (inst->best_sol[yT_pos(i, inst)] > 0.5) // succ[i] != node ??
            {
                if (!flag)
                {
                    flag = 1;
                    printf(" %d", i);
                }
                else
                    printf(", %d", i);
            }
        }
    }

    {
        int flag = 0;
        printf("\nyD:");
        for (int i = 0; i < inst->dimension; i++)
        {
            if (inst->best_sol[yD_pos(i, inst)] > 0.5) // succ[i] != node ??
            {
                if (!flag)
                {
                    flag = 1;
                    printf(" %d", i);
                }
                else
                    printf(", %d", i);
            }
        }
    }

    {
        int flag = 0;
        printf("\nyC:");
        for (int i = 0; i < inst->dimension; i++)
        {
            if (inst->best_sol[yC_pos(i, inst)] > 0.5) // succ[i] != node ??
            {
                if (!flag)
                {
                    flag = 1;
                    printf(" %d", i);
                }
                else
                    printf(", %d", i);
            }
        }
    }

    printf("\n\nDRONE LEGS:\n");
    printf("  [i]    [j]    [k]      [travel time]        [min]            [max] \n");
    for (int j = 1; j < inst->dimension - 1; j++)
    {
        if (inst->best_sol[yD_pos(j, inst) < 0.5])
            continue;
        int max_i = -1;
        int max_k = -1;
        double max = 0.0;
        for (int i = 0; i < inst->dimension; i++)
        {
            for (int k = 0; k < inst->dimension; k++)
            {
                if (inst->best_sol[z_pos(i, j, k, inst)] > 0.5)
                {
                    max_i = i;
                    max_k = k;
                    max = inst->best_sol[z_pos(i, j, k, inst)];
                }
            }
        }
        if (max_i == -1)
            continue;
        printf("  %2d --> %2d --> %2d | %16.6f ", max_i, j, max_k, inst->best_sol[a_pos(max_k, inst)] - inst->best_sol[a_pos(max_i, inst)]);
        printf("%16.6f %16.6f\n", inst->min_time_drone[max_i][j][max_k], inst->max_time_drone[max_i][j][max_k]);
        double time_ijk = inst->best_sol[a_pos(max_k, inst)] - inst->best_sol[a_pos(max_i, inst)];
        if (time_ijk < inst->min_feas_time_drone[max_i][j][max_k] || time_ijk > inst->max_feas_time_drone[max_i][j][max_k]){
            printf("*** drone leg infeasible *** \n");
        }
    }

    printf("\nArrival Times:\n");
    for (int i = 0; i < inst->dimension; i++)
    {
        printf(" %d = %f \n", i, inst->best_sol[a_pos(i, inst)]);
    }

    // printf("\nu:\n");
    // for (int i = 0; i < inst->dimension; i++)
    // {
    //     printf(" %d = %f \n", i, inst->best_sol[u_pos(i, inst)]);
    // }

    // Free and close CPLEX model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    return 0;
}
