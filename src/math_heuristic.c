int math_solver(instance *inst)
{
    // Open CPLEX model
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "MATH TSP");

    // get timestamp
    inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_start) : getTimeStamp(&inst->timestamp_start);

    build_model(env, lp, inst);
    inst->cols = CPXgetnumcols(env, lp);
    inst->best_sol = (double *)calloc(inst->cols, sizeof(double));

    char path[1000];
    if (generate_path(path, "output", "model", math_model_name[inst->model_type], inst->param.name, inst->param.seed,
                      "log"))
        print_error("Unable to generate path");

    // CPLEX's parameter setting
    CPXsetlogfilename(env, path, "w");                               // Save log
    if (CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->param.seed)) // Set seed
        print_error("CPX_PARAM_RANDOMSEED error");

    if (CPXsetintparam(env, CPXPARAM_Parallel, CPX_PARALLEL_OPPORTUNISTIC)) // Set opportunistic mode
        print_error("CPXPARAM_Parallel error");

    if (CPXsetintparam(env, CPX_PARAM_CLONELOG, -1)) // CPLEX does not clone log files. (off)
        print_error("CPXPARAM_Output_CloneLog error");

    // CPLEX's precision setting
    if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) // very important if big-M is present
        print_error("CPX_PARAM_EPINT error");
    if (CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9))
        print_error("CPX_PARAM_EPRHS error");
    if (CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-5)) // abort Cplex when relative gap below this value
        print_error("CPX_PARAM_EPGAP error");

    CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION;
    if (CPXcallbacksetfunc(env, lp, contextid, callback_driver, inst))
        print_error("CPXcallbacksetfunc() error");

    // node limit = 0
    if (CPXsetintparam(env, CPX_PARAM_NODELIM, 0))
        print_error("CPX_PARAM_NODELIM error");
    // we don't want to spend all the available time for just one run!
    double time_limit_first_it = inst->time_limit / 20;
    if (inst->param.ticks)
    {
        if (CPXsetdblparam(env, CPX_PARAM_DETTILIM, time_limit_first_it))
            print_error("CPX_PARAM_DETTILIM error");
    }
    else
    {
        if (CPXsetdblparam(env, CPX_PARAM_TILIM, time_limit_first_it))
            print_error("CPX_PARAM_TILIM error");
    }
    if (!inst->param.ticks)
        printf("First iteration: time limit = %f seconds\n", time_limit_first_it);
    else
        printf("First iteration: time limit = %f ticks\n", time_limit_first_it);

    // initial solution
    if (CPXmipopt(env, lp))
        print_error("CPXmipopt() error");

    // solution status of the problem
    int lpstat = CPXgetstat(env, lp);
    printf("\tCPLEX status: %d\n", lpstat);

    if (lpstat == 108)
        print_error("Time limit exceeded; no integer solution");

    // Use the optimal solution found by CPLEX
    if (CPXgetx(env, lp, inst->best_sol, 0, inst->cols - 1))
        print_error("CPXgetx() error");

    CPXgetobjval(env, lp, &inst->z_best);      // Best objective value
    CPXgetbestobjval(env, lp, &inst->best_lb); // Best lower bound
    printf("\tIncumbent: %f\n\n", inst->z_best);
    gather_solution(inst, inst->best_sol, 0);
    save_and_plot_solution(inst, 0);

    printf("\nSOLUTION -----------------------------------------------\n");
    printf("\nRUNNING : %s\n", math_model_full_name[inst->model_type]);

    // hard-fixing heuristic
    if (inst->model_type == 0)
    {
        if (inst->param.verbose >= NORMAL)
            printf("Initial incumbent: %f\n\n", inst->z_best);
        if (CPXsetintparam(env, CPX_PARAM_NODELIM, 2100000000))
            print_error("CPX_PARAM_NODELIM error");
        hard_fixing_heuristic(env, lp, inst, (int)inst->time_limit / 20, 0.8);
    }
    // soft-fixing heuristic
    else if (inst->model_type == 1)
    {
        if (inst->param.verbose >= NORMAL)
            printf("Initial incumbent: %f\n", inst->z_best);
        // reset the node limit to default
        if (CPXsetintparam(env, CPX_PARAM_NODELIM, 2100000000))
            print_error("CPX_PARAM_NODELIM error");
        soft_fixing_heuristic(env, lp, inst, (int)inst->time_limit / 20);
    }

    // if (inst->n_edges != inst->dimension)
    //     print_error("not a tour.");

    printf("\nObjective value: %lf\n", inst->z_best);
    printf("Lower bound: %lf\n", inst->best_lb);

    // get timestamp
    inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_finish) : getTimeStamp(&inst->timestamp_finish);

    // Plot optimal solution
    gather_solution(inst, inst->best_sol, 0);
    save_and_plot_solution(inst, -1);

    // Free and close CPLEX model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    return 0;
}

void hard_fixing_heuristic(CPXENVptr env, CPXLPptr lp, instance *inst, int time_limit_iter, double fix_ratio)
{
    int iter = 1;
    double gap;
    CPXgetdblparam(env, CPX_PARAM_EPGAP, &gap);
    while (1)
    {
        // update time left
        double ts_current;
        inst->param.ticks ? CPXgetdettime(env, &ts_current) : getTimeStamp(&ts_current);
        double time_left = inst->time_limit - (ts_current - inst->timestamp_start);
        if (time_left < inst->param.time_threshold)
            return;
        if (time_left < time_limit_iter)
            CPXsetdblparam(env, CPX_PARAM_TILIM, time_left);
        if (inst->param.verbose >= DEBUG)
            printf("*** time left = %f\n", time_left);
        if (iter % 100 == 0)
        {
            printf("\nIteration %d, time left %f:\n", iter + 1, time_left);
        }
        // allocate two arrays with size ncols
        int *indices = (int *)malloc(inst->cols * sizeof(int));
        double *values = (double *)malloc(inst->cols * sizeof(double));
        char senses[inst->cols]; // we only need to change the lower bound of our variables
        int nedges = inst->dimension * (inst->dimension - 1) / 2;
        for (int k = 0; k < nedges; k++)
        {
            indices[k] = k;
            values[k] = 0.0;
            senses[k] = 'L';
            // random selection a subset of arcs of the best solution found so far
            //printf("*** best_sol[%d] = %f\n",k,inst->best_sol[k]);
            if (inst->best_sol[k] > 0.5)
                if ((rand() % 100) + 1 < fix_ratio * 100)
                    values[k] = 1.0;
        }

        // change the lower bounds
        if (CPXchgbds(env, lp, nedges, indices, senses, values))
            print_error("CPXchgbds error");

        // solve with the new constraints
        if (CPXmipopt(env, lp))
            print_error("CPXmipopt() error");

        // retrieve the incumbent of the current solution and the best known lower bound
        double current_incumbent, lb_incumbent;
        CPXgetobjval(env, lp, &current_incumbent);
        CPXgetbestobjval(env, lp, &lb_incumbent);
        if (inst->param.verbose >= DEBUG)
            printf("*** current_incumbent = %f\n", current_incumbent);

        // check if the current solution is better than the best so far
        if (current_incumbent < inst->z_best)
        {
            // update best incumbent
            inst->z_best = current_incumbent;
            inst->best_lb = lb_incumbent;
            // update arcs' selection
            int status = CPXgetx(env, lp, inst->best_sol, 0, inst->cols - 1);
            if (status)
                print_error_status("Failed to obtain the values in hard_fixing_heuristic method", status);
            if (inst->param.verbose >= NORMAL)
            {
                if (!inst->param.ticks)
                    printf("New incumbent: %f, time left = %f seconds\n", inst->z_best, time_left);
                else
                    printf("New incumbent: %f, time left = %f ticks\n", inst->z_best, time_left);
            }
            gather_solution(inst, inst->best_sol, 0);
            save_and_plot_solution(inst, iter);

            int beg = 0;
            if (CPXaddmipstarts(env, lp, 1, nedges, &beg, indices, values, CPX_MIPSTART_AUTO, NULL))
                print_error("CPXaddmipstarts error");
        }

        iter++;
        free(indices);
        free(values);
    }
}

void soft_fixing_heuristic(CPXENVptr env, CPXLPptr lp, instance *inst, int time_limit_iter)
{
    int iter = 1;
    int k = 2;
    double gap;
    CPXgetdblparam(env, CPX_PARAM_EPGAP, &gap);
    while (1)
    {
        // update time left
        double ts_current;
        inst->param.ticks ? CPXgetdettime(env, &ts_current) : getTimeStamp(&ts_current);
        double time_left = inst->time_limit - (ts_current - inst->timestamp_start);
        if (time_left < inst->param.time_threshold)
            return;
        if (time_left < time_limit_iter)
            CPXsetdblparam(env, CPX_PARAM_TILIM, time_left);
        if (inst->param.verbose >= DEBUG)
            printf("*** time left = %f\n", time_left);

        double rhs = inst->dimension - k;
        char sense = 'G';
        char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
        rname[0] = (char *)calloc(100, sizeof(char));
        sprintf(rname[0], "soft_fixing(%d)", iter);

        int row = CPXgetnumrows(env, lp); // get the maximum number of row inside the model
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("CPXnewrows error");

        int nedges = inst->dimension * (inst->dimension - 1) / 2;
        int rows[inst->dimension];
        int *indices = (int *)malloc(inst->dimension * sizeof(int));
        double *coeffs = (double *)malloc(inst->dimension * sizeof(double));
        int j = 0; // counter: 0 -> nnodes-1
        for (int i = 0; i < nedges; i++)
        {
            if (inst->best_sol[i] > 0.5)
            {
                rows[j] = row;
                indices[j] = i;
                coeffs[j] = 1.0;
                j++;
            }
        }
        if (CPXchgcoeflist(env, lp, inst->dimension, rows, indices, coeffs))
            print_error("wrong CPXchgcoeflist");

        // solve with the new constraint
        if (CPXmipopt(env, lp))
            print_error("CPXmipopt() error");

        // retrieve the incumbent of the current solution and the best known lower bound
        double current_incumbent, lb_incumbent;
        CPXgetobjval(env, lp, &current_incumbent);
        CPXgetbestobjval(env, lp, &lb_incumbent);
        if (inst->param.verbose >= DEBUG)
            printf("*** k = %d, current_incumbent = %f\n", k, current_incumbent);
        // check if the current solution is better than the best so far
        if (current_incumbent < inst->z_best)
        {
            // update best incumbent
            inst->z_best = current_incumbent;
            inst->best_lb = lb_incumbent;
            // update best sol
            int status = CPXgetx(env, lp, inst->best_sol, 0, inst->cols - 1);
            if (status)
                print_error_status("Failed to obtain the values in soft_fixing_heuristic method", status);
            if (inst->param.verbose >= NORMAL)
            {
                if (!inst->param.ticks)
                    printf("New incumbent: %f, time left = %f seconds\n", inst->z_best, time_left);
                else
                    printf("New incumbent: %f, time left = %f ticks\n", inst->z_best, time_left);
            }
            gather_solution(inst, inst->best_sol, 0);
            save_and_plot_solution(inst, iter);

            int beg = 0;
            if (CPXaddmipstarts(env, lp, 1, inst->dimension, &beg, indices, coeffs, CPX_MIPSTART_AUTO, NULL))
                print_error("CPXaddmipstarts error");
        }
        else
        {
            if (k < inst->dimension && k < 20)
                k++;
            else
            {
                // reached the max neighborhood size without improving! STOP
                printf("The procedure has been stopped before the time limit because it was reached the max neighborhood size without improving!\n");
                return;
            }
        }
        if (CPXdelrows(env, lp, row, row))
            print_error("CPXdelrows error");
        iter++;
        free(rname[0]);
        free(rname);
        free(indices);
        free(coeffs);
    }
}