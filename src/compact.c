int optimal_solver(instance *inst)
{
    // Open CPLEX model
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");
    // get timestamp

    inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_start) : getTimeStamp(&inst->timestamp_start);

    build_model(env, lp, inst);
    inst->cols = CPXgetnumcols(env, lp);
    inst->best_sol = (double *)calloc(inst->cols, sizeof(double));

    char path[1000];
    if (generate_path(path, "output", "model", optimal_model_name[inst->model_type], inst->param.name, inst->param.seed,
                      "log"))
        print_error("Unable to generate path");

    // CPLEX's parameter setting
    CPXsetlogfilename(env, path, "w");                               // Save log
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

    // CPLEX's precision setting
    if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) // very important if big-M is present
        print_error("CPX_PARAM_EPINT error");
    if (CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9))
        print_error("CPX_PARAM_EPRHS error");
    if (CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-5)) // abort Cplex when relative gap below this value
        print_error("CPX_PARAM_EPGAP error");

    if (inst->model_type == 11)
    { // callback method
        CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION;
        if (CPXcallbacksetfunc(env, lp, contextid, callback_driver, inst))
            print_error("CPXcallbacksetfunc() error");
    }

    if (CPXmipopt(env, lp))
        print_error("CPXmipopt() error");

    // solution status of the problem
    int lpstat = CPXgetstat(env, lp);
    printf("CPLEX status: %d\n", lpstat);
    if (lpstat == 108)
        print_error("Time limit exceeded; no integer solution");

    // Use the optimal solution found by CPLEX
    if (CPXgetx(env, lp, inst->best_sol, 0, inst->cols - 1))
        print_error("CPXgetx() error");

    CPXgetobjval(env, lp, &inst->z_best);      // Best objective value
    CPXgetbestobjval(env, lp, &inst->best_lb); // Best lower bound

    printf("\nSOLUTION -----------------------------------------------\n");
    printf("\nRUNNING : %s\n", optimal_model_full_name[inst->model_type]);

    if (inst->model_type == 0 || inst->model_type == 9 || inst->model_type == 10 ||
        inst->model_type == 11)
    { // undirected graph

        if (inst->model_type == 9 || inst->model_type == 10)
        { // Benders
            benders(env, lp, inst);
            // no need to gather the solution because benders update inst->succ directly

            CPXgetobjval(env, lp, &inst->z_best);      // Best objective value
            CPXgetbestobjval(env, lp, &inst->best_lb); // Best lower bound
        }
        else
        {
            if (CPXgetx(env, lp, inst->best_sol, 0, inst->cols - 1))
                print_error("CPXgetx() error");
            gather_solution(inst, inst->best_sol, 0);
        }
    }
    else
    { // directed graph
        gather_solution(inst, inst->best_sol, 1);
    }

    printf("\nObjective value: %lf\n", inst->z_best);
    printf("Lower bound: %lf\n", inst->best_lb);

    // get timestamp
    inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_finish) : getTimeStamp(&inst->timestamp_finish);

    // Plot optimal solution
    save_and_plot_solution(inst, -1);

    // Free and close CPLEX model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    return 0;
}

void MTZ_static(CPXENVptr env, CPXLPptr lp, instance *inst)
{

    basic_model_directed(env, lp, inst);

    char binary = 'B';  // B => binary variable flag
    char integer = 'I'; // I => integer variable flag

    // cname: columns' names (column = variable)
    char **cname = (char **)calloc(1, sizeof(char *)); // array of strings to store the column names
    cname[0] = (char *)calloc(100, sizeof(char));

    // rname: rows' names (row = constraint)
    char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
    rname[0] = (char *)calloc(100, sizeof(char));

    double M = inst->dimension - 1;

    // Add u-variables one for each node ( u_0 = 0 )
    for (int i = 0; i < inst->dimension; i++)
    {
        sprintf(cname[0], "u(%d)", i + 1);
        double obj = 0.0;
        double lb = 0.0;
        double ub = M;
        if (i == 0)
            ub = 0.0;
        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname))
            print_error(" wrong CPXnewcols on u var.s");
        if (CPXgetnumcols(env, lp) - 1 != upos(i, inst))
            print_error("[position_d] wrong position for u var.s");
    }

    // Add static MTZ constraints: 1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1, for each arc (i,j) not touching node 0
    double rhs = M - 1;
    char sense = 'L'; // L stands for less than or equal
    for (int i = 1; i < inst->dimension; i++)
    {
        for (int j = 1; j < inst->dimension; j++)
        {
            if (i == j)
                continue;
            int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
            sprintf(rname[0], "u_consistency for arc (%d,%d)", i + 1, j + 1);
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
                print_error("wrong CPXnewrows [degree]");
            if (CPXchgcoef(env, lp, row, upos(i, inst), 1.0)) // 1.0 * u_i
                print_error("wrong CPXchgcoef [degree]");
            if (CPXchgcoef(env, lp, row, upos(j, inst), -1.0)) // - 1.0 * u_j
                print_error("wrong CPXchgcoef [degree]");
            if (CPXchgcoef(env, lp, row, xpos_dir(i, j, inst), M)) // M * x_ij
                print_error("wrong CPXchgcoef [degree]");
        }
    }
    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}

void MTZ_static_mod(CPXENVptr env, CPXLPptr lp, instance *inst)
{
    basic_model_directed(env, lp, inst);

    char binary = 'B';  // B => binary variable flag
    char integer = 'I'; // I => integer variable flag

    // cname: columns' names (column = variable)
    char **cname = (char **)calloc(1, sizeof(char *)); // array of strings to store the column names
    cname[0] = (char *)calloc(100, sizeof(char));

    // rname: rows' names (row = constraint)
    char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
    rname[0] = (char *)calloc(100, sizeof(char));

    double M = inst->dimension - 1;
    // Add u-variables one for each node ( u_0 = 0 )
    for (int i = 0; i < inst->dimension; i++)
    {
        sprintf(cname[0], "u(%d)", i + 1);
        double obj = 0.0;
        double lb = 0.0;
        double ub = M;
        if (i == 0)
            ub = 0;
        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname))
            print_error(" wrong CPXnewcols on u var.s");
        if (CPXgetnumcols(env, lp) - 1 != upos(i, inst))
            print_error("[position_d] wrong position for u var.s");
    }

    // Add static MTZ constraints: 1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1, for each arc (i,j) not touching node 0
    double rhs = M - 1;
    char sense = 'L';                         // L stands for less than or equal
    for (int i = 0; i < inst->dimension; i++) // *** mod: including i=0 ***
    {
        for (int j = 1; j < inst->dimension; j++)
        {
            if (i == j)
                continue;
            int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
            sprintf(rname[0], "u_consistency for arc (%d,%d)", i + 1, j + 1);
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
                print_error("wrong CPXnewrows [degree]");
            if (CPXchgcoef(env, lp, row, upos(i, inst), 1.0)) // 1.0 * u_i
                print_error("wrong CPXchgcoef [degree]");
            if (CPXchgcoef(env, lp, row, upos(j, inst), -1.0)) // - 1.0 * u_j
                print_error("wrong CPXchgcoef [degree]");
            if (CPXchgcoef(env, lp, row, xpos_dir(i, j, inst), M)) // M * x_ij
                print_error("wrong CPXchgcoef [degree]");
        }
    }
    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}

void MTZ_lazy(CPXENVptr env, CPXLPptr lp, instance *inst)
{
    basic_model_directed(env, lp, inst);
    char binary = 'B';  // B => binary variable flag
    char integer = 'I'; // I => integer variable flag

    // cname: columns' names (column = variable)
    char **cname = (char **)calloc(1, sizeof(char *)); // array of strings to store the column names
    cname[0] = (char *)calloc(100, sizeof(char));

    // rname: rows' names (row = constraint)
    char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
    rname[0] = (char *)calloc(100, sizeof(char));

    // Add u-variables one for each node ( u_0 = 0 )
    for (int i = 0; i < inst->dimension; i++)
    {
        sprintf(cname[0], "u(%d)", i + 1);
        double obj = 0.0;
        double lb = 0.0;
        double ub = inst->dimension - 1;
        if (i == 0)
            ub = 0;
        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname))
            print_error(" wrong CPXnewcols on u var.s");
        if (CPXgetnumcols(env, lp) - 1 != upos(i, inst))
            print_error("[position_d] wrong position for u var.s");
    }

    int izero = 0;
    int index[3];
    double value[3];

    // add lazy constraints  1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1, for each arc (i,j) not touching node 0
    double big_M = inst->dimension - 1.0;
    double rhs = big_M - 1.0;
    char sense = 'L';
    int nnz = 3;
    for (int i = 1; i < inst->dimension; i++) // excluding node 0
    {
        for (int j = 1; j < inst->dimension; j++) // excluding node 0
        {
            if (i == j)
                continue;
            sprintf(rname[0], "u_consistency for arc (%d,%d)", i + 1, j + 1);
            index[0] = upos(i, inst);
            value[0] = 1.0;
            index[1] = upos(j, inst);
            value[1] = -1.0;
            index[2] = xpos_dir(i, j, inst);
            value[2] = big_M;
            if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, rname))
                print_error("wrong CPXlazyconstraints() for u-consistency");
        }
    }

    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}

void MTZ_lazy_sec(CPXENVptr env, CPXLPptr lp, instance *inst)
{
    basic_model_directed(env, lp, inst);
    char binary = 'B';  // B => binary variable flag
    char integer = 'I'; // I => integer variable flag

    // cname: columns' names (column = variable)
    char **cname = (char **)calloc(1, sizeof(char *)); // array of strings to store the column names
    cname[0] = (char *)calloc(100, sizeof(char));

    // rname: rows' names (row = constraint)
    char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
    rname[0] = (char *)calloc(100, sizeof(char));

    // Add u-variables one for each node ( u_0 = 0 )
    double big_M = inst->dimension - 1.0;
    for (int i = 0; i < inst->dimension; i++)
    {
        sprintf(cname[0], "u(%d)", i + 1);
        double obj = 0.0;
        double lb = 0.0;
        double ub = big_M;
        if (i == 0)
            ub = 0;
        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname))
            print_error(" wrong CPXnewcols on u var.s");
        if (CPXgetnumcols(env, lp) - 1 != upos(i, inst))
            print_error("[position_d] wrong position for u var.s");
    }

    int izero = 0;
    int index[3];
    double value[3];

    // add lazy constraints  1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1, for each arc (i,j) not touching node 0
    double rhs = big_M - 1.0;
    char sense = 'L';
    int nnz = 3;
    for (int i = 1; i < inst->dimension; i++) // excluding node 0
    {
        for (int j = 1; j < inst->dimension; j++) // excluding node 0
        {
            if (i == j)
                continue;
            sprintf(rname[0], "u_consistency for arc (%d,%d)", i + 1, j + 1);
            index[0] = upos(i, inst);
            value[0] = 1.0;
            index[1] = upos(j, inst);
            value[1] = -1.0;
            index[2] = xpos_dir(i, j, inst);
            value[2] = big_M;
            if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, rname))
                print_error("wrong CPXlazyconstraints() for u-consistency");
        }
    }

    // Add static 2-SEC contraints: x(i, j) + x(j, i) <= 1 for every i < j
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = i + 1; j < inst->dimension; j++)
        {

            int lastrow = CPXgetnumrows(env, lp);
            double rhs = 1.0;
            char sense = 'L';

            sprintf(cname[0], "2-SEC(%d, %d)", i + 1, j + 1);
            int *beg = (int *)calloc(2, sizeof(int));
            int *ind = (int *)calloc(2, sizeof(int));
            double *val = (double *)calloc(2, sizeof(double));

            int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
                print_error("wrong CPXnewrows [degree]");
            if (CPXchgcoef(env, lp, row, xpos_dir(i, j, inst), 1.0)) // 1.0 * x_ij
                print_error("wrong CPXchgcoef [degree]");
            if (CPXchgcoef(env, lp, row, xpos_dir(j, i, inst), 1.0)) // 1.0 * x_ji
                print_error("wrong CPXchgcoef [degree]");
            free(beg);
            free(ind);
        }
    }

    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}

void GG(CPXENVptr env, CPXLPptr lp, instance *inst)
{
    basic_model_directed(env, lp, inst);

    char binary = 'B';  // B => binary variable flag
    char integer = 'I'; // I => integer variable flag

    // cname: columns' names (column = variable)
    char **cname = (char **)calloc(1, sizeof(char *)); // array of strings to store the column names
    cname[0] = (char *)calloc(100, sizeof(char));

    // rname: rows' names (row = constraint)
    char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
    rname[0] = (char *)calloc(100, sizeof(char));

    // Add y-variables one for each arc (i,j) with i!=j and i,j > 0
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            sprintf(cname[0], "y(%d,%d)", i + 1, j + 1);
            double obj = 0.0;
            double lb = 0.0;
            double ub = inst->dimension - 2;
            if (i == 0)
                ub = inst->dimension - 1;
            if (i == j || j == 0)
                ub = 0;
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname))
                print_error(" wrong CPXnewcols on y var.s");
            if (CPXgetnumcols(env, lp) - 1 != ypos(i, j, inst))
                print_error("[position_d] wrong position for y var.s");
        }
    }

    // Add in-flow out-flow differential constraints for nodes h > 0: sum_i (1.0 * y_ih) + sum_i (- 1.0 * y_hi) = 1
    for (int h = 1; h < inst->dimension; h++) // exludes node 0
    {
        double rhs = 1.0;
        char sense = 'E';                 // E stands for equal
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        sprintf(rname[0], "in_flow out_flow node (%d)", h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [degree]");
        for (int i = 0; i < inst->dimension; i++)
        {
            if (h == i)
                continue;
            if (CPXchgcoef(env, lp, row, ypos(i, h, inst), 1.0)) // 1.0 * y_ih
                print_error("wrong CPXchgcoef [degree]");
            if (CPXchgcoef(env, lp, row, ypos(h, i, inst), -1.0)) // - 1.0 * y_hi
                print_error("wrong CPXchgcoef [degree]");
        }
    }

    // Add out-flow constraints for node 0: 1 * y_0j - (nnodes -1) * x_0j = 0, for each arc (0,j)
    for (int j = 1; j < inst->dimension; j++) // exclude  arc (0,0)
    {
        double rhs = 0;
        char sense = 'E';                 // E stands for equal
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        sprintf(rname[0], "out_flow(1,%d)", j + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [degree]");
        if (CPXchgcoef(env, lp, row, ypos(0, j, inst), 1.0)) // 1.0 * y_0j
            print_error("wrong CPXchgcoef [degree]");

        if (CPXchgcoef(env, lp, row, xpos_dir(0, j, inst), 1 - inst->dimension)) // - (nnodes -1) * x_0j
            print_error("wrong CPXchgcoef [degree]");
    }

    // linking constraints: y_ij <= (n-2)x_ij for each i,j > 0 and i!=j
    for (int i = 1; i < inst->dimension; i++) // exludes node 0
    {
        for (int j = 1; j < inst->dimension; j++)
        {
            if (i == j)
                continue;
            double rhs = 0.0;
            char sense = 'L';                 // E stands for equal
            int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
            sprintf(rname[0], "linking constraints (%d,%d)", i + 1, j + 1);
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
                print_error("wrong CPXnewrows [degree]");

            if (CPXchgcoef(env, lp, row, ypos(i, j, inst), 1.0)) // 1.0 * y_ij
                print_error("wrong CPXchgcoef [degree]");
            if (CPXchgcoef(env, lp, row, xpos_dir(i, j, inst), 2 - inst->dimension)) // - (2 - nnodes) * x_ij
                print_error("wrong CPXchgcoef [degree]");
        }
    }

    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}

void GG_lazy(CPXENVptr env, CPXLPptr lp, instance *inst)
{
    basic_model_directed(env, lp, inst);

    char binary = 'B';  // B => binary variable flag
    char integer = 'I'; // I => integer variable flag

    // cname: columns' names (column = variable)
    char **cname = (char **)calloc(1, sizeof(char *)); // array of strings to store the column names
    cname[0] = (char *)calloc(100, sizeof(char));

    // rname: rows' names (row = constraint)
    char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
    rname[0] = (char *)calloc(100, sizeof(char));

    // Add y-variables one for each arc (i,j) with i!=j and i,j > 0
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            sprintf(cname[0], "y(%d,%d)", i + 1, j + 1);
            double obj = 0.0;
            double lb = 0.0;
            double ub = inst->dimension - 2;
            if (i == 0)
                ub = inst->dimension - 1;
            if (i == j || j == 0)
                ub = 0;
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname))
                print_error(" wrong CPXnewcols on y var.s");
            if (CPXgetnumcols(env, lp) - 1 != ypos(i, j, inst))
                print_error("[position_d] wrong position for y var.s");
        }
    }

    // Add in-flow out-flow differential constraints for nodes h > 0: sum_i (1.0 * y_ih) + sum_i (- 1.0 * y_hi) = 1
    for (int h = 1; h < inst->dimension; h++) // exludes node 0
    {
        double rhs = 1.0;
        char sense = 'E';                 // E stands for equal
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        sprintf(rname[0], "in_flow out_flow node (%d)", h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [degree]");
        for (int i = 0; i < inst->dimension; i++)
        {
            if (h == i)
                continue;
            if (CPXchgcoef(env, lp, row, ypos(i, h, inst), 1.0)) // 1.0 * y_ih
                print_error("wrong CPXchgcoef [degree]");
            if (CPXchgcoef(env, lp, row, ypos(h, i, inst), -1.0)) // - 1.0 * y_hi
                print_error("wrong CPXchgcoef [degree]");
        }
    }

    // Add out-flow constraints for node 0: 1 * y_0j - (nnodes -1) * x_0j = 0, for each arc (0,j)
    for (int j = 1; j < inst->dimension; j++) // exclude  arc (0,0)
    {
        double rhs = 0;
        char sense = 'E';                 // E stands for equal
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        sprintf(rname[0], "out_flow(1,%d)", j + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [degree]");
        if (CPXchgcoef(env, lp, row, ypos(0, j, inst), 1.0)) // 1.0 * y_0j
            print_error("wrong CPXchgcoef [degree]");

        if (CPXchgcoef(env, lp, row, xpos_dir(0, j, inst), 1 - inst->dimension)) // - (nnodes -1) * x_0j
            print_error("wrong CPXchgcoef [degree]");
    }

    int izero = 0;
    int index[2];
    double value[2];

    // add lazy linking constraints  y_ij <= (n-2)x_ij for each i,j > 0 and i!=j
    double rhs = 0.0;
    char sense = 'L';
    int nnz = 2;
    for (int i = 1; i < inst->dimension; i++) // excluding node 0
    {
        for (int j = 1; j < inst->dimension; j++) // excluding node 0
        {
            if (i == j)
                continue;
            sprintf(rname[0], "lazy linking constraints (%d,%d)", i + 1, j + 1);
            index[0] = ypos(i, j, inst);
            value[0] = 1.0;
            index[1] = xpos_dir(i, j, inst);
            value[1] = 2 - inst->dimension;
            if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, rname))
                print_error("wrong CPXlazyconstraints() for u-consistency");
        }
    }

    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}

void GG_lazy_sec(CPXENVptr env, CPXLPptr lp, instance *inst)
{
    basic_model_directed(env, lp, inst);

    char binary = 'B';  // B => binary variable flag
    char integer = 'I'; // I => integer variable flag

    // cname: columns' names (column = variable)
    char **cname = (char **)calloc(1, sizeof(char *)); // array of strings to store the column names
    cname[0] = (char *)calloc(100, sizeof(char));

    // rname: rows' names (row = constraint)
    char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
    rname[0] = (char *)calloc(100, sizeof(char));

    // Add y-variables one for each arc (i,j) with i!=j and i,j > 0
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            sprintf(cname[0], "y(%d,%d)", i + 1, j + 1);
            double obj = 0.0;
            double lb = 0.0;
            double ub = inst->dimension - 2;
            if (i == 0)
                ub = inst->dimension - 1;
            if (i == j || j == 0)
                ub = 0;
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname))
                print_error(" wrong CPXnewcols on y var.s");
            if (CPXgetnumcols(env, lp) - 1 != ypos(i, j, inst))
                print_error("[position_d] wrong position for y var.s");
        }
    }

    // Add in-flow out-flow differential constraints for nodes h > 0: sum_i (1.0 * y_ih) + sum_i (- 1.0 * y_hi) = 1
    for (int h = 1; h < inst->dimension; h++) // exludes node 0
    {
        double rhs = 1.0;
        char sense = 'E';                 // E stands for equal
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        sprintf(rname[0], "in_flow out_flow node (%d)", h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [degree]");
        for (int i = 0; i < inst->dimension; i++)
        {
            if (h == i)
                continue;
            if (CPXchgcoef(env, lp, row, ypos(i, h, inst), 1.0)) // 1.0 * y_ih
                print_error("wrong CPXchgcoef [degree]");
            if (CPXchgcoef(env, lp, row, ypos(h, i, inst), -1.0)) // - 1.0 * y_hi
                print_error("wrong CPXchgcoef [degree]");
        }
    }

    // Add out-flow constraints for node 0: 1 * y_0j - (nnodes -1) * x_0j = 0, for each arc (0,j)
    for (int j = 1; j < inst->dimension; j++) // exclude  arc (0,0)
    {
        double rhs = 0;
        char sense = 'E';                 // E stands for equal
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        sprintf(rname[0], "out_flow(1,%d)", j + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [degree]");
        if (CPXchgcoef(env, lp, row, ypos(0, j, inst), 1.0)) // 1.0 * y_0j
            print_error("wrong CPXchgcoef [degree]");

        if (CPXchgcoef(env, lp, row, xpos_dir(0, j, inst), 1 - inst->dimension)) // - (nnodes -1) * x_0j
            print_error("wrong CPXchgcoef [degree]");
    }

    int izero = 0;
    int index[2];
    double value[2];

    // add lazy linking constraints  y_ij <= (n-2)x_ij for each i,j > 0 and i!=j
    double rhs = 0.0;
    char sense = 'L';
    int nnz = 2;
    for (int i = 1; i < inst->dimension; i++) // excluding node 0
    {
        for (int j = 1; j < inst->dimension; j++) // excluding node 0
        {
            if (i == j)
                continue;
            sprintf(rname[0], "lazy linking constraints (%d,%d)", i + 1, j + 1);
            index[0] = ypos(i, j, inst);
            value[0] = 1.0;
            index[1] = xpos_dir(i, j, inst);
            value[1] = 2 - inst->dimension;
            if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, rname))
                print_error("wrong CPXlazyconstraints() for u-consistency");
        }
    }

    // Add static 2-SEC contraints: x(i, j) + x(j, i) <= 1 for every i < j
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = i + 1; j < inst->dimension; j++)
        {

            int lastrow = CPXgetnumrows(env, lp);
            double rhs = 1.0;
            char sense = 'L';

            sprintf(cname[0], "2-SEC(%d, %d)", i + 1, j + 1);
            int *beg = (int *)calloc(2, sizeof(int));
            int *ind = (int *)calloc(2, sizeof(int));
            double *val = (double *)calloc(2, sizeof(double));

            int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
                print_error("wrong CPXnewrows [degree]");
            if (CPXchgcoef(env, lp, row, xpos_dir(i, j, inst), 1.0)) // 1.0 * x_ij
                print_error("wrong CPXchgcoef [degree]");
            if (CPXchgcoef(env, lp, row, xpos_dir(j, i, inst), 1.0)) // 1.0 * x_ji
                print_error("wrong CPXchgcoef [degree]");
            free(beg);
            free(ind);
        }
    }

    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}

void GG_original(CPXENVptr env, CPXLPptr lp, instance *inst)
{
    basic_model_directed(env, lp, inst);

    char binary = 'B';  // B => binary variable flag
    char integer = 'I'; // I => integer variable flag

    // cname: columns' names (column = variable)
    char **cname = (char **)calloc(1, sizeof(char *)); // array of strings to store the column names
    cname[0] = (char *)calloc(100, sizeof(char));

    // rname: rows' names (row = constraint)
    char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
    rname[0] = (char *)calloc(100, sizeof(char));

    // Add y-variables one for each arc (i,j) with i!=j and i,j > 0
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            sprintf(cname[0], "y(%d,%d)", i + 1, j + 1);
            double obj = 0.0;
            double lb = 0.0;
            double ub = inst->dimension - 2;
            if (i == 0)
                ub = inst->dimension - 1;
            if (i == j || j == 0)
                ub = 0;
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname))
                print_error(" wrong CPXnewcols on y var.s");
            if (CPXgetnumcols(env, lp) - 1 != ypos(i, j, inst))
                print_error("[position_d] wrong position for y var.s");
        }
    }

    // Add in-flow out-flow differential constraints for nodes h > 0: sum_i (1.0 * y_ih) + sum_i (- 1.0 * y_hi) = 1
    for (int h = 1; h < inst->dimension; h++) // exludes node 0
    {
        double rhs = 1.0;
        char sense = 'E';                 // E stands for equal
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        sprintf(rname[0], "in_flow out_flow node (%d)", h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [degree]");
        for (int i = 0; i < inst->dimension; i++)
        {
            if (h == i)
                continue;
            if (CPXchgcoef(env, lp, row, ypos(i, h, inst), 1.0)) // 1.0 * y_ih
                print_error("wrong CPXchgcoef [degree]");
            if (CPXchgcoef(env, lp, row, ypos(h, i, inst), -1.0)) // - 1.0 * y_hi
                print_error("wrong CPXchgcoef [degree]");
        }
    }

    // Add out-flow constraint for node 0: sum(y_0j) = n-1
    double rhs = inst->dimension - 1;
    char sense = 'E'; // E stands for equal
    sprintf(rname[0], "out_flow(0)");
    int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
    if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
        print_error("wrong CPXnewrows [degree]");
    for (int j = 1; j < inst->dimension; j++) // exclude  arc (0,0)
    {
        if (CPXchgcoef(env, lp, row, ypos(0, j, inst), 1.0)) // 1.0 * y_0j
            print_error("wrong CPXchgcoef [degree]");
    }

    // original linking constraints: y_ij <= (n-1)x_ij for each i != j
    for (int i = 0; i < inst->dimension; i++) // exludes node 0
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            if (i == j)
                continue;
            double rhs = 0.0;
            char sense = 'L';
            int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
            sprintf(rname[0], "linking constraints (%d,%d)", i + 1, j + 1);
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
                print_error("wrong CPXnewrows [degree]");

            if (CPXchgcoef(env, lp, row, ypos(i, j, inst), 1.0)) // 1.0 * y_ij
                print_error("wrong CPXchgcoef [degree]");
            if (i == 0)
            {
                if (CPXchgcoef(env, lp, row, xpos_dir(i, j, inst), 1 - inst->dimension)) // - (1 - nnodes) * x_ij
                    print_error("wrong CPXchgcoef [degree]");
            }
            else if (CPXchgcoef(env, lp, row, xpos_dir(i, j, inst), 2 - inst->dimension)) // - (2 - nnodes) * x_ij
                print_error("wrong CPXchgcoef [degree]");
        }
    }

    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}

void benders(CPXENVptr env, CPXLPptr lp, instance *inst)
{
    // Application of Benders method
    int done = 0;
    int it = 0; // iteration number
    while (!done)
    {
        // update time left
        double ts_current;
        inst->param.ticks ? CPXgetdettime(env, &ts_current) : getTimeStamp(&ts_current);
        double time_left = inst->time_limit - (ts_current - inst->timestamp_start);
        if (time_left <= 0.5 + inst->param.ticks * 500)
            return;
        CPXsetdblparam(env, CPX_PARAM_TILIM, time_left);

        double *xstar = (double *)calloc(inst->cols, sizeof(double));

        int status = CPXgetx(env, lp, xstar, 0, inst->cols - 1);
        if (status)
        {
            print_error_status("Failed to obtain the values in LOOP method", status);
        }

        int *comp = (int *)calloc(inst->dimension, sizeof(int));
        int c = 0; // number of connected components
        int *length_comp;

        // Retrieve the number of connected components so far
        if (inst->model_type == 9)
            findConnectedComponents(xstar, inst, inst->succ, comp, &c, &length_comp);
        else if (inst->model_type == 10)
            findConnectedComponents_kruskal(xstar, inst, inst->succ, comp, &c, &length_comp);
        char str[10] = "TIME";
        if (inst->param.ticks)
            sprintf(str, "TICKS");
        printf("\nITERATION: %d\tCONNECTED COMPONENTS FOUND: %d \t%s LEFT: %f\n", ++it, c, str, time_left);

        if (c == 1)
        {
            // If exactly one component is found, end the loop and exit
            if (inst->param.verbose >= NORMAL)
            {
                for (int n = 0; n < c; n++)
                {
                    printf("\t -- COMPONENT %d : %d NODES\n", n + 1, length_comp[n]);
                }
            }
            done = 1;
        }
        else
        {
            // If more than one component plot and check the partial solution
            if (inst->param.verbose >= NORMAL)
            {
                for (int n = 0; n < c; n++)
                {
                    printf("\t -- COMPONENT %d : %d NODES\n", n + 1, length_comp[n]);
                }
            }

            save_and_plot_solution(inst, it);

            // Add a SEC for each component found

            // rname: rows' names (row = constraint)
            char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
            rname[0] = (char *)calloc(100, sizeof(char));

            for (int mycomp = 0; mycomp < c; mycomp++)
            {

                double rhs = length_comp[mycomp] - 1.0;
                char sense = 'L';
                int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
                sprintf(rname[0], "SEC(%d,%d)", mycomp, it);
                if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
                    print_error("wrong CPXnewrows");

                for (int i = 0; i < inst->dimension; i++)
                {
                    if (inst->model_type == 9)
                    {
                        if (comp[i] != mycomp)
                            continue;
                        for (int j = i + 1; j < inst->dimension; j++)
                        {
                            if (comp[j] != mycomp)
                                continue;
                            // add (i,j) to SEC
                            if (CPXchgcoef(env, lp, row, xpos(i, j, inst), 1.0)) // 1.0 * x_ij
                                print_error("wrong CPXchgcoef");
                        }
                    }
                    else if (inst->model_type == 10)
                    {
                        if (comp[i] != inst->succ[mycomp])
                            continue;
                        for (int j = i + 1; j < inst->dimension; j++)
                        {
                            if (comp[j] != inst->succ[mycomp])
                                continue;
                            // add (i,j) to SEC
                            if (CPXchgcoef(env, lp, row, xpos(i, j, inst), 1.0)) // 1.0 * x_ij
                                print_error("wrong CPXchgcoef");
                        }
                    }
                }
            }

            // solve with the new constraints
            if (CPXmipopt(env, lp))
                print_error("CPXmipopt() error");
            free(rname[0]);
            free(rname);
        }

        free(comp);
        free(xstar);
        free(length_comp);
    }
}
