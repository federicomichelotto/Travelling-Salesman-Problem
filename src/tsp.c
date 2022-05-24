#include "../include/utils.h"

int tsp_solver(instance_TSP *inst)
{
    printf("Computing the optimal TSP solution...\n");
    // Open CPLEX model
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");
    // get timestamp

    inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_start) : getTimeStamp(&inst->timestamp_start);

    // build the GG model with lazy constraints and static 2-degree constraints
    GG_lazy_2sec(env, lp, inst);
    inst->cols = CPXgetnumcols(env, lp);
    inst->best_sol = (double *)calloc(inst->cols, sizeof(double));

    // CPLEX's parameter setting
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

    // directed graph
    gather_solution_TSP(inst, inst->best_sol, 1);

    printf("\tOPTIMAL TSP SOLUTION: 0");
    int node = 0;
    for (int i = 0; i < inst->dimension; i++)
    {
        node = inst->succ[node];
        printf(" -> %d", node);
    }

    printf("\n\tObjective value: %lf\n", inst->z_best);
    printf("\tLower bound: %lf\n", inst->best_lb);

    // get timestamp
    inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_finish) : getTimeStamp(&inst->timestamp_finish);

    // Plot optimal solution
    save_and_plot_solution_TSP(inst);

    // Free and close CPLEX model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    return inst->z_best;
}

void GG_lazy_2sec(CPXENVptr env, CPXLPptr lp, instance_TSP *inst)
{
    basic_model_directed_TSP(env, lp, inst);

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

void basic_model_directed_TSP(CPXENVptr env, CPXLPptr lp, instance_TSP *inst)
{
    char binary = 'B';  // B => binary variable flag
    char integer = 'I'; // I => integer variable flag

    // cname: columns' names (column = variable)
    char **cname = (char **)calloc(1, sizeof(char *)); // array of strings to store the column names
    cname[0] = (char *)calloc(100, sizeof(char));

    // rname: rows' names (row = constraint)
    char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
    rname[0] = (char *)calloc(100, sizeof(char));

    // Add binary variables x(i,j) for each (i,j)
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
            double obj = inst->weights[i][j]; // cost == time
            double lb = 0.0;
            double ub = 1.0;
            if (i == j)
                ub = 0.0;
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
                print_error(" wrong CPXnewcols on x var.s");
            if (CPXgetnumcols(env, lp) - 1 != xpos_dir(i, j, inst))
                print_error("[position_d] wrong position for x var.s");
        }
    }

    // Add the in-degree constraints
    for (int h = 0; h < inst->dimension; h++)
    {
        int row = CPXgetnumrows(env, lp); // get the maximum number of row inside the model
        double rhs = 1.0;
        char sense = 'E'; // E stands for equality constraint
        sprintf(rname[0], "in_degree(%d)", h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [degree]");
        for (int i = 0; i < inst->dimension; i++)
        {
            if (CPXchgcoef(env, lp, row, xpos_dir(i, h, inst), 1.0))
                print_error("wrong CPXchgcoef [degree]");
        }
    }

    // Add the out-degree constraints
    for (int h = 0; h < inst->dimension; h++)
    {
        int row = CPXgetnumrows(env, lp); // get the maximum number of row inside the model
        double rhs = 1.0;
        char sense = 'E'; // E stands for equality constraint
        sprintf(rname[0], "out_degree(%d)", h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [degree]");
        for (int i = 0; i < inst->dimension; i++)
        {
            if (CPXchgcoef(env, lp, row, xpos_dir(h, i, inst), 1.0))
                print_error("wrong CPXchgcoef [degree]");
        }
    }
    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}

double gather_solution_TSP(instance_TSP *inst, const double *xstar, int type)
{
    int n_edges = 0;
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            if (xstar[xpos_dir(i, j, inst)] > 0.5)
            {
                inst->succ[i] = j;
                if (++n_edges > inst->dimension)
                    print_error("more arcs than nodes, not a hamiltonian tour.");
            }
        }
    }
}

int xpos_dir(int i, int j, instance_TSP *inst)
{

    if (i < 0 || j < 0)
        print_error("Negative indexes are not valid!");
    if (i > inst->dimension || j > inst->dimension)
        print_error("Indexes exceeding the dimension are not valid!");
    return i * inst->dimension + j;
}

int upos(int i, instance_TSP *inst)
{
    if (i < 0)
        print_error("Negative index is not valid!");
    return xpos_dir(inst->dimension - 1, inst->dimension - 1, inst) + 1 + i;
}

int ypos(int i, int j, instance_TSP *inst)
{
    if (i < 0)
        print_error("Negative index is not valid!");
    return xpos_dir(inst->dimension - 1, inst->dimension - 1, inst) + 1 + i * inst->dimension + j;
}

int generate_TSP_instance(instance *inst, instance_TSP *inst_tsp)
{
    inst_tsp->dimension = inst->dimension - 1;
    inst_tsp->time_limit = inst->time_limit * 0.1; // assign 1/10 of the time available for the FSTSP problem
    inst_tsp->succ = (int *)calloc(inst_tsp->dimension, sizeof(int));
    inst_tsp->weights = inst->truck_times;
    inst_tsp->z_best = -1.0;
    inst_tsp->best_lb = -1.0;
    inst_tsp->timestamp_start = 0.0;
    inst_tsp->timestamp_finish = 0.0;
    inst_tsp->gnuplotPipe = inst->gnuplotPipe;

    inst_tsp->nodes = (node_TSP *)calloc(inst_tsp->dimension, sizeof(node_TSP));
    // copy nodes
    for (int i = 0; i < inst_tsp->dimension; i++)
    {
        inst_tsp->nodes[i].id = inst->nodes[i].id;
        inst_tsp->nodes[i].x = inst->nodes[i].x;
        inst_tsp->nodes[i].y = inst->nodes[i].y;
    }

    strcpy(inst_tsp->instance_path, inst->instance_path);

    inst_tsp->param.ticks = inst->param.ticks;
    inst_tsp->param.seed = inst->param.seed;
    inst_tsp->param.verbose = inst->param.verbose;
    inst_tsp->param.interactive = inst->param.interactive;
    inst_tsp->param.saveplots = inst->param.saveplots;
}

void post_tsp_sol(CPXENVptr env, CPXLPptr lp, instance_TSP *inst_tsp, instance *inst)
{
    //int nedges = inst_tsp->dimension * inst_tsp->dimension;
    int array_size = 3 * inst_tsp->dimension + 1;
    int *indices = (int *)malloc(array_size * sizeof(int));
    double *coeffs = (double *)malloc(array_size * sizeof(double));
    int count = 0; // counter: 0 -> nnodes-1
    // for (int i = 0; i < nedges; i++)
    // {
    //     if (inst_tsp->best_sol[i] > 0.5)
    //     {
    //         printf(" %d->%d\n", i / inst_tsp->dimension, i % inst_tsp->dimension);
    //         indices[j] = i;
    //         coeffs[j] = 1.0;
    //         j++;
    //     }
    // }

    // truck and drone variables
    for (int i = 0; i < inst_tsp->dimension; i++)
    {
        for (int j = 0; j < inst_tsp->dimension; j++)
        {
            if (inst_tsp->best_sol[xpos_dir(i, j, inst_tsp)] > 0.5)
            {
                //printf(" %d->%d\n", i, j);
                if (j == 0)
                {
                    indices[count] = xTruck_pos(i, inst->dimension - 1, inst);
                    coeffs[count++] = 1.0;
                    indices[count] = xDrone_pos(i, inst->dimension - 1, inst);
                    coeffs[count++] = 1.0;
                }
                else
                {
                    indices[count] = xTruck_pos(i, j, inst);
                    coeffs[count++] = 1.0;
                    indices[count] = xDrone_pos(i, j, inst);
                    coeffs[count++] = 1.0;
                }
            }
        }
    }
    //printf("count= %d\n", count);
    if (count != 2 * inst_tsp->dimension)
        print_error("Error in post_tsp_sol");

    // set all nodes as combined customers
    for (int i = 0; i < inst->dimension; i++)
    {
        indices[count] = yC_pos(i, inst);
        coeffs[count] = 1.0;
        count++;
    }
    //printf("count= %d\n", count);
    if (count != 3 * inst_tsp->dimension + 1)
        print_error("Error in post_tsp_sol");

    int beg = 0;
    int effort_level = CPX_MIPSTART_AUTO;
    if (CPXaddmipstarts(env, lp, 1, array_size, &beg, indices, coeffs, &effort_level, NULL))
        print_error("CPXaddmipstarts error");
}

void save_and_plot_solution_TSP(instance_TSP *inst)
{
    if (inst->param.interactive) // plot solution
    {
        // write solution to file

        char data_points_filename[200];
        sprintf(data_points_filename, "data_temp/data_points");
        FILE *data_points = fopen(data_points_filename, "w");

        char data_filename[100];
        sprintf(data_filename, "data_temp/tsp_opt");
        FILE *temp = fopen(data_filename, "w");
        for (int i = 0; i < inst->dimension; i++)
        {
            //Write the coordinates of the two nodes inside a temporary file
            fprintf(temp, "%lf %lf \n%lf %lf \n\n", inst->nodes[i].x, inst->nodes[i].y, inst->nodes[inst->succ[i]].x, inst->nodes[inst->succ[i]].y);
            // double '\n' to create a new edge (block of coordinates)
            fprintf(data_points, "%lf %lf %d \n", inst->nodes[i].x, inst->nodes[i].y, inst->nodes[i].id);
        }

        char plot_str[1000];

        sprintf(plot_str, "plot '%s' with linespoints pt 7 lc rgb 'blue' lw 1,\
                         '%s' using 1:2:3 with labels offset (1,1) font 'Arial'",
                data_filename, data_points_filename);
        char *commandsForGnuplot[] = {"set title 'Optimal TSP solution'",
                                      "set term wxt noraise",
                                      "unset key",
                                      "set autoscale",
                                      "set ylabel 'Y'",
                                      "set xlabel 'X'",
                                      plot_str};

        int commands = sizeof(commandsForGnuplot) / sizeof(commandsForGnuplot[0]);
        for (int i = 0; i < commands; i++)
            fprintf(inst->gnuplotPipe, "%s \n", commandsForGnuplot[i]);
        fflush(inst->gnuplotPipe); // execute the commands

        fclose(temp);
        fclose(data_points);
    }
}