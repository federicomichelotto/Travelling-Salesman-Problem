#include "../include/utils.h"

double dist(int i, int j, instance *inst)
{

    double distance = INFINITY;

    if (strncmp(inst->param.weight_type, "GEO", 3) == 0)
    {
        double deg, min;
        deg = (int)inst->nodes[i].x;
        min = inst->nodes[i].x - deg;
        double lat_i = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

        deg = (int)inst->nodes[i].y;
        min = inst->nodes[i].y - deg;
        double long_i = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

        deg = (int)inst->nodes[j].x;
        min = inst->nodes[j].x - deg;
        double lat_j = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

        deg = (int)inst->nodes[j].y;
        min = inst->nodes[j].y - deg;
        double long_j = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

        double RRR = 6378.388;

        double q1 = cos(long_i - long_j);
        double q2 = cos(lat_i - lat_j);
        double q3 = cos(lat_i + lat_j);

        distance = (int)(RRR * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    }
    else if (strncmp(inst->param.weight_type, "EUC_2D", 6) == 0)
    {
        double dx = inst->nodes[i].x - inst->nodes[j].x;
        double dy = inst->nodes[i].y - inst->nodes[j].y;
        if (!inst->integer_costs)
            return sqrt(dx * dx + dy * dy);
        int dis = sqrt(dx * dx + dy * dy) + 0.499999999;
        distance = dis + 0.0;
    }
    else if (strncmp(inst->param.weight_type, "ATT", 3) == 0)
    {
        double dx = inst->nodes[i].x - inst->nodes[j].x;
        double dy = inst->nodes[i].y - inst->nodes[j].y;
        double dis1 = sqrt((dx * dx + dy * dy) / 10.0);
        int dis2 = (int)(dis1 + 0.5);
        if (dis2 < dis1)
            distance = dis2 + 1;
        else
            distance = dis2;
    }

    return distance;
}

int xpos(int i, int j, instance *inst)
{

    if (i == j)
        print_error("Same indices are not valid!");
    if (i < 0 || j < 0)
        print_error("Negative indexes are not valid!");
    if (i > inst->dimension || j > inst->dimension)
        print_error("Indexes exceeding the dimension are not valid!");
    if (i > j)
        return xpos(j, i, inst);
    else
        return i * inst->dimension + j - (i + 1) * (i + 2) / 2;
}

int xpos_dir(int i, int j, instance *inst)
{

    if (i < 0 || j < 0)
        print_error("Negative indexes are not valid!");
    if (i > inst->dimension || j > inst->dimension)
        print_error("Indexes exceeding the dimension are not valid!");
    return i * inst->dimension + j;
}

int upos(int i, instance *inst)
{

    if (i < 0)
        print_error("Negative index is not valid!");
    return xpos_dir(inst->dimension - 1, inst->dimension - 1, inst) + 1 + i;
}

int ypos(int i, int j, instance *inst)
{

    if (i < 0)
        print_error("Negative index is not valid!");
    return xpos_dir(inst->dimension - 1, inst->dimension - 1, inst) + 1 + i * inst->dimension + j;
}

double gather_solution(instance *inst, const double *xstar, int type)
{
    int n_edges = 0;
    if (type == 0) // undirected graph
    {
        for (int i = 0; i < inst->dimension; i++)
        {
            inst->succ[i] = -1;
        }
        int init = 0;
        int node = 0;
        for (int i = 0; i < inst->dimension; i++)
        {
            if (i == node || inst->succ[i] == node)
                continue;
            if (xstar[xpos(node, i, inst)] > 0.5)
            {
                inst->succ[node] = i;
                if (++n_edges == inst->dimension)
                    break;
                node = i;
                if (i == init)
                {
                    // look for a new connected component
                    for (int j = 0; j < inst->dimension; j++)
                    {
                        if (inst->succ[j] != -1)
                        { // node not visited yet
                            node = j;
                            break;
                        }
                    }
                }
                i = -1; // restart the loop from i=0
            }
        }
    }
    else if (type == 1) // directed graph
    {
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
}

void findConnectedComponents(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp,
                             int **length_comp)
{
    // build succ[] and comp[] wrt xstar()...

    // initialization
    *ncomp = 0;
    for (int i = 0; i < inst->dimension; i++)
    {
        succ[i] = -1;
        comp[i] = -1;
    }

    for (int start = 0; start < inst->dimension; start++)
    {
        if (comp[start] >= 0)
            continue; // node "start" was already visited, just skip it

        // a new component is found
        (*ncomp)++;
        if ((*ncomp) == 1)
            (*length_comp) = (int *)calloc(1, sizeof(int));
        else
            (*length_comp) = (int *)realloc(*length_comp, (*ncomp) * sizeof(int));
        int i = start;
        int length = 1;
        int done = 0;
        while (!done) // go and visit the current component
        {
            comp[i] = (*ncomp) - 1;
            done = 1;
            for (int j = 0; j < inst->dimension; j++)
            {
                if (i == j)
                    continue;
                if (xstar[xpos(i, j, inst)] > 0.5 &&
                    comp[j] == -1) // the edge [i,j] is selected in xstar and j was not visited before
                {
                    succ[i] = j;
                    length++;
                    i = j;
                    done = 0;
                    break;
                }
            }
        }
        succ[i] = start;                       // last arc to close the cycle
        (*length_comp)[(*ncomp) - 1] = length; // save length of the cycle

        // go to the next component...
    }
}

// Kruskal algorithm to find connected components
void findConnectedComponents_kruskal(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp,
                                     int **length_comp)
{

    // Some initialization
    for (int i = 0; i < inst->dimension; i++)
    {
        comp[i] = i;
        succ[i] = -1;
    }

    // There are no components so far
    *ncomp = 0;

    // Found connected components
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = i + 1; j < inst->dimension; j++)
        {
            if (xstar[xpos(i, j, inst)] == 1)
            {
                int c1 = comp[i];
                int c2 = comp[j];
                if (c1 != c2)
                {
                    for (int k = 0; k < inst->dimension; k++)
                        if (comp[k] == c2)
                            comp[k] = c1;
                }
            }
        }
    }

    // Count how many components are present
    for (int i = 0; i < inst->dimension; i++)
    {

        int inside = 0;

        for (int j = 0; j < (*ncomp) + 1; j++)
        {
            if (succ[j] == comp[i])
            {
                inside = 1;
                break;
            }
        }

        if (!inside)
        {
            succ[(*ncomp)] = comp[i];
            (*ncomp)++;
        }
    }

    (*length_comp) = (int *)calloc((*ncomp), sizeof(int));

    for (int k = 0; k < (*ncomp); k++)
    {
        int length = 0;
        for (int h = 0; h < inst->dimension; h++)
        {

            if (comp[h] != succ[k])
                continue;
            else
                length++;
        }

        (*length_comp)[k] = length;
    }
}

static int CPXPUBLIC callback_driver(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle)
{
    instance *inst = (instance *)userhandle;
    if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
        return callback_candidate(context, contextid, userhandle);
    if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
        return callback_relaxation(context, contextid, userhandle);
    print_error("contextid unknownn in my_callback");
    return 1;
}

static int CPXPUBLIC callback_candidate(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle)
{
    instance *inst = (instance *)userhandle;
    double *xstar = (double *)malloc(inst->cols * sizeof(double));
    double objval = CPX_INFBOUND;
    /*
    if (inst->param.verbose >= NORMAL)
    {
        printf("*** callback #%d (candidate)\n", ++(inst->param.callback_counter));
    }
    */
    if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->cols - 1, &objval))
        print_error("CPXcallbackgetcandidatepoint error");

    int *comp = (int *)calloc(inst->dimension, sizeof(int));
    int *succ = (int *)calloc(inst->dimension, sizeof(int));
    int ncomp = 0; // number of connected components
    int *length_comp;

    // Retrieve the connected components of the current solution
    findConnectedComponents(xstar, inst, succ, comp, &ncomp, &length_comp);

    if (ncomp > 1)
    {
        // add one cut for each connected component
        for (int mycomp = 0; mycomp < ncomp; mycomp++)
        {
            int nnz = 0;
            int izero = 0;
            char sense = 'L';
            double rhs = length_comp[mycomp] - 1.0; // in order to have |S|-1 in the end
            int *index = (int *)calloc(inst->cols, sizeof(int));
            double *value = (double *)calloc(inst->cols, sizeof(double));

            for (int i = 0; i < inst->dimension; i++)
            {
                if (comp[i] != mycomp)
                    continue;
                for (int j = i + 1; j < inst->dimension; j++)
                {
                    if (comp[j] != mycomp)
                        continue;
                    index[nnz] = xpos(i, j, inst);
                    value[nnz++] = 1.0;
                }
            }
            if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value))
                print_error("CPXcallbackrejectcandidate() error"); // reject the solution and adds one cut
            /*
            if (inst->param.verbose >= DEBUG)
            {
                printf("### added a SEC constraint \n");
            }
            */
            free(index);
            free(value);
        }
    }
    free(comp);
    free(succ);
    free(length_comp);
    free(xstar);
    return 0;
}

static int CPXPUBLIC callback_relaxation(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle)
{
    double ticks = 0;
    CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_DETTIME, &ticks);
    if (!((int)ticks % 10))
        return 0;
    instance *inst = (instance *)userhandle;
    double *xstar = (double *)malloc(inst->cols * sizeof(double));
    double objval = CPX_INFBOUND;
    double const eps = 0.1;
    /*
    if (inst->param.verbose >= NORMAL)
    {
        printf("*** callback #%d (relaxation) \n", ++(inst->param.callback_counter));
    }
    */
    if (CPXcallbackgetrelaxationpoint(context, xstar, 0, inst->cols - 1, &objval))
        print_error("CPXcallbackgetrelaxationpoint error");

    int ncomp;
    int *comp = (int *)calloc(inst->dimension, sizeof(int));
    int *length_comp = (int *)calloc(inst->dimension, sizeof(int));
    // list of edges in "node format"
    int elist[2 * inst->cols]; // [0,1, 0,2, 0,3, ...]

    int loader = 0;
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = i + 1; j < inst->dimension; j++)
        {
            // CHECK xstar > 0 (eps)
            // CCxstar
            elist[loader++] = i;
            elist[loader++] = j;
        }
    }

    if (CCcut_connect_components(inst->dimension, inst->cols, elist, xstar, &ncomp, &length_comp, &comp))
        print_error("CCcut_connect_components error");

    doit_fn_input in;
    in.context = context;
    in.inst = inst;

    if (ncomp == 1)
    {
        if (CCcut_violated_cuts(inst->dimension, inst->cols, elist, xstar, 2 - eps, doit_fn_concorde, (void *)&in))
            print_error("CCcut_violated_cuts error");
    }

    free(comp);
    free(length_comp);
    free(xstar);
    return 0;
}

// double cutval = value of the cut
// int cutcount = number of nodes in the cut (rhs + 1 ?)
// int ∗cut = the array of the members of the cut (indeces of the nodes in the cut?)
int doit_fn_concorde(double cutval, int cutcount, int *cut, void *in)
{
    doit_fn_input *input = (doit_fn_input *)in;
    double rhs = cutcount - 1.0;
    int nnz = 0;
    char sense = 'L';
    int purgeable = CPX_USECUT_FILTER; // Let CPLEX decide whether to keep the cut or not
    int local = 0;
    int izero = 0;
    /*
    if (input->inst->param.verbose >= DEBUG)
    {
        printf("#### cutval: %f ", cutval);
        printf("#### cutcount: %d \n", cutcount);
        for (int i = 0; i < cutcount; i++)
            printf("cut[%d] = %d\n", i, cut[i]);
    }
    */
    double *value = (double *)calloc(cutcount * (cutcount - 1) / 2, sizeof(double));
    int *index = (int *)calloc(cutcount * (cutcount - 1) / 2, sizeof(int));

    // CHECK THIS
    for (int i = 0; i < cutcount; i++)
    {
        for (int j = 0; j < cutcount; j++)
        {
            if (cut[i] < cut[j])
            {
                index[nnz] = xpos(cut[i], cut[j], input->inst);
                value[nnz++] = 1.0;
            }
        }
    }
    // TODO: CHECK nnz = cutcount*(cutcount-1)/2

    if (CPXcallbackaddusercuts(input->context, 1, nnz, &rhs, &sense, &izero, index, value, &purgeable, &local))
        print_error("CPXcallbackaddusercuts() error"); // add user cut
    /*
    if (input->inst->param.verbose >= DEBUG)
    {
        printf("### added a user cut \n");
    }
    */
    free(index);
    free(value);
    return 0;
}

int optimal_solver(instance *inst)
{
    // Open CPLEX model
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");
    // get timestamp
    inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_start) : CPXgettime(env, &inst->timestamp_start);

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
    inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_finish) : CPXgettime(env, &inst->timestamp_finish);

    // Free and close CPLEX model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    return 0;
}

int math_solver(instance *inst)
{
    // Open CPLEX model
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "MATH TSP");

    // get timestamp
    inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_start) : CPXgettime(env, &inst->timestamp_start);
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
    if (inst->param.ticks)
    {
        if (CPXsetdblparam(env, CPX_PARAM_DETTILIM, inst->time_limit / 20))
            print_error("CPX_PARAM_DETTILIM error");
    }
    else
    {
        if (CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit / 20))
            print_error("CPX_PARAM_TILIM error");
    }

    // initial solution
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

    gather_solution(inst, inst->best_sol, 0);
    save_and_plot_solution(inst, 0);

    printf("\nSOLUTION -----------------------------------------------\n");
    printf("\nRUNNING : %s\n", math_model_full_name[inst->model_type]);

    // hard-fixing heuristic
    if (inst->model_type == 0)
    {
        if (inst->param.verbose >= NORMAL)
            printf("Initial incumbent: %f\n", inst->z_best);
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
    inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_finish) : CPXgettime(env, &inst->timestamp_finish);

    // Free and close CPLEX model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    return 0;
}

int heuristic_solver(instance *inst)
{
    // get timestamp
    struct timespec timestamp;
    if (clock_gettime(CLOCK_REALTIME, &timestamp) == -1)
        print_error("Error clock_gettime");
    inst->timestamp_start = timestamp.tv_sec + timestamp.tv_nsec * pow(10, -9);

    double min_obj = DBL_MAX;
    double obj_i = 0;
    double obj_opt = 0;

    printf("\nSOLUTION -----------------------------------------------\n");
    printf("\nRUNNING : %s\n", heuristic_model_full_name[inst->model_type]);
    int *succ_i = (int *)calloc(inst->dimension, sizeof(int));

    switch (inst->model_type)
    {
    case 0: // Nearest Neighbours
        for (int i = 0; i < inst->dimension; i++)
        {
            if (inst->param.grasp == 1) obj_i = nearest_neighbours(inst, i, succ_i, inst->param.grasp_choices);
            else obj_i = nearest_neighbours(inst, i, succ_i, 1);

            if (obj_i < min_obj)
            {
                min_obj = obj_i;
                inst->z_best = obj_i;
                for (int j = 0; j < inst->dimension; j++)
                    inst->succ[j] = succ_i[j];
                save_and_plot_solution(inst, i + 1);
            }
        }
        printf("2-opt: %d crossings found\n", two_opt_v2(inst));
        save_and_plot_solution(inst, inst->dimension);
        printf("Best objective value: %f\n", min_obj);
        printf("Best objective value (optimized by 2-opt): %f\n", inst->z_best);
        break;
    case 1: // Extra Mileage
        for (int i = 0; i < inst->dimension; i++)
        {
            obj_i = extra_mileage(inst, i, succ_i);
            if (obj_i < min_obj)
            {
                min_obj = obj_i;
                inst->z_best = obj_i;
                for (int j = 0; j < inst->dimension; j++)
                    inst->succ[j] = succ_i[j];
                save_and_plot_solution(inst, i + 1);
            }
        }
        printf("2-opt: %d crossings found\n", two_opt_v2(inst));
        save_and_plot_solution(inst, inst->dimension);
        printf("Best objective value: %f\n", min_obj);
        printf("Best objective value (optimized by 2-opt): %f\n", inst->z_best);
        break;
    case 2: // Extra Mileage with furthest starting nodes
        min_obj = extra_mileage_furthest_starting_nodes(inst, succ_i);
        inst->z_best = min_obj;
        for (int j = 0; j < inst->dimension; j++)
            inst->succ[j] = succ_i[j];
        save_and_plot_solution(inst, 1);

        printf("2-opt: %d crossings found\n", two_opt_v2(inst));
        save_and_plot_solution(inst, 2);
        printf("Best objective value: %f\n", min_obj);
        printf("Best objective value (optimized by 2-opt): %f\n", inst->z_best);
        break;
    case 3: // Tabu Search

        if (inst->param.grasp == 1) min_obj = nearest_neighbours(inst, 0, inst->succ, inst->param.grasp_choices);
        else min_obj = nearest_neighbours(inst, 0, inst->succ, 1);

        inst->z_best = min_obj;
        tabu_search(inst);
        printf("Best objective value: %f\n", min_obj);
        printf("Best objective value (optimized by tabu search): %f\n", inst->z_best);
        break;
    default:
        fprintf(stderr, "ERROR: Model type %d not available.\n", inst->model_type);
        break;
    }

    free(succ_i);

    // get timestamp
    if (clock_gettime(CLOCK_REALTIME, &timestamp) == -1)
        print_error("Error clock_gettime");
    inst->timestamp_finish = timestamp.tv_sec + timestamp.tv_nsec * pow(10, -9);
    return 0;
}

void build_model(CPXENVptr env, CPXLPptr lp, instance *inst)
{
    char path[1000];

    if (inst->param.solver == 0)
    { // optimal solver involved
        switch (inst->model_type)
        {
        case 0:  // basic model (no SEC)
        case 9:  // benders model (SEC)
        case 10: // benders model (SEC) - kruskal
        case 11: // callback model (SEC)
            basic_model_no_sec(env, lp, inst);
            break;
        case 1: // MTZ with static constraints
            MTZ_static(env, lp, inst);
            break;
        case 2: // MTZ (mod) with static constraints
            MTZ_static_mod(env, lp, inst);
            break;
        case 3: // MTZ with lazy constraints
            MTZ_lazy(env, lp, inst);
            break;
        case 4: // MTZ with lazy constraints and sub-tour elimination constraints of degree 2
            MTZ_lazy_sec(env, lp, inst);
            break;
        case 5: // GG
            GG(env, lp, inst);
            break;
        case 6: // GG with lazy constraints
            GG_lazy(env, lp, inst);
            break;
        case 7: // GG with lazy constraints and sub-tour elimination constraints of degree 2
            GG_lazy_sec(env, lp, inst);
            break;
        case 8: // GG orignal formulation
            GG_original(env, lp, inst);
            break;
        default:
            fprintf(stderr, "ERROR: Model type %d not available.\n", inst->model_type);
            break;
        }

        if (generate_path(path, "output", "model", optimal_model_name[inst->model_type], inst->param.name,
                          inst->param.seed, "lp"))
            print_error("Unable to generate path");
    }
    else if (inst->param.solver == 1)
    {

        switch (inst->model_type)
        {
        case 0: // hard fixing heuristic
        case 1: // soft fixing heuristic
            basic_model_no_sec(env, lp, inst);
            break;
        default:
            fprintf(stderr, "ERROR: Model type %d not available.\n", inst->model_type);
            break;
        }

        if (generate_path(path, "output", "model", math_model_name[inst->model_type], inst->param.name,
                          inst->param.seed, "lp"))
            print_error("Unable to generate path");
    }

    CPXwriteprob(env, lp, path, NULL);
}

void basic_model_no_sec(CPXENVptr env, CPXLPptr lp, instance *inst)
{
    char binary = 'B'; // B => binary variable flag
    // cname: columns' names (column = variable)
    char **cname = (char **)calloc(1, sizeof(char *)); // array of strings to store the column names
    cname[0] = (char *)calloc(100, sizeof(char));

    // rname: rows' names (row = constraint)
    char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
    rname[0] = (char *)calloc(100, sizeof(char));
    // Add binary variables x(i,j) for i < j
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = i + 1; j < inst->dimension; j++)
        {
            sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
            double obj = dist(i, j, inst); // cost == distance
            double lb = 0.0;
            double ub = 1.0;
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
                print_error(" wrong CPXnewcols on x var.s");
            if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst))
                print_error(" wrong position for x var.s");
        }
    }

    // Add the degree constraints
    for (int h = 0; h < inst->dimension; h++)
    {
        int row = CPXgetnumrows(env, lp); // get the maximum number of row inside the model
        double rhs = 2.0;
        char sense = 'E'; // E stands for equality constraint

        sprintf(rname[0], "degree(%d)", h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [degree]");
        for (int i = 0; i < inst->dimension; i++)
        {
            if (i == h)
                continue;
            if (CPXchgcoef(env, lp, row, xpos(i, h, inst), 1.0))
                print_error("[position_u] wrong CPXchgcoef [degree]");
        }
    }
    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}

void basic_model_directed(CPXENVptr env, CPXLPptr lp, instance *inst)
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
            double obj = dist(i, j, inst); // cost == distance
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
        inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_finish) : CPXgettime(env, &inst->timestamp_finish);
        double time_left = inst->time_limit - (inst->timestamp_finish - inst->timestamp_start);
        if (time_left <= 0.5)
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

        printf("\nITERATION: %d\tCONNECTED COMPONENTS FOUND: %d \tTIME LEFT: %f\n", ++it, c, time_left);

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

void hard_fixing_heuristic(CPXENVptr env, CPXLPptr lp, instance *inst, int time_limit_iter, double fix_ratio)
{
    int iter = 1;
    while (1)
    {
        // update time left
        inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_finish) : CPXgettime(env, &inst->timestamp_finish);
        double time_left = inst->time_limit - (inst->timestamp_finish - inst->timestamp_start);
        if (time_left <= 0.5)
            return;
        if (time_left <= time_limit_iter)
            CPXsetdblparam(env, CPX_PARAM_TILIM, time_left);
        if (inst->param.verbose >= DEBUG)
            printf("*** time left = %f\n", time_left);

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

        // retrieve the incumbent of the current solution
        double current_incumbent;
        CPXgetobjval(env, lp, &current_incumbent);
        if (inst->param.verbose >= DEBUG)
            printf("*** current_incumbent = %f\n", current_incumbent);
        // check if the current solution is better than the best so far
        if (current_incumbent < inst->z_best)
        {
            // update best incumbent
            inst->z_best = current_incumbent;
            // update arcs' selection
            int status = CPXgetx(env, lp, inst->best_sol, 0, inst->cols - 1);
            if (status)
                print_error_status("Failed to obtain the values in hard_fixing_heuristic method", status);
            if (inst->param.verbose >= NORMAL)
                printf("New incumbent: %f\n", inst->z_best);
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
    while (1)
    {
        // update time left
        inst->param.ticks ? CPXgetdettime(env, &inst->timestamp_finish) : CPXgettime(env, &inst->timestamp_finish);
        double time_left = inst->time_limit - (inst->timestamp_finish - inst->timestamp_start);
        if (time_left <= 0.5)
            return;
        if (time_left <= time_limit_iter)
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

        // retrieve the incumbent of the current solution
        double current_incumbent;
        CPXgetobjval(env, lp, &current_incumbent);
        if (inst->param.verbose >= DEBUG)
            printf("*** k = %d, current_incumbent = %f\n", k, current_incumbent);
        // check if the current solution is better than the best so far
        if (current_incumbent < inst->z_best)
        {
            // update best incumbent
            inst->z_best = current_incumbent;
            // update best sol
            int status = CPXgetx(env, lp, inst->best_sol, 0, inst->cols - 1);
            if (status)
                print_error_status("Failed to obtain the values in soft_fixing_heuristic method", status);
            if (inst->param.verbose >= NORMAL)
                printf("New incumbent: %f\n", inst->z_best);

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

double nearest_neighbours(instance *inst, int starting_node, int *succ, int options)
{
    double obj = 0.0;

    int *selected = (int *)calloc(inst->dimension, sizeof(int));

    // Successor array initialization
    for (int i = 0; i < inst->dimension; i++) succ[i] = -1;

    selected[starting_node] = 1;
    int current = starting_node; // Index of the current node

    if (options > 1) printf("GRASP approach selected, available option for each node %d", options);

    // Build the circuit adding inst->dimension - 1 edges
    for (int count = 0; count < inst->dimension - 1; count++) {

        // Check available nodes (remember to ignore the current one)
        if (inst->dimension - count - 1 < options) {
            printf("Too many choices (%d) w.r.t. available nodes (%d). Available choices are now (%d)\n", options, inst->dimension - count - 1,inst->dimension - count-1);
            options = inst->dimension - count - 1;
            int a = -1;
            for (int i = 0; i < inst->dimension; i++) if(succ[i] == -1) a++;
            printf("Nodes not yet selected are %d\n", a);
        }

        double min_dist[options];         // Minimum distances
        int min_node[options];            // Closest nodes indices

        for (int i = 0; i < options; ++i) {
            min_dist[i] = CPX_INFBOUND;
            min_node[i] = -1;
        }

        int *chosen = (int *)calloc(inst->dimension, sizeof(int));
        int k = 0;

        // While a choice can be made, select the closest node
        while (k < options){
            // select the closest node w.r.t. to the current node
            for (int i = 0; i < inst->dimension; i++)
            {
                if (i != current && selected[i] == 0 && chosen[i] == 0) // i has not been selected yet
                {
                    double distance = dist(current, i, inst);
                    if (distance < min_dist[k])
                    {
                        min_dist[k] = distance;
                        min_node[k] = i;
                    }
                }
            }

            chosen[min_node[k]] = 1;  // i-th node was chosen as closest
            k++;
        }

        for (int j = 0; j < options; ++j) {
            printf("%d : %f\n", min_node[j], min_dist[j]);
        }
        printf("------------\n");

        // Minimum node random selection
        srand(time(0));
        int h = rand() % options;

        if (min_node[h] == -1)
            print_error("min_node[h] == -1");

        // Add edge current node - random node at minimum distance
        succ[current] = min_node[h];
        selected[min_node[h]] = 1;

        // Update current node
        current = min_node[h];

        // Update current incumbent
        obj += min_dist[h];
    }

    // Close circuit
    succ[current] = starting_node;

    // Update incumbent
    obj += dist(current, starting_node, inst);
    printf("Best objective value for starting node %d: %f\n", starting_node + 1, obj);

    return obj;
}

double extra_mileage(instance *inst, int starting_node, int *succ)
{

    double obj = 0; // Objective value
    double *succ_dist = (double *)calloc(inst->dimension,
                                         sizeof(double)); // array of the distance between the node i and succ[i]

    for (int i = 0; i < inst->dimension; i++)
    {
        succ[i] = -1;
    }

    // Find node at maximum distance from first node
    double max = 0;
    double distance;
    int idx = starting_node;
    for (int i = 0; i < inst->dimension; i++)
    {
        if (i != starting_node)
        {
            distance = dist(starting_node, i, inst);
            if (max < distance)
            {
                max = distance;
                idx = i;
            }
        }
    }

    succ[starting_node] = idx;
    succ[idx] = starting_node;
    succ_dist[starting_node] = distance;
    succ_dist[idx] = distance;
    obj += 2 * succ_dist[starting_node]; // starting_node->idx->starting_node

    // add inst->dimension - 2 nodes
    for (int iter = 0; iter < inst->dimension - 2; iter++)
    {

        double min_value = CPX_INFBOUND; // Insertion of node selected_node in the circuit
        int idx_edge;                    // (idx_edge, succ[idx_edge]): edge edge that corresponds to the smallest extra mileage to add the node selected_node to the circuit
        int selected_node;               // node to add in the circuit

        // for each uncovered node and for each edge in the circuit, compute the extra mileage
        for (int k = 0; k < inst->dimension; k++)
        {
            for (int i = 0; i < inst->dimension; i++)
            {
                if (succ[k] == -1 && succ[i] != -1)
                {
                    // k: node not in the circuit (it has no successor)
                    // (i, succ[i]): edge
                    double extra_mileage = dist(i, k, inst) + dist(succ[i], k, inst) - succ_dist[i];
                    if (extra_mileage < min_value)
                    {
                        min_value = extra_mileage;
                        idx_edge = i; // (i, succ[i])
                        selected_node = k;
                    }
                }
            }
        }
        succ[selected_node] = succ[idx_edge];
        succ_dist[selected_node] = dist(selected_node, succ[idx_edge], inst);
        succ[idx_edge] = selected_node;
        succ_dist[idx_edge] = dist(idx_edge, selected_node, inst);

        obj += min_value; // Update objective value
    }

    printf("Objective value for starting node %d: %f\n", starting_node + 1, obj);
    free(succ_dist);
    return obj;
}

// TO CHECK BETTER
double extra_mileage_furthest_starting_nodes(instance *inst, int *succ)
{
    double obj = 0;                                                        // Objective value
    double *succ_dist = (double *)calloc(inst->dimension, sizeof(double)); // array of the distance between the node i and succ[i]

    for (int i = 0; i < inst->dimension; i++)
    {
        succ[i] = -1;
    }

    // Find the two nodes at maximum distance
    double max = 0;
    double distance;
    int node1, node2;
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            if (i != j)
            {
                distance = dist(i, j, inst);
                if (distance > max)
                {
                    max = distance;
                    node1 = i;
                    node2 = j;
                }
            }
        }
    }

    succ[node1] = node2;
    succ[node2] = node1;
    succ_dist[node1] = distance;
    succ_dist[node2] = distance;
    obj += 2 * distance; // starting_node1->node2->node1

    // add inst->dimension - 2 nodes
    for (int iter = 0; iter < inst->dimension - 2; iter++)
    {
        double min_value = CPX_INFBOUND; // Insertion of node selected_node in the circuit
        int idx_edge;                    // (idx_edge, succ[idx_edge]): edge edge that corresponds to the smallest extra mileage to add the node selected_node to the circuit
        int selected_node;               // node to add in the circuit

        // for each uncovered node and for each edge in the circuit, compute the extra mileage
        for (int k = 0; k < inst->dimension; k++)
        {
            for (int i = 0; i < inst->dimension; i++)
            {
                if (succ[k] == -1 && succ[i] != -1)
                {
                    // k: node not in the circuit (it has no successor)
                    // (i, succ[i]): edge
                    double extra_mileage = dist(i, k, inst) + dist(succ[i], k, inst) - succ_dist[i];
                    if (extra_mileage < min_value)
                    {
                        min_value = extra_mileage;
                        idx_edge = i; // (i, succ[i])
                        selected_node = k;
                    }
                }
            }
        }
        succ[selected_node] = succ[idx_edge];
        succ_dist[selected_node] = dist(selected_node, succ[idx_edge], inst);
        succ[idx_edge] = selected_node;
        succ_dist[idx_edge] = dist(idx_edge, selected_node, inst);

        obj += min_value; // Update objective value
    }
    free(succ_dist);
    return obj;
}

int two_opt(instance *inst, int maxMoves)
{
    print_message("Inside 2-opt function");
    int optimal = 0;
    int moves = 0;
    double incumbent = inst->z_best;

    // For each couple of edges (a,b) and (c,d) so that they are not subsequent,
    // Compute d(a,c) + d(b,d)
    // if d(a,b) + d(c,d) > d(a,c) + d(b,d) --> change and rebuild solution and restart
    // else continue
    //
    // If all couples where considered and no changes were needed --> optimal == 1;
    while (!optimal && (moves < maxMoves))
    {
        for (int i = 0; i < inst->dimension; i++)
        {
            for (int j = 0; j < inst->dimension; j++)
            {
                if (i == j)
                    continue;
                if (inst->succ[i] == j || i == inst->succ[j])
                    continue;

                optimal = 1;
                int a = i;
                int b = inst->succ[a];
                int c = j;
                int d = inst->succ[c];
                double originalDist = dist(a, b, inst) + dist(c, d, inst);
                double newDist = dist(a, c, inst) + dist(b, d, inst);

                if (newDist < originalDist)
                { // crossing
                    printf("Edge#1 : (%d,%d)\tEdge#2 : (%d,%d)\n", a + 1, b + 1, c + 1, d + 1);
                    printf("Old distance : %f\n", originalDist);
                    printf("New distance : %f\n", newDist);
                    printf("---------------------------------------------------------------\n");

                    // update inc
                    incumbent = incumbent - originalDist + newDist;

                    // reverse tour
                    if (reverse_successors(inst->succ, inst->dimension, b, c))
                        print_error("Error in reverse_successors");
                    inst->succ[a] = c;
                    inst->succ[b] = d;

                    optimal = 0;
                    moves++;
                    break;
                }
            }
            if (!optimal)
                break;
        }
    }

    inst->z_best = incumbent;
    return !optimal;
}

// move applied on the most negative delta
int two_opt_v2(instance *inst)
{
    print_message("Inside 2-opt_v2 function");
    printf("Initial incument = %f\n", inst->z_best);
    // For each couple of edges (a,c) and (b,d) so that they are not subsequent,
    // If the given solution does not have any crossing: return 0, else return #crossing found;
    int iter = 1;
    while (1)
    {
        double min_delta = DBL_MAX;
        int a, b;
        for (int i = 0; i < inst->dimension; i++)
        {
            for (int j = 0; j < inst->dimension; j++)
            {
                // look for two nodes that are not connected
                if (i == j)
                    continue;
                if (inst->succ[i] == j || i == inst->succ[j])
                    continue;

                double delta = dist(i, j, inst) + dist(inst->succ[i], inst->succ[j], inst) - dist(i, inst->succ[i], inst) - dist(j, inst->succ[j], inst);
                if (delta < min_delta)
                {
                    min_delta = delta;
                    a = i;
                    b = j;
                }
            }
        }
        if (min_delta == DBL_MAX)
            print_error("Error in tabu search delta computation");
        if (min_delta >= -1e-5)
            break;

        int c = inst->succ[a];
        int d = inst->succ[b];


        // reverse the path from c to b
        if (reverse_successors(inst->succ, inst->dimension, c, b))
            print_error("Error in reverse_successors");

        // make the move
        inst->succ[a] = b;
        inst->succ[c] = d;
        // update incumbent
        inst->z_best += min_delta;
        iter++;
        printf("new incumbent = %f (delta = %f)\n", inst->z_best, min_delta);
        save_and_plot_solution(inst, inst->dimension);
    }
    if (iter == 1)
        return 0;
    return iter - 1;
}

int reverse_successors(int *succ, int size, int start, int end)
{
    if (start >= size || end >= size)
        print_error("Error in reverse_successors");
    int node = start;
    int next = succ[node];
    int counter = 0;
    while (node != end)
    {
        int temp = succ[next];
        succ[next] = node;
        node = next;
        next = temp;
        if (counter++ >= size)
            return 1;
    }
    return 0;
}

void tabu_search(instance *inst)
{
    int *tabu_list = (int *)calloc(inst->dimension, sizeof(int));
    // tabu_list[node] = 0: not a tabu node
    // tabu_list[node] > 0: tabu node
    int *temp_succ = (int *)calloc(inst->dimension, sizeof(int));
    // init temp_succ
    for (int i = 0; i < inst->dimension; i++)
        temp_succ[i] = inst->succ[i];
    double incumbent = inst->z_best;
    double best_inc = incumbent;
    int iter = 1;
    int best_iter = 1;
    int tenure = inst->dimension / 10;
    printf("INTENSIFICATION PHASE\n");
    while (iter < 500)
    {
        if (!(iter % 100))
        {
            if (tenure == 20)
            {
                tenure = inst->dimension / 10;
                printf("INTENSIFICATION PHASE\n");
            }
            else
            {
                tenure = 20;
                printf("DIVERSIFICATION PHASE\n");
            }
        }
        double min_delta = DBL_MAX;
        int a, b;

        for (int i = 0; i < inst->dimension; i++)
        {
            for (int j = 0; j < inst->dimension; j++)
            {
                if (i == j)
                    continue;
                if (temp_succ[i] == j || temp_succ[j] == i) // two contiguous edges
                    continue;
                int th = iter - tenure;
                if (th < 0)
                    th = 0;
                if (tabu_list[i] > th || tabu_list[temp_succ[i]] > th || tabu_list[j] > th || tabu_list[temp_succ[j]] > th)
                    continue;

                double delta = dist(i, j, inst) + dist(temp_succ[i], temp_succ[j], inst) - dist(i, temp_succ[i], inst) - dist(j, temp_succ[j], inst);

                if (delta <= min_delta)
                {
                    min_delta = delta;
                    a = i;
                    b = j;
                }
            }
        }

        if (min_delta == DBL_MAX)
            print_error("Error in tabu search delta computation");

        if (incumbent < best_inc)
        {
            best_iter = iter;
            inst->z_best = incumbent;
            // save best solution
            for (int j = 0; j < inst->dimension; j++)
                inst->succ[j] = temp_succ[j];
            save_and_plot_solution(inst, iter);
        }

        printf("incumbent = %f, delta = %f\n", incumbent, min_delta);

        int c = temp_succ[a];
        int d = temp_succ[b];

        if (min_delta > 0)
        {
            printf("*** worsening move ***\n");
            tabu_list[d] = iter;
        }

        // reverse the path from c to b
        if (reverse_successors(temp_succ, inst->dimension, c, b))
            print_error("Error in reverse_successors");

        // make the move
        temp_succ[a] = b;
        temp_succ[c] = d;
        // update incumbent
        incumbent = incumbent + min_delta;
        iter++;
    }
    printf("best incumbent = %f, best iter = %d\n", inst->z_best, best_iter);
    free(temp_succ);
    free(tabu_list);
}

// TODO : FINISH GENETIC ALGORITHM

double genetic(instance *inst, int size, int epochs) {

    int **population = (int **) calloc(size, sizeof(int **));
    double *fitness = (double *) calloc(size, sizeof(double *));
    int *individual = (int *) calloc(inst->dimension, sizeof(int));

    double *epoch_avg_fitness = (double *) calloc(epochs, sizeof(double *));
    double *epoch_champion = (double *) calloc(epochs, sizeof(double *));

    // Initialize population
    for (int i = 0; i < size; i++) {

        // Generate random individuals
        if (i < size * 0.8) fitness[i] = random_individuals(inst, i, individual, 1);
        else fitness[i] = random_individuals(inst, i, individual, 0);

        population[i] = individual;
    }

    // For each epoch
    for (int k = 0; k < epochs; k++) {
        epoch_avg_fitness[k] = compute_average_fitness(fitness, size);
        epoch_champion[k] = compute_champion(fitness, size);

        printf("%f\n", epoch_avg_fitness[k]);
        printf("%f\n", epoch_champion[k]);
    }

    free(population);

    return 0;

}

double random_individuals(instance *inst, int seed, int *individual, int optimize) {

    if (inst->param.grasp == 1) inst->z_best  = nearest_neighbours(inst, seed, inst->succ, inst->param.grasp_choices);
    else inst->z_best  = nearest_neighbours(inst, seed, inst->succ, 1);

    // Willing to optimize current chromosome
    if (optimize == 1) two_opt_v2(inst);

    double fitness = inst->z_best;

    convert_to_list(inst->succ, individual, inst->dimension);

    return fitness;

}

void convert_to_successor(int *list, int *succ, int n){
    // From node list to successor
    succ[0] = list[0];
    for (int k = 1; k < n; k++) {
        succ[list[k - 1]] = list[k];
    }
}

void convert_to_list(int *succ, int *list, int n){
    // From successor to node list
    list[0] = succ[0];
    int next = succ[0];
    for (int i = 1; i < n; i++) {
        list[i] = succ[next];
        next = succ[next];
    }

}

double compute_average_fitness(const double *fitness, int size) {

    double sum = 0;

    for (int i = 0; i < size; i++) sum += fitness[i];

    return sum / size;

}

double compute_champion(const double *fitness, int size) {

    int index = 0;
    double min = CPX_INFBOUND;

    for (int i = 1; i < size; i++) if (fitness[i] < min) index = i;

    return index;

}

void mutation(int *individual){}

void crossover(int *parent1, int *parent2, int *offspring){}

void repair_chromosome(int *offspring){}

void evolve(){}