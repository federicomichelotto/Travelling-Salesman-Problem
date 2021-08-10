#include "../include/utils.h"
#include "compact.c"
#include "math_heuristic.c"
#include "meta_heuristic.c"

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
            // reject the solution and adds one cut
            if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value))
                print_error("CPXcallbackrejectcandidate() error");
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
    else
    {
        // the candidate solution has not connected components but could have crossings... let's apply 2-opt
        double delta = two_opt_v2(inst, succ, 0);
        if (delta < 0)
        {
            objval += delta;
            // succ -> xstar
            int nnz = 0;
            int izero = 0;
            int *index = (int *)calloc(inst->cols, sizeof(int));
            double *xstar_succ = (double *)calloc(inst->cols, sizeof(double));

            for (int i = 0; i < inst->dimension; i++)
            {
                for (int j = i + 1; j < inst->dimension; j++)
                {
                    index[nnz] = xpos(i, j, inst);
                    if (j == succ[i] || i == succ[j])
                    {
                        xstar_succ[nnz++] = 1.0;
                    }
                    else
                        xstar_succ[nnz++] = 0.0;
                }
            }
            // sanity check
            if (nnz != inst->cols)
                print_error("Error in applying 2-opt in callback_candidate");

            if (CPXcallbackpostheursoln(context, nnz, index, xstar_succ, objval, CPXCALLBACKSOLUTION_CHECKFEAS))
                print_error("CPXcallbackpostheursoln() error");
            free(xstar_succ);
            free(index);
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
