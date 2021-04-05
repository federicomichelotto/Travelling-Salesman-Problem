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

void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp) // build succ() and comp() wrt xstar()...
{

    *ncomp = 0;
    for ( int i = 0; i < inst->dimension; i++ )
    {
        succ[i] = -1;
        comp[i] = -1;
    }

    for ( int start = 0; start < inst->dimension; start++ )
    {
        if ( comp[start] >= 0 ) continue;  // node "start" was already visited, just skip it

        // a new component is found
        (*ncomp)++;
        int i = start;
        int done = 0;
        while ( !done )  // go and visit the current component
        {
            comp[i] = *ncomp;
            done = 1;
            for ( int j = 0; j < inst->dimension; j++ )
            {
                if ( i != j && xstar[xpos(i,j,inst)] > 0.5 && comp[j] == -1 ) // the edge [i,j] is selected in xstar and j was not visited before
                {
                    succ[i] = j;
                    i = j;
                    done = 0;
                    break;
                }
            }
        }
        succ[i] = start;  // last arc to close the cycle

        // go to the next component...
    }
}

// Kruskal algorithm to find connected components
int findConnectedComponents(instance *inst, int* components, int* successors, const double *xstar) {

    int counter = 0;

    // Initialize components and successors
    for (int i = 0; i < inst->dimension; i++) {
        components[i] = i;
        successors[i] = -1;
    }

    // Found connected components
    for (int i = 0; i < inst->dimension; i++) {
        for (int j = i+1; j < inst->dimension; j++) {
            int pos = xpos(i,j,inst);
            if (xstar[pos] == 1) {
                if (components[i] != components[j]) {
                    int c1 = components[i];
                    int c2 = components[j];
                    for (int v = 0; v < inst->dimension; v++) if (components[v] == c2) components[v] = c1;
                }
            }
        }
    }

    // Count how many components are
    for (int i = 0; i < inst->dimension; i++) {

        int inside = 0;

        for (int j = 0; j < counter +1 ; j++) {
            if (successors[j] == components[i]) {
                inside = 1;
                break;
            }
        }

        if (!inside) {
            successors[counter] = components[i];
            counter++;
        }

    }

    return  counter;

}

int TSPopt(instance *inst)
{

    // Open CPLEX model
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

    build_model(env, lp, inst);

    char path[1000];
    generate_path(path, "output", "model", model_name[inst->model_type], inst->param.name, inst->param.seed, "log");

    // CPLEX's parameter setting
    CPXsetlogfilename(env, path, "w");                           // Save log
    CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->param.seed); // Set seed
    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit);
    CPXsetintparam(env, CPXPARAM_Parallel, CPX_PARALLEL_OPPORTUNISTIC); // Set opportunistic mode

    // CPLEX's precision setting
    CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);
    CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);

    if (CPXmipopt(env, lp))
        print_error("CPXmipopt() error");

    // Use the optimal solution found by CPLEX
    int cols = CPXgetnumcols(env, lp);
    double *xstar = (double *)calloc(cols, sizeof(double));
    if (CPXgetx(env, lp, xstar, 0, cols - 1))
        print_error("CPXgetx() error");

    inst->n_edges = 0;

    if (inst->model_type == 0 || inst->model_type == 6) // undirected graph
    {
        for (int i = 0; i < inst->dimension; i++)
        {
            for (int j = i + 1; j < inst->dimension; j++)
            {
                if (xstar[xpos(i, j, inst)] > 0.5)
                {
                    printf("  ... x(%3d,%3d) = 1\n", i + 1, j + 1);
                    inst->edges[inst->n_edges].dist = dist(i, j, inst);
                    inst->edges[inst->n_edges].prev = i;
                    inst->edges[inst->n_edges].next = j;
                    if (++inst->n_edges > inst->dimension)
                        print_error("more edges than nodes, not a hamiltonian tour.");
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < inst->dimension; i++)
        {
            for (int j = 0; j < inst->dimension; j++)
            {
                if (xstar[xpos_dir(i, j, inst)] > 0.5)
                {
                    printf("  ... x(%3d,%3d) = 1\n", i + 1, j + 1);
                    inst->edges[inst->n_edges].dist = dist(i, j, inst);
                    inst->edges[inst->n_edges].prev = i;
                    inst->edges[inst->n_edges].next = j;
                    if (++inst->n_edges > inst->dimension)
                        print_error("more edges than nodes, not a hamiltonian tour.");
                }
            }
        }
    }

    if (inst->n_edges != inst->dimension)
        print_error("not a tour.");
    CPXgetobjval(env, lp, &inst->z_best);      // Best objective value
    CPXgetbestobjval(env, lp, &inst->best_lb); // Best lower bound

    printf("Best objective value: %lf\n", inst->z_best);
    printf("Best lower bound: %lf\n", inst->best_lb);

    free(xstar);

    // Free and close CPLEX model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return 0;
}

void build_model(CPXENVptr env, CPXLPptr lp, instance *inst)
{
    if (inst->model_type == 0) // basic model (no SEC) for undirected graphs
    {
        basic_model_no_sec(env, lp, inst);
    }
    else if (inst->model_type == 1) // TMZ with static constraints
    {
        TMZ_static(env, lp, inst);
    }
    else if (inst->model_type == 2) // TMZ with static constraints
    {
//        TMZ_static_mod(env, lp, inst);
    }
    else if (inst->model_type == 3) // TMZ with lazy constraints
    {
        TMZ_lazy(env, lp, inst);
    }
    else if (inst->model_type == 4) // TMZ with lazy constraints and subtour elimination
    {
        TMZ_lazy_sec(env, lp, inst);
    }
    else if (inst->model_type == 5) // GG
    {
        GG(env, lp, inst);
    }
    else if (inst->model_type == 6) // Bender
    {
        bender(env, lp, inst);
    }
    else
    {
        printf("ERROR: Model type %d not available.\n", inst->model_type);
        print_error("Model type.");
    }

    char path[1000];
    generate_path(path, "output", "model", model_name[inst->model_type], inst->param.name, inst->param.seed, "lp");
    // path : "../output/model_[type]_[name].lp"

    CPXwriteprob(env, lp, path, NULL);
}

int xpos(int i, int j, instance *inst)
{

    if (i == j)
        print_error("Same indexes are not valid!");
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

void TMZ_static(CPXENVptr env, CPXLPptr lp, instance *inst)
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

    // Add static TMZ constraints: 1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1, for each arc (i,j) not touching node 0
    double M = inst->dimension - 1;
    double rhs = M - 1;
    char sense = 'L';                         // L stands for less than or equal
    for (int i = 1; i < inst->dimension; i++) // *** faster including i=0 ? yes -> TMZ_static_mod ***
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

//void TMZ_static_mod(CPXENVptr env, CPXLPptr lp, instance *inst)
//{
//    char binary = 'B';  // B => binary variable flag
//    char integer = 'I'; // I => integer variable flag
//
//    // cname: columns' names (column = variable)
//    char **cname = (char **)calloc(1, sizeof(char *)); // array of strings to store the column names
//    cname[0] = (char *)calloc(100, sizeof(char));
//
//    // rname: rows' names (row = constraint)
//    char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
//    rname[0] = (char *)calloc(100, sizeof(char));
//    // Add binary variables x(i,j) for each (i,j)
//    for (int i = 0; i < inst->dimension; i++)
//    {
//        for (int j = 0; j < inst->dimension; j++)
//        {
//            sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
//            double obj = dist(i, j, inst); // cost == distance
//            double lb = 0.0;
//            double ub = 1.0;
//            if (i == j)
//                ub = 0.0;
//            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
//                print_error(" wrong CPXnewcols on x var.s");
//            if (CPXgetnumcols(env, lp) - 1 != xpos_dir(i, j, inst))
//                print_error("[position_d] wrong position for x var.s");
//        }
//    }
//
//    // Add u-variables one for each node ( u_0 = 0 )
//    for (int i = 0; i < inst->dimension; i++)
//    {
//        sprintf(cname[0], "u(%d)", i + 1);
//        double obj = 0.0;
//        double lb = 0.0;
//        double ub = inst->dimension - 1;
//        if (i == 0)
//            ub = 0;
//        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname))
//            print_error(" wrong CPXnewcols on u var.s");
//        if (CPXgetnumcols(env, lp) - 1 != upos(i, inst))
//            print_error("[position_d] wrong position for u var.s");
//    }
//
//    // Add the in-degree constraints
//    for (int h = 0; h < inst->dimension; h++)
//    {
//
//        int row = CPXgetnumrows(env, lp); // get the maximum number of row inside the model
//        double rhs = 1.0;
//        char sense = 'E'; // E stands for equality constraint
//        sprintf(rname[0], "in_degree(%d)", h + 1);
//        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
//            print_error("wrong CPXnewrows [degree]");
//        for (int i = 0; i < inst->dimension; i++)
//        {
//            if (CPXchgcoef(env, lp, row, xpos_dir(i, h, inst), 1.0))
//                print_error("wrong CPXchgcoef [degree]");
//        }
//    }
//
//    // Add the out-degree constraints
//    for (int h = 0; h < inst->dimension; h++)
//    {
//
//        int row = CPXgetnumrows(env, lp); // get the maximum number of row inside the model
//        double rhs = 1.0;
//        char sense = 'E'; // E stands for equality constraint
//        sprintf(rname[0], "out_degree(%d)", h + 1);
//        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
//            print_error("wrong CPXnewrows [degree]");
//        for (int i = 0; i < inst->dimension; i++)
//        {
//            if (CPXchgcoef(env, lp, row, xpos_dir(h, i, inst), 1.0))
//                print_error("wrong CPXchgcoef [degree]");
//        }
//    }
//
//    // Add static TMZ constraints: 1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1, for each arc (i,j) not touching node 0
//    double M = inst->dimension - 1;
//    double rhs = M - 1;
//    char sense = 'L';                         // L stands for less than or equal
//    for (int i = 0; i < inst->dimension; i++) // *** faster including i=0 ? yes ***
//    {
//        for (int j = 1; j < inst->dimension; j++)
//        {
//            if (i == j)
//                continue;
//            int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
//            sprintf(rname[0], "u_consistency for arc (%d,%d)", i + 1, j + 1);
//            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
//                print_error("wrong CPXnewrows [degree]");
//            if (CPXchgcoef(env, lp, row, upos(i, inst), 1.0)) // 1.0 * u_i
//                print_error("wrong CPXchgcoef [degree]");
//            if (CPXchgcoef(env, lp, row, upos(j, inst), -1.0)) // - 1.0 * u_j
//                print_error("wrong CPXchgcoef [degree]");
//            if (CPXchgcoef(env, lp, row, xpos_dir(i, j, inst), M)) // M * x_ij
//                print_error("wrong CPXchgcoef [degree]");
//        }
//    }
//    free(cname[0]);
//    free(cname);
//    free(rname[0]);
//    free(rname);
//}

void TMZ_lazy(CPXENVptr env, CPXLPptr lp, instance *inst)
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

void TMZ_lazy_sec(CPXENVptr env, CPXLPptr lp, instance *inst)
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
        }
    }
    /*
    for (int i = 0; i < inst->dimension; i++)
    { // x(i, j) + x(j, i) <= 1 for every i < j
        for (int j = i + 1; j < inst->dimension; j++)
        {

            int lastrow = CPXgetnumrows(env, lp);
            double rhs = 1.0;
            char sense = 'L';

            sprintf(cname[0], "no_binary_loops_LAZY(%d, %d)", i + 1, j + 1);
            int *beg = (int *)calloc(2, sizeof(int));
            int *ind = (int *)calloc(2, sizeof(int));
            double *val = (double *)calloc(2, sizeof(double));

            ind[0] = xpos_dir(i, j, inst);
            ind[1] = xpos_dir(j, i, inst);
            val[0] = 1.0;
            val[1] = 1.0;
            beg[0] = 0;
            beg[1] = inst->dimension;
            if (CPXaddlazyconstraints(env, lp, 1, 2, &rhs, &sense, beg, ind, val, cname))
            {
                print_error("WRONG LAZY [y2]");
            }
        }
    }
*/
    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}

void GG(CPXENVptr env, CPXLPptr lp, instance *inst)
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
            if (CPXchgcoef(env, lp, row, xpos_dir(i, j, inst), 2 - inst->dimension)) // (2 - nnodes) * y_ij
                print_error("wrong CPXchgcoef [degree]");
        }
    }

    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}

void bender(CPXENVptr env, CPXLPptr lp, instance *inst) {
    char binary = 'B'; // B => binary variable flag
    // cname: columns' names (column = variable)
    char **cname = (char **) calloc(1, sizeof(char *)); // array of strings to store the column names
    cname[0] = (char *) calloc(100, sizeof(char));

    // rname: rows' names (row = constraint)
    char **rname = (char **) calloc(1, sizeof(char *)); // array of strings to store the row names
    rname[0] = (char *) calloc(100, sizeof(char));
    // Add binary variables x(i,j) for i < j
    for (int i = 0; i < inst->dimension; i++) {
        for (int j = i + 1; j < inst->dimension; j++) {
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
    for (int h = 0; h < inst->dimension; h++) {

        int row = CPXgetnumrows(env, lp); // get the maximum number of row inside the model
        double rhs = 2.0;
        char sense = 'E'; // E stands for equality constraint

        sprintf(rname[0], "degree(%d)", h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [degree]");
        for (int i = 0; i < inst->dimension; i++) {
            if (i == h)
                continue;
            if (CPXchgcoef(env, lp, row, xpos(i, h, inst), 1.0))
                print_error("[position_u] wrong CPXchgcoef [degree]");
        }
    }

    // Application of Benders method

    int done = 0;
    int i = 0;      // number of iteration
    int c = 0;      // number of connected components

    while (!done) {

        CPXmipopt(env, lp);

        int ncols = CPXgetnumcols(env, lp);
        double *xstar = (double *) calloc(ncols, sizeof(double));

        int status = CPXgetx(env, lp, xstar, 0, ncols - 1);
        if (status) { print_error_status("Failed to obtain the values in LOOP method", status); }

        int *comp = (int *) calloc(inst->dimension, sizeof(int));
        int *succ = (int *) calloc(inst->dimension, sizeof(int));

        // Retrieve the number of connected components so far
        c = findConnectedComponents(inst, comp, succ, xstar);
//        build_sol(xstar,inst, succ, comp, &c);

        printf("ITERATION: %d\tCONNECTED COMPONENTS FOUND: %d\n", i++, c);

        if (c == 1) {
            // If exactly one component is found end the loop
            done = 1;
        } else {
            // If more than one component is found add SEC to each component
            char sense = 'L';
            for (int k = 0; k < c; k++) {

                int nnz = 0;
                double rhs = -1.0;

                int *index = (int *) calloc(ncols, sizeof(int));
                for (int h = 0; h < inst->dimension; h++) {

                    if (comp[h] != succ[k]) continue;
                    else rhs++;

                    for (int j = h + 1; j < inst->dimension; j++) {
                        if (comp[j] == succ[k]) index[nnz++] = xpos(h, j, inst);
                    }
                }

                sprintf(cname[0], "SEC(%d)", k + 1);
                int lastrow = CPXgetnumrows(env, lp);

                status = CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname);
                if (status) print_error_status(" wrong CPXnewrows SEC", status);

                for (int z = 0; z < nnz; z++) {
                    status = CPXchgcoef(env, lp, lastrow, index[z], 1.0);
                    if (status) print_error_status(" wrong CPXchgcoef in SEC", status);
                }
                free(index);
            }
        }
        free(comp);
        free(succ);
        free(xstar);
    }

    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}


void print_time_csv()
{
}
