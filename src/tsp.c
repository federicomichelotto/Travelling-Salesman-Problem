#include "../include/utils.h"

double dist(int i, int j, instance *inst) {

    double distance = INFINITY;

    if (strncmp(inst->param.weight_type, "GEO", 3) == 0) {
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

        distance = (int)(RRR * acos(0.5*((1.0 + q1)*q2 - (1.0 - q1)*q3)) + 1.0);

    } else if (strncmp(inst->param.weight_type, "EUC_2D", 6) == 0) {
        double dx = inst->nodes[i].x - inst->nodes[j].x;
        double dy = inst->nodes[i].y - inst->nodes[j].y;
        if (!inst->integer_costs)
            return sqrt(dx * dx + dy * dy);

        int dis = sqrt(dx * dx + dy * dy) + 0.499999999;
        distance = dis + 0.0;
    } else if (strncmp(inst->param.weight_type, "ATT", 3) == 0) {
        double dx = inst->nodes[i].x - inst->nodes[j].x;
        double dy = inst->nodes[i].y - inst->nodes[j].y;

        int dis = sqrt((dx * dx + dy * dy) / 10.0) + 0.499999999;
        distance = dis + 0.0;
    }

    return distance;
}

int TSPopt(instance *inst) {

    // Open CPLEX model
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

    build_model(env, lp, inst);

    char path[1000];
    generate_path(path, "output", "log", inst->param.name, "txt");

    // CPLEX's parameter setting
    CPXsetlogfilename(env, path, "w"); // Save log
    CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 1234);        // Use different seed
    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit);

    if (CPXmipopt(env, lp))
        print_error("CPXmipopt() error");

    // Use the optimal solution found by CPLEX
    int cols = CPXgetnumcols(env, lp);
    double *xstar = (double *) calloc(cols, sizeof(double));
    if (CPXgetx(env, lp, xstar, 0, cols - 1))
        print_error("CPXgetx() error");

    inst->n_edges = 0;
    for (int i = 0; i < inst->dimension; i++) {
        for (int j = i + 1; j < inst->dimension; j++) {
            if (xstar[position(i, j, inst)] > 0.5) {
                printf("  ... x(%3d,%3d) = 1\n", i + 1, j + 1);
                inst->edges[inst->n_edges].dist = dist(i,j,inst);
                inst->edges[inst->n_edges].prev = i;
                inst->edges[inst->n_edges].next = j;
                if (++inst->n_edges > inst->dimension) print_error("more edges than nodes, not a tour.");
            }
        }
    }

    if (inst->n_edges != inst->dimension)
        print_error("not a tour.");

    free(xstar);

    // Free and close CPLEX model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return 0;
}

void build_model(CPXENVptr env, CPXLPptr lp, instance *inst) {

    double zero = 0.0;
    char binary = 'B'; // B => binary variables flag

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
            if (CPXgetnumcols(env, lp) - 1 != position(i, j, inst))
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
            print_error(" wrong CPXnewrows [degree]");
        for (int i = 0; i < inst->dimension; i++) {
            if (i == h)
                continue;
            if (CPXchgcoef(env, lp, row, position(i, h, inst), 1.0))
                print_error(" wrong CPXchgcoef [degree]");
        }
    }

    char path[1000];
    generate_path(path, "output", "model", inst->param.name, "lp");
    CPXwriteprob(env, lp, path, "RLP");

    free(cname[0]);
    free(cname);
}

int position(int i, int j, instance *inst) {

    if (i == j)
        print_error("Same indexes are not valid!");
    if (i < 0 || j < 0)
        print_error("Negative indexes are not valid!");
    if (i > inst->dimension || j > inst->dimension)
        print_error("Indexes exceeding the dimension are not valid!");

    if (i > j)
        return position(j, i, inst);
    else
        return i * inst->dimension + j - (i + 1) * (i + 2) / 2;
}
