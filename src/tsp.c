#include "../include/utils.h"

double dist(int i, int j, instance *inst) {
    double dx = inst->x[i] - inst->x[j];
    double dy = inst->y[i] - inst->y[j];
    if (!inst->integer_costs) return sqrt(dx * dx + dy * dy);
    int dis = sqrt(dx * dx + dy * dy) + 0.499999999;
    return dis + 0.0;
}

int TSPopt(instance *inst) {

    // Open CPLEX model
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

    build_model(env, lp, inst);

    // CPLEX's parameter setting
    CPXsetlogfilename(env, "../output/model_opt.txt", "w");    // Save log
    CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 1234);            // Use different seed
    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit);

    if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");

    // Use the optimal solution found by CPLEX
    int cols = CPXgetnumcols(env, lp);
    double *xstar = (double *) calloc(cols, sizeof(double));
    if (CPXgetx(env, lp, xstar, 0, cols - 1)) print_error("CPXgetx() error");

    for (int i = 0; i < inst->nodes; i++) {
        for (int j = i + 1; j < inst->nodes; j++) {
            if (xstar[position(i, j, inst)] > 0.5) {
                printf("  ... x(%3d,%3d) = 1\n", i + 1, j + 1);
            }
        }
    }
    free(xstar);

    // Free and close CPLEX model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return 0;
}

void build_model(CPXENVptr env, CPXLPptr lp, instance *inst) {

    double zero = 0.0;
    char binary = 'B';

    char **cname = (char **) calloc(1, sizeof(char *));     // set an array of strings to store the column names
    cname[0] = (char *) calloc(100, sizeof(char));

    // Add binary variables x(i,j) for i < j
    for (int i = 0; i < inst->nodes; i++) {
        for (int j = i + 1; j < inst->nodes; j++) {
            sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
            double obj = dist(i, j, inst); // cost == distance
            double lb = 0.0;
            double ub = 1.0;
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error(" wrong CPXnewcols on x var.s");
            if (CPXgetnumcols(env, lp) - 1 != position(i, j, inst)) print_error(" wrong position for x var.s");
        }
    }

    // Add the degree constraints
    for (int h = 0; h < inst->nodes; h++) {

        int row = CPXgetnumrows(env, lp);   // get the maximum number of row inside the model
        double rhs = 2.0;
        char sense = 'E';                   // E stands for equality constraint

        sprintf(cname[0], "degree(%d)", h + 1);

        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(" wrong CPXnewrows [degree]");
        for (int i = 0; i < inst->nodes; i++) {
            if (i == h) continue;
            if (CPXchgcoef(env, lp, row, position(i, h, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
        }
    }

//    if (verbose >= QUIET) {
//        printf("vgfaa");
        CPXwriteprob(env, lp, "../output/tsp_model.lp", "RLP" );
//    }

    free(cname[0]);
    free(cname);

}

int position(int i, int j, instance *inst) {

    if (i == j) print_error("Same indexes are not valid!");
    if (i < 0 || j < 0) print_error("Negative indexes are not valid!");
    if (i > inst->nodes || j > inst->nodes) print_error("Indexes exceeding the dimension are not valid!");

    if (i > j) return position(j, i, inst);
    else return i * inst->nodes + j - (i + 1) * (i + 2) / 2;
}
