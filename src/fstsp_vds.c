#include "../include/utils.h"
#include "compact.c"

double dist(int i, int j, instance *inst)
{
    return dist_GEO(i, j, inst);
}

double dist_GEO(int i, int j, instance *inst)
{
    double radius_of_earth = 6378100.0; // [meters]

    // deg -> rad
    double lat1 = inst->nodes[i].x * M_PI / 180;
    double lat2 = inst->nodes[j].x * M_PI / 180;
    double long1 = inst->nodes[i].y * M_PI / 180;
    double long2 = inst->nodes[j].y * M_PI / 180;

    double distance = 2 * radius_of_earth * asin(sqrt(pow(sin((lat2 - lat1) / 2), 2) + cos(lat1) * cos(lat2) * pow(sin((long2 - long1) / 2), 2)));
    return distance;
}

double dist_EUC_2D(int i, int j, instance *inst)
{
    double dx = inst->nodes[i].x - inst->nodes[j].x;
    double dy = inst->nodes[i].y - inst->nodes[j].y;
    return sqrt(dx * dx + dy * dy);
}

// CPLEX vector will be organized as follow:
//      xTruck_pos | xDrone_pos | yT_pos | yD_pos | yC_pos

// retrieve the position of the decision var. x^T(i,j) in the cplex vector
int xTruck_pos(int i, int j, instance *inst)
{

    if (i < 0 || j < 0)
        print_error("Error in xTruck_pos: negative indexes are not valid!");
    if (i >= inst->dimension || j >= inst->dimension)
        print_error("Error in xTruck_pos: indexes exceeding the dimension are not valid!");
    return i * inst->dimension + j;
}

// retrieve the position of the decision var. x^D(i,j) in the cplex vector
int xDrone_pos(int i, int j, instance *inst)
{
    if (i < 0 || j < 0)
        print_error("Error in xDrone_pos: negative indexes are not valid!");
    if (i > inst->dimension || j > inst->dimension)
        print_error("Error in xDrone_pos: indexes exceeding the dimension are not valid!");
    return xTruck_pos(inst->dimension - 1, inst->dimension - 1, inst) + 1 + i * inst->dimension + j;
}

// retrieve the position of the decision var. Y^T(i) in the cplex vector
int yT_pos(int i, instance *inst)
{
    if (i < 0 || i >= inst->dimension)
        print_error("Error in yT_pos: i must be in the range [0, inst->dimension]");
    return xDrone_pos(inst->dimension - 1, inst->dimension - 1, inst) + 1 + i;
}

// retrieve the position of the decision var. Y^D(i) in the cplex vector
int yD_pos(int i, instance *inst)
{
    if (i < 0 || i >= inst->dimension)
        print_error("Error in yD_pos: i must be in the range [0, inst->dimension]");
    return yT_pos(inst->dimension - 1, inst) + 1 + i;
}

// retrieve the position of the decision var. Y^C(i) in the cplex vector
int yC_pos(int i, instance *inst)
{
    if (i < 0 || i >= inst->dimension)
        print_error("Error in yC_pos: i must be in the range [0, inst->dimension]");
    return yD_pos(inst->dimension - 1, inst) + 1 + i;
}

// retrieve the position of the decision var. a_i in the cplex vector
// arrival time of the truck or the drone (or both) at node i ∈ V
int a_pos(int i, instance *inst)
{
    if (i < 0 || i >= inst->dimension)
        print_error("Error in a_pos: i must be in the range [0, inst->dimension]");
    return yC_pos(inst->dimension - 1, inst) + 1 + i;
}

// retrieve the position of the decision var. z(i,j,k) in the cplex vector
int z_pos(int i, int j, int k, instance *inst)
{
    if (i < 0 || j < 0 || k < 0)
        print_error("Error in z_pos: negative indexes are not valid!");
    if (i >= inst->dimension || j >= inst->dimension || k >= inst->dimension)
        print_error("Error in z_pos: indexes exceeding the dimension are not valid!");
    return a_pos(inst->dimension - 1, inst) + 1 + i * inst->dimension * inst->dimension + j * inst->dimension + k;
}

// redundant decision variable u_i
// u_i = position of node i in the path
// retrieve the position of the decision var. u_i in the cplex vector
int u_pos(int i, instance *inst)
{
    if (i < 0 || i >= inst->dimension)
        print_error("Error in a_pos: i must be in the range [0, inst->dimension]");
    return z_pos(inst->dimension - 1, inst->dimension - 1, inst->dimension - 1, inst) + 1 + i;
}

void convert_seq2succ(instance *inst, int *seq, int *succ)
{
    for (int i = 0; i < inst->dimension; i++)
        succ[i] = -1;

    for (int i = 0; i < inst->dimension; i++)
    {
        succ[seq[i]] = seq[i + 1];
        if (seq[i+1] == inst->dimension - 1)
            break;
    }
}

//// to fix
double gather_solution(instance *inst, const double *xstar)
{

    // TRUCK PATH
    int n_edges = 0;
    for (int i = 0; i < inst->dimension - 1; i++)
    {
        if (inst->best_sol[yD_pos(i, inst)] > 0.5)
        { // i is a drone customer
            inst->succ_truck[i] = -1;
            continue;
        }
        for (int j = 1; j < inst->dimension; j++)
        {
            if (xstar[xTruck_pos(i, j, inst)] > 0.5)
            {
                inst->succ_truck[i] = j;
                if (++n_edges > inst->dimension - 1)
                    printf("WARNING: Truck path is not feasible");
            }
        }
    }

    // DRONE PATH
    n_edges = 0;
    for (int i = 0; i < inst->dimension - 1; i++)
    {
        if (inst->best_sol[yT_pos(i, inst)] > 0.5)
        { // i is a truck customer
            inst->succ_drone[i] = -1;
            continue;
        }
        for (int j = 1; j < inst->dimension; j++)
        {
            if (xstar[xDrone_pos(i, j, inst)] > 0.5)
            {
                inst->succ_drone[i] = j;
                if (++n_edges > inst->dimension - 1)
                    printf("WARNING: Truck path is not feasible");
            }
        }
    }
}

void build_model(CPXENVptr env, CPXLPptr lp, instance *inst, double bigM)
{
    char binary = 'B';     // B => binary variable flag
    char integer = 'I';    // I => integer variable flag
    char continuous = 'C'; // C => continuous variable flag

    // cname: columns' names (column = variable)
    char **cname = (char **)calloc(1, sizeof(char *)); // array of strings to store the column names
    cname[0] = (char *)calloc(100, sizeof(char));

    // rname: rows' names (row = constraint)
    char **rname = (char **)calloc(1, sizeof(char *)); // array of strings to store the row names
    rname[0] = (char *)calloc(100, sizeof(char));

    // ************************ VARIABLES ************************ //

    // Add a binary variable 'x_T' for each edge (i,j)
    // x_T(i,j) = 1 if truck crosses edge (i,j)
    // x_T(i,j) = 0 otherwise
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            sprintf(cname[0], "x_T(%d,%d)", i, j);
            double obj = 0.0;
            double lb = 0.0;
            double ub = 1.0;
            if (i == j)
                ub = 0.0;
            if (j == 0) // no ingoing truck edges to 0 (starting depot) are admissible
                ub = 0.0;
            if (i == (inst->dimension - 1)) // no outgoing truck edges from (N-1) (returning depot) are admissible
                ub = 0.0;
            if (inst->truck_times[i][j] < 1e-6) // time[i][j] = 0 => edge (i,j) is not valid (traversable)
                ub = 0.0;
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
                print_error(" wrong CPXnewcols on x_T variables");
            if (CPXgetnumcols(env, lp) - 1 != xTruck_pos(i, j, inst))
                print_error("wrong position for x_T variables");
        }
    }

    // Add a binary variable 'x_D' for each edge (i,j)
    // x_D(i,j) = 1 if drone crosses edge (i,j)
    // x_D(i,j) = 0 otherwise
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            sprintf(cname[0], "x_D(%d,%d)", i, j);
            double obj = 0.0;
            double lb = 0.0;
            double ub = 1.0;
            if (i == j)
                ub = 0.0;
            if (j == 0) // no ingoing drone edges to 0 (starting depot) are admissible
                ub = 0.0;
            if (i == (inst->dimension - 1)) // no outgoing drone edges from (N-1) (returning depot) are admissible
                ub = 0.0;
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
                print_error("wrong CPXnewcols on x_D variables");
            if (CPXgetnumcols(env, lp) - 1 != xDrone_pos(i, j, inst))
                print_error("wrong position for x_D variables");
        }
    }

    // Add a binary variable 'y_T' for each node
    // y_T(i) = 1 if node i is a truck (only) customer
    for (int i = 0; i < inst->dimension; i++)
    {
        sprintf(cname[0], "y_T(%d)", i);
        double obj = 0.0;
        double lb = 0.0;
        double ub = 1.0;
        if (i == 0 || i == inst->dimension - 1)
            ub = 0.0;
        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
            print_error("wrong CPXnewcols on y_T variables");
        if (CPXgetnumcols(env, lp) - 1 != yT_pos(i, inst))
            print_error("wrong position for y_T variables");
    }

    // Add a binary variable 'y_D' for each node
    // y_D(i) = 1 if node i is a drone (only) customer
    for (int i = 0; i < inst->dimension; i++)
    {
        sprintf(cname[0], "y_D(%d)", i);
        double obj = 0.0;
        double lb = 0.0;
        double ub = 1.0;
        if (i == 0 || i == inst->dimension - 1)
            ub = 0.0;
        if (inst->nodes[i].weight > WEIGHT_LIMIT)
            ub = 0.0;
        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
            print_error("wrong CPXnewcols on y_D variables");
        if (CPXgetnumcols(env, lp) - 1 != yD_pos(i, inst))
            print_error("wrong position for y_D variables");
    }

    // Add a binary variable 'y_C' for each node
    // y_C(i) = 1 if node i is a combined customer
    for (int i = 0; i < inst->dimension; i++)
    {
        sprintf(cname[0], "y_C(%d)", i);
        double obj = 0.0;
        double lb = 0.0;
        double ub = 1.0;
        if (i == 0 || i == inst->dimension - 1)
            lb = 1.0;
        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
            print_error("wrong CPXnewcols on y_C variables");
        if (CPXgetnumcols(env, lp) - 1 != yC_pos(i, inst))
            print_error("wrong position for y_C variables");
    }

    // Add a continuous variable 'a' for each node
    // *** a(i) = arrival time at node i of the truck or the drone (or both) ***
    // a(i) = departure time from node i

    for (int i = 0; i < inst->dimension; i++)
    {
        sprintf(cname[0], "a(%d)", i);
        double obj = 0.0;
        double lb = 0.0;
        double ub = DBL_MAX;
        if (bigM > 0)
            ub = bigM;
        if (i == 0)
            ub = 0.0;
        if (i == inst->dimension - 1)
            obj = 1.0;

        if (CPXnewcols(env, lp, 1, &obj, &lb, NULL, &continuous, cname))
            print_error("wrong CPXnewcols on 'a' variables");

        if (CPXgetnumcols(env, lp) - 1 != a_pos(i, inst))
            print_error("wrong position for 'a' variables");
    }

    // Add a binary variable 'z' for each drone leg i-->j-->k
    // z(i,j,k) = 1 if drone leg i-->j-->k is selected
    // z(i,j,k) = 0 otherwise
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            for (int k = 0; k < inst->dimension; k++)
            {

                sprintf(cname[0], "z(%d,%d,%d)", i, j, k);
                double obj = 0.0;
                double lb = 0.0;
                double ub = 1.0;
                if (i == k)
                    ub = 0.0;
                if (j == i || j == k) // the intermediate node must be different from the starting and ending node
                    ub = 0.0;
                if (j == 0 || j == inst->dimension - 1) // the intermediate node cannot be a depot
                    ub = 0.0;
                if (inst->min_time_drone[i][j][k] < 1e-6) // time[i][j][k] = 0 => drone leg i-->j-->k is not valid (traversable)
                    ub = 0.0;
                if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
                    print_error("wrong CPXnewcols on z variables");
                if (CPXgetnumcols(env, lp) - 1 != z_pos(i, j, k, inst))
                    print_error("wrong position for z variables");
            }
        }
    }

    // Add a decision variable 'u' for each node
    // *** u(i) = poisition of node i in the path ***
    // for (int i = 0; i < inst->dimension; i++)
    // {
    //     sprintf(cname[0], "u(%d)", i);
    //     double obj = 0.0;
    //     double lb = 0.0;
    //     double ub = inst->dimension - 1;
    //     if (i == 0)
    //         ub = 0.0;

    //     if (CPXnewcols(env, lp, 1, &obj, &lb, NULL, &integer, cname))
    //         print_error("wrong CPXnewcols on 'u' variables");

    //     if (CPXgetnumcols(env, lp) - 1 != u_pos(i, inst))
    //         print_error("wrong position for 'u' variables");
    // }

    // bigM definition
    double M = bigM;
    if (M == 0)
    {
        for (int i = 0; i < inst->dimension - 1; i++)
        {
            for (int j = 0; j < inst->dimension - 1; j++)
            {
                if (i == j)
                    continue;
                if (inst->truck_times[i][j] > M)
                    M = inst->truck_times[i][j];
            }
        }
        M = M * (inst->dimension - 1);
    }

    // ************************ CONSTRAINTS ************************ //

    // nodes classification:
    // 0 = starting depot
    // [1, N-1] = customers
    // N = returning depot

    // ************ CUSTOMER CONSTRAINTS ************ //

    // Add "y" constraints
    // each node node must be a [TRUCK|DRONE|COMBINED] customer
    // *
    for (int i = 1; i < inst->dimension - 1; i++)
    {
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        double rhs = 1.0;
        char sense = 'E'; // E stands for equality constraint
        sprintf(rname[0], "y(%d)", i);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [y]");
        if (CPXchgcoef(env, lp, row, yT_pos(i, inst), 1.0))
            print_error("wrong CPXchgcoef [y]");
        if (CPXchgcoef(env, lp, row, yD_pos(i, inst), 1.0))
            print_error("wrong CPXchgcoef [y]");
        if (CPXchgcoef(env, lp, row, yC_pos(i, inst), 1.0))
            print_error("wrong CPXchgcoef [y]");
    }

    // ************ TRUCK CONSTRAINTS ************ //

    // Add the "in-flow = out-flow TRUCK" constraint for each customer node
    // *
    for (int i = 1; i < inst->dimension - 1; i++)
    {
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        double rhs = 0.0;
        char sense = 'E'; // E stands for equality constraint
        sprintf(rname[0], "in_flow = out_flow TRUCK(%d)", i);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [in_flow = out_flow TRUCK]");
        // in-flow
        for (int j = 0; j < inst->dimension; j++)
        {
            if (CPXchgcoef(env, lp, row, xTruck_pos(j, i, inst), 1.0))
                print_error("wrong CPXchgcoef [in-flow = out-flow TRUCK]");
        }
        // out-flow
        for (int j = 0; j < inst->dimension; j++)
        {
            if (CPXchgcoef(env, lp, row, xTruck_pos(i, j, inst), -1.0))
                print_error("wrong CPXchgcoef [in-flow = out-flow TRUCK]");
        }
    }
    // Add "truck_initial_flow" constraint
    // *
    {
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        double rhs = 1.0;
        char sense = 'E'; // E stands for equality constraint
        sprintf(rname[0], "truck_initial_flow");
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [truck_initial_flow]");
        // truck out-flow starting depot
        for (int i = 1; i < inst->dimension; i++)
        {
            if (CPXchgcoef(env, lp, row, xTruck_pos(0, i, inst), 1.0))
                print_error("wrong CPXchgcoef [truck_initial_flow]");
        }
    }

    // Add "truck_ending_flow" constraint
    // *
    {
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        double rhs = 1.0;
        char sense = 'E'; // E stands for equality constraint
        sprintf(rname[0], "truck_ending_flow");
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [truck_ending_flow]");
        // truck in-flow node final depot
        for (int i = 0; i < inst->dimension - 1; i++)
        {
            if (CPXchgcoef(env, lp, row, xTruck_pos(i, inst->dimension - 1, inst), 1.0))
                print_error("wrong CPXchgcoef [truck_ending_flow]");
        }
    }

    // Add the "truck_out_flow" constraint for each customer node
    // *
    for (int i = 1; i < inst->dimension - 1; i++)
    {
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        double rhs = 0.0;
        char sense = 'E'; // E stands for equality constraint
        sprintf(rname[0], "truck_out_flow(%d)", i);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [truck_out_flow]");
        // out-flow
        for (int j = 0; j < inst->dimension; j++)
        {
            if (CPXchgcoef(env, lp, row, xTruck_pos(i, j, inst), 1.0))
                print_error("wrong CPXchgcoef [truck_out_flow]");
        }
        if (CPXchgcoef(env, lp, row, yT_pos(i, inst), -1.0))
            print_error("wrong CPXchgcoef [truck_out_flow]");
        if (CPXchgcoef(env, lp, row, yC_pos(i, inst), -1.0))
            print_error("wrong CPXchgcoef [truck_out_flow]");
    }

    // ************ DRONE CONSTRAINTS ************ //

    // Add the "in-flow = out-flow DRONE" constraints for each customer
    // *
    for (int i = 1; i < inst->dimension - 1; i++)
    {
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        double rhs = 0.0;
        char sense = 'E'; // E stands for equality constraint
        sprintf(rname[0], "in_flow = out_flow DRONE(%d)", i);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [in-flow = out-flow DRONE]");
        // in-flow
        for (int j = 0; j < inst->dimension; j++)
        {
            if (CPXchgcoef(env, lp, row, xDrone_pos(j, i, inst), 1.0))
                print_error("wrong CPXchgcoef [in-flow = out-flow DRONE]");
        }
        // out-flow
        for (int j = 0; j < inst->dimension; j++)
        {
            if (CPXchgcoef(env, lp, row, xDrone_pos(i, j, inst), -1.0))
                print_error("wrong CPXchgcoef [in-flow = out-flow DRONE]");
        }
    }

    // Add "drone_initial_flow" constraint
    // *
    {
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        double rhs = 1.0;
        char sense = 'E'; // E stands for equality constraint
        sprintf(rname[0], "drone_initial_flow");
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [drone_initial_flow]");
        // drone out-flow from the starting depot
        for (int i = 1; i < inst->dimension; i++)
        {
            if (CPXchgcoef(env, lp, row, xDrone_pos(0, i, inst), 1.0))
                print_error("wrong CPXchgcoef [drone_initial_flow]");
        }
    }

    // Add "drone_ending_flow" constraint
    // *
    {
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        double rhs = 1.0;
        char sense = 'E'; // E stands for equality constraint
        sprintf(rname[0], "drone_ending_flow");
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [drone_ending_flow]");
        // drone in-flow in the final depot
        for (int i = 0; i < inst->dimension - 1; i++)
        {
            if (CPXchgcoef(env, lp, row, xDrone_pos(i, inst->dimension - 1, inst), 1.0))
                print_error("wrong CPXchgcoef [drone_ending_flow]");
        }
    }

    // Add the "drone_out_flow" constraint for each customer node
    // *
    for (int i = 1; i < inst->dimension - 1; i++)
    {
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        double rhs = 0.0;
        char sense = 'E'; // E stands for equality constraint
        sprintf(rname[0], "drone_out_flow(%d)", i);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [drone_out_flow]");
        // out-flow
        for (int j = 0; j < inst->dimension; j++)
        {
            if (CPXchgcoef(env, lp, row, xDrone_pos(i, j, inst), 1.0))
                print_error("wrong CPXchgcoef [drone_out_flow]");
        }
        // -rhs
        if (CPXchgcoef(env, lp, row, yD_pos(i, inst), -1.0))
            print_error("wrong CPXchgcoef [drone_out_flow]");
        if (CPXchgcoef(env, lp, row, yC_pos(i, inst), -1.0))
            print_error("wrong CPXchgcoef [drone_out_flow]");
    }

    // Add the drone-leg selection constraint (y)
    // *
    for (int j = 1; j < inst->dimension - 1; j++)
    {
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        double rhs = 0.0;
        char sense = 'E'; // "="
        sprintf(rname[0], "drone_leg_selection_y(%d)", j);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [drone_leg_selection_y]");
        if (CPXchgcoef(env, lp, row, yD_pos(j, inst), 1.0))
            print_error("wrong CPXchgcoef [drone_leg_selection_y]");
        for (int i = 0; i < inst->dimension - 1; i++)
        {
            for (int k = 1; k < inst->dimension; k++)
            {
                if (CPXchgcoef(env, lp, row, z_pos(i, j, k, inst), -1.0))
                    print_error("wrong CPXchgcoef [drone_leg_selection_y]");
            }
        }
    }

    // Add the drone-leg selection constraint
    // *
    for (int i = 0; i < inst->dimension - 1; i++)
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            if (i == j)
                continue;
            int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
            double rhs = 0.0;
            char sense = 'G'; // "<="
            sprintf(rname[0], "drone_leg_selection(%d,%d)", i, j);
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
                print_error("wrong CPXnewrows [drone_leg_selection]");
            if (CPXchgcoef(env, lp, row, xDrone_pos(i, j, inst), 1.0))
                print_error("wrong CPXchgcoef [drone_leg_selection]");
            for (int k = 0; k < inst->dimension; k++)
            {
                if (k == i || k == j)
                    continue;
                if (k > 0)
                {
                    if (CPXchgcoef(env, lp, row, z_pos(i, j, k, inst), -1.0))
                        print_error("wrong CPXchgcoef [drone_leg_selection]");
                }
                if (k < inst->dimension)
                {
                    if (CPXchgcoef(env, lp, row, z_pos(k, i, j, inst), -1.0))
                        print_error("wrong CPXchgcoef [drone_leg_selection]");
                }
            }
        }
    }

    // // if drone traverse edge (i,j), edge (i,j) must be traversed by the truck, OR it must be part of a drone leg
    // // *
    // for (int i = 0; i < inst->dimension - 1; i++)
    // {
    //     for (int j = 0; j < inst->dimension; j++)
    //     {
    //         if (i == j)
    //             continue;
    //         int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
    //         double rhs = 0.0;
    //         char sense = 'L'; // "<="
    //         sprintf(rname[0], "drone_leg_selection(%d,%d)", i, j);
    //         if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
    //             print_error("wrong CPXnewrows [drone_leg_selection]");
    //         if (CPXchgcoef(env, lp, row, xDrone_pos(i, j, inst), 1.0))
    //             print_error("wrong CPXchgcoef [drone_leg_selection]");
    //         for (int k = 0; k < inst->dimension; k++)
    //         {
    //             if (k == i || k == j)
    //                 continue;
    //             if (k > 0)
    //             {
    //                 if (CPXchgcoef(env, lp, row, z_pos(i, j, k, inst), -1.0))
    //                     print_error("wrong CPXchgcoef [drone_leg_selection]");
    //             }
    //             if (k < inst->dimension)
    //             {
    //                 if (CPXchgcoef(env, lp, row, z_pos(k, i, j, inst), -1.0))
    //                     print_error("wrong CPXchgcoef [drone_leg_selection]");
    //             }
    //         }
    //         if (CPXchgcoef(env, lp, row, xTruck_pos(i, j, inst), -1.0))
    //             print_error("wrong CPXchgcoef [drone_leg_selection]");
    //     }
    // }

    // ************ TIME CONSTRAINTS ************ //

    // big-M constraint:
    // set M equal to the time associated to the optimal TSP solution.
    // In this way, the only feasible solutions that violate this constraint,
    // are the ones with a makespan greater than the makespan of the optimal
    // TSP solution.
    // This is fine since a TSP solution is also a TSP-D solution.

    // M = time(TSP_best_sol) - [ min(t_D(0,i)) ] ??

    // Add time constraints

    // a(j) >= a(i) + t(i,j)
    // -> bigM
    // a(j) + M(1 - x_ij) >= a(i) + t(i,j)

    // TRUCK time constraints ORIGINAL
    // *
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
            double rhs = inst->truck_times[i][j] - M;
            char sense = 'G'; // ">="
            sprintf(rname[0], "truck_time_constraint(%d,%d)", i, j);
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
                print_error("wrong CPXnewrows [degree]");

            if (CPXchgcoef(env, lp, row, a_pos(j, inst), 1.0))
                print_error("wrong CPXchgcoef [truck_time_constraint]");
            if (CPXchgcoef(env, lp, row, a_pos(i, inst), -1.0))
                print_error("wrong CPXchgcoef [truck_time_constraint]");
            if (CPXchgcoef(env, lp, row, xTruck_pos(i, j, inst), -M))
                print_error("wrong CPXchgcoef [truck_time_constraint]");
        }
    }

    // // TRUCK time constraints EXTENDED for the drone leg i-->k-->j
    // for (int i = 0; i < inst->dimension - 1; i++)
    // {
    //     for (int j = 1; j < inst->dimension; j++)
    //     {
    //         if (i == j)
    //             continue;
    //         int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
    //         double rhs = M;
    //         char sense = 'L';
    //         sprintf(rname[0], "truck_time_constraint(%d,%d)", i, j);
    //         if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
    //             print_error("wrong CPXnewrows [degree]");

    //         if (CPXchgcoef(env, lp, row, a_pos(j, inst), -1.0))
    //             print_error("wrong CPXchgcoef [truck_time_constraint]");
    //         if (CPXchgcoef(env, lp, row, a_pos(i, inst), 1.0))
    //             print_error("wrong CPXchgcoef [truck_time_constraint]");
    //         if (CPXchgcoef(env, lp, row, xTruck_pos(i, j, inst), inst->truck_times[i][j] + M))
    //             print_error("wrong CPXchgcoef [truck_time_constraint]");
    //         for (int k = 1; k < inst->dimension - 1; k++)
    //         {
    //             if (k == i || k == j)
    //                 continue;
    //             double coeff = inst->min_time_drone[i][k][j] - inst->truck_times[i][j];
    //             if (coeff > 0)
    //             {
    //                 if (CPXchgcoef(env, lp, row, z_pos(i, k, j, inst), coeff))
    //                     print_error("wrong CPXchgcoef [truck_time_constraint]");
    //             }
    //         }
    //     }
    // }

    // remove single drone edges
    // xD(i,j) <= xT(i,j) + yD(i) + yD(j)
    for (int i = 0; i < inst->dimension - 1; i++)
    {
        for (int j = 1; j < inst->dimension; j++)
        {
            int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
            double rhs = 0.0;
            char sense = 'L'; // "<="
            sprintf(rname[0], "single_drone_edge_removal(%d,%d)", i, j);
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
                print_error("wrong CPXnewrows [single_drone_edge_removal]");

            if (CPXchgcoef(env, lp, row, xDrone_pos(i, j, inst), 1.0))
                print_error("wrong CPXchgcoef [single_drone_edge_removal]");
            if (CPXchgcoef(env, lp, row, xTruck_pos(i, j, inst), -1.0))
                print_error("wrong CPXchgcoef [single_drone_edge_removal]");
            if (CPXchgcoef(env, lp, row, yD_pos(i, inst), -1.0))
                print_error("wrong CPXchgcoef [single_drone_edge_removal]");
            if (CPXchgcoef(env, lp, row, yD_pos(j, inst), -1.0))
                print_error("wrong CPXchgcoef [single_drone_edge_removal]");
        }
    }

    // DRONE min time constraints
    //drone leg: i-->j-->k
    //*
    for (int i = 0; i < inst->dimension - 1; i++)
    {
        for (int k = 1; k < inst->dimension; k++)
        {
            if (i == k)
                continue;

            int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
            double rhs = M;
            char sense = 'L'; // "<="
            sprintf(rname[0], "drone_min_time_constraint(%d,%d)", i, k);
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
                print_error("wrong CPXnewrows [drone_min_time_constraint]");

            if (CPXchgcoef(env, lp, row, a_pos(i, inst), 1.0))
                print_error("wrong CPXchgcoef [drone_min_time_constraint]");
            if (CPXchgcoef(env, lp, row, a_pos(k, inst), -1.0))
                print_error("wrong CPXchgcoef [drone_min_time_constraint]");

            for (int j = 1; j < inst->dimension - 1; j++) // intermediate node cannot be a depot
            {
                if (j == i || j == k)
                    continue;
                if (CPXchgcoef(env, lp, row, z_pos(i, j, k, inst), inst->min_time_drone[i][j][k] + M))
                    print_error("wrong CPXchgcoef [drone_min_time_constraint]");
            }
        }
    }

    // DRONE max time constraints
    // drone leg: i-->j-->k
    // *
    for (int i = 0; i < inst->dimension - 1; i++)
    {
        for (int k = 1; k < inst->dimension; k++)
        {
            if (i == k)
                continue;

            int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
            double rhs = -M;
            char sense = 'G'; // ">="
            sprintf(rname[0], "drone_max_time_constraint(%d,%d)", i, k);
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
                print_error("wrong CPXnewrows [drone_max_time_constraint]");

            if (CPXchgcoef(env, lp, row, a_pos(i, inst), 1.0))
                print_error("wrong CPXchgcoef [drone_max_time_constraint]");
            if (CPXchgcoef(env, lp, row, a_pos(k, inst), -1.0))
                print_error("wrong CPXchgcoef [drone_max_time_constraint]");

            for (int j = 1; j < inst->dimension - 1; j++) // intermediate node cannot be a depot
            {
                if (j == i || j == k)
                    continue;
                if (CPXchgcoef(env, lp, row, z_pos(i, j, k, inst), inst->max_time_drone[i][j][k] - M))
                    print_error("wrong CPXchgcoef [drone_max_time_constraint]");
            }
        }
    }

    // eligibility of x_D(i,j) or x_D(j,i)
    // "if and only if the truck visits at least one of the
    // two customers i and j, thus ensuring that each drone
    // leg consists of a single drone customer"
    // x_D(i,j) + x_D(j,i) <= yC(i) + yC(j)
    // *
    // for (int i = 1; i < inst->dimension - 2; i++)
    // {
    //     for (int j = i + 1; j < inst->dimension - 1; j++)
    //     {
    //         int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
    //         double rhs = 0.0;
    //         char sense = 'L'; // "<="
    //         sprintf(rname[0], "drone_travel(%d,%d)", i, j);
    //         if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
    //             print_error("wrong CPXnewrows [drone_visits]");

    //         if (CPXchgcoef(env, lp, row, xDrone_pos(i, j, inst), 1.0))
    //             print_error("wrong CPXchgcoef [drone_travel]");
    //         if (CPXchgcoef(env, lp, row, xDrone_pos(j, i, inst), 1.0))
    //             print_error("wrong CPXchgcoef [drone_travel]");
    //         if (CPXchgcoef(env, lp, row, yC_pos(i, inst), -1.0))
    //             print_error("wrong CPXchgcoef [drone_travel]");
    //         if (CPXchgcoef(env, lp, row, yC_pos(j, inst), -1.0))
    //             print_error("wrong CPXchgcoef [drone_travel]");
    //     }
    // }

    // // yD_j = 1 => a_j = 0
    // for (int j = 1; j < inst->dimension - 1; j++)
    // {
    //     int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
    //     double rhs = M;
    //     char sense = 'L'; // "<="
    //     sprintf(rname[0], "yD_a_j(%d)", j);
    //     if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
    //         print_error("wrong CPXnewrows [drone_min_time_constraint]");

    //     if (CPXchgcoef(env, lp, row, a_pos(j, inst), 1.0))
    //         print_error("wrong CPXchgcoef [drone_min_time_constraint]");
    //     if (CPXchgcoef(env, lp, row, yD_pos(j, inst), M))
    //         print_error("wrong CPXchgcoef [drone_min_time_constraint]");
    // }

    // **** makespan lower bounds ****

    // a(0') >= truck travel time
    {
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        double rhs = 0.0;
        char sense = 'G'; // ">="
        sprintf(rname[0], "makespan_lb_truck");
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [makespan_lb_truck]");
        if (CPXchgcoef(env, lp, row, a_pos(inst->dimension - 1, inst), 1.0))
            print_error("wrong CPXchgcoef [makespan_lb_truck]");
        for (int i = 0; i < inst->dimension; i++)
        {
            for (int j = 0; j < inst->dimension; j++)
            {
                if (CPXchgcoef(env, lp, row, xTruck_pos(i, j, inst), -inst->truck_times[i][j]))
                    print_error("wrong CPXchgcoef [makespan_lb_truck]");
            }
        }
    }

    // a(0') >= drone travel time ORIGINAL
    {
        int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
        double rhs = 0.0;
        char sense = 'G'; // ">="
        sprintf(rname[0], "makespan_lb_drone");
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [makespan_lb_drone]");
        if (CPXchgcoef(env, lp, row, a_pos(inst->dimension - 1, inst), 1.0))
            print_error("wrong CPXchgcoef [makespan_lb_drone]");
        for (int i = 0; i < inst->dimension; i++)
        {
            for (int j = 0; j < inst->dimension; j++)
            {
                for (int k = 0; k < inst->dimension; k++)
                {
                    if (CPXchgcoef(env, lp, row, z_pos(i, j, k, inst), -inst->min_time_drone[i][j][k]))
                        print_error("wrong CPXchgcoef [makespan_lb_drone]");
                }
            }
        }
    }

    // // a(0') <= drone travel time EXTENDED
    // {
    //     int row = CPXgetnumrows(env, lp); // get the number of rows inside the model
    //     double rhs = 0.0;
    //     char sense = 'L';
    //     sprintf(rname[0], "makespan_lb_drone");
    //     if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
    //         print_error("wrong CPXnewrows [makespan_lb_drone]");
    //     if (CPXchgcoef(env, lp, row, a_pos(inst->dimension - 1, inst), -1.0))
    //         print_error("wrong CPXchgcoef [makespan_lb_drone]");
    //     for (int i = 0; i < inst->dimension; i++)
    //     {
    //         for (int j = 0; j < inst->dimension; j++)
    //         {
    //             for (int k = 0; k < inst->dimension; k++)
    //             {
    //                 if (CPXchgcoef(env, lp, row, z_pos(i, j, k, inst), inst->min_time_drone[i][j][k]))
    //                     print_error("wrong CPXchgcoef [makespan_lb_drone]");
    //                 if (CPXchgcoef(env, lp, row, z_pos(i, j, k, inst), - inst->truck_times[i][j] - inst->truck_times[j][k]))
    //                     print_error("wrong CPXchgcoef [makespan_lb_drone]");
    //             }
    //             if (CPXchgcoef(env, lp, row, xDrone_pos(i, j, inst), inst->truck_times[i][j]))
    //                 print_error("wrong CPXchgcoef [makespan_lb_drone]");
    //         }
    //     }
    // }

    // ********** OPTIONAL CONSTRAINTS **********

    // sum { yD(i) } <= ceil(n_customers / 2)
    {
        int row = CPXgetnumrows(env, lp);
        double rhs = ceil((inst->dimension - 2) / 2.0);
        char sense = 'L'; // "<="
        sprintf(rname[0], "sum_yD");
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
            print_error("wrong CPXnewrows [sum_yD]");

        for (int i = 1; i < inst->dimension - 1; i++)
        {
            if (CPXchgcoef(env, lp, row, yD_pos(i, inst), 1.0))
                print_error("wrong CPXchgcoef [sum_yD]");
        }
    }

    // {
    //     int row = CPXgetnumrows(env, lp);
    //     double rhs = inst->dimension - 1;
    //     char sense = 'E';
    //     sprintf(rname[0], "u_returning_depot");
    //     if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
    //         print_error("wrong CPXnewrows [u_returning_depot]");

    //     if (CPXchgcoef(env, lp, row, u_pos(inst->dimension - 1, inst), 1.0))
    //         print_error("wrong CPXchgcoef [u_returning_depot]");
    //     for (int i = 1; i < inst->dimension - 1; i++)
    //     {
    //         if (CPXchgcoef(env, lp, row, yT_pos(i, inst), 1.0))
    //             print_error("wrong CPXchgcoef [u_returning_depot]");
    //     }
    // }

    // u_j >= u_i +1 if (x_ij = 1)
    // for (int i = 1; i < inst->dimension - 1; i++)
    // {
    //     for (int j = 1; j < inst->dimension; j++)
    //     {
    //         int row = CPXgetnumrows(env, lp);
    //         double rhs = (inst->dimension - 1.0) - 1.0;
    //         char sense = 'L'; // "<="
    //         sprintf(rname[0], "u_drone(%d,%d)", i, j);
    //         if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
    //             print_error("wrong CPXnewrows [u]");
    //         if (CPXchgcoef(env, lp, row, u_pos(i, inst), 1.0))
    //             print_error("wrong CPXchgcoef [u]");
    //         if (CPXchgcoef(env, lp, row, u_pos(j, inst), -1.0))
    //             print_error("wrong CPXchgcoef [u]");
    //         if (CPXchgcoef(env, lp, row, xDrone_pos(i, j, inst), inst->dimension - 1.0))
    //             print_error("wrong CPXchgcoef [u]");

    //         row = CPXgetnumrows(env, lp);
    //         rhs = (inst->dimension - 1.0) - 1.0;
    //         sense = 'L'; // "<="
    //         sprintf(rname[0], "u_truck(%d,%d)", i, j);
    //         if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
    //             print_error("wrong CPXnewrows [u]");
    //         if (CPXchgcoef(env, lp, row, u_pos(i, inst), 1.0))
    //             print_error("wrong CPXchgcoef [u]");
    //         if (CPXchgcoef(env, lp, row, u_pos(j, inst), -1.0))
    //             print_error("wrong CPXchgcoef [u]");
    //         if (CPXchgcoef(env, lp, row, xTruck_pos(i, j, inst), inst->dimension - 1.0))
    //             print_error("wrong CPXchgcoef [u]");
    //     }
    // }

    // //    5 -->  3 --> 11       710.248510
    // {
    //     int row = CPXgetnumrows(env, lp);
    //     double rhs = -2*710.248510;
    //     char sense = 'G';
    //     sprintf(rname[0], "drone_leg_5_3_11");
    //     if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
    //         print_error("wrong CPXnewrows [drone_leg_5_3_11]");
    //     if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, rname))
    //         print_error("wrong CPXnewrows [drone_leg_5_3_11]");
    //     if (CPXchgcoef(env, lp, row, a_pos(11, inst), 1.0))
    //         print_error("wrong CPXchgcoef [drone_leg_5_3_11]");
    //     if (CPXchgcoef(env, lp, row, a_pos(5, inst), -1.0))
    //         print_error("wrong CPXchgcoef [drone_leg_5_3_11]");

    //     if (CPXchgcoef(env, lp, row, xDrone_pos(5,3, inst), -710.248510))
    //         print_error("wrong CPXchgcoef [drone_leg_5_3_11]");
    //     if (CPXchgcoef(env, lp, row, xDrone_pos(3,11, inst), -710.248510))
    //         print_error("wrong CPXchgcoef [drone_leg_5_3_11]");
    //     if (CPXchgcoef(env, lp, row, xTruck_pos(5,11, inst), -710.248510))
    //         print_error("wrong CPXchgcoef [drone_leg_5_3_11]");
    // }

    free(cname[0]);
    free(cname);
    free(rname[0]);
    free(rname);
}

// static int CPXPUBLIC callback_driver(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle)
// {
//     if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
//         return callback_candidate(context, contextid, userhandle);
//     if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
//         return callback_relaxation(context, contextid, userhandle);
//     print_error("contextid unknownn in my_callback");
//     return 1;
// }

// static int CPXPUBLIC callback_candidate(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle)
// {
//     instance *inst = (instance *)userhandle;
//     double *xstar = (double *)malloc(inst->cols * sizeof(double));
//     double objval = CPX_INFBOUND;
//     if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->cols - 1, &objval))
//         print_error("CPXcallbackgetcandidatepoint error");

//     int *comp = (int *)calloc(inst->dimension, sizeof(int));
//     int *succ = (int *)calloc(inst->dimension, sizeof(int));
//     int ncomp = 0; // number of connected components
//     int *length_comp;

//     // Retrieve the connected components of the current solution
//     findConnectedComponents(xstar, inst, succ, comp, &ncomp, &length_comp);

//     if (ncomp > 1)
//     {
//         // add one cut for each connected component
//         for (int mycomp = 0; mycomp < ncomp; mycomp++)
//         {
//             int nnz = 0;
//             int izero = 0;
//             char sense = 'L';
//             double rhs = length_comp[mycomp] - 1.0; // |S|-1
//             int *index = (int *)calloc(inst->cols, sizeof(int));
//             double *value = (double *)calloc(inst->cols, sizeof(double));

//             for (int i = 0; i < inst->dimension; i++)
//             {
//                 if (comp[i] != mycomp)
//                     continue;
//                 for (int j = i + 1; j < inst->dimension; j++)
//                 {
//                     if (comp[j] != mycomp)
//                         continue;
//                     index[nnz] = xpos(i, j, inst);
//                     value[nnz++] = 1.0;
//                 }
//             }
//             // reject the solution and adds one cut
//             if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value))
//                 print_error("CPXcallbackrejectcandidate() error");
//             free(index);
//             free(value);
//         }
//     }
//     else if (inst->param.opt)
//     {
//         // the candidate solution has not connected components but could have crossings... let's apply 2-opt
//         double delta = two_opt(inst, succ, 0);
//         if (delta < 0)
//         {
//             objval += delta;
//             // succ -> xstar
//             int nnz = 0;
//             int izero = 0;
//             int *index = (int *)calloc(inst->cols, sizeof(int));
//             double *xstar_succ = (double *)calloc(inst->cols, sizeof(double));

//             for (int i = 0; i < inst->dimension; i++)
//             {
//                 for (int j = i + 1; j < inst->dimension; j++)
//                 {
//                     index[nnz] = xpos(i, j, inst);
//                     if (j == succ[i] || i == succ[j])
//                     {
//                         xstar_succ[nnz++] = 1.0;
//                     }
//                     else
//                         xstar_succ[nnz++] = 0.0;
//                 }
//             }
//             // sanity check
//             if (nnz != inst->cols)
//                 print_error("Error in applying 2-opt in callback_candidate");

//             if (CPXcallbackpostheursoln(context, nnz, index, xstar_succ, objval, CPXCALLBACKSOLUTION_CHECKFEAS))
//                 print_error("CPXcallbackpostheursoln() error");
//             if (inst->param.verbose >= DEBUG)
//                 printf("[callback_candidate] Posted a new solution to CPLEX with incumbent: %f\n", objval);
//             free(xstar_succ);
//             free(index);
//         }
//     }
//     free(comp);
//     free(succ);
//     free(length_comp);
//     free(xstar);
//     return 0;
// }

// static int CPXPUBLIC callback_relaxation(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle)
// {
//     //printf("*** Inside callback_relaxation\n");
//     int node_depth;
//     CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODEDEPTH, &node_depth);
//     if (!node_depth)
//         return 0;
//     instance *inst = (instance *)userhandle;
//     double *xstar = (double *)malloc(inst->cols * sizeof(double));
//     double objval = CPX_INFBOUND;
//     double const eps = 0.1;

//     if (CPXcallbackgetrelaxationpoint(context, xstar, 0, inst->cols - 1, &objval))
//         print_error("CPXcallbackgetrelaxationpoint error");

//     int ncomp;
//     int *comp = (int *)calloc(inst->dimension, sizeof(int));
//     int *length_comp = (int *)calloc(inst->dimension, sizeof(int));
//     // list of edges in "node format" [ i_1, j_1 , i_2, j_2, i_3, j_3, ...]
//     int *elist = malloc(2 * inst->cols * sizeof(int));

//     int loader = 0;
//     int ecount = 0; // edge count
//     double *x = malloc(inst->cols * sizeof(double));
//     for (int i = 0; i < inst->dimension; i++)
//     {
//         for (int j = i + 1; j < inst->dimension; j++)
//         {
//             if (xstar[xpos(i, j, inst)] > 0.001) // just save the selected (also partially) edges
//             {
//                 elist[loader++] = i;
//                 elist[loader++] = j;
//                 x[ecount++] = xstar[xpos(i, j, inst)];
//             }
//         }
//     }

//     // realloc
//     x = realloc(x, ecount * sizeof(double));
//     elist = realloc(elist, 2 * ecount * sizeof(int));
//     //printf("*** ecount = %d\n", ecount);

//     if (CCcut_connect_components(inst->dimension, ecount, elist, x, &ncomp, &length_comp, &comp))
//         print_error("CCcut_connect_components error");

//     //printf("*** ncomp = %d\n", ncomp);
//     if (ncomp == 1) // CCcut_violated_cuts() assumes graph is connected
//     {
//         doit_fn_input in;
//         in.ecount = ecount;
//         in.elist = elist;
//         in.x = x;
//         in.context = context;
//         in.inst = inst;
//         if (CCcut_violated_cuts(inst->dimension, ecount, elist, x, 2 - eps, doit_fn_concorde, (void *)&in))
//             print_error("CCcut_violated_cuts error");
//     }
//     else if (ncomp > 1)
//     {
//         // add one cut for each connected component
//         for (int mycomp = 0; mycomp < ncomp; mycomp++)
//         {
//             int nnz = 0;
//             int izero = 0;
//             char sense = 'L';
//             double rhs = length_comp[mycomp] - 1.0; // |S|-1
//             int *index = (int *)calloc(inst->cols, sizeof(int));
//             double *value = (double *)calloc(inst->cols, sizeof(double));

//             for (int i = 0; i < inst->dimension; i++)
//             {
//                 if (comp[i] != mycomp)
//                     continue;
//                 for (int j = i + 1; j < inst->dimension; j++)
//                 {
//                     if (comp[j] != mycomp)
//                         continue;
//                     index[nnz] = xpos(i, j, inst);
//                     value[nnz++] = 1.0;
//                 }
//             }
//             // reject the solution and adds one cut
//             int purgeable = CPX_USECUT_FILTER; // Let CPLEX decide whether to keep the cut or not
//             int local = 0;                     // Global cut
//             if (CPXcallbackaddusercuts(context, 1, nnz, &rhs, &sense, &izero, index, value, &purgeable, &local))
//                 print_error("CPXcallbackaddusercuts() error in callback_relaxation()");
//             free(index);
//             free(value);
//         }
//     }
//     else
//     {
//         print_error("Invalip ncomp in callback_relaxation()");
//     }

//     free(comp);
//     free(length_comp);
//     free(xstar);
//     return 0;
// }

// // double cutval = value of the cut
// // int cutcount = number of nodes in the cut
// // int ∗cut = the array of the members of the cut (indices of the nodes in the cut) =>
// // => the nodes of one of the two connected components (look cut_st.c of the concorde library)
// int doit_fn_concorde(double cutval, int cutcount, int *cut, void *in)
// {

//     doit_fn_input *input = (doit_fn_input *)in;

//     // DEBUG
//     // printf("*** Inside doit_fn_concorde\n");
//     // printf("#### cutval: %f ", cutval);
//     // printf("#### cutcount: %d \n", cutcount);
//     // for (int i = 0; i < cutcount; i++)
//     //     printf("%d,", cut[i]);
//     // printf("\n");
//     // for (int i = 0; i < input->ecount; i++)
//     // {
//     //     printf("x[%d,%d] = %f\n", input->elist[2 * i], input->elist[2 * i + 1], input->x[i]);
//     // }

//     // SEC for the current connected component
//     double *value = (double *)calloc(cutcount * (cutcount - 1) / 2, sizeof(double));
//     int *index = (int *)calloc(cutcount * (cutcount - 1) / 2, sizeof(int));

//     double rhs = cutcount - 1.0;
//     int nnz = 0;
//     char sense = 'L';
//     int purgeable = CPX_USECUT_FILTER; // Let CPLEX decide whether to keep the cut or not
//     int local = 0;                     // Global cut
//     int izero = 0;

//     for (int i = 0; i < cutcount; i++)
//     {
//         for (int j = 0; j < cutcount; j++)
//         {
//             if (i == j)
//                 continue;
//             if (cut[i] < cut[j])
//             {
//                 index[nnz] = xpos(cut[i], cut[j], input->inst);
//                 value[nnz++] = 1.0;
//             }
//         }
//     }
//     // SANITY CHECK
//     if (nnz != cutcount * (cutcount - 1) / 2)
//         print_error("nnz != cutcount*(cutcount-1)/2");

//     if (CPXcallbackaddusercuts(input->context, 1, nnz, &rhs, &sense, &izero, index, value, &purgeable, &local))
//         print_error("CPXcallbackaddusercuts() error"); // add user cut

//     free(index);
//     free(value);
//     return 0;
// }

// void findConnectedComponents(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp,
//                              int **length_comp)
// {
//     // build succ[] and comp[] wrt xstar()...

//     // initialization
//     *ncomp = 0;
//     for (int i = 0; i < inst->dimension; i++)
//     {
//         succ[i] = -1;
//         comp[i] = -1;
//     }

//     for (int start = 0; start < inst->dimension; start++)
//     {
//         if (comp[start] >= 0)
//             continue; // node "start" was already visited, just skip it

//         // a new component is found
//         (*ncomp)++;
//         if ((*ncomp) == 1)
//             (*length_comp) = (int *)calloc(1, sizeof(int));
//         else
//             (*length_comp) = (int *)realloc(*length_comp, (*ncomp) * sizeof(int));
//         int i = start;
//         int length = 1;
//         int done = 0;
//         while (!done) // go and visit the current component
//         {
//             comp[i] = (*ncomp) - 1;
//             done = 1;
//             for (int j = 0; j < inst->dimension; j++)
//             {
//                 if (i == j)
//                     continue;
//                 if (xstar[xpos(i, j, inst)] > 0.5 &&
//                     comp[j] == -1) // the edge [i,j] is selected in xstar and j was not visited before
//                 {
//                     succ[i] = j;
//                     length++;
//                     i = j;
//                     done = 0;
//                     break;
//                 }
//             }
//         }
//         succ[i] = start;                       // last arc to close the cycle
//         (*length_comp)[(*ncomp) - 1] = length; // save length of the cycle

//         // go to the next component...
//     }
// }

// retrieve the max drone leg time for any landing node given the takeoff node and the customer served by the drone
void retrieve_max_i_drone_leg_times(instance *inst)
{
    // initialize the array max_i_drone_leg_times
    inst->max_i_drone_leg_times = (double *)calloc(inst->dimension, sizeof(double));
    // retrieve the max for each (i,j)
    for (int i = 0; i < inst->dimension; i++)
    {
        inst->max_i_drone_leg_times[i] = 0.0;
        for (int j = 1; j < inst->dimension - 1; j++)
        {
            if (inst->nodes[j].weight > WEIGHT_LIMIT)
                continue;
            for (int k = 1; k < inst->dimension; k++)
            {
                if (inst->max_time_drone[i][j][k] > inst->max_i_drone_leg_times[i])
                    inst->max_i_drone_leg_times[i] = inst->max_time_drone[i][j][k];
            }
        }
    }
}

// retrieve the max drone leg time for any landing node given the takeoff node and the customer served by the drone
void retrieve_max_ij_drone_leg_times(instance *inst)
{
    // initialize the array max_i_drone_leg_times
    inst->max_ij_drone_leg_times = (double **)calloc(inst->dimension, sizeof(double *));

    for (int i = 0; i < inst->dimension; i++)
    {
        inst->max_ij_drone_leg_times[i] = (double *)calloc(inst->dimension, sizeof(double));
    }

    // retrieve the max for each (i,j)
    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = 1; j < inst->dimension - 1; j++)
        {
            inst->max_ij_drone_leg_times[i][j] = 0.0;
            if (inst->nodes[j].weight > WEIGHT_LIMIT)
                continue;
            for (int k = 1; k < inst->dimension; k++)
            {
                if (inst->max_time_drone[i][j][k] > inst->max_ij_drone_leg_times[i][j])
                    inst->max_ij_drone_leg_times[i][j] = inst->max_time_drone[i][j][k];
            }
        }
    }
}
