#include "../include/utils.h"

// seq: rules of precedence
// if two nodes i,j are visited by the same vehicle
// then i must appear before j in the sequence seq
// double compute_sequence_min_time(int *seq, instance *inst)
// {
//     // array of minimum arrival times T
//     // T[i] = minimum arrival time at node i is visited by the truck
//     // (since the last node must be visited by the truck, this assumption is ok)
//     double T[inst->dimension];

//     // initialize the arrival times T[i] equal to the arrival times related to the case in which the customers are served only by the truck
//     T[0] = 0.0;
//     for (int i = 0; i < inst->dimension - 1; i++)
//     {
//         T[i + 1] = T[i] + inst->truck_times[seq[i]][seq[i + 1]];
//     }

//     T[0] = 0.0;
//     for (int i = 0; i < inst->dimension - 1; i++)
//     {
//         // move 1: move the truck with the drone onboard by one position
//         if (T[i] + inst->truck_times[seq[i]][seq[i + 1]] < T[i + 1])
//         {
//             T[i + 1] = T[i] + inst->truck_times[seq[i]][seq[i + 1]];
//         }

//         // move 2: launch the drone from the current node i, to visit node j, and rejoin the truck with the drone at node k
//         //printf("\ni = %d (%d)\n", i, seq[i]);
//         double truck_t = 0.0;
//         for (int k = i + 2; k < inst->dimension; k++)
//         {
//             if (k == i + 2)
//                 truck_t += inst->truck_times[seq[i]][seq[k]];
//             else
//                 truck_t += inst->truck_times[seq[k - 1]][seq[k]];

//             if (truck_t > inst->max_i_drone_leg_times[seq[i]])
//                 break; // does not exist a drone leg that start from node i that can last up to the time truck_t
//             double truck_tj = truck_t;
//             //printf("\n\tk = %d(%d): ", k, seq[k]);
//             for (int j = i + 1; j < k; j++)
//             {
//                 if (inst->nodes[seq[j]].weight > WEIGHT_LIMIT) // node j cannot be a drone customer
//                     continue;
//                 // adapt the travel time of the truck path, since now the customer served by the drone is j
//                 if (j > i + 1)
//                 {
//                     truck_tj = truck_tj - inst->truck_times[seq[j - 2]][seq[j]] - inst->truck_times[seq[j]][seq[j + 1]];
//                     truck_tj = truck_tj + inst->truck_times[seq[j - 2]][seq[j - 1]] + inst->truck_times[seq[j - 1]][seq[j + 1]];
//                 }
//                 // check if the drone leg i-->j-->k is feasible
//                 if (truck_tj > inst->max_time_drone[seq[i]][seq[j]][seq[k]]) // the drone leg i-->j-->k is infeasible
//                     continue;
//                 // // drone leg is feasible
//                 //printf("%d(%d) ", j, seq[j]);
//                 double tmp = inst->min_time_drone[seq[i]][seq[j]][seq[k]];
//                 if (tmp < truck_tj)
//                     tmp = truck_tj;
//                 if (T[i] + tmp < T[k])
//                     T[k] = T[i] + tmp;
//             }
//         }
//     }
//     // for (int i = 0; i < inst->dimension; i++)
//     //     printf("T[%d] = %f\n", seq[i], T[i]);
//     return T[inst->dimension - 1];
// }

double compute_sequence_min_time(int *seq, instance *inst, int *truck_seq, int *drone_seq)
{
    // array of minimum arrival times T
    // T[i] = minimum arrival time at node i, when i is visited by the truck
    // (since the last node must be visited by the truck, this assumption is ok)
    double T[inst->dimension];

    // move1[k] = 1 => optimal move that ends in k is of type 1 (truck moves forward with the drone onboard)
    // move1[k] = 0 => optimal move that ends in k is of type 2 (drone leg which ends in k)
    int move1[inst->dimension];
    // move2_drone_cust[k] = j if j is the drone customer visited by the drone leg associated to the landing node k
    // move2_takeoff_node[k] = i if i is the takeoff node from which it starts the drone leg associated to the landing node k
    int move2_drone_cust[inst->dimension];
    int move2_takeoff_node[inst->dimension];

    // initialize arrival times T[i] equal to the arrival times related to the truck only case
    T[0] = 0.0;
    for (int i = 0; i < inst->dimension - 1; i++)
    {
        T[i + 1] = T[i] + inst->truck_times[seq[i]][seq[i + 1]];
        move1[i + 1] = 1;
    }

    for (int i = 0; i < inst->dimension - 1; i++)
    {
        // move 1: move the truck with the drone onboard by one position
        if (T[i] + inst->truck_times[seq[i]][seq[i + 1]] < T[i + 1])
        {
            T[i + 1] = T[i] + inst->truck_times[seq[i]][seq[i + 1]];
            move1[i + 1] = 1;
        }

        // move 2: launch the drone from the current node i, to visit node j, and rejoin the truck with the drone at node k
        //** printf("\ni = %d (%d)\n", i, seq[i]);

        for (int j = i + 1; j < inst->dimension - 1; j++)
        {
            if (inst->nodes[seq[j]].weight > WEIGHT_LIMIT) // node j cannot be a drone customer
                continue;
            //** printf("\n\tj = %d(%d) : ", j, seq[j]);
            double truck_t = 0.0;
            for (int k = j + 1; k < inst->dimension; k++)
            {
                // update truck travel time
                if (k == j + 1)
                {
                    for (int w = i + 1; w < j; w++)
                    {
                        truck_t += inst->truck_times[seq[w - 1]][seq[w]];
                    }
                    truck_t += inst->truck_times[seq[j - 1]][seq[k]];
                }
                else
                    truck_t += inst->truck_times[seq[k - 1]][seq[k]];

                // if the truck_t > max_k {drone travel time i-->j-->k}   STOP
                if (truck_t > inst->max_ij_drone_leg_times[seq[i]][seq[j]])
                {
                    // does not exist a landing node k associated to a feasible drone leg i-->j-->k
                    break;
                }
                // check if the drone leg i-->j-->k is feasible
                if (truck_t > inst->max_time_drone[seq[i]][seq[j]][seq[k]]) // the drone leg i-->j-->k is infeasible
                {
                    continue;
                }

                // drone leg is feasible
                //** printf("%d(%d) ", k, seq[k]);

                double tmp = inst->min_time_drone[seq[i]][seq[j]][seq[k]];

                if (tmp < truck_t)
                    tmp = truck_t;

                if (T[i] + tmp < T[k])
                {
                    T[k] = T[i] + tmp;
                    move1[k] = 0;
                    move2_drone_cust[k] = j;
                    move2_takeoff_node[k] = i;
                }
            }
        }
    }
    // for (int i = 0; i < inst->dimension; i++)
    //     printf("T[%d] = %f\n", seq[i], T[i]);
    // store truck and drone sequences

    int truck_idx = inst->dimension - 1;
    int drone_idx = inst->dimension - 1;
    truck_seq[truck_idx--] = inst->dimension - 1;
    drone_seq[drone_idx--] = inst->dimension - 1;

    for (int k = inst->dimension - 1; k > 0; k--)
    {
        if (move1[k] == 1)
        {
            truck_seq[truck_idx--] = seq[k - 1];
            drone_seq[drone_idx--] = seq[k - 1];
        }
        else
        {
            drone_seq[drone_idx--] = seq[move2_drone_cust[k]];
            drone_seq[drone_idx--] = seq[move2_takeoff_node[k]];
            int i = k - 1;
            while (i > move2_takeoff_node[k])
            {
                if (i == move2_drone_cust[k])
                {
                    i--;
                    continue;
                }
                truck_seq[truck_idx--] = seq[i];
                i--;
            }
            truck_seq[truck_idx--] = seq[move2_takeoff_node[k]];
            // go directly to k = move2_takeoff_node[k]
            k = move2_takeoff_node[k] + 1;
        }
    }

    // shift the sequences such that 0 is at the beginning
    memmove(truck_seq, truck_seq + truck_idx + 1, (inst->dimension - truck_idx) * sizeof(int));
    memmove(drone_seq, drone_seq + drone_idx + 1, (inst->dimension - drone_idx) * sizeof(int));
    return T[inst->dimension - 1];
}

void generate_random_sequence(int *rand_seq, instance *inst)
{
    // fix starting and ending node
    rand_seq[0] = 0;
    rand_seq[inst->dimension - 1] = inst->dimension - 1;

    int numbers[inst->dimension - 2];
    for (int i = 0; i < inst->dimension - 2; i++)
    {
        numbers[i] = i + 1;
    }

    for (int i = 0; i < inst->dimension - 2; i++)
    {
        int idx = rand() % (inst->dimension - 2 - i);
        rand_seq[i + 1] = numbers[idx];
        numbers[idx] = numbers[inst->dimension - 3 - i];
    }
}

// revert the operation with the greatest time gain (if it exists)
// let us call C the current operation that we are evaluating
// case 1: the operation C is not preceded or succeeded by other operations
// case 2: the operation C is not preceded by another operation, but it is succeded by another operation
// case 3: the operation C is preceded by another operation, but it is not succeded by another operation
// case 4: the operation C is both preceded and succeded by another operation
double opt_drone_legs_direction(int *truck_succ, int *drone_succ, instance *inst)
{
    int last_takeoff_node = -1;
    // a drone leg that starts from node 0 must not be inverted
    int i = drone_succ[0];
    int prev_truck_node = 0;
    // check if from node 0 it starts a drone leg
    if (drone_succ[0] != truck_succ[0])
    {
        last_takeoff_node = 0;
        i = drone_succ[drone_succ[0]];
        while (truck_succ[prev_truck_node] != i)
        {
            prev_truck_node = truck_succ[prev_truck_node];
        }
    }

    double delta_max = 0;
    // a drone leg that ends in the final node must not be inverted
    while (drone_succ[i] != inst->dimension - 1 && drone_succ[drone_succ[i]] != inst->dimension - 1)
    {
        if (drone_succ[i] != truck_succ[i]) // drone leg found
        {
            // printf("DRONE LEG FOUND: %d --> %d --> %d\n", i, drone_succ[i], drone_succ[drone_succ[i]]);
            // printf("{prev_truck_node = %d}\n", prev_truck_node);
            // drone leg: i --> j --> k
            // we want to invert it: k --> j --> i
            int j = drone_succ[i];
            int k = drone_succ[j];

            int new_prev_truck_node; // next value to assign to the variable prev_truck_node at the end of the if
            // compute the truck travel time to traverse the path k -> i and the path i -> k
            int tmp = i;
            double truck_t_ki = 0.0;
            double truck_t_ik = 0.0;
            while (tmp != k)
            {
                truck_t_ik += inst->truck_times[tmp][truck_succ[tmp]];
                truck_t_ki += inst->truck_times[truck_succ[tmp]][tmp];
                new_prev_truck_node = tmp;
                tmp = truck_succ[tmp];
            }
            // check if the drone leg k --> j --> i is feasible
            //printf("\ttruck_t_ki = %f\n", truck_t_ki);
            //printf("\tmax_time_drone[k][j][i] = %f\n", inst->max_time_drone[k][j][i]);
            if (truck_t_ki > inst->max_time_drone[k][j][i])
            {
                //printf("*** drone leg k --> j --> i is infeasible ***\n");
                prev_truck_node = new_prev_truck_node;
                last_takeoff_node = i;
                i = k;
                continue;
            }
            // compute the minimum feasible travel times related to the drone legs i --> j --> k and k --> j --> i
            double t_ijk = inst->min_time_drone[i][j][k];
            double t_kji = inst->min_time_drone[k][j][i];
            if (t_ijk < truck_t_ik)
                t_ijk = truck_t_ik;
            if (t_kji < truck_t_ki)
                t_kji = truck_t_ki;
            // printf("\tt_ijk = %f\n", t_ijk);
            // printf("\tt_kji = %f\n", t_kji);
            // printf("{last_takeoff_node = %d}\n", last_takeoff_node);
            if (last_takeoff_node == -1 || drone_succ[drone_succ[last_takeoff_node]] != i) // case 1|2
            {
                if (drone_succ[k] == truck_succ[k]) // case 1
                {
                    //printf("\tCASE 1\n");
                    double delta = -inst->truck_times[prev_truck_node][i] - t_ijk - inst->truck_times[k][truck_succ[k]] + inst->truck_times[prev_truck_node][k] + t_kji + inst->truck_times[i][truck_succ[k]];
                    //printf("\tdelta = %f\n", delta);
                    if (delta < 0)
                    {
                        // it is convenient to revert the drone leg
                        // printf("*** it is convenient to reverse the drone leg %d --> %d --> %d\n", i, j, k);
                        // printf("\tdelta = %f\n", delta);
                        if (delta < delta_max)
                        {
                            delta_max = delta;
                        }
                    }
                }
                else // case 2
                {
                    //printf("\tCASE 2\n");
                    // from k it begins the drone leg: k --> u --> v
                    int u = drone_succ[k];
                    int v = drone_succ[u];
                    if (inst->min_time_drone[i][u][v] < 1e-6)
                    {
                        //printf("*** drone leg i --> u --> v is infeasible ***\n");
                        prev_truck_node = new_prev_truck_node;
                        last_takeoff_node = i;
                        i = k;
                        continue;
                    }
                    double truck_t_kv = 0.0;
                    int tmp = k;
                    while (tmp != v)
                    {
                        truck_t_kv += inst->truck_times[tmp][truck_succ[tmp]];

                        tmp = truck_succ[tmp];
                    }
                    double truck_t_iv = truck_t_kv - inst->truck_times[k][truck_succ[k]] + inst->truck_times[i][truck_succ[k]];
                    // check if the drone leg can be inverted
                    //printf("\ttruck_t_iv = %f\n", truck_t_iv);
                    //printf("\tmax_time_drone[i][u][v] = %f\n", inst->max_time_drone[i][u][v]);
                    if (truck_t_iv > inst->max_time_drone[i][u][v])
                    {
                        //printf("*** drone leg i --> u --> v is infeasible ***\n");
                        prev_truck_node = new_prev_truck_node;
                        last_takeoff_node = i;
                        i = k;
                        continue;
                    }
                    // compute the minimum feasible travel times related to the drone legs k --> u --> v and i --> u --> v
                    double t_kuv = inst->min_time_drone[k][u][v];
                    double t_iuv = inst->min_time_drone[i][u][v];
                    if (t_kuv < truck_t_kv)
                        t_kuv = truck_t_kv;
                    if (t_iuv < truck_t_iv)
                        t_iuv = truck_t_iv;
                    //printf("\tt_kuv = %f\n", t_kuv);
                    //printf("\tt_iuv = %f\n", t_iuv);
                    double delta = -inst->truck_times[prev_truck_node][i] - t_ijk - t_kuv + inst->truck_times[prev_truck_node][k] + t_kji + t_iuv;
                    //printf("\tdelta = %f\n", delta);
                    if (delta < 0)
                    {
                        // it is convenient to revert the drone leg
                        //printf("*** it is convenient to reverse the drone leg %d --> %d --> %d\n", i, j, k);
                        //printf("\tdelta = %f\n", delta);
                        if (delta < delta_max)
                        {
                            delta_max = delta;
                        }
                    }
                }
            }
            else // case 3|4
            {
                // a drone lands at node i
                // lets call this drone leg: m --> n --> i
                int m = last_takeoff_node;
                int n = drone_succ[m];
                if (inst->min_time_drone[m][n][k] < 1e-6)
                {
                    //printf("*** drone leg m --> n --> k is infeasible ***\n");
                    prev_truck_node = new_prev_truck_node;
                    last_takeoff_node = i;
                    i = k;
                    continue;
                }

                double truck_t_mi = 0.0;
                int tmp = m;
                while (tmp != i)
                {
                    truck_t_mi += inst->truck_times[tmp][truck_succ[tmp]];
                    tmp = truck_succ[tmp];
                }
                double truck_t_mk = truck_t_mi - inst->truck_times[prev_truck_node][i] + inst->truck_times[prev_truck_node][k];
                // check if the drone leg can be inverted
                // printf("\ttruck_t_mi = %f\n", truck_t_mi);
                // printf("\ttruck_t_mk = %f\n", truck_t_mk);
                // printf("\tmax_time_drone[m][n][i] = %f\n", inst->max_time_drone[m][n][i]);
                if (truck_t_mk > inst->max_time_drone[m][n][k])
                {
                    //printf("*** drone leg m --> n --> k is infeasible ***\n");
                    prev_truck_node = new_prev_truck_node;
                    last_takeoff_node = i;
                    i = k;
                    continue;
                }
                // compute the minimum feasible travel times related to the drone legs m --> n --> i and m --> n --> k
                double t_mni = inst->min_time_drone[m][n][i];
                double t_mnk = inst->min_time_drone[m][n][k];
                if (t_mni < truck_t_mi)
                    t_mni = truck_t_mi;
                if (t_mnk < truck_t_mk)
                    t_mnk = truck_t_mk;
                // printf("\tt_mni = %f\n", t_mni);
                // printf("\tt_mnk = %f\n", t_mnk);
                if (drone_succ[k] == truck_succ[k]) // case 3
                {
                    //printf("\tCASE 3\n");
                    double delta = -t_mni - t_ijk - inst->truck_times[k][drone_succ[k]] + t_mnk + t_kji + inst->truck_times[i][drone_succ[k]];
                    //printf("\tdelta = %f\n", delta);
                    if (delta < 0)
                    {
                        // it is convenient to revert the drone leg
                        // printf("*** it is convenient to reverse the drone leg %d --> %d --> %d\n", i, j, k);
                        // printf("\tdelta = %f\n", delta);
                        if (delta < delta_max)
                        {
                            delta_max = delta;
                        }
                    }
                }
                else // case 4
                {
                    //printf("\tCASE 4\n");
                    // from k it begins the drone leg: k --> u --> v
                    int u = drone_succ[k];
                    int v = drone_succ[u];
                    if (inst->min_time_drone[i][u][v] < 1e-6)
                    {
                        //printf("*** drone leg i --> u --> v is infeasible ***\n");
                        prev_truck_node = new_prev_truck_node;
                        last_takeoff_node = i;
                        i = k;
                        continue;
                    }
                    double truck_t_kv = 0.0;
                    int tmp = k;
                    while (tmp != v)
                    {
                        truck_t_kv += inst->truck_times[tmp][truck_succ[tmp]];
                        tmp = truck_succ[tmp];
                    }
                    double truck_t_iv = truck_t_kv - inst->truck_times[k][truck_succ[k]] + inst->truck_times[i][truck_succ[k]];
                    // printf("\ttruck_t_kv = %f\n", truck_t_kv);
                    // printf("\ttruck_t_iv = %f\n", truck_t_iv);
                    // printf("\tmax_time_drone[i][u][v] = %f\n", inst->max_time_drone[i][u][v]);
                    // check if the drone leg can be inverted
                    if (truck_t_iv > inst->max_time_drone[i][u][v])
                    {
                        //printf("*** drone leg i --> u --> v is infeasible ***\n");
                        prev_truck_node = new_prev_truck_node;
                        last_takeoff_node = i;
                        i = k;
                        continue;
                    }
                    // compute the minimum feasible travel times related to the drone legs k --> u --> v and i --> u --> v
                    double t_kuv = inst->min_time_drone[k][u][v];
                    double t_iuv = inst->min_time_drone[i][u][v];
                    if (t_kuv < truck_t_kv)
                        t_kuv = truck_t_kv;
                    if (t_iuv < truck_t_iv)
                        t_iuv = truck_t_iv;
                    // printf("\tt_kuv = %f\n", t_kuv);
                    // printf("\tt_iuv = %f\n", t_iuv);
                    double delta = -t_mni - t_ijk - t_kuv + t_mnk + t_kji + t_iuv;
                    //printf("\tdelta = %f\n", delta);
                    if (delta < 0)
                    {
                        // it is convenient to revert the drone leg
                        // printf("*** it is convenient to reverse the drone leg %d --> %d --> %d\n", i, j, k);
                        // printf("\tdelta = %f\n", delta);
                        if (delta < delta_max)
                        {
                            delta_max = delta;
                        }
                    }
                }
            }
            // set the last node served by the truck
            prev_truck_node = new_prev_truck_node;
            // set the last node from which the drone took off
            last_takeoff_node = i;
            // advance i
            i = k;
        }
        else
        {
            // advance the node
            prev_truck_node = i;
            i = drone_succ[i];
        }
    }
    return delta_max;
}

// revert the operation with the greatest time gain (if it exists)
// let us call C the current operation that we are evaluating
// case 1: the operation C is not preceded or succeeded by other operations
// case 2: the operation C is not preceded by another operation, but it is succeded by another operation
// case 3: the operation C is preceded by another operation, but it is not succeded by another operation
// case 4: the operation C is both preceded and succeded by another operation
double reverse_drone_leg(int *truck_seq, int *drone_seq, instance *inst)
{
    // a drone leg that starts from node 0 must not be inverted
    // check if from node 0 it starts a drone leg
    int t_i = 0;                                  // index node i in the truck sequence
    int d_i = 0;                                  // index node i in the drone sequence
    if (drone_seq[d_i + 1] != truck_seq[t_i + 1]) // drone leg that starts from 0
    {
        t_i++;
        d_i += 2;
        while (truck_seq[t_i] != drone_seq[d_i])
        {
            t_i++;
        }
    }
    else
    {
        t_i++;
        d_i++;
    }
    // at this point: truck_seq[t_i] = drone_seq[d_i]

    // truck and drone indices of the takeoff and landing nodes related to the drone leg with the greatest gain
    int best_t_i = -1;
    int best_t_k = -1;
    int best_d_i = -1;
    int best_d_k = -1;

    double delta_max = 0;
    // a drone leg that ends in the final node must not be inverted
    while (truck_seq[t_i + 1] != inst->dimension - 1 && drone_seq[d_i + 2] != inst->dimension - 1)
    {
        if (truck_seq[t_i] != drone_seq[d_i])
            print_error("*** truck_seq[t_i] != drone_seq[d_i] ***\n");
        // printf("drone_seq[d_i] = %d\n", drone_seq[d_i]);
        // printf("truck_seq[t_i] = %d\n", truck_seq[t_i]);
        if (truck_seq[t_i + 1] != drone_seq[d_i + 1]) // drone leg found
        {
            // printf("DRONE LEG FOUND: %d --> %d --> %d\n", drone_seq[d_i], drone_seq[d_i + 1], drone_seq[d_i + 2]);
            // printf("{prev_truck_node = %d}\n", prev_truck_node);
            // drone leg: i --> j --> k
            // we want to invert it: k --> j --> i
            int d_j = d_i + 1;
            int d_k = d_i + 2;

            // compute the truck travel time to traverse the path k -> i and the path i -> k
            int tmp = t_i;
            double truck_t_ki = 0.0;
            double truck_t_ik = 0.0;
            while (truck_seq[tmp] != drone_seq[d_k])
            {
                truck_t_ik += inst->truck_times[truck_seq[tmp]][truck_seq[tmp + 1]];
                truck_t_ki += inst->truck_times[truck_seq[tmp + 1]][truck_seq[tmp]];
                tmp++;
            }
            int t_k = tmp; // truck index associated to the rendezvouz node k

            // check if the drone leg k --> j --> i is feasible
            //printf("\ttruck_t_ki = %f\n", truck_t_ki);
            //printf("\tmax_time_drone[k][j][i] = %f\n", inst->max_time_drone[k][j][i]);
            if (truck_t_ki > inst->max_time_drone[drone_seq[d_k]][drone_seq[d_j]][drone_seq[d_i]])
            {
                //printf("*** drone leg k --> j --> i is infeasible ***\n");
                // jump to the rendezvouz node
                t_i = t_k;
                d_i = d_k;
                continue;
            }
            // compute the minimum feasible travel times related to the drone legs i --> j --> k and k --> j --> i
            double t_ijk = inst->min_time_drone[drone_seq[d_i]][drone_seq[d_j]][drone_seq[d_k]];
            double t_kji = inst->min_time_drone[drone_seq[d_k]][drone_seq[d_j]][drone_seq[d_i]];
            if (t_ijk < truck_t_ik)
                t_ijk = truck_t_ik;
            if (t_kji < truck_t_ki)
                t_kji = truck_t_ki;
            // printf("\tt_ijk = %f\n", t_ijk);
            // printf("\tt_kji = %f\n", t_kji);
            // printf("{last_takeoff_node = %d}\n", last_takeoff_node);
            if (drone_seq[d_i - 1] == truck_seq[t_i - 1]) // case 1|2
            {
                if (drone_seq[d_k + 1] == truck_seq[t_k + 1]) // case 1
                {
                    // printf("\tCASE 1\n");
                    // delta = time gained reversing the drone leg
                    double delta = inst->truck_times[truck_seq[t_i - 1]][truck_seq[t_i]] + t_ijk + inst->truck_times[truck_seq[t_k]][truck_seq[t_k + 1]];
                    delta = delta - inst->truck_times[truck_seq[t_i - 1]][truck_seq[t_k]] - t_kji - inst->truck_times[truck_seq[t_i]][truck_seq[t_k + 1]];
                    //printf("\tdelta = %f\n", delta);
                    if (delta > 0)
                    {
                        // it is convenient to revert the drone leg
                        // printf("*** it is convenient to reverse the drone leg %d --> %d --> %d\n", drone_seq[d_i], drone_seq[d_j], drone_seq[d_k]);
                        // printf("\tdelta = %f\n", delta);
                        if (delta > delta_max)
                        {
                            delta_max = delta;
                            best_t_i = t_i;
                            best_t_k = t_k;
                            best_d_i = d_i;
                            best_d_k = d_k;
                        }
                    }
                }
                else // case 2
                {
                    // printf("\tCASE 2\n");
                    // from k it begins the drone leg: k --> u --> v
                    int d_u = d_k + 1;
                    int d_v = d_k + 2;
                    if (inst->min_time_drone[drone_seq[d_i]][drone_seq[d_u]][drone_seq[d_v]] < 1e-6)
                    {
                        //printf("*** drone leg i --> u --> v is infeasible ***\n");
                        // jump to the rendezvouz node
                        t_i = t_k;
                        d_i = d_k;
                        continue;
                    }
                    // compute the truck path length to go from k to v
                    double truck_t_kv = 0.0;
                    int tmp = t_k;
                    while (truck_seq[tmp] != drone_seq[d_v])
                    {
                        truck_t_kv += inst->truck_times[truck_seq[tmp]][truck_seq[tmp + 1]];
                        tmp++;
                    }
                    double truck_t_iv = truck_t_kv - inst->truck_times[truck_seq[t_k]][truck_seq[t_k + 1]] + inst->truck_times[truck_seq[t_i]][truck_seq[t_k + 1]];
                    // check if the drone leg can be inverted
                    //printf("\ttruck_t_iv = %f\n", truck_t_iv);
                    //printf("\tmax_time_drone[i][u][v] = %f\n", inst->max_time_drone[i][u][v]);
                    if (truck_t_iv > inst->max_time_drone[drone_seq[d_i]][drone_seq[d_u]][drone_seq[d_v]])
                    {
                        //printf("*** drone leg i --> u --> v is infeasible ***\n");
                        // jump to the rendezvouz node
                        t_i = t_k;
                        d_i = d_k;
                        continue;
                    }
                    // compute the minimum feasible travel times related to the drone legs k --> u --> v and i --> u --> v
                    double t_kuv = inst->min_time_drone[drone_seq[d_k]][drone_seq[d_u]][drone_seq[d_v]];
                    double t_iuv = inst->min_time_drone[drone_seq[d_i]][drone_seq[d_u]][drone_seq[d_v]];
                    if (t_kuv < truck_t_kv)
                        t_kuv = truck_t_kv;
                    if (t_iuv < truck_t_iv)
                        t_iuv = truck_t_iv;
                    //printf("\tt_kuv = %f\n", t_kuv);
                    //printf("\tt_iuv = %f\n", t_iuv);
                    //double delta = -inst->truck_times[prev_truck_node][i] - t_ijk - t_kuv + inst->truck_times[prev_truck_node][k] + t_kji + t_iuv;
                    double delta = inst->truck_times[truck_seq[t_i - 1]][truck_seq[t_i]] + t_ijk + t_kuv;
                    delta = delta - inst->truck_times[truck_seq[t_i - 1]][truck_seq[t_k]] - t_kji - t_iuv;
                    //printf("\tdelta = %f\n", delta);
                    if (delta > 0)
                    {
                        // it is convenient to revert the drone leg
                        // printf("*** it is convenient to reverse the drone leg %d --> %d --> %d\n", drone_seq[d_i], drone_seq[d_j], drone_seq[d_k]);
                        // printf("\tdelta = %f\n", delta);
                        if (delta > delta_max)
                        {
                            delta_max = delta;
                            best_t_i = t_i;
                            best_t_k = t_k;
                            best_d_i = d_i;
                            best_d_k = d_k;
                        }
                    }
                }
            }
            else // case 3|4
            {
                // a drone lands at node i
                // lets call this drone leg: m --> n --> i
                int d_n = d_i - 1;
                int d_m = d_i - 2;
                if (inst->min_time_drone[drone_seq[d_m]][drone_seq[d_n]][drone_seq[d_i]] < 1e-6)
                {
                    //printf("*** drone leg m --> n --> k is infeasible ***\n");
                    // jump to the rendezvouz node
                    t_i = t_k;
                    d_i = d_k;
                    continue;
                }

                double truck_t_mi = 0.0;
                int tmp = t_i;
                while (truck_seq[tmp] != drone_seq[d_m])
                {
                    truck_t_mi += inst->truck_times[truck_seq[tmp - 1]][truck_seq[tmp]];
                    tmp--;
                }
                double truck_t_mk = truck_t_mi - inst->truck_times[truck_seq[t_i - 1]][truck_seq[t_i]] + inst->truck_times[truck_seq[t_i - 1]][truck_seq[t_k]];
                // check if the drone leg can be inverted
                // printf("\ttruck_t_mi = %f\n", truck_t_mi);
                // printf("\ttruck_t_mk = %f\n", truck_t_mk);
                // printf("\tmax_time_drone[m][n][i] = %f\n", inst->max_time_drone[m][n][i]);
                if (truck_t_mk > inst->max_time_drone[drone_seq[d_m]][drone_seq[d_n]][drone_seq[d_k]])
                {
                    //printf("*** drone leg m --> n --> k is infeasible ***\n");
                    // jump to the rendezvouz node
                    t_i = t_k;
                    d_i = d_k;
                    continue;
                }
                // compute the minimum feasible travel times related to the drone legs m --> n --> i and m --> n --> k
                double t_mni = inst->min_time_drone[drone_seq[d_m]][drone_seq[d_n]][drone_seq[d_i]];
                double t_mnk = inst->min_time_drone[drone_seq[d_m]][drone_seq[d_n]][drone_seq[d_k]];
                if (t_mni < truck_t_mi)
                    t_mni = truck_t_mi;
                if (t_mnk < truck_t_mk)
                    t_mnk = truck_t_mk;
                // printf("\tt_mni = %f\n", t_mni);
                // printf("\tt_mnk = %f\n", t_mnk);
                if (drone_seq[d_k + 1] == truck_seq[t_k + 1]) // case 3
                {
                    // printf("\tCASE 3\n");
                    //double delta = -t_mni - t_ijk - inst->truck_times[k][drone_succ[k]] + t_mnk + t_kji + inst->truck_times[i][drone_succ[k]];
                    double delta = t_mni + t_ijk + inst->truck_times[truck_seq[t_k]][truck_seq[t_k + 1]];
                    delta = delta - t_mnk - t_kji - inst->truck_times[truck_seq[t_i]][truck_seq[t_k + 1]];
                    //printf("\tdelta = %f\n", delta);
                    if (delta > 0)
                    {
                        // it is convenient to revert the drone leg
                        // printf("*** it is convenient to reverse the drone leg %d --> %d --> %d\n", drone_seq[d_i], drone_seq[d_j], drone_seq[d_k]);
                        // printf("\tdelta = %f\n", delta);
                        if (delta > delta_max)
                        {
                            delta_max = delta;
                            best_t_i = t_i;
                            best_t_k = t_k;
                            best_d_i = d_i;
                            best_d_k = d_k;
                        }
                    }
                }
                else // case 4
                {
                    // printf("\tCASE 4\n");
                    // from k it begins the drone leg: k --> u --> v
                    int d_u = d_k + 1;
                    int d_v = d_k + 2;
                    if (inst->min_time_drone[drone_seq[d_i]][drone_seq[d_u]][drone_seq[d_v]] < 1e-6)
                    {
                        //printf("*** drone leg i --> u --> v is infeasible ***\n");
                        // jump to the rendezvouz node
                        t_i = t_k;
                        d_i = d_k;
                        continue;
                    }
                    double truck_t_kv = 0.0;
                    int tmp = t_k;
                    while (truck_seq[tmp] != drone_seq[d_v])
                    {
                        truck_t_kv += inst->truck_times[truck_seq[tmp]][truck_seq[tmp + 1]];
                        tmp++;
                    }
                    double truck_t_iv = truck_t_kv - inst->truck_times[truck_seq[t_k]][truck_seq[t_k + 1]] + inst->truck_times[truck_seq[t_i]][truck_seq[t_k + 1]];
                    // printf("\ttruck_t_kv = %f\n", truck_t_kv);
                    // printf("\ttruck_t_iv = %f\n", truck_t_iv);
                    // printf("\tmax_time_drone[i][u][v] = %f\n", inst->max_time_drone[i][u][v]);
                    // check if the drone leg can be inverted
                    if (truck_t_iv > inst->max_time_drone[drone_seq[d_i]][drone_seq[d_u]][drone_seq[d_v]])
                    {
                        //printf("*** drone leg i --> u --> v is infeasible ***\n");
                        // jump to the rendezvouz node
                        t_i = t_k;
                        d_i = d_k;
                        continue;
                    }
                    // compute the minimum feasible travel times related to the drone legs k --> u --> v and i --> u --> v
                    double t_kuv = inst->min_time_drone[drone_seq[d_k]][drone_seq[d_u]][drone_seq[d_v]];
                    double t_iuv = inst->min_time_drone[drone_seq[d_i]][drone_seq[d_u]][drone_seq[d_v]];
                    if (t_kuv < truck_t_kv)
                        t_kuv = truck_t_kv;
                    if (t_iuv < truck_t_iv)
                        t_iuv = truck_t_iv;
                    // printf("\tt_kuv = %f\n", t_kuv);
                    // printf("\tt_iuv = %f\n", t_iuv);
                    double delta = t_mni + t_ijk + t_kuv;
                    delta = delta - t_mnk - t_kji - t_iuv;
                    //printf("\tdelta = %f\n", delta);
                    if (delta > 0)
                    {
                        // it is convenient to revert the drone leg
                        // printf("*** it is convenient to reverse the drone leg %d --> %d --> %d\n", drone_seq[d_i], drone_seq[d_j], drone_seq[d_k]);
                        // printf("\tdelta = %f\n", delta);
                        if (delta > delta_max)
                        {
                            delta_max = delta;
                            best_t_i = t_i;
                            best_t_k = t_k;
                            best_d_i = d_i;
                            best_d_k = d_k;
                        }
                    }
                }
            }
            // jump to the rendezvouz node
            t_i = t_k;
            d_i = d_k;
        }
        else
        {
            // advance the node
            d_i++;
            t_i++;
        }
    }
    if (delta_max > 0) // reverse the most convenient drone leg found (if any)
    {
        // reverse the drone leg
        int tmp = drone_seq[best_d_i];
        drone_seq[best_d_i] = drone_seq[best_d_k];
        drone_seq[best_d_k] = tmp;
        // reverse the truck leg
        int left = best_t_i;
        int right = best_t_k;
        while (left < right)
        {
            int tmp = truck_seq[left];
            truck_seq[left] = truck_seq[right];
            truck_seq[right] = tmp;
            left++;
            right--;
        }
    }
    return delta_max;
}

void nearest_neighbours(instance *inst, int *seq, int options)
{
    if (options < 1)
        options = 1;
    // increase the number of options by 1, in order to handle an extra case

    // boolean array : selected[i] = 1 if node i has been already selected
    int *selected = (int *)calloc(inst->dimension, sizeof(int));
    int current = 0; // Index of the current node
    seq[0] = current;
    selected[current] = 1;

    seq[inst->dimension - 1] = inst->dimension - 1;
    selected[inst->dimension - 1] = 1;

    // flag = 1 if in the last iteration the current node has been updated (0 otherwise)
    int flag = 1;

    // Build the circuit adding inst->dimension - 2 nodes (the first and the last node are fixed)
    for (int count = 1; count < inst->dimension - 1; count++)
    {

        // Check the number of nodes not yet selected
        // options must be smaller than or equal to number of available nodes
        if (inst->dimension - 1 - count < options)
        {
            options = inst->dimension - 1 - count;
        }

        double min_time[options]; // Minimum truck travel times
        int min_node[options];    // Closest nodes indices

        for (int i = 0; i < options; ++i)
        {
            min_time[i] = DBL_MAX;
            min_node[i] = -1;
        }

        // select the closest (options) nodes w.r.t. to the current one
        for (int i = 1; i < inst->dimension - 1; i++)
        {
            if (selected[i] == 0) // node i has not been selected yet
            {
                double time = inst->truck_times[current][i];
                if (time < min_time[options - 1])
                {
                    int k;
                    for (k = options - 2; k >= 0; k--)
                    { // if options == 1 => it does not enter in the loop => k+1 = 0: OK
                        if (time >= min_time[k])
                            break; // node i is the (k+1)-th nearest node (currently)
                    }

                    // the selected node is the (k+1)-th closest node
                    // right-shift the elements that are farthest by one
                    for (int j = options - 2; j > k; j--)
                    { // if options == 1 => k=-1, j=-1 => it does not enter in the loop
                        min_time[j + 1] = min_time[j];
                        min_node[j + 1] = min_node[j];
                    }
                    min_time[k + 1] = time;
                    min_node[k + 1] = i;
                }
            }
        }

        // Minimum node random selection
        int h = rand() % options;
        seq[count] = min_node[h];
        selected[min_node[h]] = 1;

        // Update the current node with prob 0.5
        // If in the last iteration the node has not been updated, it must be updated in this iteration
        if (flag && rand() % 2)
        {
            flag = 0;
        }
        else
        {
            current = min_node[h];
            flag = 1;
        }
    }

    free(selected);
}

double compute_score(int *truck_seq, int *drone_seq, instance *inst)
{
    if (truck_seq[0] != 0 || drone_seq[0] != 0)
        return -1;

    double obj = 0.0;
    int t_i = 0;
    int d_i = 0;
    while (truck_seq[t_i] != inst->dimension - 1)
    {
        if (truck_seq[t_i + 1] != drone_seq[d_i + 1]) // drone leg found
        {
            int d_j = d_i + 1;
            int d_k = d_i + 2;

            // compute the truck travel time to traverse the path i -> k
            int tmp = t_i;
            double truck_t_ik = 0.0;
            while (truck_seq[tmp] != drone_seq[d_k])
            {
                truck_t_ik += inst->truck_times[truck_seq[tmp]][truck_seq[tmp + 1]];
                tmp++;
            }
            int t_k = tmp; // truck index associated to the rendezvouz node k

            // compute the minimum feasible travel times related to the drone legs i --> j --> k and k --> j --> i
            if (truck_t_ik > inst->max_time_drone[drone_seq[d_i]][drone_seq[d_j]][drone_seq[d_k]])
            {
                printf("truck path exceeds max drone travel time %d->%d->%d\n", drone_seq[d_i], drone_seq[d_j], drone_seq[d_k]);
                return -1;
            }
            double t_ijk = inst->min_time_drone[drone_seq[d_i]][drone_seq[d_j]][drone_seq[d_k]];
            if (t_ijk < truck_t_ik)
                t_ijk = truck_t_ik;
            obj += t_ijk;
            t_i = t_k;
            d_i = d_k;
        }
        else
        {
            obj += inst->truck_times[truck_seq[t_i]][truck_seq[t_i + 1]];
            t_i++;
            d_i++;
        }
    }
    return obj;
}

void partial_fstsp(instance *inst, int *seq, int seq_dim, instance *new_inst)
{
    initialize_instance(new_inst);
    new_inst->dimension = seq_dim;
    new_inst->nodes = (node *)calloc(seq_dim, sizeof(node));
    // initialize truck travel times matrix
    allocate_mem_truck_times_and_dists(new_inst);
    // initialize drone travel times matrix
    allocate_mem_drone_min_max_times(new_inst);

    // copy nodes info
    for (int i = 0; i < seq_dim; i++)
    {
        new_inst->nodes[i].id = inst->nodes[seq[i]].id;
        new_inst->nodes[i].x = inst->nodes[seq[i]].x;
        new_inst->nodes[i].y = inst->nodes[seq[i]].y;
        new_inst->nodes[i].weight = inst->nodes[seq[i]].weight;
    }
    // copy truck travel times for the new instance
    for (int i = 0; i < seq_dim; i++)
    {
        for (int j = 0; j < seq_dim; j++)
        {
            if (i == j)
                continue;
            new_inst->truck_times[i][j] = inst->truck_times[seq[i]][seq[j]];
            new_inst->truck_dists[i][j] = inst->truck_dists[seq[i]][seq[j]];
        }
    }
    // copy min max drone times for the new instance
    for (int i = 0; i < seq_dim; i++)
    {
        for (int j = 0; j < seq_dim; j++)
        {
            if (i == j)
                continue;
            for (int k = 0; k < seq_dim; k++)
            {
                if (i == k || j == k)
                    continue;
                new_inst->min_time_drone[i][j][k] = inst->min_time_drone[seq[i]][seq[j]][seq[k]];
                new_inst->max_time_drone[i][j][k] = inst->max_time_drone[seq[i]][seq[j]][seq[k]];
            }
        }
    }
}