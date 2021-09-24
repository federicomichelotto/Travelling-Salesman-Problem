int heuristic_solver(instance *inst)
{
    // get timestamp
    getTimeStamp(&inst->timestamp_start);

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
            obj_i = nearest_neighbours(inst, i, succ_i, inst->param.grasp_choices);
            printf("Best objective value for node %d: %f\n", i + 1, obj_i);

            if (obj_i < min_obj)
            {
                min_obj = obj_i;
                inst->z_best = obj_i;
                for (int j = 0; j < inst->dimension; j++)
                    inst->succ[j] = succ_i[j];
                save_and_plot_solution(inst, i + 1);
            }
        }
        if (inst->param.opt)
            inst->z_best += two_opt_v2(inst, inst->succ, 0);

        save_and_plot_solution(inst, -1);

        printf("\nBest objective value: %f\n", min_obj);
        if (inst->param.opt)
            printf("Best objective value (optimized by 2-opt): %f\n", inst->z_best);
        break;
    case 1: // GRASP Nearest Neighbours + 2-opt (random starting node)
        min_obj = nearest_neighbours(inst, rand() % inst->dimension, succ_i, inst->param.grasp_choices);
        for (int j = 0; j < inst->dimension; j++)
            inst->succ[j] = succ_i[j];
        inst->z_best = min_obj;
        inst->z_best += two_opt_v2(inst, inst->succ, 0);

        save_and_plot_solution(inst, -1);
        
        printf("\nBest objective value: %f\n", min_obj);
        if (inst->param.opt)
            printf("Best objective value (optimized by 2-opt): %f\n", inst->z_best);
        break;
    case 2: // Extra Mileage
        for (int i = 0; i < inst->dimension; i++)
        {
            obj_i = extra_mileage(inst, succ_i, i);
            if (obj_i < min_obj)
            {
                min_obj = obj_i;
                inst->z_best = obj_i;
                for (int j = 0; j < inst->dimension; j++)
                    inst->succ[j] = succ_i[j];
                save_and_plot_solution(inst, i + 1);
            }
        }
        if (inst->param.opt)
            inst->z_best += two_opt_v2(inst, inst->succ, 0);

        save_and_plot_solution(inst, -1);

        printf("\nBest objective value: %f\n", min_obj);
        if (inst->param.opt)
            printf("Best objective value (optimized by 2-opt): %f\n", inst->z_best);
        break;
    case 3: // Extra Mileage with furthest starting nodes
        min_obj = extra_mileage_furthest_starting_nodes(inst, succ_i);
        inst->z_best = min_obj;
        for (int j = 0; j < inst->dimension; j++)
            inst->succ[j] = succ_i[j];
        if (inst->param.opt)
            inst->z_best += two_opt_v2(inst, inst->succ, 0);

        save_and_plot_solution(inst, -1);

        printf("\nBest objective value: %f\n", min_obj);
        if (inst->param.opt)
            printf("Best objective value (optimized by 2-opt): %f\n", inst->z_best);
        break;

    default:
        fprintf(stderr, "ERROR: Model type %d not available.\n", inst->model_type);
        break;
    }

    free(succ_i);

    // get timestamp
    getTimeStamp(&inst->timestamp_finish);

    return 0;
}

double nearest_neighbours(instance *inst, int starting_node, int *succ, int options)
{
    double obj = 0.0;

    int *selected = (int *)calloc(inst->dimension, sizeof(int));

    // Successor array initialization
    for (int i = 0; i < inst->dimension; i++)
        succ[i] = -1;

    selected[starting_node] = 1;
    int current = starting_node; // Index of the current node

    if (inst->param.verbose >= DEBUG && options > 1)
        printf("\nGRASP approach selected, available option for each node %d\n", options);

    // Build the circuit adding inst->dimension - 1 edges
    for (int count = 0; count < inst->dimension - 1; count++)
    {

        // Check available nodes (remember to ignore the current one)
        if (inst->dimension - count - 1 < options)
        {
            options = inst->dimension - count - 1;
        }

        double min_dist[options]; // Minimum distances
        int min_node[options];    // Closest nodes indices

        for (int i = 0; i < options; ++i)
        {
            min_dist[i] = DBL_MAX;
            min_node[i] = -1;
        }
       
        // select the closest (options) nodes w.r.t. to the current node
        for (int i = 0; i < inst->dimension; i++)
        {
            if (selected[i] == 0) // i has not been selected yet
            {
                double distance = dist(current, i, inst);
                if (distance < min_dist[options-1])
                {
                    int k;
                    for (k = options - 2; k >= 0 ; k--){ // if options == 1 => it does not enter in the loop => k+1 = 0: OK
                        if (distance >= min_dist[k])
                            break; // node i is the (k+1)-th nearest node (currently)
                    }
                    
                    // shift elements right by one 
                    for (int j = options-2; j>k; j--){ // if options == 1 => k=-1, j=-1 => it does not enter in the loop 
                        min_dist[j+1] = min_dist[j];
                        min_node[j+1] = min_node[j];
                    }
                    min_dist[k+1] = distance;
                    min_node[k+1] = i;
                }
            }
        }
        
        // Minimum node random selection
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

    free(selected);

    return obj;
}

// node: node to add in the circuit
double add_node_extra_mileage(instance *inst, int *succ, int node)
{
    double min_extra_mileage = DBL_MAX;
    int min_i = -1;
    // compute the extra mileage to add the node "node"
    for (int i = 0; i < inst->dimension; i++)
    {
        if (succ[i] != -1) // node i is part of the circuit
        {
            double extra_mileage = dist(i, node, inst) + dist(node, succ[i], inst) - dist(i, succ[i], inst);
            if (extra_mileage < min_extra_mileage)
            {
                min_extra_mileage = extra_mileage;
                min_i = i; // represent the edge (i, succ[i])
            }
        }
    }
    if (min_i == -1)
        print_error("Error add_node_extra_mileage.");
    succ[node] = succ[min_i];
    succ[min_i] = node;
    return min_extra_mileage;
}

double extra_mileage(instance *inst, int *succ, int starting_node)
{

    double obj = 0;                                                        // Objective value
    double *succ_dist = (double *)calloc(inst->dimension, sizeof(double)); // array of the distance between the node i and succ[i]

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

        double min_extra_mileage = DBL_MAX; // Insertion of node selected_node in the circuit
        int min_i = -1;                     // (min_i, succ[min_i]): edge edge that corresponds to the smallest extra mileage to add the node selected_node to the circuit
        int selected_node;                  // node to add in the circuit

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
                    if (extra_mileage < min_extra_mileage)
                    {
                        min_extra_mileage = extra_mileage;
                        min_i = i; // (i, succ[i])
                        selected_node = k;
                    }
                }
            }
        }

        if (min_i == -1)
            print_error("Error add_node_extra_mileage.");

        succ[selected_node] = succ[min_i];
        succ_dist[selected_node] = dist(selected_node, succ[min_i], inst);
        succ[min_i] = selected_node;
        succ_dist[min_i] = dist(min_i, selected_node, inst);

        obj += min_extra_mileage; // Update objective value
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
        double min_value = DBL_MAX; // Insertion of node selected_node in the circuit
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

double two_opt(instance *inst, int *succ, int maxMoves)
{
    // if (inst->param.verbose >= DEBUG)
    //     print_message("Inside 2-opt function");
    int optimal = 0;
    int moves = 0;
    int iter = 0;
    double total_delta = 0.0;
    // if (inst->param.verbose >= DEBUG)
    //     printf("Initial incumbent = %f\n", inst->z_best);

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
                if (!inst->param.ticks)
                {
                    double ts_current;
                    getTimeStamp(&ts_current);
                    double time_left = inst->time_limit - (ts_current - inst->timestamp_start);
                    if (time_left < 1.0)
                        return total_delta;
                }
                if (i == j)
                    continue;
                if (succ[i] == j || i == succ[j])
                    continue;
                if (succ[i] == -1 || succ[j] == -1)
                    continue;

                optimal = 1;

                int a = i;
                int b = succ[a];
                int c = j;
                int d = succ[c];
                double originalDist = dist(a, b, inst) + dist(c, d, inst);
                double newDist = dist(a, c, inst) + dist(b, d, inst);

                if (newDist < originalDist)
                { // crossing

                    // update inc
                    double delta = newDist - originalDist;
                    total_delta += delta;
                    //                     if (inst->param.verbose >= DEBUG)
                    //                         printf("%d° iteration - total delta = %f (delta = %f)\n", iter + 1, total_delta, delta);

                    // reverse tour
                    if (reverse_successors(succ, inst->dimension, b, c))
                        print_error("Error in reverse_successors");
                    succ[a] = c;
                    succ[b] = d;

                    optimal = 0;
                    moves++;
                    break;
                }
            }
            if (!optimal)
                break;
        }
        iter++;
        save_and_plot_solution_general(inst, succ, iter);
    }

    return total_delta;
}

// move applied on the most negative delta
double two_opt_v2(instance *inst, int *succ, int maxMoves)
{
    //print_message("Inside 2-opt_v2 function");
    //printf("Initial incumbent = %f\n", inst->z_best);
    // For each couple of edges (a,c) and (b,d) so that they are not subsequent,
    // If the given solution does not have any crossing: return 0, else return #crossing found;
    int iter = 1;
    double total_delta = 0.0;
    while (!maxMoves || iter <= maxMoves)
    {
        double min_delta = DBL_MAX;
        int a, b;
        for (int i = 0; i < inst->dimension; i++)
        {
            for (int j = 0; j < inst->dimension; j++)
            {
                if (!inst->param.ticks)
                {
                    double ts_current;
                    getTimeStamp(&ts_current);
                    double time_left = inst->time_limit - (ts_current - inst->timestamp_start);
                    if (time_left < 1.0)
                        return total_delta;
                }
                // look for two nodes that are not connected
                if (i == j)
                    continue;
                if (succ[i] == j || i == succ[j])
                    continue;
                if (succ[i] == -1 || succ[j] == -1)
                    continue;

                double delta = dist(i, j, inst) + dist(succ[i], succ[j], inst) - dist(i, succ[i], inst) - dist(j, succ[j], inst);
                if (delta < min_delta)
                {
                    min_delta = delta;
                    a = i;
                    b = j;
                }
            }
        }
        if (min_delta == DBL_MAX)
            print_error("Error min_delta in two_opt_v2");
        if (min_delta >= -1e-5)
            break;

        int c = succ[a];
        int d = succ[b];

        // reverse the path from c to b
        if (reverse_successors(succ, inst->dimension, c, b))
            print_error("Error in reverse_successors in two_opt_v2");

        // make the move
        succ[a] = b;
        succ[c] = d;
        // update incumbent
        total_delta += min_delta;
        //printf("%d° iteration - new incumbent = %f (delta = %f)\n", iter, inst->z_best, min_delta);
        iter++;
        save_and_plot_solution_general(inst, succ, iter);
    }

    return total_delta;
}

int reverse_successors(int *succ, int size, int start, int end)
{
    if (start >= size || end >= size)
        print_error("Error in reverse_successors : index out of bound");
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
        {
            print_error("Error in reverse_successors : counter > size");
            return 1;
        }
    }

    return 0;
}