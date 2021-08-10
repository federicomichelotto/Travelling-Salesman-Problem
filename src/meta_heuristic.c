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
            if (inst->param.grasp == 1)
                obj_i = nearest_neighbours(inst, i, succ_i, inst->param.grasp_choices);
            else
                obj_i = nearest_neighbours(inst, i, succ_i, 1);

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

        inst->z_best += two_opt_v2(inst, inst->succ, 0);
        save_and_plot_solution(inst, inst->dimension);

        printf("Best objective value: %f\n", min_obj);
        printf("Best objective value (optimized by 2-opt): %f\n", inst->z_best);
        break;
    case 1: // Extra Mileage
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

        inst->z_best += two_opt_v2(inst, inst->succ, 0);
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

        inst->z_best += two_opt_v2(inst, inst->succ, 0);
        save_and_plot_solution(inst, 2);

        printf("Best objective value: %f\n", min_obj);
        printf("Best objective value (optimized by 2-opt): %f\n", inst->z_best);
        break;
    case 3: // Tabu Search

        if (inst->param.grasp == 1)
            min_obj = nearest_neighbours(inst, 0, inst->succ, inst->param.grasp_choices);
        else
            min_obj = nearest_neighbours(inst, 0, inst->succ, 1);
        save_and_plot_solution(inst, 1);
        inst->z_best = min_obj;
        tabu_search(inst);

        save_and_plot_solution(inst, -1);

        printf("Best objective value: %f\n", min_obj);
        printf("Best objective value (optimized by tabu search): %f\n", inst->z_best);
        break;

    case 4: // Genetic algorithm
        genetic_v2(inst);
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
        printf("GRASP approach selected, available option for each node %d\n", options);

    // Build the circuit adding inst->dimension - 1 edges
    for (int count = 0; count < inst->dimension - 1; count++)
    {

        // Check available nodes (remember to ignore the current one)
        if (inst->dimension - count - 1 < options)
        {
            if (inst->param.verbose >= DEBUG)
                printf("Too many choices (%d) w.r.t. available nodes (%d). Available choices are now (%d)\n", options, inst->dimension - count - 1, inst->dimension - count - 1);

            options = inst->dimension - count - 1;
            int not_selected = -1;
            for (int i = 0; i < inst->dimension; i++)
                if (succ[i] == -1)
                    not_selected++;
            if (inst->param.verbose >= DEBUG)
                printf("Nodes not yet selected are %d\n", not_selected);
        }

        double min_dist[options]; // Minimum distances
        int min_node[options];    // Closest nodes indices

        for (int i = 0; i < options; ++i)
        {
            min_dist[i] = CPX_INFBOUND;
            min_node[i] = -1;
        }

        int *chosen = (int *)calloc(inst->dimension, sizeof(int));
        int k = 0;

        // While a choice can be made, select the closest node
        while (k < options)
        {
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

            chosen[min_node[k]] = 1; // i-th node was chosen as closest
            k++;
        }

        //        for (int j = 0; j < options; ++j) {
        //            printf("%d : %f\n", min_node[j], min_dist[j]);
        //        }
        //        printf("------------\n");

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

        free(chosen);
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

double two_opt(instance *inst, int *succ, int maxMoves)
{
    // if (inst->param.verbose >= DEBUG)
    //     print_message("Inside 2-opt function");
    int optimal = 0;
    int moves = 0;
    int iter = 0;
    double incumbent = inst->z_best;
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
                    incumbent = incumbent + delta;
                    // if (inst->param.verbose >= DEBUG)
                    //     printf("%d° iteration - new incumbent = %f (delta = %f)\n", iter + 1, incumbent, delta);

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

    return incumbent;
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

void tabu_search(instance *inst)
{
    int *tabu_list = (int *)calloc(inst->dimension, sizeof(int));
    // tabu_list[node] = 0: not a tabu node
    // tabu_list[node] > 0: tabu node
    int *temp_succ = (int *)calloc(inst->dimension, sizeof(int));
    // init temp_succ
    for (int i = 0; i < inst->dimension; i++)
        temp_succ[i] = inst->succ[i];
    double incumbent = inst->z_best; // current incumbent
    int iter = 1;
    int best_iter = 1;
    int tenure = inst->dimension / 10;
    printf("Initial phase: INTENSIFICATION PHASE\n");

    double time_left;

    // get timestamp
    struct timespec timestamp;
    if (clock_gettime(CLOCK_REALTIME, &timestamp) == -1)
        print_error("Error clock_gettime");
    inst->timestamp_finish = timestamp.tv_sec + timestamp.tv_nsec * pow(10, -9);
    // update time left
    time_left = inst->time_limit - (inst->timestamp_finish - inst->timestamp_start);

    while (time_left > 0.5)
    {
        if (!(iter % 100))
        {
            if (tenure == 20)
            {
                tenure = inst->dimension / 10;
                printf("*** INTENSIFICATION PHASE activated, time left = %f seconds ***\n", time_left);
            }
            else
            {
                tenure = 20;
                printf("*** DIVERSIFICATION PHASE activated, time left = %f seconds ***\n", time_left);
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

        int c = temp_succ[a];
        int d = temp_succ[b];

        // reverse the path from c to b
        if (reverse_successors(temp_succ, inst->dimension, c, b))
            print_error("Error in reverse_successors");

        // make the move
        temp_succ[a] = b;
        temp_succ[c] = d;
        // update incumbent
        incumbent += min_delta;

        if (min_delta > 0)
            tabu_list[d] = iter; // add the node in the tabu list
        if (incumbent < inst->z_best)
        {
            printf("new incumbent found = %f\n", incumbent);
            best_iter = iter;
            inst->z_best = incumbent;
            // save best solution
            for (int j = 0; j < inst->dimension; j++)
                inst->succ[j] = temp_succ[j];
            save_and_plot_solution(inst, iter);
        }
        iter++;

        // get timestamp
        if (clock_gettime(CLOCK_REALTIME, &timestamp) == -1)
            print_error("Error clock_gettime");
        inst->timestamp_finish = timestamp.tv_sec + timestamp.tv_nsec * pow(10, -9);
        // update time limit
        time_left = inst->time_limit - (inst->timestamp_finish - inst->timestamp_start);
    }
    free(temp_succ);
    free(tabu_list);
}

// TODO : FINISH GENETIC ALGORITHM

int genetic_v2(instance *inst)
{

    int size, children_size;
    printf("How many individual should have your population: ");
    scanf("%d", &size);

    printf("How many children should conceive your population each epoch: ");
    scanf("%d", &children_size);

    int areas = 5;
    int size_final = children_size * 2;
    int ind_per_area = (int)size_final / areas;
    if (ind_per_area > size)
        ind_per_area = size;
    population *final_population = (population *)calloc(size_final, sizeof(population));

    for (int i = 0; i < size_final; i++)
        final_population[i].chromosome = (int *)calloc(inst->dimension, sizeof(int));

    for (int area = 0; area < areas; area++)
    {

        population *individuals = (population *)calloc(size, sizeof(population));
        printf("\nAREA %d:\nGENERATION STEP: ", area);
        printf("%d individuals will be generated (# = %d individuals generated): \n", size, (int)ceil(size / 10));

        int start_population = rand() % inst->dimension;
        // Population
        for (int i = 0; i < size; i++)
        {
            // if (i < size * 0.8)
            // random_individual(inst, &individuals[i], start_population + i, 0);
            // //else
            //random_individual_2(inst, &individuals[i], start_population + i, 1);
            if (i % 10 == 0)
                random_individual(inst, &individuals[i], start_population + i, 1);
            else
                random_individual(inst, &individuals[i], start_population + i, 0);
            printf("Individual %d generated: fitness = %f\n", i + 1, individuals[i].fitness);
        }

        // get timestamp
        struct timespec timestamp;
        double time_left = DBL_MAX;
        // reset time_start
        if (clock_gettime(CLOCK_REALTIME, &timestamp) == -1)
            print_error("Error clock_gettime");
        inst->timestamp_start = timestamp.tv_sec + timestamp.tv_nsec * pow(10, -9);

        int no_improvement = 0;
        int epochs = 1;

        // Epochs
        double *average = (double *)calloc(1, sizeof(double));
        population *champion = (population *)calloc(1, sizeof(population));
        champion[epochs - 1].chromosome = (int *)calloc(inst->dimension, sizeof(int));

        int tournment_selection_size = 8;

        // Since there's time left
        while (time_left > 0.5)
        {
            printf("\nEPOCH %d:\n", epochs);
            if (epochs % 20 == 0)
            {
                tournment_selection_size = 3 + (rand() % 6);
                printf("*** new tournment_selection_size = %d\n", tournment_selection_size);
            }
            // // optimize
            // for (int i = 0; i < size; i++)
            // {
            //     if (rand() % 10 == 0)
            //         individuals[i].fitness += two_opt(inst, individuals[i].chromosome, 20);
            // }
            // Offspring generation
            for (int k = 0; k < children_size; k++)
            {
                general_alg(inst, individuals, tournment_selection_size, size, children_size);
            }

            epoch_champion_and_average(inst, individuals, size, &champion[epochs - 1], &average[epochs - 1]); // put the champion in individuals[0]
            population champion_child;
            champion_child.chromosome = (int *)calloc(inst->dimension, sizeof(int));

            one_point_crossover(inst, &individuals[0], &individuals[1 + rand() % (size - 2)], &champion_child);
            if (champion_child.fitness < individuals[size - 1].fitness)
            {
                printf("*** worst individual replaced: %f -> %f\n", individuals[size - 1].fitness, champion_child.fitness);
                free(individuals[size - 1].chromosome);
                individuals[size - 1].chromosome = champion_child.chromosome;
                individuals[size - 1].fitness = champion_child.fitness;
            }
            save_and_plot_solution_general(inst, champion[epochs - 1].chromosome, epochs - 1);

            printf("\t- Champion's fitness : %f\n", champion[epochs - 1].fitness);
            printf("\t- Average fitness : %f\n", average[epochs - 1]);

            if (epochs != 1 && champion[epochs - 1].fitness == champion[epochs - 2].fitness)
                no_improvement++;
            else
                no_improvement = 0;

            if (no_improvement > 20)
            {
                print_message("\nExit, no improvement were found in the last 20 epochs.");
                break;
            }

            if (clock_gettime(CLOCK_REALTIME, &timestamp) == -1)
                print_error("Error clock_gettime");
            inst->timestamp_finish = timestamp.tv_sec + timestamp.tv_nsec * pow(10, -9);
            // update time limit
            time_left = inst->time_limit - (inst->timestamp_finish - inst->timestamp_start);

            if (time_left > 0.5)
            {
                // OPTIMIZE the champion with one two 2-opt move
                // IF THE CHAMPION CANNOT BE OPTIMIZED ANYMORE (no_improvement > 0)
                // or there is a new champion, try to apply a two-opt move on another individual
                if (no_improvement == 0)
                {
                    //individuals[0].fitness += two_opt_v2(inst, individuals[0].chromosome, 1);
                }
                else
                {
                    printf("*** no_improvement *** \n");
                    //int index = 1 + rand() % (size - 1);
                    //individuals[index].fitness += two_opt_v2(inst, individuals[index].chromosome, 1);
                }
                epochs++;
                // realloc and alloc
                average = (double *)realloc(average, epochs * sizeof(double));
                champion = (population *)realloc(champion, epochs * sizeof(population));
                champion[epochs - 1].chromosome = (int *)calloc(inst->dimension, sizeof(int));
            }
        }

        printf("\n\nEPOCH SUMMARY");
        printf("\n%-10s| %-15s| %-15s", "Epoch", "Average", "Champion");
        for (int e = 0; e < epochs; e++)
        {
            printf("\n%-10d| %-15.4f| %-15.4f", e + 1, average[e], champion[e].fitness);
            //printf("* e = %d\n", e);
        }

        // rank individuals
        //rank(inst, individuals, size);
        // copy individuals in the final population
        int unit_step = (int)size / ind_per_area;
        unit_step = 1;

        // copy first individual
        for (int j = 0; j < inst->dimension; j++)
        {
            final_population[(area * ind_per_area)].chromosome[j] = individuals[0].chromosome[j];
        }
        final_population[(area * ind_per_area)].fitness = individuals[0].fitness;
        int count = 1;

        for (int i = 1; i < size; i++)
        {
            if (count == ind_per_area)
                break;
            // avoid duplicates (if it is possible)
            if (ind_per_area - count < size - i)
            {
                if (final_population[(area * ind_per_area) + count - 1].fitness == individuals[i].fitness)
                    continue;
            }
            // copy individual
            for (int j = 0; j < inst->dimension; j++)
            {
                final_population[(area * ind_per_area) + count].chromosome[j] = individuals[i].chromosome[j];
            }
            final_population[(area * ind_per_area) + count].fitness = individuals[i].fitness;
            count++;
        }

        free(average);
        free(champion);
        free(individuals);
    }
    printf("\n\nCHAMPIONS BY AREA SUMMARY\n");
    for (int i = 0; i < areas; i++)
        printf("\n%d) %-15.4f", i, final_population[i * ind_per_area].fitness);

    /// *** POPULATION champions *** ///

    population *individuals = final_population;

    printf("\n\nFINAL POPULATION:\n");
    for (int i = 0; i < size_final; i++)
        printf("\n%d) %-15.4f", i, final_population[i].fitness);

    // get timestamp
    struct timespec timestamp;
    double time_left = DBL_MAX;
    // reset time_start
    if (clock_gettime(CLOCK_REALTIME, &timestamp) == -1)
        print_error("Error clock_gettime");
    inst->timestamp_start = timestamp.tv_sec + timestamp.tv_nsec * pow(10, -9);

    int no_improvement = 0;
    int epochs = 1;

    // Epochs
    double *average = (double *)calloc(1, sizeof(double));
    population *champion = (population *)calloc(1, sizeof(population));
    champion[epochs - 1].chromosome = (int *)calloc(inst->dimension, sizeof(int));
    int tournment_selection_size = 3;

    epoch_champion_and_average(inst, individuals, size_final, &champion[0], &average[0]); // put the champion in individuals[0]
    // Since there's time left
    while (time_left > 0.5)
    {
        printf("\nEPOCH %d:\n", epochs);
        if (epochs % 20 == 0)
        {
            tournment_selection_size = 3 + (rand() % 6);
            printf("*** new tournment_selection_size = %d\n", tournment_selection_size);
        }
        // optimize
        // for (int i = 0; i < size; i++)
        // {
        //     if (rand() % 10 == 0)
        //         individuals[i].fitness += two_opt(inst, individuals[i].chromosome, 20);
        // }

        // Offspring generation
        for (int k = 0; k < children_size; k++)
        {
            general_alg(inst, individuals, tournment_selection_size, size_final, children_size);
        }
        epoch_champion_and_average(inst, individuals, size_final, &champion[epochs - 1], &average[epochs - 1]); // put the champion in individuals[0]
        save_and_plot_solution_general(inst, champion[epochs - 1].chromosome, epochs - 1);

        printf("\t- Champion's fitness* : %f\n", champion[epochs - 1].fitness);
        printf("\t- Average fitness* : %f\n", average[epochs - 1]);

        if (epochs != 1 && champion[epochs - 1].fitness == champion[epochs - 2].fitness)
            no_improvement++;
        else
            no_improvement = 0;

        if (no_improvement > 20)
        {
            print_message("\nNo improvement were found in the last 20 epochs*.");
            break;
        }

        if (clock_gettime(CLOCK_REALTIME, &timestamp) == -1)
            print_error("Error clock_gettime");
        inst->timestamp_finish = timestamp.tv_sec + timestamp.tv_nsec * pow(10, -9);
        // update time limit
        time_left = inst->time_limit - (inst->timestamp_finish - inst->timestamp_start);

        if (time_left > 0.5)
        {
            // OPTIMIZE the champion with one two 2-opt move
            // IF THE CHAMPION CANNOT BE OPTIMIZED ANYMORE (no_improvement > 0)
            // or there is a new champion, try to apply a two-opt move on another individual
            // if (no_improvement == 0)
            // {
            //     individuals[0].fitness += two_opt_v2(inst, individuals[0].chromosome, 1);
            // }
            // else
            // {
            //     printf("*** no_improvement *** \n");
            //     int index = 1 + rand() % (size - 1);
            //     individuals[index].fitness += two_opt_v2(inst, individuals[index].chromosome, 1);
            // }

            epochs++;
            // realloc and alloc
            average = (double *)realloc(average, epochs * sizeof(double));
            champion = (population *)realloc(champion, epochs * sizeof(population));
            champion[epochs - 1].chromosome = (int *)calloc(inst->dimension, sizeof(int));
        }
    }

    printf("\n\nEPOCH SUMMARY*");
    printf("\n%-10s| %-15s| %-15s", "Epoch*", "Average*", "Champion*");
    for (int e = 0; e < epochs; e++)
        printf("\n%-10d| %-15.4f| %-15.4f", e + 1, average[e], champion[e].fitness);

    printf("\nOptimizing the champion...\n");
    individuals[0].fitness += two_opt_v2(inst, individuals[0].chromosome, 0);
    printf("Champion's fitness after optimization: %f\n", individuals[0].fitness);
    return 0;
}

void random_individual(instance *inst, population *individual, int seed, int optimize)
{

    individual->fitness = 0.0;
    individual->chromosome = (int *)calloc(inst->dimension, sizeof(int));

    //if (inst->param.grasp == 1)
    individual->fitness = nearest_neighbours(inst, seed % inst->dimension, individual->chromosome, 1 + rand() % 3);

    //else
    //    individual->fitness = nearest_neighbours(inst, seed % inst->dimension, individual->chromosome, 1);

    // Willing to optimize current chromosome
    // if (optimize == 1)
    // {
    //     individual->fitness += two_opt_v2(inst, individual->chromosome, 5);
    //     printf("\n\t\t- %4d° individual has fitness : %f (OPT) (RI V1)", seed + 1, individual->fitness);
    // }
    // else
    // {
    //     printf("\n\t\t- %4d° individual has fitness : %f (RI V1)", seed + 1, individual->fitness);
    // }
    //if (seed % 2)
    if (optimize == 1)
    {
        printf("fitness(pre two_opt) = %f, ", individual->fitness);
        fflush(stdout);
        individual->fitness += two_opt_v2(inst, individual->chromosome, 30);
        save_and_plot_solution_general(inst, individual->chromosome, seed);
    }
}

void random_individual_2(instance *inst, population *individual, int seed, int optimize)
{

    seed = seed % inst->dimension;

    individual->fitness = 0.0;
    individual->chromosome = (int *)calloc(inst->dimension, sizeof(int));

    // Bucket initialization
    int *bucket = (int *)calloc(inst->dimension, sizeof(int));
    for (int i = 0; i < inst->dimension; i++)
        bucket[i] = i;

    if (seed != inst->dimension - 1)
        bucket[seed] = bucket[inst->dimension - 1];

    int current = seed;

    for (int iter = 1; iter < inst->dimension; iter++)
    {
        int r = rand() % (inst->dimension - iter);

        // Fill chromosome with choice
        individual->chromosome[current] = bucket[r];
        // Reduce bucket by swapping selected element with last
        //        printf("\nBucket [%d] (r -> [%d], bucket[r] -> [%d]): ", iter, r, bucket[r]);
        //        for (int k = 0; k < inst->dimension-iter; k++) printf("%d ", bucket[k]);
        if (r != inst->dimension - iter - 1)
            bucket[r] = bucket[inst->dimension - iter - 1];
        individual->fitness += dist(current, individual->chromosome[current], inst);
        current = individual->chromosome[current];
    }

    individual->chromosome[current] = seed;
    individual->fitness += dist(current, seed, inst);

    //    printf("\nChromosome : ");
    //    for (int k = 0; k < inst->dimension; k++) printf("%d ", individual->chromosome[k]);
    //    printf("\nFitness : %f \n", individual->fitness);

    int maxMoves = 50;

    if (optimize == 1)
    {
        individual->fitness += two_opt(inst, individual->chromosome, maxMoves);
    }
}

void refine_population(instance *inst, population *individuals, int size)
{

    printf("\n\t- Survivor refinement ... ");

    for (int i = 0; i < size; i++)
        individuals[i].fitness += two_opt_v2(inst, individuals[i].chromosome, 5);
}

// void survivor_selection_A(instance *inst, population *individuals, population *offsprings, int individuals_size, int offsprings_size)
// {
//     for (int i = 0; i < offsprings_size; i++)
//     {
//         int index = 1 + (rand() % (individuals_size - 1));
//         //printf("\nindividuals[index].fitness = %f, offsprings[i].fitness = %f\n", individuals[index].fitness, offsprings[i].fitness);                                                          // preserve the champion of the previous epoch
//         int delta = (int)((individuals[index].fitness - offsprings[i].fitness) * 100 / individuals[index].fitness); // delta > 0 => offspring[i] is better than individuals[index]
//         //printf("delta = %d\n", delta);
//         if (delta == 0)
//         {
//             // 50% probabiliy of surviving
//             if (rand() % 2)
//             {
//                 // add the offspring in the population
//                 free(individuals[index].chromosome);
//                 individuals[index].chromosome = offsprings[i].chromosome;
//                 individuals[index].fitness = offsprings[i].fitness;
//             }
//         }
//         else if (delta > 0)
//         {
//             // limit delta to be smaller than 90
//             if (delta > 90)
//                 delta = 90;
//             int r = 1 + ((delta - 1) / 10); // 1 <= r <= 9
//             //printf("mod %d\n", r);
//             if (rand() % (1 + r))
//             {
//                 // add the offspring in the population
//                 free(individuals[index].chromosome);
//                 individuals[index].chromosome = offsprings[i].chromosome;
//                 individuals[index].fitness = offsprings[i].fitness;
//             }
//         }
//         else
//         {
//             // limit delta to be greater than -90
//             if (delta < -90)
//                 delta = -90;
//             int r = 1 + ((-delta - 1) / 10); // 1 <= r <= 9
//             //printf("mod %d\n", r);
//             if ((rand() % (1 + r)) == 0)
//             {
//                 // add the offspring in the population
//                 free(individuals[index].chromosome);
//                 individuals[index].chromosome = offsprings[i].chromosome;
//                 individuals[index].fitness = offsprings[i].fitness;
//             }
//         }
//     }
// }

void survivor_selection_A(instance *inst, population *individuals, population *offsprings, int individuals_size, int offsprings_size)
{
    for (int i = 0; i < offsprings_size; i++)
    {
        int index = 1 + (rand() % (individuals_size - 1)); // preserve the champion of the previous epoch
        if (individuals[index].fitness > offsprings[i].fitness)
        //if (rand() % 2)
        {
            // add the offspring in the population
            free(individuals[index].chromosome);
            individuals[index].chromosome = offsprings[i].chromosome;
            individuals[index].fitness = offsprings[i].fitness;
        }
    }
}

void survivor_selection_B(instance *inst, population *individuals, population *offsprings, int individuals_size, int offsprings_size)
{

    rank(inst, individuals, individuals_size);
    rank(inst, offsprings, offsprings_size);

    //      All individuals are splitted in three classes
    //      - high_ranked : 30% individuals
    //      - mid_ranked  : 60% individuals
    //      - low_ranked  : 10% individuals

    int upper_bound = floor(individuals_size * 0.3);
    int lower_bound = floor(upper_bound + individuals_size * 0.6);

    //    printf("\nlower %d", lower_bound);
    //    printf("\nupper %d", upper_bound);

    int *selected = (int *)calloc(lower_bound - upper_bound, sizeof(int));

    // Offspring chance to survive
    for (int i = 0; i < offsprings_size; i++)
    {

        // Challenging mid-rank individual
        int r = (rand() % (lower_bound - upper_bound)) + upper_bound;
        while (selected[r - upper_bound] == 1)
            r = (rand() % (lower_bound - upper_bound)) + upper_bound;
        selected[r - upper_bound] = 1;

        //        printf("\nr %d", r);
        //        printf("\nr (adjusted) %d", r - upper_bound);
        //        printf("\n");
        //        for (int j = 0; j < lower_bound - upper_bound; ++j) {
        //            printf("%d ", selected[j]);
        //        }

        if (individuals[r].fitness > offsprings[i].fitness)
        {
            if ((rand() % 100) < 75)
            {
                individuals[r].chromosome = individuals[i].chromosome;
                individuals[r].fitness = offsprings[i].fitness;
            }
        }
        else
        {
            if ((rand() % 100) < 25)
            {
                individuals[r].chromosome = individuals[i].chromosome;
                individuals[r].fitness = offsprings[i].fitness;
            }
        }
    }
    free(selected);

    printf("finished\n");
}

// save the champion (and put it on top) and return the average
void epoch_champion_and_average(instance *inst, population *individuals, int size, population *champion, double *average)
{
    double sum = individuals[0].fitness;
    // Find champion among all individual
    champion->fitness = individuals[0].fitness;
    int champion_index = 0;
    int worst_index = 0;
    double worst_fitness = individuals[0].fitness;
    for (int i = 1; i < size; i++)
    {
        sum += individuals[i].fitness;

        // look for the champion
        if (individuals[i].fitness < champion->fitness)
        {
            champion->fitness = individuals[i].fitness;
            champion_index = i;
        }

        // look for the worst individual
        if (individuals[i].fitness > worst_fitness)
        {
            worst_fitness = individuals[i].fitness;
            worst_index = i;
        }
    }

    // PUT THE CHAMPION ON TOP OF THE POPULATION'S LIST
    // swap memory pointers in the array individuals
    int *chromosome_temp = individuals[champion_index].chromosome;
    individuals[champion_index].chromosome = individuals[0].chromosome;
    individuals[0].chromosome = chromosome_temp;
    // swap fitness
    individuals[champion_index].fitness = individuals[0].fitness;
    individuals[0].fitness = champion->fitness;

    // PUT THE WORST INDIVIDUAL IN THE LAST POSITION
    // swap memory pointers in the array individuals
    chromosome_temp = individuals[worst_index].chromosome;
    individuals[worst_index].chromosome = individuals[size - 1].chromosome;
    individuals[size - 1].chromosome = chromosome_temp;
    // swap fitness
    individuals[worst_index].fitness = individuals[size - 1].fitness;
    individuals[size - 1].fitness = worst_fitness;

    // SAVE THE CHAMPION
    // hard copy the champion's chromosome
    for (int i = 0; i < inst->dimension; i++)
    {
        champion->chromosome[i] = individuals[0].chromosome[i];
    }
    champion->fitness = individuals[0].fitness;

    *average = sum / size;
}

void rank(instance *inst, population *individuals, int size)
{

    int i, j, min_idx;

    for (i = 0; i < size - 1; i++)
    {

        // Find individual with minimum fitness
        min_idx = i;
        for (j = i + 1; j < size; j++)
            if (individuals[j].fitness < individuals[min_idx].fitness)
                min_idx = j;

        // Swap individual with minimum fitness with the first individual in the unordered array
        double fitness_best = individuals[min_idx].fitness;
        individuals[min_idx].fitness = individuals[i].fitness;
        individuals[i].fitness = fitness_best;

        int *chromosome_best = individuals[min_idx].chromosome;
        individuals[min_idx].chromosome = individuals[i].chromosome;
        individuals[i].chromosome = chromosome_best;
    }
}

void swap_genes(instance *inst, population *individual, int n)
{
    double delta = 0.0;
    for (int i = 0; i < n; i++)
    {
        // pick 2 random genes
        int a = rand() % inst->dimension;
        int b = rand() % inst->dimension;

        if (a == b)
        {
            b = (b + 1) % inst->dimension;
        }
        // Contiguous node were selected
        if (individual->chromosome[a] == b || individual->chromosome[b] == a)
        {
            return;
        }

        // Make the move
        int c = individual->chromosome[a];
        int d = individual->chromosome[b];

        // Compute the delta
        delta = dist(a, b, inst) + dist(c, d, inst) - dist(a, c, inst) - dist(b, d, inst);

        // Reverse the path from c to b
        if (reverse_successors(individual->chromosome, inst->dimension, c, b))
            print_error("Error in reverse_successors");

        individual->chromosome[a] = b;
        individual->chromosome[c] = d;

        individual->fitness += delta;
    }
}

void roulette_wheel_selection(population *individuals, int size, int *selection)
{

    double t_sum = 0.0;

    for (int i = 0; i < size; i++)
        t_sum += individuals[i].fitness;

    //    printf("\nTotal %f\n", t_sum);

    // Select two parent
    for (int parent = 0; parent < 2; parent++)
    {

        double p_sum = 0.0;
        double point = ((double)rand() / (RAND_MAX));

        //        printf("Wheel point : %f\n", point);

        for (int i = 0; i < size; ++i)
        {

            // Compute pies wheel value (remember lower fitness means better solution)
            p_sum += (individuals[i].fitness / t_sum);
            //            printf("%d° partial sum : %f\n", i, p_sum);

            if (p_sum > point)
            {
                selection[parent] = i;
                break;
            }
        }
    }
}

// selection: array of 2 elements to be filled with two parents
void general_alg(instance *inst, population *individuals, int k, int size, int children_size)
{
    // parents selection
    int parent[2];
    for (int p = 0; p < 2; p++)
    {
        int index = rand() % size;
        // try another index if the second parent has the same fitness of the first
        if (p == 1 && individuals[parent[0]].fitness == individuals[index].fitness)
        {
            index = (index + 1) % size;
        }
        parent[p] = index;
        double winner = individuals[index].fitness;

        // Add k-1 partecipants randomly in the tournament
        for (int i = 0; i < k - 1; i++)
        {
            index = rand() % size;
            // select the best one
            if (individuals[index].fitness < winner)
            {
                // avoid to select two parents with the same fitness
                if (p == 1 && individuals[parent[0]].fitness == individuals[index].fitness)
                    continue;
                winner = individuals[index].fitness;
                parent[p] = index;
            }
        }
    }

    // ************* //
    // child generation
    population child;
    child.chromosome = (int *)calloc(inst->dimension, sizeof(int));
    for (int j = 0; j < inst->dimension; j++)
        child.chromosome[j] = -1;
    int genes_parent0 = one_point_crossover(inst, &individuals[parent[0]], &individuals[parent[1]], &child);

    // substitute the child with parent[index] if it has a better fitness
    int index = 0;
    if (genes_parent0 < (int)inst->dimension / 2)
        index = 1;
    if (child.fitness < individuals[parent[index]].fitness)
    {
        free(individuals[parent[index]].chromosome);
        individuals[parent[index]].chromosome = child.chromosome;
        individuals[parent[index]].fitness = child.fitness;
    }
    if (rand() % children_size == 0)
    {
        int r = rand() % 2;
        printf("*** optimizing the individual %d with 10 two_opt(v2) moves : %f", parent[r], individuals[parent[r]].fitness);
        fflush(stdout);
        individuals[parent[r]].fitness += two_opt_v2(inst, individuals[parent[r]].chromosome, 10);
        printf(" -> %f\n", individuals[parent[r]].fitness);
    }
}

// selection: array of 2 elements to be filled with two parents
void tournament_selection(population *individuals, int k, int size, int *selection)
{

    for (int parent = 0; parent < 2; parent++)
    {
        double *partecipant = (double *)calloc(k, sizeof(double));
        int *ids = (int *)calloc(k, sizeof(int));
        double winner = DBL_MAX;
        // Add k partecipants randomly in the tournament (not the champion (individuals[0]))
        for (int i = 0; i < k; ++i)
        {
            ids[i] = rand() % size;
            partecipant[i] = individuals[ids[i]].fitness;
            // Declare the tournament winner
            if (partecipant[i] < winner)
            {
                if (parent == 1 && selection[0] == ids[i]) // avoid the same parent
                    continue;
                winner = partecipant[i];
                selection[parent] = ids[i];
            }
        }
        free(ids);
        free(partecipant);
    }
}

void rank_selection(instance *inst, population *individuals, int size, int *selection)
{

    rank(inst, individuals, size);

    // Select two parent
    for (int parent = 0; parent < 2; parent++)
    {

        int i = rand() % size;

        selection[parent] = i;
    }
}

void random_selection(population *individuals, int size, int *selection)
{

    // Select first parent
    selection[0] = rand() % size;

    // Select second parent
    selection[1] = selection[0];
    while (selection[0] == selection[1])
        selection[1] = rand() % size;
}

// return crossover_point (#genes of parentA copied)
int one_point_crossover(instance *inst, population *parentA, population *parentB, population *offspring)
{

    // Initialization
    int *selected = (int *)calloc(inst->dimension, sizeof(int));
    int selected_edges_count = 0;
    for (int k = 0; k < inst->dimension; k++)
        offspring->chromosome[k] = -1;

    // Select a crossover point w.r.t. parent fitness
    //    int crossover_point = floor(inst->dimension / (parentA.fitness + parentB.fitness) * parentA.fitness);
    int crossover_point = (int)inst->dimension * 0.2 + (rand() % (int)inst->dimension * 0.6); // crossover_point between 1 and (inst->dimension - 1)
    // int crossover_point = (int)(inst->dimension / 2); // crossover_point between 1 and (inst->dimension - 1)
    // int sign = rand() % 2;
    // if (sign == 0)
    //     sign = -1;
    // crossover_point += sign * ((inst->dimension / 8) + rand() % (int)(inst->dimension / 4));
    int starting_point = rand() % inst->dimension;

    //    printf("\nParent A\n");
    //    for (int k = 0; k < inst->dimension; k++)
    //        printf("%d ", parentA->chromosome[k]);
    //    printf("\n------------------------------------------------\n");
    //    printf("Parent B\n");
    //    for (int k = 0; k < inst->dimension; k++)
    //        printf("%d ", parentB->chromosome[k]);
    //    printf("\n------------------------------------------------\n");
    //    printf("\ncrossover %d \nstarting %d", crossover_point, starting_point);

    // Offspring generation

    int current = starting_point;
    int next = -1;

    selected[starting_point] = 1;

    offspring->fitness = 0.0;
    // First parent sub-path
    for (int i = 0; i < crossover_point; i++)
    {
        // DEBUG
        /*
        printf("\nOffspring : ");
        for (int k = 0; k < inst->dimension; k++)
            printf("%d ", offspring.chromosome[k]);
        printf("\nSelected : ");
        for (int k = 0; k < inst->dimension; k++)
            printf("%d ", selected[k]);
        */
        next = parentA->chromosome[current];
        selected_edges_count++;
        // DEBUG
        //printf("\n[one_point_crossover]next %d, current %d", next, current);

        offspring->chromosome[current] = next;
        offspring->fitness += dist(current, next, inst);

        // Mark node as selected
        selected[next] = 1;

        // Move the current to the next index
        current = next;
    }

    //    printf("\n-------------------------------------------");

    next = parentB->chromosome[current];
    //printf("\n[one_point_crossover] current %d, next, %d", current, next);
    int counter = 0;
    while (next != starting_point)
    {
        if (counter++ > inst->dimension)
            print_error("error in while of one_point_crossover\n");
        // DEBUG

        // printf("\nOffspring : ");
        // for (int k = 0; k < inst->dimension; k++)
        //     printf("%d ", offspring->chromosome[k]);

        // printf("\nSelected : ");
        // for (int k = 0; k < inst->dimension; k++)
        //     printf("%d ", selected[k]);

        if (selected[next] == 0)
        {
            // select the node "next" as the successor of the node "current"
            offspring->chromosome[current] = next;
            offspring->fitness += dist(current, next, inst);
            selected[next] = 1;
            current = next;
            next = parentB->chromosome[current];
            selected_edges_count++;
        }
        else
        {
            // advance the "next" pointer
            //            printf(" node %d already selected", next);
            next = parentB->chromosome[next];
        }
        // DEBUG
        //        printf("\ncurrent %d, next %d", current, next);
    }

    // Closing the circuit
    offspring->chromosome[current] = starting_point;
    offspring->fitness += dist(current, starting_point, inst);

    //    printf("\nChild (partial)\n");
    //    for (int k = 0; k < inst->dimension; k++)
    //        printf("%d ", offspring.chromosome[k]);

    // ADD IN THE CIRCUIT THE REMAINING NODES
    for (int j = 0; j < inst->dimension; j++)
    {
        if (selected[j] == 0)
        {
            offspring->fitness += add_node_extra_mileage(inst, offspring->chromosome, j);
        }
    }

    //    printf("\nChild\n");
    //    for (int k = 0; k < inst->dimension; k++)
    //        printf("%d ", offspring.chromosome[k]);
    return crossover_point;
}

void epoch_average_fitness(population *individuals, double *average, int size)
{

    for (int i = 0; i < size; i++)
        *average += individuals[i].fitness / size;
}

double epoch_percent_deviation(population *individuals, int size)
{

    double average = 0.0;
    double deviation = 0.0;
    double percent_deviation = 0.0;

    epoch_average_fitness(individuals, &average, size);
    //    printf("\nAverage %f", average);

    for (int i = 0; i < size; ++i)
        deviation += fabs(individuals[i].fitness - average);
    deviation = deviation / size;

    percent_deviation = deviation / average * 100;
    //    printf("\nDeviation %f", percent_deviation);
    return percent_deviation;
}