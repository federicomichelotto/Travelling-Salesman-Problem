int meta_heuristic_solver(instance *inst)
{
    // get timestamp
    getTimeStamp(&inst->timestamp_start);

    double min_obj = DBL_MAX;
    double obj_i = 0;
    double obj_opt = 0;

    printf("\nSOLUTION -----------------------------------------------\n");
    printf("\nRUNNING : %s\n", meta_heuristic_model_full_name[inst->model_type]);

    switch (inst->model_type)
    {
    case 0: // Tabu Search
        min_obj = nearest_neighbours(inst, 0, inst->succ, inst->param.grasp_choices);

        save_and_plot_solution(inst, 1);
        inst->z_best = min_obj;
        tabu_search(inst);

        save_and_plot_solution(inst, -1);

        printf("Best objective value: %f\n", min_obj);
        printf("Best objective value (optimized by tabu search): %f\n", inst->z_best);
        break;

    case 1: // Genetic V1
        genetic(inst);
        save_and_plot_solution(inst, -1);

        printf("Best objective value (optimized by genetic algorithm): %f\n", inst->z_best);

        break;
    case 2: // Genetic V2
        genetic_v2(inst);
        break;

    default:
        fprintf(stderr, "ERROR: Model type %d not available.\n", inst->model_type);
        break;
    }

    // get timestamp
    getTimeStamp(&inst->timestamp_finish);

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

    // update time_left
    double ts_current;
    getTimeStamp(&ts_current);
    double time_left = inst->time_limit - (ts_current - inst->timestamp_start);

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

        // update time_left
        getTimeStamp(&ts_current);
        time_left = inst->time_limit - (ts_current - inst->timestamp_start);
    }
    free(temp_succ);
    free(tabu_list);
}

int genetic(instance *inst)
{

    int size = inst->param.pop_size;
    int children_size = inst->param.off_size;
    //    printf("How many individual should have your population: ");
    //    scanf("%d", &size);
    //
    //    printf("How many children should conceive your population each epoch: ");
    //    scanf("%d", &children_size);

    population *individuals = (population *)calloc(size, sizeof(population));

    printf("\nGENERATION STEP");
    printf("\n\t%d individuals will be generated (# = %d individuals generated)\n", size, (int)ceil(size / 10));

    // Population
    for (int i = 0; i < size; i++)
    {

        // if (i % (int)ceil(size / 10) == 0)
        // {
        //     printf("#");
        //     fflush(stdout);
        // }

        // Generate random individuals
        printf("\n\tINDIVIDUAL #%5d FITNESS -> ", i + 1);

        if (i < size * 0.15)
            random_individual(inst, &individuals[i], i, 1);
        else if (i >= size * 0.15 && i < size * 0.5)
            random_individual(inst, &individuals[i], i, 0);
        else
            random_individual_2(inst, &individuals[i], i, 0);
    }

    // update time_left
    double ts_current;
    getTimeStamp(&ts_current);
    double time_left = inst->time_limit - (ts_current - inst->timestamp_start);

    printf("\nEVOLUTION STEP\n");

    int no_improvement = 0;
    int epochs = 1;

    // Epochs
    double *average = (double *)calloc(1, sizeof(double));

    population *champion = (population *)calloc(1, sizeof(population));
    champion[epochs - 1].chromosome = (int *)calloc(inst->dimension, sizeof(int));

    // Since there's time left or fitness isn't too stable
    while (time_left > 0.5)
    {

        printf("\nEPOCH #%d [%.2f sec]\n", epochs - 1, time_left);

        population *offsprings = (population *)calloc(children_size, sizeof(population));
        for (int i = 0; i < children_size; i++)
        {
            offsprings[i].chromosome = (int *)calloc(inst->dimension, sizeof(int));
            for (int j = 0; j < inst->dimension; j++)
                offsprings[i].chromosome[j] = -1;
        }

        // Offspring generation
        printf("\t- Generation of %d offspring ... \n", children_size);
        for (int k = 0; k < children_size; k++)
        {
            // Parent selection
            int *parent = (int *)calloc(2, sizeof(int)); // array where the two parents are going to be stored

            switch (inst->param.par_sel) {
                case 0:
                    tournament_selection(individuals, 3, size, parent);
                    break;
                case 1:
                    roulette_wheel_selection(individuals,size,parent);
                    break;
                case 2:
                    rank_selection(inst, individuals, size, parent);
                    break;
                case 3:
                    random_selection(individuals, size, parent);
                    break;
                default:
                    print_error("Parent selection algorithm not found");
            }

            //            printf("\t\t+ Crossover #%d between parent %d and %d\n", k + 1, parent[0], parent[1]);
            one_point_crossover(inst, individuals, &offsprings[k], parent[0], parent[1], 1);

            free(parent);
        }

        printf("\n\t- Survivor selection ... \n");

        switch (inst->param.sur_sel) {
            case 0:
                survivor_selection_A(inst, individuals, offsprings, size, children_size);
                break;
            case 1:
                 survivor_selection_B(inst, individuals, offsprings, size, children_size);
                break;
            default:
                print_error("Survival selection algorithm not found");
        }

        // Mutate population
        // How much population is stressed
        double stress = (rand() % (60 - 20 + 1)) + 20;
        printf("\n\t- Mutation occur on %4.2f%% of individuals ... \n", stress);
        int mutation_size = floor(size * stress / 100);

        for (int k = 0; k < mutation_size; k++)
        {
            int m = 1 + rand() % (size - 1); // pick randomly an individual (not the champion)

            //            printf("\t\t- Individual %d :\n", m );
            swap_genes(inst, &individuals[m], 1);
        }

        if (epochs % 50 == 0)
        {
            printf("\n\t- Random individual full optimization ... \n");
            int m = 1 + rand() % (size - 1); // pick randomly an individual (not the champion)
            printf("individuals[m].fitness = %f \n", individuals[m].fitness);
            individuals[m].fitness += two_opt_v2(inst, individuals[m].chromosome, 0);
            printf("individuals[m].fitness = %f \n", individuals[m].fitness);
            no_improvement = 0;
        }

        printf("\n\t- Summary .. \n");
        epoch_champion_and_average(inst, individuals, size, &champion[epochs - 1], &average[epochs - 1]); // put the champion in individuals[0]

        // Print champion solution inside current epoch
        save_and_plot_solution_general(inst, champion[epochs - 1].chromosome, epochs - 1);

        printf("\t\t- Champion's fitness : %f\n", champion[epochs - 1].fitness);
        printf("\t\t- Average fitness : %f\n", average[epochs - 1]);

        if (epochs > 1 && champion[epochs - 1].fitness == champion[epochs - 2].fitness)
            no_improvement++;
        else
            no_improvement = 0;

        if (no_improvement == 20)
        {
            printf("\tWARNING : No improvement were found in the last %d epochs.\n", no_improvement);
            printf("\tStarting quick optimization... ");
            fflush(stdout);
            fast_population_refinement(inst, individuals, size, 20);
            printf("\tComplete\n");

            no_improvement = 0;
        }

        if (epochs > 1000)
            break;

        // update time_left
        ts_current;
        getTimeStamp(&ts_current);
        time_left = inst->time_limit - (ts_current - inst->timestamp_start);
        if (time_left > 0.5)
        {
            // A new epoch it's about to start
            epochs++;
            // realloc and alloc
            average = (double *)realloc(average, epochs * sizeof(double));
            champion = (population *)realloc(champion, epochs * sizeof(population));
            champion[epochs - 1].chromosome = (int *)calloc(inst->dimension, sizeof(int));
        }

        // Free memory
        free(offsprings);
    }

    printf("\n\t- Improving champion solution\n");
    champion[epochs - 1].fitness += two_opt_v2(inst, champion[epochs - 1].chromosome, 0);
    printf("\t\t- Fitness : %f\n", champion[epochs - 1].fitness);

    printf("\n\nEPOCH SUMMARY");
    printf("\n%-10s| %-15s| %-15s", "Epoch", "Average", "Champion");
    for (int e = 0; e < epochs; e++)
        printf("\n%-10d| %-15.4f| %-15.4f", e, average[e], champion[e].fitness);

    inst->z_best = champion[epochs - 1].fitness;

    free(average);
    free(champion);
    free(individuals);

    return 0;
}

int genetic_v2(instance *inst)
{

    int size = inst->param.pop_size;
    int children_size = inst->param.off_size;
    //    printf("How many individual should have your population: ");
    //    scanf("%d", &size);
    //
    //    printf("How many children should conceive your population each epoch: ");
    //    scanf("%d", &children_size);

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

            one_point_crossover(inst, individuals, &champion_child, 0, floor(1 + rand() % (size - 2)), 0);
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

    individual->fitness = nearest_neighbours(inst, seed % inst->dimension, individual->chromosome, inst->param.grasp_choices);

    printf("BASE [%f]", individual->fitness);
    fflush(stdout);
    if (optimize == 1)
    {
        individual->fitness += two_opt_v2(inst, individual->chromosome, 30);
        printf(" | OPTIMIZED [%f]", individual->fitness);
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

    if (optimize == 1)
    {
        individual->fitness += two_opt(inst, individual->chromosome, 50);
    }
}

void refine_population(instance *inst, population *individuals, int size)
{

    printf("\n\t- Survivor refinement ... ");

    for (int i = 0; i < size; i++)
        individuals[i].fitness += two_opt_v2(inst, individuals[i].chromosome, 5);
}

void fast_population_refinement(instance *inst, population *individuals, int size, int moves)
{

    for (int i = 0; i < size; i++)
        individuals[i].fitness += two_opt(inst, individuals[i].chromosome, moves);
}

void deep_population_refinement(instance *inst, population *individuals, int size, int moves)
{

    for (int i = 0; i < size; i++)
        individuals[i].fitness += two_opt_v2(inst, individuals[i].chromosome, moves);
}

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

        if (offsprings[i].fitness < individuals[r].fitness)
        {
            if ((rand() % 100) < 75)
            {
                free(individuals[r].chromosome);
                individuals[r].chromosome = offsprings[i].chromosome;
                individuals[r].fitness = offsprings[i].fitness;
            }
        }
        else
        {
            if ((rand() % 100) < 25)
            {
                free(individuals[r].chromosome);
                individuals[r].chromosome = offsprings[i].chromosome;
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
    int genes_parent0 = one_point_crossover(inst, individuals, &child, parent[0], parent[1], 0);

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

void roulette_wheel_selection(population *individuals, int size, int *selection)
{

    double t_sum = 0.0;

    // Compute total sum
    for (int i = 0; i < size; i++)
        t_sum += individuals[i].fitness;

    // Select two parent
    for (int i = 0; i < 2; i++)
    {

        double p_sum = 0.0;
        double pointer = ((double)rand() / (RAND_MAX));

        for (int j = 0; j < size; j++)
        {

            // Compute pies wheel value (remember lower fitness means better solution)
            p_sum += (1.0 - (individuals[j].fitness / t_sum)) / (size - 1);

            if (p_sum > pointer)
            {

                selection[i] = j;
                break;
            }
        }
    }
}

void rank_selection(instance *inst, population *individuals, int size, int *selection)
{
    // Need to rank all individuals
    rank(inst, individuals, size);

    int *c_sum = (int *)calloc(size, sizeof(int));

    c_sum[0] = 1;
    for (int i = 1; i < size; i++)
        c_sum[i] = c_sum[i - 1] + i + 1;

    //    printf("\n\nRANK SUMMARY");
    //    printf("\n%-10s| %-15s| %-15s", "IND", "Fitness", "Rank");
    //    for (int e = 0; e < size; e++)
    //        printf("\n%-10d| %-15.4f| %-15.4d", e, individuals[e].fitness, c_sum[e]);

    // Select two parent
    for (int i = 0; i < 2; i++)
    {

        int pointer = 1 + (rand() % size * (size + 1) / 2);

        // Apply 2-nd order equation solution
        int index = (-1 + sqrt(1 + 8 * pointer)) / 2.0;

        //        printf("\nRandom value : %d", random_value);
        //        printf("\nSelected index %d -> value %d", index, c_sum[index]);
        selection[i] = index;
    }

    free(c_sum);
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

int one_point_crossover(instance *inst, population *parent, population *offspring, int A, int B, int weighted)
{

    // Initialization
    int *selected = (int *)calloc(inst->dimension, sizeof(int));
    for (int k = 0; k < inst->dimension; k++)
        offspring->chromosome[k] = -1;

    int crossover_point = -1;

    if (weighted)
    {

        // Select a crossover point w.r.t. parent fitness

        if (parent[A].fitness > parent[B].fitness)
        {
            int temp = A;
            A = B;
            B = temp;
        }

        crossover_point = floor(inst->dimension / (parent[A].fitness + parent[B].fitness) * parent[A].fitness);
    }
    else
    {

        // Select a crossover point randomly between 1 and (inst->dimension - 1)

        crossover_point = 1 + (rand() % (inst->dimension - 1));
        // crossover_point = (int)inst->dimension * 0.2 + (rand() % (int)inst->dimension * 0.6); // crossover_point between 1 and (inst->dimension - 1)
        // crossover_point = (int)(inst->dimension / 2); // crossover_point between 1 and (inst->dimension - 1)
        // int sign = rand() % 2;
        // if (sign == 0)
        //     sign = -1;
        // crossover_point += sign * ((inst->dimension / 8) + rand() % (int)(inst->dimension / 4));
    }

    int starting_point = rand() % inst->dimension;

    //        printf("\nParent A\n");
    //        for (int k = 0; k < inst->dimension; k++)
    //            printf("%d ", parent[A].chromosome[k]);
    //        printf("\n------------------------------------------------\n");
    //        printf("Parent B\n");
    //        for (int k = 0; k < inst->dimension; k++)
    //            printf("%d ", parent[B].chromosome[k]);
    //        printf("\n------------------------------------------------\n");
    //        printf("\ncrossover %d \nstarting %d", crossover_point, starting_point);

    // Offspring generation

    int current = starting_point;
    int next = -1;

    selected[starting_point] = 1;
    offspring->fitness = 0.0;

    // First parent sub-path
    for (int i = 0; i < crossover_point; i++)
    {
        // DEBUG

        //        printf("\nOffspring : ");
        //        for (int k = 0; k < inst->dimension; k++)
        //            printf("%d ", offspring->chromosome[k]);
        //        printf("\nSelected : ");
        //        for (int k = 0; k < inst->dimension; k++)
        //        printf("%d ", selected[k]);

        next = parent[A].chromosome[current];
        // DEBUG
        // printf("\nnext %d, temp %d, current %d", next, temp, current);

        offspring->chromosome[current] = next;
        offspring->fitness += dist(current, next, inst);

        // Mark node as selected
        selected[next] = 1;

        // Move the current to the next index
        current = next;
    }

    //    printf("\n-------------------------------------------");

    next = parent[B].chromosome[current];
    //    printf("\ncurrent %d, next, %d", current, next);

    while (next != starting_point)
    {
        // DEBUG

        //        printf("\nOffspring : ");
        //        for (int k = 0; k < inst->dimension; k++)
        //            printf("%d ", offspring->chromosome[k]);
        //        printf("\nSelected : ");
        //        for (int k = 0; k < inst->dimension; k++)
        //            printf("%d ", selected[k]);

        if (selected[next] == 0)
        {
            // select the node "next" as the successor of the node "current"
            offspring->chromosome[current] = next;
            offspring->fitness += dist(current, next, inst);
            selected[next] = 1;
            current = next;
            next = parent[B].chromosome[current];
        }
        else
        {
            // advance the "next" pointer
            //            printf(" node %d already selected", next);
            next = parent[B].chromosome[next];
        }
        // DEBUG
        //        printf("\ncurrent %d, next %d", current, next);
    }

    // Closing the circuit
    offspring->chromosome[current] = starting_point;
    offspring->fitness += dist(current, starting_point, inst);

    //        printf("\nChild (partial)\n");
    //        for (int k = 0; k < inst->dimension; k++)
    //            printf("%d ", offspring->chromosome[k]);

    // ADD IN THE CIRCUIT THE REMAINING NODES
    for (int j = 0; j < inst->dimension; j++)
    {
        if (selected[j] == 0)
        {
            offspring->fitness += add_node_extra_mileage(inst, offspring->chromosome, j);
        }
    }

    //        printf("\nChild\n");
    //        for (int k = 0; k < inst->dimension; k++)
    //            printf("%d ", offspring->chromosome[k]);

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