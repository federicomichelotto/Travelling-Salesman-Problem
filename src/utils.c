#include "../include/utils.h"
#include <sys/stat.h>

void parse_command_line(int argc, char **argv, instance *inst)
{

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s -h for help\n", argv[0]);
        exit(1);
    }
    else if (argc == 2 && strcmp(argv[1], "-h") == 0)
    {
        print_help();
    }
    else
    {

        for (int i = 1; i < argc; i++)
        {

            if (strcmp(argv[i], "-f") == 0)
            { // Input file
                strcpy(inst->param.input_file, argv[++i]);
                continue;
            }

            if (strcmp(argv[i], "-t") == 0)
            { // Time limit (seconds)
                inst->time_limit = strtof(argv[++i], NULL);
                continue;
            }

            if (strcmp(argv[i], "-m") == 0)
            { // Model type
                inst->model_type = strtol(argv[++i], NULL, 10);
                continue;
            }

            if (strcmp(argv[i], "-s") == 0)
            { // seed
                inst->param.seed = atoi(argv[++i]);
                continue;
            }

            if (strcmp(argv[i], "-r") == 0)
            { // run
                inst->param.run = atoi(argv[++i]);
                continue;
            }

            if (strcmp(argv[i], "-v") == 0)
            { // Verbosity
                switch (strtol(argv[++i], NULL, 10))
                {
                case QUIET:
                    inst->param.verbose = QUIET;
                    break;
                case NORMAL:
                    inst->param.verbose = NORMAL;
                    break;
                case VERBOSE:
                    inst->param.verbose = VERBOSE;
                    break;
                case NERD:
                    inst->param.verbose = NERD;
                    break;
                case DEBUG:
                    inst->param.verbose = DEBUG;
                    break;
                default:
                    inst->param.verbose = NORMAL;
                    break;
                }
                continue;
            }

            if (strcmp(argv[i], "--ticks") == 0)
            {
                inst->param.ticks = 1;
                continue;
            }
            if (strcmp(argv[i], "--interactive") == 0)
            {
                inst->param.interactive = 1;
                continue;
            }

            if (strcmp(argv[i], "--saveplots") == 0)
            {
                inst->param.saveplots = 1;
                continue;
            }

            if (strcmp(argv[i], "--math") == 0)
            {
                inst->param.solver = 1;
                continue;
            }

            if (strcmp(argv[i], "--heur") == 0)
            {
                inst->param.solver = 2;
                continue;
            }

            if (strcmp(argv[i], "--meta") == 0)
            {
                inst->param.solver = 3;
                continue;
            }

            if (strcmp(argv[i], "--popsize") == 0)
            { // population
                inst->param.pop_size = atoi(argv[++i]);
                continue;
            }

            if (strcmp(argv[i], "--offsize") == 0)
            { // population
                inst->param.off_size = atoi(argv[++i]);
                continue;
            }

            if (strcmp(argv[i], "--grasp") == 0)
            {
                inst->param.grasp = 1;
                inst->param.grasp_choices = atoi(argv[++i]);
                continue;
            }
        }

        print_command_line(inst);
    }
}

void parse_instance(instance *inst)
{

    FILE *fp = NULL;
    char *filename = inst->param.input_file;

    fp = fopen(filename, "r");
    if (fp == NULL)
        print_error("Could not open the file");

    char line[180];
    int section;

    // Read each line from the input file:
    // - compare the first token of the line to check the section
    // - split the line using strtok and consider the generated tokens
    while (fgets(line, sizeof(line), fp) != NULL)
    {

        if (inst->param.verbose > DEBUG)
            printf("%s\n", line);
        line[strcspn(line, "\n")] = 0; // removing trailing \n

        char delimiters[] = " :\n\r\t";
        char *parameter;
        char *value;

        parameter = strtok(line, delimiters);

        // Check the current section
        if (strncmp(parameter, "NODE_COORD_SECTION", 18) == 0)
        {
            section = NODE_COORD;
        }
        else if (strncmp(parameter, "EDGE_WEIGHT_SECTION", 19) == 0)
        {
            section = EDGE_WEIGHT;
        }
        else if (strncmp(parameter, "DISPLAY_DATA_SECTION", 20) == 0)
        {
            section = DISPLAY_DATA;
        }
        else
        {
            section = PARAMETERS;
        }

        if (strncmp(parameter, "EOF", 3) == 0)
        {
            if (inst->param.verbose > DEBUG)
                printf("(EOF) found, instance parsing complete\n");
            break;
        }

        if (section == PARAMETERS)
        {
            if (strncmp(parameter, "NAME", 4) == 0)
            {
                check_format(inst->param.name);
                value = strtok(NULL, delimiters);
                strcpy(inst->param.name, value);
                if (inst->param.verbose > DEBUG)
                    printf("NAME %s\n\n", inst->param.name);
            }
            else if (strncmp(parameter, "COMMENT", 7) == 0)
            {
                value = strtok(NULL, ":");
                if (strncmp(inst->param.comment, "NULL", 4) != 0)
                    strcat(inst->param.comment, value); // if more than one comment, append
                else
                    strcpy(inst->param.comment, value);
                if (inst->param.verbose > DEBUG)
                    printf("Solving instance %s with model %d\n\n", inst->param.comment, inst->model_type);
            }
            else if (strncmp(parameter, "TYPE", 4) == 0)
            {
                check_format(inst->param.type);
                value = strtok(NULL, delimiters);
                if (strncmp(value, "TSP", 3) != 0 && strncmp(value, "ATSP", 4) != 0)
                    print_error("(TYPE) only TSP and ATSP implemented so far");
                else
                {
                    strcpy(inst->param.type, value);
                    if (inst->param.verbose > DEBUG)
                        printf("TYPE %s\n\n", inst->param.type);
                }
            }
            else if (strncmp(parameter, "DIMENSION", 9) == 0)
            {
                if (inst->dimension != -1)
                    print_error("Bad input format");
                value = strtok(NULL, delimiters);
                inst->dimension = strtol(value, NULL, 10);
                if (inst->param.verbose > DEBUG)
                    printf("NODES %d\n", inst->dimension);
                inst->succ = (int *)calloc(inst->dimension, sizeof(int));
                inst->nodes = (node *)calloc(inst->dimension, sizeof(node));
            }
            else if (strncmp(parameter, "EDGE_WEIGHT_TYPE", 16) == 0)
            {
                check_format(inst->param.weight_type);
                value = strtok(NULL, delimiters);
                strcpy(inst->param.weight_type, value);
                if (strncmp(inst->param.type, "TSP", 3) == 0)
                {
                    if (strncmp(inst->param.weight_type, "EXPLICIT", 8) != 0 &&
                        strncmp(inst->param.weight_type, "CEIL_2D", 7) != 0 &&
                        strncmp(inst->param.weight_type, "EUC_2D", 6) != 0 &&
                        strncmp(inst->param.weight_type, "MAN_2D", 6) != 0 &&
                        strncmp(inst->param.weight_type, "MAX_2D", 6) != 0 &&
                        strncmp(inst->param.weight_type, "ATT", 3) != 0 &&
                        strncmp(inst->param.weight_type, "GEO", 3) != 0)
                    {
                        print_error("(EDGE_WEIGHT_TYPE) only EUC_2D, ATT, MAN_2D, MAX_2D, CEIL_2D, GEO implemented so far");
                    }
                }
                else if (strncmp(inst->param.type, "ATSP", 4) == 0)
                {
                    if (strncmp(inst->param.weight_type, "EXPLICIT", 8) != 0)
                    {
                        print_error("(EDGE_WEIGHT_TYPE) only EXPLICIT implemented so far");
                    }
                }
            }
            else if (strncmp(parameter, "EDGE_WEIGHT_FORMAT", 18) == 0)
            {
                check_format(inst->param.weight_format);
                value = strtok(NULL, delimiters);
                strcpy(inst->param.weight_format, value);
                if (strncmp(inst->param.weight_format, "FULL_MATRIX", 11) != 0 &&
                    strncmp(inst->param.weight_format, "FUNCTION", 8) != 0)
                {
                    print_error("(EDGE_WEIGHT_FORMAT) only FULL_MATRIX and FUNCTION implemented so far");
                }
            }
            else if (strncmp(parameter, "DISPLAY_DATA_TYPE", 17) == 0)
            {
                check_format(inst->param.data_type);
                value = strtok(NULL, delimiters);
                strcpy(inst->param.data_type, value);
                if (strncmp(inst->param.data_type, "COORD_DISPLAY", 13) != 0)
                    print_error("(DISPLAY_DATA_TYPE) only COORD_DISPLAY implemented so far");
            }
        }
        else if (section == NODE_COORD)
        {

            for (int i = 0; i < inst->dimension; i++)
            {

                fgets(line, sizeof(line), fp);
                value = strtok(line, delimiters);

                inst->nodes[i].id = strtol(value, NULL, 10);
                inst->nodes[i].x = strtof(strtok(NULL, delimiters), NULL);
                inst->nodes[i].y = strtof(strtok(NULL, delimiters), NULL);
                if (inst->param.verbose > DEBUG)
                    printf("NODES %d at coordinates (%15.7lf , %15.7lf)\n", inst->nodes[i].id, inst->nodes[i].x, inst->nodes[i].y);
            }
        }
        else if (section == EDGE_WEIGHT)
        {
            // TODO implementing the parsing of the EDGE WEIGHT SECTION
        }
        else if (section == DISPLAY_DATA)
        {
            // TODO implementing the parsing of the DISPLAY DATA SECTION
        }
    }

    fclose(fp);
}

void initialize_instance(instance *inst)
{
    srand(time(NULL));     // Initialize random seed for the computation

    inst->model_type = 0;
    inst->time_limit = CPX_INFBOUND;
    inst->param.seed = 1;
    inst->param.run = 0;
    inst->param.verbose = NORMAL;
    inst->param.callback_counter = 0;
    inst->param.ticks = 0;
    inst->param.solver = 0; // default
    inst->param.interactive = 0;
    inst->param.saveplots = 0;
    inst->param.grasp = 0;
    inst->param.grasp_choices = 1;  // default base case 1 possible choice

    // Genetic parameter
    inst->param.pop_size = 100;
    inst->param.off_size = 10;

    inst->dimension = -1;
    inst->nodes = NULL;
    inst->succ = NULL;

    inst->weights = NULL;
    inst->integer_costs = 0;
    inst->z_best = -1.0;

    inst->timestamp_start = 0.0;
    inst->timestamp_finish = 0.0;

    strcpy(inst->param.input_file, "NULL");
    strcpy(inst->param.name, "NULL");
    strcpy(inst->param.type, "NULL");
    strcpy(inst->param.comment, "NULL");
    strcpy(inst->param.weight_type, "NULL");
    strcpy(inst->param.weight_format, "NULL");
    strcpy(inst->param.data_type, "NULL");
    inst->gnuplotPipe = popen("gnuplot -persistent", "w");
}

void check_format(char *param)
{
    if (strncmp(param, "NULL", 4) != 0)
    {
        print_error("Bad input format");
    }
}

void free_instance(instance *inst)
{
    // close gnuplot pipe
    if (pclose(inst->gnuplotPipe) == -1)
        print_error("pclose error");

    if (inst->param.solver != 2 && inst->param.solver != 3)
        free(inst->best_sol);

    free(inst->nodes);
    free(inst->succ);
    free(inst->weights);
    // todo write the code to free the allocated memory within the instance (bottom-up approach)
}

void print_command_line(instance *inst)
{
    printf("\nPARAMETERS ---------------------------------------------\n");
    printf("-f %s\n", inst->param.input_file);
    if (inst->param.ticks == 1)
    {
        printf("-t %f ticks\n", inst->time_limit);
        printf("--ticks (ACTIVE -> %d)\n", inst->param.ticks);
    }
    else
    {
        printf("-t %f seconds\n", inst->time_limit);
    }

    switch (inst->param.solver)
    {
    case 0:
        printf("-m %d (%s)\n", inst->model_type, optimal_model_name[inst->model_type]);
        break;
    case 1:
        printf("-m %d (%s)\n", inst->model_type, math_model_name[inst->model_type]);
        printf("--math (MATH HEURISTICS SOLVER) (ACTIVE -> %d)\n", inst->param.solver);
        break;
    case 2:
        printf("-m %d (%s)\n", inst->model_type, heuristic_model_name[inst->model_type]);
        printf("--heur (HEURISTICS SOLVER) (ACTIVE -> %d)\n", inst->param.solver);
        break;
    case 3:
        printf("-m %d (%s)\n", inst->model_type, meta_heuristic_model_name[inst->model_type]);
        printf("--meta (META HEURISTICS SOLVER) (ACTIVE -> %d)\n", inst->param.solver);
        break;
    default:
        print_error("No implemented solver selected");
    }

    if (inst->param.interactive == 1)
    {
        printf("--interactive (ACTIVE -> %d)\n", inst->param.interactive);
    }

    if (inst->param.grasp == 1)
    {
        printf("--grasp (ACTIVE -> %d)\n", inst->param.grasp);
        printf("--grasp_choices (%d)\n", inst->param.grasp_choices);
    }

    if (inst->param.saveplots == 1)
    {
        printf("--saveplots (ACTIVE -> %d)\n", inst->param.saveplots);
    }

    printf("-s %d\n", inst->param.seed);
    printf("-v %d (%s)\n", inst->param.verbose, verbose_name[inst->param.verbose]);
    printf("--------------------------------------------------------\n\n");
}

void print_instance(instance *inst)
{
    printf("\nINSTANCE -----------------------------------------------\n");
    printf("Name: %s\n", inst->param.name);
    printf("Type: %s\n", inst->param.type);
    printf("Comment: %s\n", inst->param.comment);
    printf("Dimension: %d \n", inst->dimension);
    printf("Edge weight type: %s\n", inst->param.weight_type);
    printf("Edge weight format: %s\n", inst->param.weight_format);

    switch (inst->param.solver)
    {
    case 0:
        printf("Model type: %d (%s)\n", inst->model_type, optimal_model_name[inst->model_type]);
        break;
    case 1:
        printf("Model type: %d (%s)\n", inst->model_type, math_model_name[inst->model_type]);
        break;
    case 2:
        printf("Model type: %d (%s)\n", inst->model_type, heuristic_model_name[inst->model_type]);
        break;
    case 3:
        printf("Model type: %d (%s)\n", inst->model_type, meta_heuristic_model_name[inst->model_type]);
        break;
    default:
        print_error("No implemented solver selected");
    }

    printf("Input file path: %s\n", inst->param.input_file);
    printf("Data type: %s\n", inst->param.data_type);
    // TODO implement the EDGE WEIGHT SECTION
    // TODO implement the DISPLAY DATA SECTION
    if (inst->nodes != NULL)
    {
        printf("NODE COORD SECTION\n");
        for (int i = 0; i < inst->dimension; i++)
            printf("%d\t(%15.7lf, %15.7lf)\n", inst->nodes[i].id, inst->nodes[i].x, inst->nodes[i].y);
    }
    printf("--------------------------------------------------------\n\n");
}

void print_help()
{
    printf("\nHELP ---------------------------------------------------\n");
    printf("-f <path>       : used to pass the relative instance path \n");
    printf("-t <time>       : used to pass the total running time allowed in seconds\n");
    printf("--ticks         : used to set the way time is interpreted inside the solver (optional)\n");
    printf("--interactive   : used to plot all the solutions (by default none plot is displayed)\n");
    printf("--saveplots     : used to save all the solutions' plots (by default only the final solution's plot is saved)\n");
    printf("--math          : used to set the math-heuristic solver (optional)\n");
    printf("--heur          : used to set the heuristic solver (optional)\n");
    printf("--meta          : used to set the meta-heuristic solver (optional)\n");
    printf("--grasp <value> : used to set GRASP approach and possible choices (optional)\n");
    printf("-m <model>      : used to set the model type (based on the solver)\n");
    printf("-s <seed>       : used to set the seed\n");
    printf("-v <value>      : used to set the verbosity, from QUIET (0) up to DEBUG (4)\n");
    printf("--------------------------------------------------------\n\n");
    exit(1);
}

void print_error(const char *err)
{
    fprintf(stderr, "\nERROR: %s \n\n", err);
    fflush(NULL);
    exit(1);
}

void print_error_status(const char *err, int e)
{
    printf("\n\n ERROR: exit with status %d. %s. \n\n", e, err);
    fflush(NULL);
    exit(1);
}

void print_message(const char *msg)
{
    fprintf(stdout, "\nMESSAGE: %s \n\n", msg);
    fflush(NULL);
}

int save_and_plot_solution(instance *inst, int iter)
{

    if (inst->param.saveplots || inst->param.interactive || iter == -1)
    {
        // write solution to file
        FILE *temp = fopen("data.temp", "w");
        for (int i = 0; i < inst->dimension; i++)
        {
            //Write the coordinates of the two nodes inside a temporary file
            fprintf(temp, "%lf %lf \n%lf %lf \n\n", inst->nodes[i].x, inst->nodes[i].y, inst->nodes[inst->succ[i]].x, inst->nodes[inst->succ[i]].y);
            // double '\n' to create a new edge (block of coordinates)
        }
        fclose(temp);

        if (inst->param.saveplots || iter == -1) // save plot
        {
            FILE *gnuplotPipe = popen("gnuplot", "w"); // local gnuplotPipe to avoid conflicts with the global gnuplotPipe
            char *out = (char *)calloc(1000, sizeof(char));
            char *plot = (char *)calloc(1000, sizeof(char));
            char *model_name = (char *)calloc(100, sizeof(char));

            switch (inst->param.solver)
            {
            case 0:
                strcpy(model_name, optimal_model_name[inst->model_type]);
                break;
            case 1:
                strcpy(model_name, math_model_name[inst->model_type]);
                break;
            case 2:
                strcpy(model_name, heuristic_model_name[inst->model_type]);
                break;
            }

            sprintf(out, "set output '../output/plot/%s_%s_%d.jpg'", inst->param.name, model_name, iter);
            sprintf(plot, "plot 'data.temp' with linespoints pt 7 lc rgb 'blue' lw 1 notitle");
            if (iter == -1) // final solution
            {
                sprintf(out, "set output '../output/plot/%s_%s_final.jpg'", inst->param.name, model_name);
                sprintf(plot, "plot 'data.temp' with linespoints pt 7 lc rgb 'red' lw 1 notitle");
            }

            char *commandsForGnuplot[] = {
                "set title 'Solution plot'",
                "set terminal jpeg size 1024,768",
                out,
                "unset key",
                "set autoscale",
                "set ylabel 'Y'",
                "set xlabel 'X'",
                plot,
                "clear"};

            int commands = sizeof(commandsForGnuplot) / sizeof(commandsForGnuplot[0]);

            //Send commands to gnuplot one by one.
            for (int i = 0; i < commands; i++)
                fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);

            pclose(gnuplotPipe); // execute the commands
            free(plot);
            free(out);
            free(model_name);
        }

        if (inst->param.interactive) // plot solution
        {
            // here we use a global gnuplotPipe to remain on the same window
            char *commandsForGnuplot[] = {"set title 'Solution plot'",
                                          "unset key",
                                          "set autoscale",
                                          "set ylabel 'Y'",
                                          "set xlabel 'X'",
                                          "plot 'data.temp' with linespoints pt 7 lc rgb 'blue' lw 1"};

            int commands = sizeof(commandsForGnuplot) / sizeof(commandsForGnuplot[0]);
            for (int i = 0; i < commands; i++)
                fprintf(inst->gnuplotPipe, "%s \n", commandsForGnuplot[i]);
            fflush(inst->gnuplotPipe); // execute the commands
        }
    }
    return 0;
}

int save_and_plot_solution_general(instance *inst, int *succ, int iter)
{

    if (inst->param.saveplots || inst->param.interactive || iter == -1)
    {
        // write solution to file
        FILE *temp = fopen("data.temp", "w");
        for (int i = 0; i < inst->dimension; i++)
        {
            //Write the coordinates of the two nodes inside a temporary file
            fprintf(temp, "%lf %lf \n%lf %lf \n\n", inst->nodes[i].x, inst->nodes[i].y, inst->nodes[succ[i]].x, inst->nodes[succ[i]].y);
            // double '\n' to create a new edge (block of coordinates)
        }
        fclose(temp);

        if (inst->param.saveplots || iter == -1) // save plot
        {
            FILE *gnuplotPipe = popen("gnuplot", "w"); // local gnuplotPipe to avoid conflicts with the global gnuplotPipe
            char *out = (char *)calloc(1000, sizeof(char));
            char *plot = (char *)calloc(1000, sizeof(char));
            char *model_name = (char *)calloc(100, sizeof(char));

            switch (inst->param.solver)
            {
                case 0:
                    strcpy(model_name, optimal_model_name[inst->model_type]);
                    break;
                case 1:
                    strcpy(model_name, math_model_name[inst->model_type]);
                    break;
                case 2:
                    strcpy(model_name, heuristic_model_name[inst->model_type]);
                    break;
                case 3:
                    strcpy(model_name, meta_heuristic_model_name[inst->model_type]);
                    break;
            }

            sprintf(out, "set output '../output/plot/%s_%s_%d.jpg'", inst->param.name, model_name, iter);
            sprintf(plot, "plot 'data.temp' with linespoints pt 7 lc rgb 'blue' lw 1 notitle");
            if (iter == -1) // final solution
            {
                sprintf(out, "set output '../output/plot/%s_%s_final.jpg'", inst->param.name, model_name);
                sprintf(plot, "plot 'data.temp' with linespoints pt 7 lc rgb 'red' lw 1 notitle");
            }

            char *commandsForGnuplot[] = {
                    "set title 'Solution plot'",
                    "set terminal jpeg size 1024,768",
                    out,
                    "unset key",
                    "set autoscale",
                    "set ylabel 'Y'",
                    "set xlabel 'X'",
                    plot,
                    "clear"};

            int commands = sizeof(commandsForGnuplot) / sizeof(commandsForGnuplot[0]);

            //Send commands to gnuplot one by one.
            for (int i = 0; i < commands; i++)
                fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);

            pclose(gnuplotPipe); // execute the commands
            free(plot);
            free(out);
            free(model_name);
        }

        if (inst->param.interactive) // plot solution
        {
            // here we use a global gnuplotPipe to remain on the same window
            char *commandsForGnuplot[] = {"set title 'Solution plot'",
                                          "unset key",
                                          "set autoscale",
                                          "set ylabel 'Y'",
                                          "set xlabel 'X'",
                                          "plot 'data.temp' with linespoints pt 7 lc rgb 'blue' lw 1"};

            int commands = sizeof(commandsForGnuplot) / sizeof(commandsForGnuplot[0]);
            for (int i = 0; i < commands; i++)
                fprintf(inst->gnuplotPipe, "%s \n", commandsForGnuplot[i]);
            fflush(inst->gnuplotPipe); // execute the commands
        }
    }
    return 0;
}

int generate_path(char *path, char *folder, char *type, const char *model, char *filename, int seed, char *extension)
{
    snprintf(path, 1000, "../%s/%s/%s_%s_%s_%d.%s", folder, extension, filename, type, model, seed, extension);
    return 0;
}

int generate_csv_record(char *instance_name, int seed, const char *model, double z_best, double time_elapsed, int run)
{
    FILE *csv;
    char filename[100];
    printf("run: %d\n", run);
    sprintf(filename, "../output/scores_%d.csv", run);

    if (access(filename, F_OK) == 0)
    {
        print_message("scores.csv exists, adding the results");
        csv = fopen(filename, "a");
        fprintf(csv, "%s, %d, %s, %f, %f, run-%d\n", instance_name, seed, model, z_best, time_elapsed, run);
    }
    else
    {
        print_message("scores.csv does not exists yet, creating file and adding the results");
        csv = fopen(filename, "w");
        fprintf(csv, "%s, %d, %s, %f, %f, run-%d\n", instance_name, seed, model, z_best, time_elapsed, run);
    }

    if (fclose(csv))
    {
        print_error("generate_csv_record fclose error");
    }

    return 0;
}

int checkFolders()
{
    if (!IsPathExist("../output"))
    {
        if (mkdir("../output", 0777))
            print_error("error creating output folder");
        if (mkdir("../output/log", 0777))
            print_error("error creating log folder");
        if (mkdir("../output/lp", 0777))
            print_error("error creating lp folder");
        if (mkdir("../output/plot", 0777))
            print_error("error creating plot folder");
    }
    else
    {
        if (!IsPathExist("../output/log"))
            if (mkdir("../output/log", 0777))
                print_error("error creating log folder");
        if (!IsPathExist("../output/lp"))
            if (mkdir("../output/lp", 0777))
                print_error("error creating lp folder");
        if (!IsPathExist("../output/plot"))
            if (mkdir("../output/plot", 0777))
                print_error("error creating plot folder");
    }
    return 0;
}

int IsPathExist(const char *s)
{
    struct stat buffer;
    return !stat(s, &buffer);
}