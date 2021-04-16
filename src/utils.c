#include "../include/utils.h"

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
                inst->time_left = inst->time_limit;
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

                inst->nodes = (node *)calloc(inst->dimension, sizeof(node));
                inst->edges = (edge *)calloc(inst->dimension, sizeof(edge));
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

    inst->model_type = 0;
    inst->time_limit = CPX_INFBOUND;
    inst->time_left = CPX_INFBOUND;
    inst->param.seed = 1;
    inst->param.run = -1;
    inst->param.verbose = NORMAL;
    inst->param.callback_counter = 0;

    inst->dimension = -1;
    inst->nodes = NULL;
    inst->edges = NULL;
    inst->n_edges = -1;

    inst->weights = NULL;
    inst->integer_costs = 1;
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

    // todo write the code to free the allocated memory within the instance (bottom-up approach)
}

void print_command_line(instance *inst)
{
    printf("\nPARAMETERS ---------------------------------------------\n");
    printf("-f %s\n", inst->param.input_file);
    printf("-t %f seconds\n", inst->time_limit);
    printf("-m %d (%s)\n", inst->model_type, model_name[inst->model_type]);
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
    printf("Model type: %d (%s)\n", inst->model_type, model_name[inst->model_type]);
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
    printf("-f <path>  : used to pass the relative instance path \n");
    printf("-t <time>  : used to pass the total running time allowed in seconds\n");
    printf("-m <model> : used to set the model type\n");
    printf("-s <seed>  : used to set the seed\n");
    printf("-v <value> : used to set the verbosity, from QUIET (0) up to DEBUG (4)\n");
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

int plot_solution(instance *inst)
{
    char *commandsForGnuplot[] = {"set title 'Solution plot'",
                                  //"set terminal svg size 350,262",
                                  //"set output '../output/plot/test.svg", // TODO : find a way to modify the nama of the output file
                                  "unset key",
                                  "set autoscale",
                                  "set ylabel 'Y'",
                                  "set xlabel 'X",
                                  "plot 'data.temp' with linespoints pt 7 lc rgb 'blue' lw 1"};
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    FILE *temp = fopen("data.temp", "w");

    for (int i = 0; i < inst->n_edges; i++)
    {
        int prev = inst->edges[i].prev;
        int next = inst->edges[i].next;
        //Write the coordinates of the two nodes inside a temporary file
        fprintf(temp, "%lf %lf \n%lf %lf \n\n", inst->nodes[prev].x, inst->nodes[prev].y, inst->nodes[next].x,
                inst->nodes[next].y);
        // double '\n' to create a new block of coordinates
    }
    fclose(temp);

    int commands = sizeof(commandsForGnuplot) / sizeof(commandsForGnuplot[0]);
    for (int i = 0; i < commands; i++)
    {
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
    }
    if (pclose(gnuplotPipe) == -1)
    {
        print_error("pclose error");
        return -1;
    }
    return 0;
}

int plot_solution_edges(int n_edges, node *nodes, edge *edges)
{
    char *commandsForGnuplot[] = {"set title 'Solution plot'",
                                  "unset key",
                                  "set autoscale",
                                  "set ylabel 'Y'",
                                  "set xlabel 'X",
                                  "plot 'data.temp' with linespoints pt 7 lc rgb 'blue' lw 1"};
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    FILE *temp = fopen("data.temp", "w");

    for (int i = 0; i < n_edges; i++)
    {
        //Write the coordinates of the two nodes inside a temporary file
        fprintf(temp, "%lf %lf \n", nodes[edges[i].prev].x, nodes[edges[i].prev].y);   // first node
        fprintf(temp, "%lf %lf \n\n", nodes[edges[i].next].x, nodes[edges[i].next].y); // second node
        // double '\n' to create a new block of coordinates
    }
    fclose(temp);

    int commands = sizeof(commandsForGnuplot) / sizeof(commandsForGnuplot[0]);
    for (int i = 0; i < commands; i++)
    {
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
    }
    if (pclose(gnuplotPipe) == -1)
    {
        print_error("pclose error");
        return -1;
    }
    return 0;
}

int generate_path(char *path, char *folder, char *type, const char *model, char *filename, int seed, char *extension)
{
    snprintf(path, 1000, "../%s/%s/%s_%s_%s_%d.%s", folder, extension, filename, type, model, seed, extension);
    return 0;
}

int generate_csv_record(char *instance_name, int seed, int model_type, double z_best, double time_elapsed, int run)
{
    FILE *csv;
    char filename[100];
    printf("run: %d\n", run);
    sprintf(filename, "../output/scores_%d.csv", run);

    if (access(filename, F_OK) == 0)
    {
        print_message("scores.csv exists, adding the results");
        csv = fopen(filename, "a");
        fprintf(csv, "%s, %d, %s, %f, %f, run-%d\n", instance_name, seed, model_name[model_type], z_best, time_elapsed, run);
    }
    else
    {
        print_message("scores.csv does not exists yet, creating file and adding the results");
        csv = fopen(filename, "w");
        fprintf(csv, "%s, %d, %s, %f, %f, run-%d\n", instance_name, seed, model_name[model_type], z_best, time_elapsed, run);
    }

    if (fclose(csv))
    {
        print_error("generate_csv_record fclose error");
    }

    return 0;
}