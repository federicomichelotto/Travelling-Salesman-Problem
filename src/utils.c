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
                strcpy(inst->instance_path, argv[++i]);
                if (inst->instance_path[strlen(inst->instance_path) - 1] == '/')
                    inst->instance_path[strlen(inst->instance_path) - 1] = '\0';
                char *last = strrchr(inst->instance_path, '/');
                if (last == NULL)
                    strcpy(inst->instance_name, inst->instance_path);
                else
                    strcpy(inst->instance_name, last + 1);
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

            if (strcmp(argv[i], "-scores") == 0)
            { // name file score
                strcpy(inst->param.scores, argv[++i]);
                continue;
            }

            if (strcmp(argv[i], "-v") == 0)
            { // Verbosity
                i++;
                int v = strtol(argv[i], NULL, 10);
                if (v > DEBUG)
                    inst->param.verbose = DEBUG;
                else if (v < QUIET)
                    inst->param.verbose = QUIET;
                else
                {
                    switch (v)
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
                }
                continue;
            }

            if (strcmp(argv[i], "--ticks") == 0)
            { // use ticks instead of seconds
                inst->param.ticks = 1;
                continue;
            }
            if (strcmp(argv[i], "--interactive") == 0)
            { // show plots
                inst->param.interactive = 1;
                continue;
            }

            if (strcmp(argv[i], "--saveplots") == 0)
            { // save plots
                inst->param.saveplots = 1;
                continue;
            }
            if (strcmp(argv[i], "--heur") == 0)
            { // use the heuristic method instead of the exact one
                inst->param.heur = 1;
                continue;
            }
        }

        print_command_line(inst);
    }
}

void parse_instance(instance *inst)
{
    parse_locations(inst);
    allocate_mem_truck_times_and_dists(inst);
    parse_truck_travel_data(inst);
    allocate_mem_drone_min_max_times(inst);
    if (!inst->param.heur)
        parse_min_max_drone_legs_times(inst);
    else
        compute_opt_speeds(inst);
    // retrieve the maximum drone travel time for each (i,j,-)
    retrieve_max_ij_drone_leg_times(inst);
    // allocate the memory for the successor array (array in which the successor of each node is stored)
    inst->succ_truck = (int *)calloc(inst->dimension - 1, sizeof(int));
    inst->succ_drone = (int *)calloc(inst->dimension - 1, sizeof(int));
}

// for instances in https://github.com/optimatorlab/mFSTSP-VDS/tree/master/Problems
// nodeID, nodeType, latDeg, lonDeg, altMeters, parcelWtLbs
// nodeType 0 represents the depot, nodeType 1 represents a customer.
// The depot is always always assigned to nodeID 0.
// Each node has a corresponding latitude and longitude (specified in degrees).
// The altitude is always 0. Customer nodes have a corresponding non-zero parcel weight (in [pounds]).
// There is no parcel associated with the depot.
void parse_locations(instance *inst)
{
    printf("\n\nParsing the locations...\n");
    FILE *fp = NULL;
    char filename[1000];
    strcpy(filename, inst->instance_path);
    strcat(filename, "/tbl_locations.csv");

    fp = fopen(filename, "r");
    if (fp == NULL)
        print_error("Could not open the file ");

    char line[1000]; // tmp string where to store a line

    // Read each line from the input file:
    // split the line using strtok and consider the generated tokens

    // read (ignore) the first line
    if (fgets(line, sizeof(line), fp) == NULL)
        print_error("Error reading the first file of the input data.");

    inst->dimension = 1000; // initial instance dimension
    // allocate a temporary memory that will be reallocated when the number of nodes will be known or when the memory is full
    inst->nodes = (node *)calloc(inst->dimension, sizeof(node));

    int i = 0; // #nodes counter
    while (fgets(line, sizeof(line), fp) != NULL)
    {
        if (inst->param.verbose > DEBUG)
            printf("%s\n", line);
        line[strcspn(line, "\n")] = 0; // removing trailing \n

        // double the memory of the inst->nodes array
        if (i == inst->dimension)
        {
            inst->dimension *= 2;
            inst->nodes = (node *)realloc(inst->nodes, inst->dimension * sizeof(node));
        }
        char delimiters[] = ",";

        inst->nodes[i].id = strtol(strtok(line, delimiters), NULL, 10);
        int nodeType = strtol(strtok(NULL, delimiters), NULL, 10); // nodeType
        if (!i && nodeType)
            print_error("Error: node 0 must be the depot (nodeType = 0).");

        inst->nodes[i].x = strtod(strtok(NULL, delimiters), NULL);                   // latitude
        inst->nodes[i].y = strtod(strtok(NULL, delimiters), NULL);                   // longitude
        strtok(NULL, delimiters);                                                    // ignore altitude
        inst->nodes[i].weight = strtof(strtok(NULL, delimiters), NULL) * 0.45359237; // parcel weight (lbs -> kg)

        if (inst->param.verbose >= DEBUG)
            printf("\tNODES %2d at coordinates (%11.6lf , %11.6lf) \t parcel weight = %7.2lf kg \n", inst->nodes[i].id, inst->nodes[i].x, inst->nodes[i].y, inst->nodes[i].weight);
        i++;
    }

    // adding an extra node that represent the returning depot
    inst->nodes[i].id = i;
    inst->nodes[i].x = inst->nodes[0].x;
    inst->nodes[i].y = inst->nodes[0].y;
    inst->nodes[i].weight = inst->nodes[0].weight;

    inst->dimension = ++i;                                                      // correct number of nodes in the instance
    inst->nodes = (node *)realloc(inst->nodes, inst->dimension * sizeof(node)); // realloc

    fclose(fp);

    if (inst->dimension < 3)
        print_error("EXIT: the instance does not have any customer.");
}

// tbl_truck_travel_data_PG.csv contains the directed truck travel time and distance information from one node to another.
// All time and distance values were obtained by pgRouting, using OpenStreetMaps data.
// from location i, to location j, time [sec], distance [meters]
void parse_truck_travel_data(instance *inst)
{
    printf("\n\nParsing truck travel times and distances...\n");

    FILE *fp = NULL;
    char filename[1000];
    strcpy(filename, inst->instance_path);
    strcat(filename, "/tbl_truck_travel_data_PG.csv");

    fp = fopen(filename, "r");
    if (fp == NULL)
        print_error("Could not open the file");

    char line[1000]; // tmp string where to store a line
    // Read each line from the input file:
    // split the line using strtok and consider the generated tokens

    // read (ignore) the first line (header)
    if (fgets(line, sizeof(line), fp) == NULL)
        print_error("Error reading the first file of the input data.");

    int nlines = 0; // #lines counter
    while (fgets(line, sizeof(line), fp) != NULL)
    {
        line[strcspn(line, "\n")] = 0; // removing trailing \n
        char delimiters[] = ",";

        int i = strtol(strtok(line, delimiters), NULL, 10);               // node i
        int j = strtol(strtok(NULL, delimiters), NULL, 10);               // node j
        inst->truck_times[i][j] = strtod(strtok(NULL, delimiters), NULL); // travel time
        if (j != inst->dimension - 1)
            inst->truck_times[i][j] += 30;                                // delivery time
        inst->truck_dists[i][j] = strtod(strtok(NULL, delimiters), NULL); // travel distance

        nlines++;
    }
    if (nlines != (inst->dimension - 1) * (inst->dimension - 1))
        printf("WARNING: for some (i,j) truck travel times/distances are missing, by default these values are set to 0.\n");

    fclose(fp);

    // define the time/distance [customer -> returning depot] equal to time/distance [customer -> starting depot]
    for (int i = 0; i < inst->dimension - 1; i++)
    {
        inst->truck_times[i][inst->dimension - 1] = inst->truck_times[i][0];
        inst->truck_dists[i][inst->dimension - 1] = inst->truck_dists[i][0];
    }

    if (inst->param.verbose >= DEBUG)
    {
        for (int i = 0; i < inst->dimension; i++)
        {
            for (int j = 0; j < inst->dimension; j++)
            {
                printf("\ttruck_time(%3d,%3d) = %15.6lf \t\t dist(%3d,%3d) = %15.6lf \n", i, j, inst->truck_times[i][j], i, j, inst->truck_dists[i][j]);
            }
        }
    }
}

void allocate_mem_truck_times_and_dists(instance *inst)
{
    // initialize the matrix of the truck travel times
    inst->truck_dists = (double **)calloc(inst->dimension, sizeof(double *));
    inst->truck_times = (double **)calloc(inst->dimension, sizeof(double *));
    for (int i = 0; i < inst->dimension; i++)
    {
        inst->truck_dists[i] = (double *)calloc(inst->dimension, sizeof(double));
        inst->truck_times[i] = (double *)calloc(inst->dimension, sizeof(double));
    }
    // truck times/dists can be accessed as follows:
    // truck_dists[i][j]
    // truck_times[i][j]
}

void allocate_mem_drone_min_max_times(instance *inst)
{
    // initialize the 3D matrices of the min/max drone's travel times
    inst->min_time_drone = (double ***)calloc(inst->dimension, sizeof(double *));
    inst->min_feas_time_drone = (double ***)calloc(inst->dimension, sizeof(double *));

    inst->max_time_drone = (double ***)calloc(inst->dimension, sizeof(double *));
    inst->max_feas_time_drone = (double ***)calloc(inst->dimension, sizeof(double *));

    for (int i = 0; i < inst->dimension; i++)
    {
        inst->min_time_drone[i] = (double **)calloc(inst->dimension, sizeof(double *));
        inst->min_feas_time_drone[i] = (double **)calloc(inst->dimension, sizeof(double *));

        inst->max_feas_time_drone[i] = (double **)calloc(inst->dimension, sizeof(double *));
        inst->max_time_drone[i] = (double **)calloc(inst->dimension, sizeof(double *));
    }

    for (int i = 0; i < inst->dimension; i++)
    {
        for (int j = 0; j < inst->dimension; j++)
        {
            inst->min_time_drone[i][j] = (double *)calloc(inst->dimension, sizeof(double));
            inst->min_feas_time_drone[i][j] = (double *)calloc(inst->dimension, sizeof(double));

            inst->max_feas_time_drone[i][j] = (double *)calloc(inst->dimension, sizeof(double));
            inst->max_time_drone[i][j] = (double *)calloc(inst->dimension, sizeof(double));
        }
    }
    // drone min/max times can be accessed as follows:
    // min_time_drone[i][j][k]
    // max_time_drone[i][j][k]
}

// tbl_truck_travel_data_PG.csv contains the directed truck travel time and distance information from one node to another.
// All time and distance values were obtained by pgRouting, using OpenStreetMaps data.
// from location i, to location j, time [sec], distance [meters]
void parse_min_max_drone_legs_times(instance *inst)
{
    printf("\n\nParsing the min/max drone legs travel times...\n");

    FILE *fp = NULL;
    char filename[1000];
    strcpy(filename, inst->instance_path);
    strcat(filename, "/min_max.csv");

    fp = fopen(filename, "r");
    if (fp == NULL)
        print_error("Could not open the file");

    char line[1000]; // tmp string where to store a line
    // Read each line from the input file:
    // split the line using strtok and consider the generated tokens

    // read (ignore) the first line (header)
    if (fgets(line, sizeof(line), fp) == NULL)
        print_error("Error reading the first file of the input data.");

    int nlines = 0;     // #lines counter
    int flag_index = 0; // flag set to 1 if indices not valid are found
    while (fgets(line, sizeof(line), fp) != NULL)
    {
        line[strcspn(line, "\n")] = 0; // removing trailing \n
        char delimiters[] = ",";

        int i = strtol(strtok(line, delimiters), NULL, 10); // node i
        int j = strtol(strtok(NULL, delimiters), NULL, 10); // node j
        int k = strtol(strtok(NULL, delimiters), NULL, 10); // node k

        if (i < 0 || j < 0 || k < 0 || i >= inst->dimension || j >= inst->dimension || k >= inst->dimension)
        {
            flag_index = 1;
            strtok(NULL, delimiters);
            strtok(NULL, delimiters);
            strtok(NULL, delimiters);
            continue;
        }

        inst->min_time_drone[i][j][k] = strtod(strtok(NULL, delimiters), NULL);      // min
        inst->min_feas_time_drone[i][j][k] = strtod(strtok(NULL, delimiters), NULL); //  min feasible

        inst->max_feas_time_drone[i][j][k] = strtod(strtok(NULL, delimiters), NULL); // max feasible
        inst->max_time_drone[i][j][k] = strtod(strtok(NULL, delimiters), NULL);      //  max

        if (inst->param.verbose >= DEBUG)
            printf("\t[%3d,%3d,%3d] = \t %12.6lf \t %12.6lf \t %12.6lf \t %12.6lf \n", i, j, k, inst->min_time_drone[i][j][k], inst->min_feas_time_drone[i][j][k], inst->max_feas_time_drone[i][j][k], inst->max_time_drone[i][j][k]);
        nlines++;
    }

    if (flag_index)
        printf("WARNING [parse_min_max_drone_legs_times]: some index is not valid.");

    fclose(fp);
}

void initialize_instance(instance *inst)
{
    srand(time(NULL)); // Initialize a random seed for the random numbers generator
    inst->model_type = 0;
    inst->time_limit = CPX_INFBOUND;
    inst->param.seed = 1;
    inst->param.run = 1;
    inst->param.verbose = NORMAL;
    inst->param.ticks = 0;
    inst->param.interactive = 0;
    inst->param.saveplots = 0;
    inst->param.heur = 0;

    inst->dimension = -1;
    inst->nodes = NULL;
    inst->succ_truck = NULL;
    inst->succ_drone = NULL;
    inst->truck_dists = NULL;
    inst->truck_times = NULL;

    inst->z_best = -1.0;

    inst->timestamp_start = 0.0;
    inst->timestamp_finish = 0.0;

    inst->timestamp_last_plot = 0.0;
    inst->plot_counter = 0;

    strcpy(inst->instance_path, "NULL");
    strcpy(inst->instance_name, "NULL");
    strcpy(inst->param.scores, "0");
    inst->gnuplotPipe = popen("gnuplot -persistent", "w");
}

void free_instance(instance *inst)
{
    // close gnuplot pipe
    if (pclose(inst->gnuplotPipe) == -1)
        print_error("pclose error");

    free(inst->nodes);
    free(inst->succ_drone);
    free(inst->succ_truck);
}

void print_command_line(instance *inst)
{
    printf("\nPARAMETERS ---------------------------------------------\n");
    printf("-f (file path) %s\n", inst->instance_path);
    if (inst->param.ticks == 1)
    {
        printf("-t (ticks limit) %.0f ticks\n", inst->time_limit);
        printf("--ticks (ACTIVE -> %d)\n", inst->param.ticks);
    }
    else
    {
        printf("-t (time limit) %.0f seconds\n", inst->time_limit);
    }

    if (inst->param.interactive == 1)
    {
        printf("--interactive (ACTIVE -> %d)\n", inst->param.interactive);
    }

    if (inst->param.saveplots == 1)
    {
        printf("--saveplots (ACTIVE -> %d)\n", inst->param.saveplots);
    }

    printf("-s (seed) %d\n", inst->param.seed);
    printf("-v (verbosity) %d (%s)\n", inst->param.verbose, verbose_name[inst->param.verbose]);
    printf("--------------------------------------------------------\n\n");
}

void print_instance(instance *inst)
{
    // last node is the "fake" returning depot
    printf("\n\nINSTANCE -----------------------------------------------\n");
    printf("Name: %s\n", inst->instance_name);
    printf("Path: %s\n", inst->instance_path);
    printf("Dimension: %d (+ a returning depot)\n", inst->dimension - 1);

    if (inst->nodes != NULL)
    {
        printf("\nNode_id %*s Latitude %*s Longitude   Parcel_weight \n", 7, "", 5, "");
        for (int i = 0; i < inst->dimension - 1; i++)
            printf("%7d\t %15.6lf %15.6lf %12.2lf kg\n", inst->nodes[i].id, inst->nodes[i].x, inst->nodes[i].y, inst->nodes[i].weight);
    }
    printf("*%6d\t %15.6lf %15.6lf %12.2lf kg\n", inst->nodes[inst->dimension - 1].id, inst->nodes[inst->dimension - 1].x, inst->nodes[inst->dimension - 1].y, inst->nodes[inst->dimension - 1].weight);

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

void save_and_plot_solution(instance *inst, int iter)
{
    save_and_plot_solution_general(inst, inst->succ_truck, inst->succ_drone, iter);
}

void save_and_plot_solution_general(instance *inst, int *succ_truck, int *succ_drone, int iter)
{
    if (inst->param.saveplots || inst->param.interactive || iter == -1)
    {
        // write solution to file

        char data_points_filename[200];
        char truck_edges_filename[200];
        char drone_edges_filename[200];

        sprintf(data_points_filename, "data_temp/data_points_%d", inst->plot_counter);
        sprintf(truck_edges_filename, "data_temp/truck_edges_%d", inst->plot_counter);
        sprintf(drone_edges_filename, "data_temp/drone_edges_%d", inst->plot_counter);
        FILE *data_points = fopen(data_points_filename, "w");
        FILE *truck_edges = fopen(truck_edges_filename, "w");
        FILE *drone_edges = fopen(drone_edges_filename, "w");

        // prepare the data to pass to gnuplot
        for (int i = 0; i < inst->dimension - 1; i++)
        {
            fprintf(data_points, "%lf %lf %d \n", inst->nodes[i].x, inst->nodes[i].y, inst->nodes[i].id);

            if (succ_truck[i] != -1)
            {
                double delta_x = inst->nodes[succ_truck[i]].x - inst->nodes[i].x;
                double delta_y = inst->nodes[succ_truck[i]].y - inst->nodes[i].y;
                if (fabs(delta_x) > fabs(delta_y))
                {
                    if (delta_x > 0)
                        delta_x -= 0.003;
                    else
                        delta_x += 0.003;
                }
                else
                {
                    if (delta_y > 0)
                        delta_y -= 0.003;
                    else
                        delta_y += 0.003;
                }
                fprintf(truck_edges, "%lf %lf %lf %lf \n", inst->nodes[i].x, inst->nodes[i].y, delta_x, delta_y);
            }
            if (succ_drone[i] != -1)
            {
                double delta_x = inst->nodes[succ_drone[i]].x - inst->nodes[i].x;
                double delta_y = inst->nodes[succ_drone[i]].y - inst->nodes[i].y;
                if (fabs(delta_x) > fabs(delta_y))
                {
                    if (delta_x > 0)
                        delta_x -= 0.003;
                    else
                        delta_x += 0.003;
                }
                else
                {
                    if (delta_y > 0)
                        delta_y -= 0.003;
                    else
                        delta_y += 0.003;
                }
                fprintf(drone_edges, "%lf %lf %lf %lf \n", inst->nodes[i].x, inst->nodes[i].y, delta_x, delta_y);
            }
        }
        fclose(data_points);
        fclose(truck_edges);
        fclose(drone_edges);

        if (inst->param.saveplots || iter == -1) // save plot
        {
            FILE *gnuplotPipe = popen("gnuplot", "w"); // local gnuplotPipe to avoid conflicts with the global gnuplotPipe
            char *out = (char *)calloc(1000, sizeof(char));
            char *plot_str = (char *)calloc(1000, sizeof(char));

            sprintf(out, "set output '../output/%s/seed_%d/run_%d/plot/sol_%d.jpg'", inst->instance_name, inst->param.seed, inst->param.run, iter);
            sprintf(plot_str, "plot '%s' with vectors head filled lc rgb 'blue' lw 1,\
                                '%s' with vectors head filled lc rgb 'red' dashtype '-_' lw 1,\
                                '%s' using 1:2:(0.003) with circles linecolor rgb 'white' lw 1 fill solid border lc lt 0,\
                                '%s' using 1:2:3 with labels offset (0,0) font 'Arial'",
                    truck_edges_filename, drone_edges_filename, data_points_filename, data_points_filename);

            // sprintf(plot_str, "plot '%s' using 1:2:3:4 with linespoints pt 7 lc rgb 'black' lw 1,\
            //                     '%s' using 1:2:3:4 with linespoints pt 7 lc rgb 'red' lw 1,\
            //                     '%s' using 1:2:3 with labels offset (-1,-1)",
            //                     truck_edges_filename, drone_edges_filename, data_points_filename);
            char *commandsForGnuplot_drone[] = {"set title 'FSTSP solution'",
                                                "set terminal jpeg size 1024,768",
                                                out,
                                                "unset key",
                                                //"set autoscale",
                                                "set ylabel 'Y'",
                                                "set xlabel 'X'",
                                                plot_str};

            //Send commands to gnuplot one by one.
            int commands = sizeof(commandsForGnuplot_drone) / sizeof(commandsForGnuplot_drone[0]);
            for (int i = 0; i < commands; i++)
                fprintf(gnuplotPipe, "%s \n", commandsForGnuplot_drone[i]);

            pclose(gnuplotPipe); // execute the commands
            free(plot_str);
            free(out);
        }

        if (inst->param.interactive) // plot solution
        {
            char plot_str[5000];
            sprintf(plot_str, "plot '%s' with vectors head filled lc rgb 'blue' lw 1,\
                                '%s' with vectors head filled lc rgb 'red' dashtype '-_' lw 1,\
                                '%s' using 1:2:(0.003) with circles linecolor rgb 'white' lw 1 fill solid border lc lt 0,\
                                '%s' using 1:2:3 with labels offset (0,0) font 'Arial'",
                    truck_edges_filename, drone_edges_filename, data_points_filename, data_points_filename);

            // sprintf(plot_str, "plot '%s' using 1:2:3:4 with linespoints pt 7 lc rgb 'black' lw 1,\
            //                     '%s' using 1:2:3:4 with linespoints pt 7 lc rgb 'red' lw 1,\
            //                     '%s' using 1:2:3 with labels offset (-1,-1)",
            //                     truck_edges_filename, drone_edges_filename, data_points_filename);
            char *commandsForGnuplot_drone[] = {"set title 'Optimal FSTSP solution'",
                                                "set term wxt noraise",
                                                "unset key",
                                                //"set autoscale",
                                                "set ylabel 'Y'",
                                                "set xlabel 'X'",
                                                plot_str};

            int commands = sizeof(commandsForGnuplot_drone) / sizeof(commandsForGnuplot_drone[0]);
            for (int i = 0; i < commands; i++)
                fprintf(inst->gnuplotPipe, "%s \n", commandsForGnuplot_drone[i]);
            fflush(inst->gnuplotPipe); // execute the commands
        }
    }
}

int generate_csv_record(char *path, char *instance_name, int seed, int model, int run, double z_best, double time_elapsed, int ticks)
{
    FILE *csv;
    char score_path[1000];
    strcpy(score_path, path);
    int flag_new_file = 0;
    if (!IsPathExist(score_path))
        flag_new_file = 1;
    csv = fopen(score_path, "a");
    if (csv == NULL)
        print_error("generate_csv_record fopen() error");
    if (flag_new_file)
    {
        if (!ticks)
            fprintf(csv, "# name instance, model, seed, run, incumbent, time elapsed [s]\n");
        else
            fprintf(csv, "# name instance, model, seed, run, incumbent, time elapsed [ticks]\n");
    }
    fprintf(csv, "%s, %d, %d, %d, %f, %f\n", instance_name, model, seed, run, z_best, time_elapsed);

    if (fclose(csv))
        print_error("generate_csv_record fclose() error");

    return 0;
}

void createInstanceFolders(instance *inst)
{
    if (!IsPathExist("../output"))
    {
        if (mkdir("../output", 0777))
            print_error("error creating output folder");
    }
    char instance_folder_path[1000];
    sprintf(instance_folder_path, "../output/%s", inst->instance_name);
    if (!IsPathExist(instance_folder_path))
    {
        if (mkdir(instance_folder_path, 0777))
            print_error("error creating instance folder");
    }
    char seed_path[1050];
    sprintf(seed_path, "%s/seed_%d", instance_folder_path, inst->param.seed);
    if (!IsPathExist(seed_path))
    {
        if (mkdir(seed_path, 0777))
            print_error("error creating seed folder");
    }
    setRunNumber(inst);
    char run_path[1100];
    sprintf(run_path, "%s/run_%d", seed_path, inst->param.run);
    if (mkdir(run_path, 0777))
        print_error("error creating run folder");
    char plot_path[1110];
    sprintf(plot_path, "%s/plot", run_path);
    if (mkdir(plot_path, 0777))
        print_error("error creating plot folder");

    if (!IsPathExist("data_temp"))
        if (mkdir("data_temp", 0777))
            print_error("error creating temp folder");
}

int IsPathExist(const char *s)
{
    struct stat buffer;
    return !stat(s, &buffer);
}

void getTimeStamp(double *ts)
{
    struct timespec timestamp;
    if (clock_gettime(CLOCK_REALTIME, &timestamp) == -1)
        print_error("Error clock_gettime");
    *ts = timestamp.tv_sec + timestamp.tv_nsec * pow(10, -9);
}

// comparison function for qsort
int compare(const void *a, const void *b)
{
    if (*(double *)a > *(double *)b)
        return 1;
    if (*(double *)a < *(double *)b)
        return -1;
    return 0;
}

void compute_opt_speeds(instance *inst)
{
    // double *sorted_weights = (double *)calloc(inst->dimension, sizeof(double));
    // int j = 0;
    // for (int i = 0; i < inst->dimension; i++)
    // {
    //     if (inst->nodes[i].weight > 0 && inst->nodes[i].weight < WEIGHT_LIMIT)
    //         sorted_weights[j++] = inst->nodes[i].weight;
    // }
    // sorted_weights = (double *)realloc(sorted_weights, j * sizeof(double));
    // printf("%d\n", j);
    // for (int i = 0; i < j; i++)
    // {
    //     printf("%f, ", sorted_weights[i]);
    // }
    // qsort(sorted_weights, j, sizeof(double), compare);
    // printf("\n");
    // for (int i = 0; i < j; i++)
    // {
    //     printf("%f, ", sorted_weights[i]);
    // }
    // printf("\n");
    // // delete duplicates
    // int last = 0;
    // for (int i = 1; i < j; i++)
    // {
    //     if (sorted_weights[i] != sorted_weights[last])
    //     {
    //         sorted_weights[++last] = sorted_weights[i];
    //     }
    // }
    // sorted_weights = (double *)realloc(sorted_weights, (last + 1) * sizeof(double));
    // printf("\n");
    // for (int i = 0; i < last + 1; i++)
    // {
    //     printf("%f, ", sorted_weights[i]);
    // }

    double B = 0.5 * pow(10, 6); // battery capacity [J]

    // speed range for the minimum time
    double max_speed_min = 40.0;
    double min_speed_min = 13.0;

    // speed range for the maximum time
    double min_speed_max = 1.0;
    double max_speed_max = 23.0;

    double granularity = 0.5;
    int nspeeds = (max_speed_min - min_speed_max) / granularity + 1;
    // speeds
    double speeds[nspeeds];

    // initialize speeds
    // initialize times w.r.t. to each speed
    for (int i = 0; i < nspeeds; i++)
    {
        speeds[i] = min_speed_max + (i + 1) * granularity;
    }

    double *cruise_power_0 = (double *)calloc(nspeeds, sizeof(double));
    cruise_power_consumption(0, nspeeds, speeds, cruise_power_0);

    double takeoff_speed = 10.0;   // [m/s]
    double landing_speed = 5.0;    // [m/s]
    double cruise_altitude = 50.0; // [m]
    double takeoff_power_0 = vertical_power_consumption(0, takeoff_speed);
    double landing_power_0 = vertical_power_consumption(0, landing_speed);

    for (int i = 0; i < inst->dimension - 1; i++)
    {
        for (int j = 1; j < inst->dimension - 1; j++)
        {
            if (j == i)
                continue;
            for (int k = 1; k < inst->dimension; k++)
            {
                if (k == j || k == i)
                    continue;
                double w = inst->nodes[j].weight;
                if (w > WEIGHT_LIMIT)
                    continue;
                double d1 = dist(i, j, inst);
                double d2 = dist(j, k, inst);
                //printf("d1 = %f, d2 = %f\n", d1, d2);

                int opt_min_idx_s1, opt_min_idx_s2;
                int opt_max_idx_s1, opt_max_idx_s2;

                // compute the energy for each couple of speed (s1,s2)
                double *cruise_power_w = (double *)calloc(nspeeds, sizeof(double));
                cruise_power_consumption(w, nspeeds, speeds, cruise_power_w);

                // travel times
                double T1[nspeeds];
                double T2[nspeeds];
                for (int i = 0; i < nspeeds; i++)
                {
                    T1[i] = d1 * (1 / speeds[i]);
                    T2[i] = d2 * (1 / speeds[i]);
                }

                double takeoff_power_w = vertical_power_consumption(w, takeoff_speed);
                double landing_power_w = vertical_power_consumption(w, landing_speed);
                // remeaning battery capacity for the horizontal flights [J]
                double newB = B - (cruise_altitude / takeoff_speed) * (takeoff_power_0 + takeoff_power_w) - (cruise_altitude / landing_speed) * (landing_power_w + landing_power_0);
                //printf("newB = %f\n", newB);
                // compute the optimal feasible couple of speed and the corresponding travel time
                printf("i = %d, j = %d, k = %d \n", i, j, k);
                double min_tt = min_time(newB, cruise_power_w, cruise_power_0, min_speed_min, granularity, nspeeds, speeds, T1, T2, &opt_min_idx_s1, &opt_min_idx_s2);
                if (min_tt < 0)
                {
                    printf("No feasible speeds have been found. \n\n");
                }
                else
                {
                    inst->min_time_drone[i][j][k] = min_tt + 150 + 2 * (cruise_altitude / takeoff_speed) + 2 * (cruise_altitude / landing_speed);
                    printf("\t Min time: %.3f seconds \n", inst->min_time_drone[i][j][k]);
                    //printf("\t [Min] Optimal speeds: (%.1f m/s, %.1f m/s) \n", speeds[opt_min_idx_s1], speeds[opt_min_idx_s2]);
                    //printf("\t [Min] Energy = %f kJ \n\n", (T1[opt_min_idx_s1] * cruise_power_w[opt_min_idx_s1] + T2[opt_min_idx_s2] * cruise_power_0[opt_min_idx_s2]) / 1000);

                    double max_tt = max_time(newB, cruise_power_w, cruise_power_0, max_speed_max, granularity, nspeeds, speeds, T1, T2, &opt_max_idx_s1, &opt_max_idx_s2);
                    if (max_tt < 0)
                    {
                        printf("No feasible speeds have been found. \n\n");
                    }
                    else
                    {
                        inst->max_time_drone[i][j][k] = max_tt + 150 + 2 * (cruise_altitude / takeoff_speed) + 2 * (cruise_altitude / landing_speed);
                        printf("\t Max time: %.3f seconds \n", inst->max_time_drone[i][j][k]);
                        //printf("\t [Max] Optimal speeds: (%.1f m/s, %.1f m/s) \n", speeds[opt_max_idx_s1], speeds[opt_max_idx_s2]);
                        //printf("\t [Max] Energy (to go from i to k) = %f kJ \n\n", (T1[opt_max_idx_s1] * cruise_power_w[opt_max_idx_s2] + T2[opt_min_idx_s2] * cruise_power_0[opt_min_idx_s2]) / 1000);
                    }
                }
                // free cruise_power_w array
                free(cruise_power_w);
            }
        }
    }

    // struct timespec timestamp;
    // if (clock_gettime(CLOCK_REALTIME, &timestamp) == -1)
    //     print_error("Error clock_gettime");
    // double timestamp_start = timestamp.tv_sec + timestamp.tv_nsec * pow(10, -9);

    // if (clock_gettime(CLOCK_REALTIME, &timestamp) == -1)
    //     print_error("Error clock_gettime");
    // double timestamp_finish = timestamp.tv_sec + timestamp.tv_nsec * pow(10, -9);

    // double time_elapsed = timestamp_finish - timestamp_start;
    // printf("Elapsed time: %f seconds\n", time_elapsed);

    // free cruise_power_0 array
    free(cruise_power_0);
}

// compute the takeoff/landing power consumption
double vertical_power_consumption(double w, double speed)
{
    double k1 = 0.8554;
    double k2 = 0.3051;
    double c2 = 0.3177;
    double g = 9.8;
    double W = 1.5;
    double P = k1 * (W + w) * g * (speed / 2 + sqrt(speed * speed / 4 + (W + w) * g / (k2 * k2))) + c2 * pow((W + w) * g, 3 / 2);
    return P;
}

// compute the cruise power consumption to travel one second at speed s[i] with payload w
void cruise_power_consumption(double w, int nspeeds, double *speeds, double *P)
{
    double c1 = 2.8037;
    double c2 = 0.3177;
    double c4 = 0.0296;
    double c5 = 0.0279;
    double g = 9.8;
    double W = 1.5;

    double gforce_1 = (W + w) * g;
    double gforce_2 = W * g;
    double c12 = c1 + c2;
    double alpha_rad = 10.0 / 180 * M_PI;
    float c5_f = c5 * cos(alpha_rad) * cos(alpha_rad);
    double c4_f = c4 * c4;

    for (int i = 0; i < nspeeds; i++)
    {
        P[i] = c12 * pow(pow((gforce_1 - c5_f * speeds[i] * speeds[i]), 2) + c4_f * pow(speeds[i], 4), 0.75) + c4 * pow(speeds[i], 3);
    }
}

// energy_grid = matrix of energy for all speed combinatation
// B = battery capacity (e.g 500 KJ)
double min_time(double B, double *power_w, double *power_0, double min_speed, double granularity, int nspeeds, double *speeds, double *T1, double *T2, int *opt_idx_s1, int *opt_idx_s2)
{
    int opt_i = -1, opt_j = -1;
    double min_time = DBL_MAX;
    int min_speed_idx = (min_speed - 1.0) / granularity;

    for (int i = min_speed_idx; i < nspeeds; i++)
    {
        for (int j = min_speed_idx; j < nspeeds; j++)
        {
            if (T1[i] * power_w[i] + T2[j] * power_0[j] > B)
                continue;
            if (T1[i] + T2[j] < min_time)
            {
                min_time = T1[i] + T2[j];
                opt_i = i;
                opt_j = j;
            }
        }
    }
    if (opt_i == -1)
    {
        return -1.0;
    }

    *opt_idx_s1 = opt_i;
    *opt_idx_s2 = opt_j;
    return min_time;
}

double max_time(double B, double *power_w, double *power_0, double max_speed, double granularity, int nspeeds, double *speeds, double *T1, double *T2, int *opt_idx_s1, int *opt_idx_s2)
{
    double c1 = 2.8037;
    double c2 = 0.3177;
    double g = 9.8;
    double W = 1.5;

    int opt_i = -1, opt_j = -1;
    double max_time = 0.0;

    double p_h = (c1 + c2) * pow((W * g), 1.5);

    int max_speed_idx = (max_speed - 1.0) / granularity;

    for (int i = 0; i < max_speed_idx; i++)
    {
        for (int j = 0; j < max_speed_idx; j++)
        {
            double e_tmp = T1[i] * power_w[i] + T2[j] * power_0[j];
            if (T1[i] * power_w[i] + T2[j] * power_0[j] > B)
                continue;
            double t_hover = (B - e_tmp) / p_h;
            if (T1[i] + T2[j] + t_hover > max_time)
            {
                max_time = T1[i] + T2[j] + t_hover;
                opt_i = i;
                opt_j = j;
            }
        }
    }
    if (opt_i == -1)
    {
        return -1.0;
    }

    *opt_idx_s1 = opt_i;
    *opt_idx_s2 = opt_j;
    return max_time;
}

void setRunNumber(instance *inst)
{
    char path[1000];
    sprintf(path, "../output/%s/seed_%d", inst->instance_name, inst->param.seed);
    inst->param.run = 1 + countDir(path);
}

int countDir(char *dir_path)
{
    struct dirent *dp;
    DIR *fd;

    if ((fd = opendir(dir_path)) == NULL)
    {
        print_error("Error in opendir\n");
        return 0;
    }
    int count = 0;
    while ((dp = readdir(fd)) != NULL)
    {
        if (!strcmp(dp->d_name, ".") || !strcmp(dp->d_name, ".."))
            continue; /* skip self and parent */
        count++;
    }
    closedir(fd);
    return count;
}