#ifndef TSP_H
#define TSP_H

#include <cplex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Data Structures
typedef struct {
    char input_file[1000];		// Path of the file
    char name[100];				// Identifies the data file
    char type[10];				// Specifies the type of data (TSP or ATSP)
    char comment[1000];			// Comment of the problem
    char weight_type[10];		// Specifies how the edge weights (or distances) are given
    char weight_format[20];     // Specifies how the edge weights (or distances) are formatted
    char data_type[20];         // Specifies how the data are displayed
} parameter;

typedef struct {

    // Input data
    int nodes;                  // Number of nodes of the problem
    double *x;                  // x coordinate
    double *y;                  // y coordinate
    double *weights;            // weights
    _Bool integer_costs;        // = TRUE (1) for integer costs (rounded distances), = FALSE (0) otherwise

    parameter param;            // Parameters of the instance

    double time_limit;          // Specifies the maximum time allowed within the execution
    int model_type;				// Specifies the type of the model
    double z_best;              // Value of the best solution available

} instance;

// Other possible structures

typedef struct {			    // Node

    _Bool flag;			        // = TRUE (1) if inside the circuit, = FALSE (0) otherwise
    int id;                     // number of the node (e.g. 1, 2, 3, ..., n)
    double x;			        // x coordinate
    double y;                   // y coordinate

} node;

typedef struct {			    // Edge in the circuit

    _Bool flag;			        // = TRUE (1) if inside the circuit, = FALSE (0) otherwise
    double dist;			    // Weight of the edge
    node prev;                  // Starting node
    node next;                  // Ending node

} edge;

enum verbose_level {
    QUIET = 0,
    NORMAL = 1,
    VERBOSE = 2,
    NERD = 3,
    DEBUG = 4
};

enum sections {
    PARAMETERS = 0,
    NODE_COORD = 1,
    EDGE_WEIGHT = 2,
    DISPLAY_DATA = 3
};

static const char* verbose_name[] = {"QUIET", "NORMAL", "VERBOSE", "NERD", "DEBUG"};
static int verbose = NORMAL;

int TSPopt(instance *inst);

void build_model(CPXENVptr env, CPXLPptr lp, instance *inst);


double dist(int i, int j, instance *inst);
int position(int i, int j, instance *inst);


#endif //TSP_H
