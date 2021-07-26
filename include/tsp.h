#ifndef TSP_H
#define TSP_H

#include <cplex.h>
#include <concorde.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <float.h>


// Data Structures
typedef struct
{
    char input_file[1000];  // Path of the file
    char name[100];         // Identifies the data file
    char type[10];          // Specifies the type of data (TSP or ATSP)
    char comment[1000];     // Comment of the problem
    char weight_type[10];   // Specifies how the edge weights (or distances) are given
    char weight_format[20]; // Specifies how the edge weights (or distances) are formatted
    char data_type[20];     // Specifies how the data are displayed

    int seed;             // Seed given to cplex
    int run;              // Number of current run
    int verbose;          // verbosity level [0,4]
    int callback_counter; // #callback's called
    int ticks;            // flag: 0->seconds, 1->ticks
    int solver;           // 0 : optimal, 1 : math, 2 : heuristic
    int interactive;      // 0: none plot is displayed, 1: all plots are printed
    int saveplots;        // 0: save only the final plot, 1: save all plots
    int grasp;            // 0: GRASP OFF, 1 : GRASP ON
    int grasp_choices;    // # of possible random choices

} parameter;

typedef struct
{ // Node

    int id;   // number of the node (e.g. 1, 2, 3, ..., n)
    double x; // x coordinate
    double y; // y coordinate
    int flag; // Inside circuit or not

} node;

typedef struct
{ // Edge in the circuit

    double dist; // Weight of the edge
    int prev;    // Starting node
    int next;    // Ending node
    int flag;    // Inside circuit or not

} edge;

typedef struct
{

    // Input data
    int dimension; // Number of nodes of the problem
    node *nodes;   // List of nodes
    int *succ;     // array of nodes' successors

    double *weights;     // weights
    _Bool integer_costs; // = TRUE (1) for integer costs (rounded distances), = FALSE (0) otherwise

    parameter param; // Parameters of the instance

    double time_limit; // Specifies the maximum time allowed within the execution
    double timestamp_start;
    double timestamp_finish;
    int model_type;   // Specifies the type of the model
    double z_best;    // Value of the best solution available (incumbent)
    double best_lb;   // best lower bound
    double *best_sol; // Best xstar found
    int cols;         // #columns in cplex
    FILE *gnuplotPipe;

} instance;

typedef struct {
    int *chromosome;
    double fitness;
} population;

typedef struct
{
    instance *inst;
    CPXCALLBACKCONTEXTptr context;
} doit_fn_input;

// Enumerations
enum verbose_level
{
    QUIET = 0,
    NORMAL = 1,
    VERBOSE = 2,
    NERD = 3,
    DEBUG = 4
};

enum sections
{
    PARAMETERS = 0,
    NODE_COORD = 1,
    EDGE_WEIGHT = 2,
    DISPLAY_DATA = 3
};

static const char *verbose_name[] = {
    "QUIET",
    "NORMAL",
    "VERBOSE",
    "NERD",
    "DEBUG"};

static const char *optimal_model_name[] = {
    "STD",
    "MTZ",
    "MTZMOD",
    "MTZL",
    "MTZLS",
    "GG",
    "GGL",
    "GGLS",
    "GG_ORIGINAL",
    "BENDERS",
    "BENDERS (KRUSKAL)",
    "CALLBACK"};

static const char *math_model_name[] = {
    "HARD FIXING HEURISTIC",
    "SOFT FIXING HEURISTIC"};

static const char *heuristic_model_name[] = {
    "GREEDY",
    "EXTRA-MILEAGE",
    "EXTRA-MILEAGE FURTHEST STARTING NODES",
    "TABU SEARCH",
    "GENETIC"};

static const char *optimal_model_full_name[] = {
    "Basic model w/o SEC",
    "Miller-Tucker-Zemlin compact model",
    "Miller-Tucker-Zemlin modified compact model",
    "Miller-Tucker-Zemlin lazy compact model",
    "Miller-Tucker-Zemlin lazy compact model w/ SEC2",
    "Garvish-Graves compact model",
    "Garvish-Graves lazy compact model",
    "Garvish-Graves lazy compact model w/ SEC2",
    "Garvish-Graves compact model original",
    "Benders' method",
    "Benders' method (Kruskal)",
    "Callback method"};

static const char *math_model_full_name[] = {
    "Hard fixing heuristic",
    "Soft fixing heuristic"};

static const char *heuristic_model_full_name[] = {
    "Nearest Neighbors (Greedy)",
    "Extra-mileage",
    "Extra-mileage furthest starting nodes",
    "Tabu Search",
    "Genetic"};


// TSP solver
int optimal_solver(instance *inst);
int math_solver(instance *inst);
int heuristic_solver(instance *inst);

// Exact model builder
void build_model(CPXENVptr env, CPXLPptr lp, instance *inst);

// Exact model prototypes

// Basic model (undirected graphs)
// model 0: basic model (no SEC)
void basic_model_no_sec(CPXENVptr env, CPXLPptr lp, instance *inst);

void basic_model_directed(CPXENVptr env, CPXLPptr lp, instance *inst);

// model 9-10: benders model (SEC)
void benders(CPXENVptr env, CPXLPptr lp, instance *inst);

// Compact model (directed graphs)
// model 1: TMZ_static
void MTZ_static(CPXENVptr env, CPXLPptr lp, instance *inst);

// model 2: TMZ_static_mod
void MTZ_static_mod(CPXENVptr env, CPXLPptr lp, instance *inst);

// model 3: TMZ_lazy
void MTZ_lazy(CPXENVptr env, CPXLPptr lp, instance *inst);

// model 4: MTZ_lazy_sec
void MTZ_lazy_sec(CPXENVptr env, CPXLPptr lp, instance *inst);

// model 5: GG
void GG(CPXENVptr env, CPXLPptr lp, instance *inst);

// model 6: GGL
void GG_lazy(CPXENVptr env, CPXLPptr lp, instance *inst);

// model 7: GGLS
void GG_lazy_sec(CPXENVptr env, CPXLPptr lp, instance *inst);

// model 8: GG_ORIGINAL
void GG_original(CPXENVptr env, CPXLPptr lp, instance *inst);

// MATH HEURISTIC
// model 0: Hard fixing heuristic
void hard_fixing_heuristic(CPXENVptr env, CPXLPptr lp, instance *inst, int time_limit_iter, double fix_ratio);

// model 1: Soft fixing heuristic
void soft_fixing_heuristic(CPXENVptr env, CPXLPptr lp, instance *inst, int time_limit_iter);

// CONSTRUCTIVE HEURISTIC
double nearest_neighbours(instance *inst, int starting_node, int *succ, int options);

double extra_mileage(instance *inst, int *succ, int starting_node);
double extra_mileage_furthest_starting_nodes(instance *inst, int *succ);
int add_node_extra_mileage(instance *inst, int *succ, int node);

int nearest_insertion(instance *inst, int n, node *node_list, double random_number);

int farthest_insertion(instance *inst, int n, node *node_list, double random_number);

// REFINEMENT HEURISTIC
// to use to refine a solution, we assume that inside inst->best_sol there is a valid solution, and the selected edges are in inst->edges
int two_opt(instance *inst, int maxMoves);
int two_opt_v2(instance *inst);
int reverse_successors(int *succ, int size, int start, int end);

//META HEURISTIC
void tabu_search(instance *inst);
double genetic(instance *inst, int size, int epoch);
double random_initialization(instance *inst, int *individual, int seed, int optimize);
void rank_individuals(population *individuals, int size);

// Mutation
double swap_mutation(instance *inst, int *chromosome, int i, int j);
double scramble_mutation(instance *inst, int *chromosome, int i, int j);
double inversion_mutation(instance *inst, int *chromosome, int i, int j);

// Parent Selection
void roulette_wheel_selection(population *individuals, int size, int *selection); // Fitness Proportionate
void tournament_selection(population *individuals, int k, int size, int *selection); // Tournament Selection
void rank_selection(population *individuals, int size, int *selection); // Rank Selection
void random_selection(population *individuals, int size, int *selection); // Random Selection

// Crossover
void one_point_crossover(instance *inst, population parentA, population parentB, population offspring);
void uniform_crossover(instance *inst, population parentA, population parentB, population offspring);

double fitness_average(population *individuals, int size);
double fitness_variability(population *individual, int size);

void convert_to_successor(int *list, int *succ, int n);
void convert_to_list(int *succ, int *list, int n);

// Some useful functions
double gather_solution(instance *inst, const double *xstar, int type);

// Retrieve the distance among each node of the instance
double dist(int i, int j, instance *inst);

double compute_edges(instance *inst, const double *xstar, int type);

// Find the connected components inside a solution
void findConnectedComponents(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp, int **length_comp);

void findConnectedComponents_kruskal(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp,
                                     int **length_comp);

// Callback

static int CPXPUBLIC callback_driver(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle);

static int CPXPUBLIC callback_candidate(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle);

static int CPXPUBLIC callback_relaxation(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle);

int doit_fn_concorde(double cutval, int cutcount, int *cut, void *void_context);

// Retrieve the position of the variable
int xpos(int i, int j, instance *inst);     // position in the model for undirected graphs
int xpos_dir(int i, int j, instance *inst); // position in the model for directed graphs
int upos(int i, instance *inst);            // position in the model of i-th u-variable
int ypos(int i, int j, instance *inst);     // position in the model of y-variable for the arc (i,j)



#endif //TSP_H
