#ifndef UTILS_H
#define UTILS_H

#include "tsp.h"

void parse_command_line(int argc, char** argv, instance *inst);
void parse_instance(instance *inst);

void print_command_line(instance *inst);
void print_instance(instance *inst);

void check_format(char *param);
void initialize_instance(instance *inst);
void free_instance(instance *inst);

void print_help();
void print_error(const char *err);
void print_error_status(const char *err, int e);
void print_message(const char *msg);
int plot_optimal_solution(instance *inst);
int plot_intermediate_solution(instance *inst, int update_plot);
int plot_solution_edges(int n_edges, node *nodes, edge *edges);

int generate_path(char *path, char *folder, char *type, const char *model, char *filename, int seed, char *extension);
int generate_csv_record(char *instance_name, int seed, const char *model, double z_best, double time_elapsed, int run);

#endif //UTILS_H
