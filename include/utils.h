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
void print_message(const char *msg);
int plot_solution(instance *inst);
int plot_solution_edges(int n_edges, node *nodes, edge *edges);

void generate_path(char *path, char *folder, char *type, const char *model, char *filename, int seed, char *extension);
void generate_csv_record(char *instance_name, int seed, int model_type, double z_best, long int time_sec, long int time_usec);

#endif //UTILS_H
