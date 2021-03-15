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
int plot_solution(int num_nodes, double *x, double *y, int num_edges, int *u, int *v);

#endif //UTILS_H
