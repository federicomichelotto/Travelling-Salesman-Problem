#ifndef UTILS_H
#define UTILS_H

#include "tsp.h"

void parse_command_line(int argc, char **argv, instance *inst);
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

void save_and_plot_solution(instance *inst, int iter);
void save_and_plot_solution_general(instance *inst, int *succ, int iter);

int generate_path(char *path, char *folder, char *type, const char *model, char *filename, int seed, char *extension);
int generate_csv_record(char *instance_name, int seed, const char *model, double z_best, double time_elapsed, int run);

int checkFolders();
int IsPathExist(const char *s);
void getTimeStamp(double *ts);

#endif //UTILS_H
