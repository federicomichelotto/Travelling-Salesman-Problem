#ifndef UTILS_H
#define UTILS_H

#include "fstsp_vds.h"

void parse_command_line(int argc, char **argv, instance *inst);
void parse_instance(instance *inst);
void parse_locations(instance *inst);
void parse_truck_travel_data(instance *inst);
void allocate_mem_truck_times_and_dists(instance *inst);
void allocate_mem_drone_min_max_times(instance *inst);
void parse_min_max_drone_legs_times(instance *inst);

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
void save_and_plot_solution_general(instance *inst, int *truck_succ, int *drone_succ, int iter);

int generate_path(char *path, char *folder, char *type, const char *model, char *filename, int seed, char *extension);
int generate_csv_record(char *path, char *instance_name, int seed, int model, int run, double z_best, double time_elapsed, int ticks);

void createInstanceFolders(instance *inst);
int IsPathExist(const char *s);
void getTimeStamp(double *ts);

void compute_opt_speeds(instance *inst);
void cruise_power_consumption(double w, int nspeed, double *speeds, double *P);
double vertical_power_consumption(double w, double speed);
double min_time(double B, double *power_w, double *power_0, double min_speed, double granularity, int nspeed, double *speeds, double *T1, double *T2, int *opt_idx_s1, int *opt_idx_s2);
double max_time(double B, double *power_w, double *power_0, double max_speed, double granularity, int nspeed, double *speeds, double *T1, double *T2, int *opt_idx_s1, int *opt_idx_s2);

void setRunNumber(instance *inst);
int countDir(char *dir_path);

#endif //UTILS_H
