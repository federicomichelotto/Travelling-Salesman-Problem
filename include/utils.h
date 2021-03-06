#ifndef UTILS_H
#define UTILS_H

#include "tsp.h"

void parse_command_line(int argc, char** argv, instance *inst);
void parse_instance(instance *inst);

void print_command_line(instance *inst);
void print_instance(instance *inst);

void initialize_instance(instance *inst);
void free_instance(instance *inst);

void print_help();
void print_error(const char *err);


#endif //UTILS_H
