#ifndef GENETIC_H
#define GENETIC_H

#include <cplex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <float.h>

double compute_sequence_min_time(int *seq, instance *inst, int *truck_succ, int *drone_succ);
void generate_random_sequence(int *rand_seq, instance *inst);
double opt_drone_legs_direction(int *truck_succ, int *drone_succ, instance *inst);
double reverse_drone_leg(int *truck_seq, int *drone_seq, instance *inst);
void nearest_neighbours(instance *inst, int *seq, int options);

double compute_score(int *truck_seq, int *drone_seq, instance *inst);

#endif //GENETIC_H