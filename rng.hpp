#ifndef RNG_H
#define RNG_H

#include <random>
#include <chrono>

extern std::default_random_engine generator;

void init_rng();
double rng();

#endif