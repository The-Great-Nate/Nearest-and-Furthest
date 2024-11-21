/////////////////////////////////////////////
// This code has been modified from SCIF20002 Workshop 16
/////////////////////////////////////////////
#include "rng.hpp"

// Initialise the random number generator:
std::default_random_engine generator;

/**
 * Random number generator initialiser for generating random points
 */
void init_rng()
{
    // Initialise a clock object
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();

    // Obtain a seed from the timer and apply it
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    generator.seed(seed); // Apply the seed

}    

/**
 * Generates random number from initialised generator.
 * - Is seperate as generator should only be initialised once ideally.
 */
double rng()
{
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    return distribution(generator);
}