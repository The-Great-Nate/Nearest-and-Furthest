#include <iostream>
#include <fstream>
#include <chrono>
#include "Vector2D.cpp"
#include "rng.cpp"
#include <omp.h>

double distance_sq, nearest_sq, nearest_wrap_sq, nearest, furthest_sq, furthest_wrap_sq, furthest;

/**
 * Generate random points using rng()
 * @param[in] std::Vector of Vector2D points
 * @param[in] int of size of vector input from command line
 */
void make_points(std::vector<Vector2D> &points, int size)
{
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (int i = 0; i < size; ++i)
        {
            points.push_back(Vector2D(rng(), rng())); // Add point to points vector based on rng
        }
    }
}

/**
 * Finds nearest and furthest distances in vector of points in standard geometry
 * @param[in] std::Vector of Vector2D points
 * @param[in] int of size of vector input from command line
 * @param[in] std::ofstream file object storing the nearest distances
 * @param[in] std::ofstream file object storing the furthest distances
 */
void find_standard_distances(std::vector<Vector2D> &points, int size, std::ofstream &near_obj, std::ofstream &far_obj)
{
    // Vector that stores nearest and furthest distances to points. Only used for outputting average
    std::vector<double> nearests, furthests;
    nearests.reserve(size * 2);
    furthests.reserve(size * 2);

// Iterating over every point
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (int i = 0; i < size; ++i)
        {
            // Initialise nearest and furthest distances for each point
            distance_sq = 0;
            nearest_sq = std::numeric_limits<double>::infinity();
            furthest_sq = 0;

            // Iterating over adjacent points
            for (int j = 0; j < size; ++j)
            {
                // If statement to prevent points from tracking distances with themselves.
                if (i != j)
                {
                    // Find distance between 2 points
                    double dx = std::abs(points[j].GetX() - points[i].GetX());
                    double dy = std::abs(points[j].GetY() - points[i].GetY());

                    // Calculate the final distance with pythagoras
                    distance_sq = dx * dx + dy * dy;

                    // Check if distance squared is smaller than current stored nearest distance
                    if (distance_sq <= nearest_sq)
                    {
                        nearest_sq = distance_sq;
                        points[i].SetWrapNeighbour(&points[j]);
                    }

                    // Check if distance squared is bigger than current stored furthest distance
                    if (distance_sq >= furthest_sq)
                    {
                        furthest_sq = distance_sq;
                        points[i].SetWrapFaraway(&points[j]);
                    }
                }
            }

            // Store nearest and furthest distances into respective vectors of distances
            nearest = std::sqrt(nearest_sq);
            furthest = std::sqrt(furthest_sq);
            nearests.push_back(nearest);
            furthests.push_back(furthest);

            // Store each distance in their respective files
            if (near_obj.is_open())
            {
                near_obj << nearest << "\n";
            }

            if (far_obj.is_open())
            {
                far_obj << furthest << "\n";
            }
        }
    }

    // Calculate the average furthest and nearest distance.
    double avrg_nearest = 0;
    double avrg_furthest = 0;

#pragma omp parallel
    {
// Summing all distances
#pragma omp for schedule(static)
        for (int i = 0; i < size * 2; i++)
        {
            avrg_nearest += nearests[i];
            avrg_furthest += furthests[i];
        }
    }

    // Averaging distances by size of vectors (which is size * 2 due to double calculations of distances sadly)
    avrg_nearest /= size * 2;
    avrg_furthest /= size * 2;

    // Output averages
    std::cout << "Average Nearest Std' Distance: " << avrg_nearest << "\n"
              << "Average Furthest Std' Distance: " << avrg_furthest << std::endl;
}

/**
 * Finds nearest and furthest distances in vector of points in wrap around geometry
 * @param[in] std::Vector of Vector2D points
 * @param[in] int of size of vector input from command line
 * @param[in] std::ofstream file object storing the nearest distances
 * @param[in] std::ofstream file object storing the furthest distances
 */
void find_wrapped_distances(std::vector<Vector2D> &points, int size, std::ofstream &near_obj, std::ofstream &far_obj)
{
    // Vector that stores nearest and furthest distances to points. Only used for outputting average
    std::vector<double> nearests, furthests;
    nearests.reserve(size * 2);
    furthests.reserve(size * 2);

#pragma omp parallel
    {
#pragma omp for schedule(static)
        // Iterating over every point
        for (int i = 0; i < size; ++i)
        {
            // Initialise nearest and furthest distances for each point
            distance_sq = 0;
            nearest_sq = std::numeric_limits<double>::infinity();
            furthest_sq = 0;

            // Iterating over adjacent points
            for (int j = 0; j < size; ++j)
            {
                // If statement to prevent points from tracking distances with themselves.
                if (i != j)
                {
                    // Find distance between 2 points
                    double dx = std::abs(points[j].GetX() - points[i].GetX());
                    double dy = std::abs(points[j].GetY() - points[i].GetY());

                    // Checking if wrapping around is worth it.
                    // If x or y is over 0.5 apart, wrap it around
                    if (dx > 0.5)
                        dx = 1.0f - dx;

                    if (dy > 0.5)
                        dy = 1.0f - dy;

                    // Calculate the final distance with pythagoras
                    distance_sq = dx * dx + dy * dy;

                    // Check if distance squared is smaller than current stored nearest distance
                    if (distance_sq <= nearest_sq)
                    {
                        nearest_sq = distance_sq;
                        points[i].SetStdNeighbour(&points[j]);
                    }

                    // Check if distance squared is bigger than current stored furthest distance
                    if (distance_sq >= furthest_sq)
                    {
                        furthest_sq = distance_sq;
                        points[i].SetStdFaraway(&points[j]);
                    }
                }
            }

            // Store nearest and furthest distances into respective vectors of distances
            nearest = std::sqrt(nearest_sq);
            furthest = std::sqrt(furthest_sq);
            nearests.push_back(nearest);
            furthests.push_back(furthest);

            // Store each distance in their respective files
            if (near_obj.is_open())
            {
                near_obj << nearest << "\n";
            }

            if (far_obj.is_open())
            {
                far_obj << furthest << "\n";
            }
        }
    }

    // Calculate the average furthest and nearest distance.
    double avrg_nearest = 0;
    double avrg_furthest = 0;

#pragma omp parallel
    {
// Summing all distances
#pragma omp for schedule(static)
        for (int i = 0; i < size * 2; i++)
        {
            avrg_nearest += nearests[i];
            avrg_furthest += furthests[i];
        }
    }

    // Averaging distances by size of vectors (which is size * 2 due to double calculations of distances sadly)
    avrg_nearest /= size * 2;
    avrg_furthest /= size * 2;

    // Output averages
    std::cout << "Average Nearest Std' Distance: " << avrg_nearest << "\n"
              << "Average Furthest Std' Distance: " << avrg_furthest << std::endl;
}

// main
int main(int argc, char **argv)
{
    // checks if user provided enough arguments in command line
    if (argc < 3)
    {
        std::cerr << "Please provide the number of points & number of threads to use as a command line argument." << std::endl;
        return EXIT_FAILURE;
    }

    // Sets size of points vector based on command line argument
    int size = atoi(argv[1]);
    std::vector<Vector2D> points;
    points.reserve(size); // Reserve space for vector in memory

    // Assign no. of threads to use based on command line argument
    int numThreads = atoi(argv[2]);
    omp_set_num_threads(numThreads);

    // Start generating points and benchmark performance
    auto start_generating = std::chrono::high_resolution_clock::now();
    init_rng();                // Initialise rng
    make_points(points, size); // Generate the points
    auto stop_generating = std::chrono::high_resolution_clock::now();
    auto duration_generating = std::chrono::duration_cast<std::chrono::microseconds>(stop_generating - start_generating);

    std::ofstream nearest_std_file, furthest_std_file, nearest_wrap_file, furthest_wrap_file;
    nearest_std_file.open("data\\nearest_std_distances.txt");
    furthest_std_file.open("data\\furthest_std_distances.txt");
    nearest_wrap_file.open("data\\nearest_wrapped_distances.txt");
    furthest_wrap_file.open("data\\furthest_wrapped_distances.txt");

    // for debugging :D
    /*
    for (int i = 0; i < size; ++i)
    {
        std::cout << "Point: (" << points[i].GetX() << ", " << points[i].GetY() << ")" << std::endl;
    }
    */

    // Start the clock for how long this simulation takes to run
    auto start_nearest = std::chrono::high_resolution_clock::now();
    find_standard_distances(points, size, nearest_std_file, furthest_std_file);
    find_wrapped_distances(points, size, nearest_wrap_file, furthest_wrap_file);

    // Stop the clock for how long this simulation takes to run
    auto end_nearest = std::chrono::high_resolution_clock::now();
    auto duration_nearest = std::chrono::duration_cast<std::chrono::seconds>(end_nearest - start_nearest);
    auto minutes_nearest = std::chrono::duration_cast<std::chrono::minutes>(duration_nearest);
    auto seconds_nearest = duration_nearest - minutes_nearest;

    std::cout << "------------------------------------------------" << std::endl;

    // Output time to generate points and general simulation run time to console
    std::cout << "Generation Runtime " << duration_generating.count() << " microseconds \n";
    std::cout << "Simulation Runtime " << minutes_nearest.count() << ":" << seconds_nearest.count() << " (MM:ss) \n";

    std::ofstream gen_times("runtimes\\generation.txt", std::ios::app);
    std::ofstream static_times("runtimes\\simulation_static.txt", std::ios::app);
    // std::ofstream parallel_times("runtimes\\simulation_parallel.txt", std::ios::app);

    gen_times << "threads_used\tduration(microseconds)" << "\n"
              << numThreads << "\t" << duration_generating.count();

    static_times << "threads_used\tduration(seconds)" << "\n"
                 << numThreads << "\t" << duration_nearest.count();

    // Close files
    nearest_std_file.close();
    furthest_std_file.close();

    return EXIT_SUCCESS;
}