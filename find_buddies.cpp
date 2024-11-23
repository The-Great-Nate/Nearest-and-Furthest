#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include "Vector2D.cpp"
#include "rng.cpp"
#include <omp.h>

/**
 * Check if inputted argument is a file_path (string) or requested number of points to generate (int).
 * @param[in] std::string command line argument
 */
bool is_path(const std::string &str)
{
    for (char c : str)
    {
        if (!std::isdigit(c))
        {
            return false;
        }
    }
    return true;
}

/**
 * Read in generated points using two column csv file of points
 * @param[in] std::Vector of Vector2D points
 * @param[in] std::ifstream* pointer to input file stream, initialised to nullptr. Used for pre_generated points
 */
std::vector<Vector2D> make_points(std::string location)
{
    std::vector<Vector2D> points;

    // Open the file for reading
    std::ifstream file(location);

    // Check if file was successfully opened
    if (!file.is_open())
    {
        // Output error message if file can't be opened
        std::cerr << "Error opening file!" << std::endl;
    }

    // Create vector of lines in csv & output string for splitting (via std::getline)
    std::vector<std::string> lines;
    std::string line;

    // Read the file sequentially and store lines in a vector
    while (std::getline(file, line))
    {
        // Add split line to lines vector
        lines.push_back(line);
    }

    // Close file
    file.close();

    // Set size variable to size of lines vector. Used for reserving memory in subsequent vectors
    int size = lines.size();

// Parallelise the reading of CSV data using OpenMP
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (int i = 0; i < size; ++i)
        {
            // Set line to std::stringstream object to be split
            std::stringstream to_split(lines[i]);

            // Initialise column variables as strings
            std::string col1, col2;

            // Check if splitting operation was successful
            if (std::getline(to_split, col1, ',') && std::getline(to_split, col2, ','))
            {
                // Convert the first and second column values to double
                double x = std::stod(col1);
                double y = std::stod(col2);

                // Check if the values are within accepted range (0, 1)
                if (x >= 0.0 && x <= 1.0 && y >= 0.0 && y <= 1.0)
                {
// If valid, assign them to points
#pragma omp critical
                    {
                        points.push_back(Vector2D(x, y));
                    }
                }
                else
                {
// Handle invalid values
#pragma omp critical
                    {
                        // Output position of invalid value in csv
                        std::cerr << "Invalid values at line " << i + 1 << ": "
                                  << "x = " << x << ", y = " << y << "\n"
                                  << "Only values between 0 and 1 inclusive allowed!" << std::endl;
                    }
                }
            }
        }
    }
    return points;
}

/**
 * Generate random points using rng()
 * @param[in] std::Vector of Vector2D points
 * @param[in] int of size of vector input from command line
 */
std::vector<Vector2D> make_points(int size)
{
    std::vector<Vector2D> points; // Vector of points
#pragma omp parallel
    {
        std::vector<Vector2D> local_points;  // Local vector for each thread
#pragma omp parallel for schedule(static)
        for (int i = 0; i < size; ++i)
        {
            local_points.push_back(Vector2D(rng(), rng())); // Add point to points vector based on rng
        }

        #pragma omp critical
        {
            // Once thread completes work, append its local vector to the main vector
            points.insert(points.end(), local_points.begin(), local_points.end());
        }
    }
    return points;
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
     
    std::vector<double> nearests, furthests; // Vector that stores nearest and furthest distances to points. Only used for outputting average
    nearests.reserve(size * 2);
    furthests.reserve(size * 2);
#pragma omp parallel
    {
        // Vector that stores nearest and furthest distances to points. Only used for outputting average
        std::vector<double> local_nearests, local_furthests;
        int local_size = size / omp_get_num_threads() + 10; // Reserve space for the points handled by this thread (add 10 in case of rounding problems with odd no. threads)
        local_nearests.reserve(local_size);
        local_furthests.reserve(local_size);

#pragma omp for schedule(static)
        // Iterating over every point
        for (int i = 0; i < size; ++i)
        {
            // Initialise nearest and furthest distances for each point
            double distance_sq = 0.0;
            double nearest_sq = std::numeric_limits<double>::infinity();
            double furthest_sq = 0.0;
            double nearest = 0.0;
            double furthest = 0.0;

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
            local_nearests.push_back(nearest);
            local_furthests.push_back(furthest);

// Write distances safely into files.
#pragma omp critical
            {
                // Write to files
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

// Merge local vectors into global vectors after the parallel region
#pragma omp critical
        {
            nearests.insert(nearests.end(), local_nearests.begin(), local_nearests.end());
            furthests.insert(furthests.end(), local_furthests.begin(), local_furthests.end());
        }
    }

    // Calculate the average furthest and nearest distance.
    double avrg_nearest = 0;
    double avrg_furthest = 0;

#pragma omp parallel
    {
// Summing all distances
#pragma omp parallel for reduction(+ : avrg_nearest, avrg_furthest) schedule(static) // Privatise avrg_nearest & avrg_furthest and then combine private copies with + operation
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
    std::cout << "Average Nearest Wraparound' Distance: " << avrg_nearest << "\n"
              << "Average Furthest Wraparound' Distance: " << avrg_furthest << std::endl;
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
    std::vector<double> nearests, furthests; // Vector that stores nearest and furthest distances to points. Only used for outputting average
    nearests.reserve(size * 2);
    furthests.reserve(size * 2);
#pragma omp parallel
    {
        // Vector that stores nearest and furthest distances to points. Only used for outputting average
        std::vector<double> local_nearests, local_furthests;
        int local_size = size / omp_get_num_threads() + 10; // Reserve space for the points handled by this thread (add 10 in case of rounding problems with odd no. threads)
        local_nearests.reserve(local_size);
        local_furthests.reserve(local_size);

#pragma omp for schedule(static)
        // Iterating over every point
        for (int i = 0; i < size; ++i)
        {
            // Initialise nearest and furthest distances for each point
            double distance_sq = 0;
            double nearest_sq = std::numeric_limits<double>::infinity();
            double furthest_sq = 0;
            double nearest = 0.0;
            double furthest = 0.0;

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
            local_nearests.push_back(nearest);
            local_furthests.push_back(furthest);

// Write distances safely into files.
#pragma omp critical
            {
                // Write to files
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

// Merge local vectors into global vectors after the parallel region
#pragma omp critical
        {
            nearests.insert(nearests.end(), local_nearests.begin(), local_nearests.end());
            furthests.insert(furthests.end(), local_furthests.begin(), local_furthests.end());
        }
    }

    // Calculate the average furthest and nearest distance.
    double avrg_nearest = 0;
    double avrg_furthest = 0;

#pragma omp parallel
    {
// Summing all distances
#pragma omp parallel for reduction(+ : avrg_nearest, avrg_furthest) schedule(static)
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
        std::cerr << "Please provide the (number of points || file path) & number of threads to use as a command line argument." << std::endl;
        return EXIT_FAILURE;
    }

    // Initialise system parameters. Dependant on user cmd line input
    int size;
    std::vector<Vector2D> points;
    std::chrono::microseconds duration_generating;
    std::string gen_method;
      
    // Check if first cmd argument is file path or integer
    if (is_path(std::to_string(atoi(argv[1]))))
    { // If arg is int (not path)
        // Sets size of points vector based on command line argument
        size = atoi(argv[1]);
        // Start generating points and benchmark performance
        auto start_generating = std::chrono::high_resolution_clock::now();
        init_rng();                // Initialise rng
        points = make_points(size); // Generate the points
        auto stop_generating = std::chrono::high_resolution_clock::now();
        auto duration_generating = std::chrono::duration_cast<std::chrono::microseconds>(stop_generating - start_generating);
        gen_method = "random";
        
    }
    else
    {
        // Start reading file of points and benchmark performance
        auto start_generating = std::chrono::high_resolution_clock::now();
        points = make_points(std::to_string(atoi(argv[1])));
        auto stop_generating = std::chrono::high_resolution_clock::now();
        auto duration_generating = std::chrono::duration_cast<std::chrono::microseconds>(stop_generating - start_generating);
        gen_method = "file of size: " + std::to_string(points.size());
    }

    // Assign no. of threads to use based on command line argument
    int numThreads = atoi(argv[2]);
    omp_set_num_threads(numThreads);

    std::ofstream nearest_std_file, furthest_std_file, nearest_wrap_file, furthest_wrap_file;
    nearest_std_file.open("data\\nearest_std_distances.txt");
    furthest_std_file.open("data\\furthest_std_distances.txt");
    nearest_wrap_file.open("data\\nearest_wrapped_distances.txt");
    furthest_wrap_file.open("data\\furthest_wrapped_distances.txt");

    // for debugging :D (output points)
    /*
    for (int i = 0; i < size; ++i)
    {
        std::cout << "Point: (" << points[i].GetX() << ", " << points[i].GetY() << ")" << std::endl;
    }
    */

    std::cout << "----------  Thread Count: " << numThreads << "  ----------" << std::endl;
    // Start benchmark performance
    auto start_nearest = std::chrono::high_resolution_clock::now();

    // Run the distance calculations
    find_standard_distances(points, size, nearest_std_file, furthest_std_file);
    find_wrapped_distances(points, size, nearest_wrap_file, furthest_wrap_file);

    // Stop benchmark performance
    auto end_nearest = std::chrono::high_resolution_clock::now();
    auto duration_nearest = std::chrono::duration_cast<std::chrono::seconds>(end_nearest - start_nearest);
    auto minutes_nearest = std::chrono::duration_cast<std::chrono::minutes>(duration_nearest);
    auto seconds_nearest = duration_nearest - minutes_nearest;

    std::cout << "------------------------------------------" << std::endl;

    // Output time to generate points and general simulation run time to console
    std::cout << "Generation Runtime " << duration_generating.count() << " microseconds \n";
    std::cout << "Simulation Runtime " << minutes_nearest.count() << ":" << seconds_nearest.count() << " (MM:ss) \n";

    std::ofstream gen_times("runtimes\\generation.txt", std::ios::app);
    std::ofstream static_times("runtimes\\simulation_static.txt", std::ios::app);
    // std::ofstream parallel_times("runtimes\\simulation_parallel.txt", std::ios::app);

    gen_times << numThreads << "\t" << duration_generating.count() << std::endl;

    static_times << numThreads << "\t" << duration_nearest.count() << std::endl;

    // Close files
    nearest_std_file.close();
    furthest_std_file.close();
    gen_times.close();
    static_times.close();

    return EXIT_SUCCESS;
}