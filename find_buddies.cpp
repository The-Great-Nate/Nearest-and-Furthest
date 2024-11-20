#include <iostream>
#include <fstream>
#include <chrono>
#include "Vector2D.cpp"
#include "rng.cpp"

double distance_sq, nearest_sq, nearest_wrap_sq, nearest, furthest_sq, furthest_wrap_sq, furthest;

// Generate random points on each slot in points vector
void make_points(std::vector<Vector2D> &points, int size)
{
    for (int i = 0; i < size; ++i)
    {
        points.push_back(Vector2D(rng(), rng())); // Add point to points vector based on rng
    }
}

void find_standard_nearest(std::vector<Vector2D> &points, int size, std::ofstream &near_obj, std::ofstream &far_obj)
{
    std::vector<double> nearests, furthests;
    nearests.reserve(size * 2);
    furthests.reserve(size * 2);

    for (int i = 0; i < size; ++i)
    {
        distance_sq = 0;
        nearest_sq = std::numeric_limits<double>::infinity();
        furthest_sq = 0;
        for (int j = 0; j < size; ++j)
        {
            if (i != j)
            {
                // Find distance between 2 points
                double dx = std::abs(points[j].GetX() - points[i].GetX());
                double dy = std::abs(points[j].GetY() - points[i].GetY());

                // Calculate the final distance with pythagoras
                distance_sq = dx * dx + dy * dy;

                if (distance_sq <= nearest_sq)
                {
                    nearest_sq = distance_sq;
                    points[i].SetStdNeighbour(&points[j]);
                }

                if (distance_sq >= furthest_sq)
                {
                    furthest_sq = distance_sq;
                    points[i].SetStdFaraway(&points[j]);
                }
            }
        }

        nearest = std::sqrt(nearest_sq);
        furthest = std::sqrt(furthest_sq);
        nearests.push_back(nearest);
        furthests.push_back(furthest);

        if (near_obj.is_open())
        {
            near_obj << nearest << "\n";
        }

        if (far_obj.is_open())
        {
            far_obj << furthest << "\n";
        }
    }

    double avrg_nearest = 0;
    double avrg_furthest = 0;
    for (int i = 0; i < size * 2; i++)
    {
        avrg_nearest += nearests[i];
        avrg_furthest += furthests[i];
    }
    avrg_nearest /= size * 2;
    avrg_furthest /= size * 2;

    std::cout << "Average Nearest Std' Distance: " << avrg_nearest << "\n"
              << "Average Furthest Std' Distance: " << avrg_furthest << std::endl;
}

void find_wrapped_nearest(std::vector<Vector2D> &points, int size)
{
    for (int i = 0; i < size; ++i)
    {
        distance_sq = 0;
        nearest_sq = std::numeric_limits<double>::infinity();
        nearest_wrap_sq = std::numeric_limits<double>::infinity();
        for (int j = 0; j < size; ++j)
        {
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

                if (distance_sq <= nearest_sq)
                {
                    nearest_sq = distance_sq;
                    points[i].SetStdNeighbour(&points[j]);
                    points[i].SetWrapNeighbour(&points[j]);
                }
            }
        }
    }
}

// main
int main(int argc, char **argv)
{
    // checks if user provided enough arguments in command line
    if (argc < 2)
    {
        std::cerr << "Please provide the number of points as a command line argument." << std::endl;
        return EXIT_FAILURE;
    }

    // Sets size of points vector based on command line argument
    int size = atoi(argv[1]);
    std::vector<Vector2D> points;
    points.reserve(size); // Reserve space for vector in memory

    init_rng();                // Initialise rng
    make_points(points, size); // Generate the points

    std::ofstream nearest_std_file, furthest_std_file, nearest_wrap_file, furthest_wrap_file;
    nearest_std_file.open("data\\nearest_distances.txt");
    furthest_std_file.open("data\\furthest_distances.txt");

    // for debugging :D
    /*
    for (int i = 0; i < size; ++i)
    {
        std::cout << "Point: (" << points[i].GetX() << ", " << points[i].GetY() << ")" << std::endl;
    }
    */

    // Start the clock for how long this simulation takes to run
    auto start_nearest = std::chrono::high_resolution_clock::now();
    find_standard_nearest(points, size, nearest_std_file, furthest_std_file);

    // Stop the clock for how long this simulation takes to run
    auto end_nearest = std::chrono::high_resolution_clock::now();
    auto duration_nearest = std::chrono::duration_cast<std::chrono::seconds>(end_nearest - start_nearest);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration_nearest);
    auto seconds = duration_nearest - minutes;

    std::cout << "Simulation Runtime " << minutes.count() << ":" << seconds.count() << " (MM:ss) \n";
    nearest_std_file << "duration=" << minutes.count() << ":" << seconds.count() << "\n";
    furthest_std_file << "duration=" << minutes.count() << ":" << seconds.count() << "\n";

    nearest_std_file.close();
    furthest_std_file.close();

    return EXIT_SUCCESS;
}