#include <iostream>
#include "Vector2D.cpp"
#include "rng.cpp"

double distance_sq, closest_sq, closest_wrap_sq, closest;

// Generate random points on each slot in points vector
void make_points(std::vector<Vector2D> &points, int size)
{
    for (int i = 0; i < size; ++i)
    {
        points.push_back(Vector2D(rng(), rng())); // Add point to points vector based on rng
    }
}

void find_nearest_neighbour(std::vector<Vector2D> &points, int size)
{
    for (int i = 0; i < size; ++i)
    {
        distance_sq = 0;
        closest_sq = std::numeric_limits<double>::infinity();
        closest_wrap_sq = std::numeric_limits<double>::infinity();
        for (int j = 0; j < size; ++j)
        {
            if (i != j)
            {
                // Find distance between 2 points
                double dx = std::abs(points[j].GetX() - points[i].GetX());
                double dy = std::abs(points[j].GetY() - points[i].GetY());
                double wrap_dx, wrap_dy
                // Checking if wrapping around is worth it.
                // If x or y is over 0.5 apart, wrap it around
                if (dx > 0.5)
                    dx = 1.0f - dx;

                if (dy > 0.5)
                    dy = 1.0f - dy;

                // Calculate the final distance with pythagoras
                distance_sq = dx * dx + dy * dy;

                if (distance_sq <= closest_sq)
                {
                    closest_sq = distance_sq;
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

    // for debugging :D
    for (int i = 0; i < size; ++i)
    {
        std::cout << "Point: (" << points[i].GetX() << ", " << points[i].GetY() << ")" << std::endl;
    }
    return EXIT_SUCCESS;
}