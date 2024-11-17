#include <iostream>
#include "Vector2D.cpp"
#include "rng.cpp"

// Generate random points on each slot in points vector
void make_points(std::vector<Vector2D> &points, int size)
{
    for (int i = 0; i < size; ++i)
    {
        points.push_back(Vector2D(rng(), rng())); // Add point to points vector based on rng
    }
}

// main
int main(int argc, char **argv)
{
    // checks if user provided enough arguments in command line
    if (argc < 2) {
        std::cerr << "Please provide the number of points as a command line argument." << std::endl;
        return EXIT_FAILURE;
    }

    // Sets size of points vector based on command line argument
    int size = atoi(argv[1]);
    std::vector<Vector2D> points;
    points.reserve(size); // Reserve space for vector in memory

    init_rng(); // Initialise rng
    make_points(points, size); // Generate the points

    // for debugging :D
    for (int i = 0; i < size; ++i)
    {
        std::cout << "Point: (" << points[i].GetX() << ", " << points[i].GetY() << ")" << std::endl;
    }
    return EXIT_SUCCESS;
}