#include <iostream>
#include <cmath>

class Vector2D
{
public:
    Vector2D();                             // Default Constructor
    Vector2D(double a, double b);             // Parameterised Constructor
    double GetX() const;                     // Getter for x_
    double GetY() const;                     // Getter for y_
    void SetStdNeighbour(Vector2D *buddy);  // Setter for nearest std neighbour
    void SetWrapNeighbour(Vector2D *buddy); // Setter for nearest wrap neighbour

private:
    double x_; // X Co Ordinate
    double y_; // Y Co Ordinate
    Vector2D *std_neighbour_; // Pointer to the nearest standard neighbour - initially nothing hence nullptr in constructors
    Vector2D *wrap_neighbour_; // Pointer to the nearest wrapped neighbour - initially nothing hence nullptr in constructors
};

// Default Constructor
Vector2D::Vector2D() : x_(0.0),
                       y_(0.0),
                       std_neighbour_(nullptr),
                       wrap_neighbour_(nullptr)
{
}

// Parameterised Constructor
Vector2D::Vector2D(double a, double b) :
                                       x_(a),
                                       y_(b),
                                       std_neighbour_(nullptr),
                                       wrap_neighbour_(nullptr)
{
}

// Getter for x_
double Vector2D::GetX() const
{
    return x_;
}

// Getter for y_
double Vector2D::GetY() const
{
    return y_;
}

// Setter for Standard Neighbour
void Vector2D::SetStdNeighbour(Vector2D *buddy)
{
    std_neighbour_ = buddy;
}

// Setter for Wrapped Neighbour
void Vector2D::SetWrapNeighbour(Vector2D *buddy)
{
    wrap_neighbour_ = buddy;
}
