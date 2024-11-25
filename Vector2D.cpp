#include <iostream>
#include <cmath>

class Vector2D
{
public:
    Vector2D();                             // Default Constructor
    Vector2D(double a, double b);           // Parameterised Constructor
    double getX() const;                    // Getter for x_
    double getY() const;                    // Getter for y_
    void setStdNeighbour(Vector2D *buddy);  // Setter for nearest neighbour in standard geometry
    void setWrapNeighbour(Vector2D *buddy); // Setter for nearest neighbour in wraparound geometry
    void setStdFaraway(Vector2D *buddy);    // Setter for furthest faraway in standard geometry
    void setWrapFaraway(Vector2D *buddy);   // Setter for furthest faraway in wraparound geometry
    double getStdNeighbourDistance() const;  // Returns distance between self and nearest neighbour in standard geometry
    double getStdFarawayDistance() const;  // Returns distance between self and furthest faraway in standard geometry
    double getWrapNeighbourDistance() const; // Returns distance between self and nearest neighbour in wraparound geometry
    double getWrapFarawayDistance() const; // Returns distance between self and furthest faraway in wraparound geometry

private:
    double x_;                 // X Co Ordinate
    double y_;                 // Y Co Ordinate
    Vector2D *std_neighbour_;  // Pointer to the nearest standard neighbour - initially nothing hence nullptr in constructors
    Vector2D *wrap_neighbour_; // Pointer to the nearest wrapped neighbour - initially nothing hence nullptr in constructors
    Vector2D *std_faraway_;    // Pointer to the furthest standard faraway - initially nothing hence nullptr in constructors
    Vector2D *wrap_faraway_;   // Pointer to the furthest wrapped faraway - initially nothing hence nullptr in constructors
};

// Default Constructor
Vector2D::Vector2D() : x_(0.0),
                       y_(0.0),
                       std_neighbour_(nullptr),
                       wrap_neighbour_(nullptr),
                       std_faraway_(nullptr),
                       wrap_faraway_(nullptr)
{
}

// Parameterised Constructor
Vector2D::Vector2D(double a, double b) : x_(a),
                                         y_(b),
                                         std_neighbour_(nullptr),
                                         wrap_neighbour_(nullptr),
                                         std_faraway_(nullptr),
                                         wrap_faraway_(nullptr)
{
}

// Getter for x_
double Vector2D::getX() const
{
    return x_;
}

// Getter for y_
double Vector2D::getY() const
{
    return y_;
}

// Setter for Standard Neighbour
void Vector2D::setStdNeighbour(Vector2D *buddy)
{
    std_neighbour_ = buddy;
}

// Setter for Wrapped Neighbour
void Vector2D::setWrapNeighbour(Vector2D *buddy)
{
    wrap_neighbour_ = buddy;
}

// Setter for Standard Faraway
void Vector2D::setStdFaraway(Vector2D *buddy)
{
    std_faraway_ = buddy;
}

// Setter for Wrapped Faraway
void Vector2D::setWrapFaraway(Vector2D *buddy)
{
    wrap_faraway_ = buddy;
}

double Vector2D::getStdNeighbourDistance() const
{
    double dx = std::abs(wrap_neighbour_->getX() - x_);
    double dy = std::abs(wrap_neighbour_->getY() - y_);
    return sqrt(dx * dx + dy * dy);
}

double Vector2D::getStdFarawayDistance() const
{
    double dx = std::abs(std_faraway_->getX() - x_);
    double dy = std::abs(std_faraway_->getY() - y_);
    return sqrt(dx * dx + dy * dy);
}

double Vector2D::getWrapNeighbourDistance() const
{
    if (std_neighbour_ == wrap_neighbour_)
        return getStdNeighbourDistance();
    else
    {
        double dx = std::abs(wrap_neighbour_->getX() - x_);
        double dy = std::abs(wrap_neighbour_->getY() - y_);
        if (dx > 0.5)
        {
            dx = 1.0f - dx;
        }

        if (dy > 0.5)
        {
            dy = 1.0f - dy;
        }
        // Calculate the final distance with Pythagoras
        return sqrt(dx * dx + dy * dy);
    }
}

double Vector2D::getWrapFarawayDistance() const
{
    if (std_neighbour_ == wrap_neighbour_)
        return getStdFarawayDistance();
    else
    {
        double dx = std::abs(wrap_faraway_->getX() - x_);
        double dy = std::abs(wrap_faraway_->getY() - y_);
        if (dx > 0.5)
        {
            dx = 1.0f - dx;
        }

        if (dy > 0.5)
        {
            dy = 1.0f - dy;
        }
        // Calculate the final distance with Pythagoras
        return sqrt(dx * dx + dy * dy);
    }
}
