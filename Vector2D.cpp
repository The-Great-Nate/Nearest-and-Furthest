#include <iostream>
#include <cmath>

class Vector2D
{
public:
    Vector2D();                             // Default Constructor
    Vector2D(double a, double b);           // Parameterised Constructor
    double GetX() const;                    // Getter for x_
    double GetY() const;                    // Getter for y_
    void SetStdNeighbour(Vector2D *buddy);  // Setter for nearest neighbour in standard geometry
    void SetWrapNeighbour(Vector2D *buddy); // Setter for nearest neighbour in wraparound geometry
    void SetStdFaraway(Vector2D *buddy);    // Setter for furthest faraway in standard geometry
    void SetWrapFaraway(Vector2D *buddy);   // Setter for furthest faraway in wraparound geometry
    double GetStdNeighourDistance() const;  // Returns distance between self and nearest neighbour in standard geometry
    double GetStdFarawayDistance() const;  // Returns distance between self and furthest faraway in standard geometry
    double GetWrapNeighourDistance() const; // Returns distance between self and nearest neighbour in wraparound geometry
    double GetWrapFarawayDistance() const; // Returns distance between self and furthest faraway in wraparound geometry

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

// Setter for Standard Faraway
void Vector2D::SetStdFaraway(Vector2D *buddy)
{
    std_faraway_ = buddy;
}

// Setter for Wrapped Faraway
void Vector2D::SetWrapFaraway(Vector2D *buddy)
{
    wrap_faraway_ = buddy;
}

double Vector2D::GetStdNeighourDistance() const
{
    double dx = std::abs(wrap_neighbour_->GetX() - x_);
    double dy = std::abs(wrap_neighbour_->GetY() - y_);
    return sqrt(dx * dx + dy * dy);
}

double Vector2D::GetStdFarawayDistance() const
{
    double dx = std::abs(std_faraway_->GetX() - x_);
    double dy = std::abs(std_faraway_->GetY() - y_);
    return sqrt(dx * dx + dy * dy);
}

double Vector2D::GetWrapNeighourDistance() const
{
    if (std_neighbour_ == wrap_neighbour_)
        return GetStdNeighourDistance();
    else
    {
        double dx = std::abs(wrap_neighbour_->GetX() - x_);
        double dy = std::abs(wrap_neighbour_->GetY() - y_);
        if (dx > 0.5)
        {
            dx = 1.0f - dx;
        }

        if (dy > 0.5)
        {
            dy = 1.0f - dy;
        }
        // Calculate the final distance with pythagoras
        return sqrt(dx * dx + dy * dy);
    }
}

double Vector2D::GetWrapFarawayDistance() const
{
    if (std_neighbour_ == wrap_neighbour_)
        return GetStdFarawayDistance();
    else
    {
        double dx = std::abs(wrap_faraway_->GetX() - x_);
        double dy = std::abs(wrap_faraway_->GetY() - y_);
        if (dx > 0.5)
        {
            dx = 1.0f - dx;
        }

        if (dy > 0.5)
        {
            dy = 1.0f - dy;
        }
        // Calculate the final distance with pythagoras
        return sqrt(dx * dx + dy * dy);
    }
}
