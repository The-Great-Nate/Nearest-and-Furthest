#ifndef VECTOR2D_H
#define VECTOR2D_H

#include <cmath>

class Vector2D
{
public:
    // Constructors
    Vector2D();                             // Default Constructor
    Vector2D(double a, double b);           // Parameterized Constructor

    // Getters
    double getX() const;                    // Getter for x_
    double getY() const;                    // Getter for y_

    // Setters for neighbours and faraways
    void setStdNeighbour(Vector2D *buddy);  // Setter for nearest neighbour in standard geometry
    void setWrapNeighbour(Vector2D *buddy); // Setter for nearest neighbour in wraparound geometry
    void setStdFaraway(Vector2D *buddy);    // Setter for furthest faraway in standard geometry
    void setWrapFaraway(Vector2D *buddy);   // Setter for furthest faraway in wraparound geometry

    // Distance getters
    double getStdNeighbourDistance() const;  // Returns distance between self and nearest neighbour in standard geometry
    double getStdFarawayDistance() const;    // Returns distance between self and furthest faraway in standard geometry
    double getWrapNeighbourDistance() const; // Returns distance between self and nearest neighbour in wraparound geometry
    double getWrapFarawayDistance() const;   // Returns distance between self and furthest faraway in wraparound geometry

private:
    double x_;                 // X Coordinate
    double y_;                 // Y Coordinate
    Vector2D *std_neighbour_;  // Pointer to the nearest standard neighbour - initially nullptr
    Vector2D *wrap_neighbour_; // Pointer to the nearest wrapped neighbour - initially nullptr
    Vector2D *std_faraway_;    // Pointer to the furthest standard faraway - initially nullptr
    Vector2D *wrap_faraway_;   // Pointer to the furthest wrapped faraway - initially nullptr
};

#endif // VECTOR2D_H
