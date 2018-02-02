#ifndef UTILS_H
#define UTILS_H

#include<algorithm>
#include<vector>
#include<unordered_map>
#include<cmath>
#include<random>
#include<iostream>
#include <algorithm>

#include<glm/vec3.hpp>

#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/Sparse>
#include<eigen3/Eigen/SparseCholesky>

//----------------------------------------------------------------------------------------------------------------------
/// @file utils.h
/// @brief The utility header, contains definitions for grid components and universal functions
/// @author Tom Hoxey
/// @version 1.0
/// @date 19/01/17 Initial version
//----------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------
/// TYPEDEFS TO ALLOW FOR EASILY SWAPPING BETWEEN VECTOR IMPLEMENTATIONS
#define USE_DOUBLE_PRECISION
//-fsingle-precision-constant
#ifdef USE_DOUBLE_PRECISION
typedef double real;
typedef Eigen::VectorXd VectorX;
typedef glm::dvec3 vec3;
typedef Eigen::MatrixXd mat;
#else
typedef float real;
typedef Eigen::VectorXf VectorX;
typedef glm::vec3 vec3;
typedef Eigen::MatrixXf mat;
#endif
typedef glm::uvec3 uvec3;
typedef glm::ivec3 ivec3;
typedef unsigned int uint;
//----------------------------------------------------------------------------------------------------------------------
enum cellType {FLUID, AIR, SOLID};

enum neighbors {RIGHT, LEFT, UP, DOWN, FORWARD, BACKWARD};
//----------------------------------------------------------------------------------------------------------------------
const vec3 rightVec = vec3(1.0,0.0,0.0);
const vec3 leftVec = vec3(-1.0,0.0,0.0);
const vec3 upVec = vec3(0.0,1.0,0.0);
const vec3 downVec = vec3(0.0,-1.0,0.0);
const vec3 forwardVec = vec3(0.0,0.0,1.0);
const vec3 backwardVec = vec3(0.0,0.0,1.0);

const uvec3 urightVec = uvec3(1,0,0);
const uvec3 uleftVec = uvec3(-1,0,0);
const uvec3 uupVec = uvec3(0,1,0);
const uvec3 udownVec = uvec3(0,-1,0);
const uvec3 uforwardVec = uvec3(0,0,1);
const uvec3 ubackwardVec = uvec3(0,0,-1);
//----------------------------------------------------------------------------------------------------------------------
namespace utils
{
//----------------------------------------------------------------------------------------------------------------------
/// @brief Returns the index of a cell depending on the size of the grid
/// @param uint _length : The amount of cells in a given axis of the grid
/// @param uvec3 _pos : The position of he cell in the grid
//----------------------------------------------------------------------------------------------------------------------
uint getIndex(uint _length, uvec3 _pos);
//----------------------------------------------------------------------------------------------------------------------
/// @brief Returns a random number between 0 and _max
/// @param real _max : The maximum size of the returned value
//----------------------------------------------------------------------------------------------------------------------
real randRange(real _max);
//----------------------------------------------------------------------------------------------------------------------
/// @brief Returns a random number between _min and _max
/// @param real _min : The minimum size of the returned value
/// @param real _max : The maximum size of the returned value
//----------------------------------------------------------------------------------------------------------------------
real randRange(real _min, real _max);
//----------------------------------------------------------------------------------------------------------------------
/// @brief Returns a value at _x between _ and _b
/// @param real _a : The point to interpolate from
/// @param real _b : The point to interpolate to
/// @param real _x : The amount to interpolate
//----------------------------------------------------------------------------------------------------------------------
real lerp(real _a, real _b, real _x);
//----------------------------------------------------------------------------------------------------------------------
/// @brief Returns the amount interpolated between a and b if x is l away
/// @param real _a : The point to interpolate from
/// @param real _b : The point to interpolate to
/// @param real _l : The length interpolated between a and b
//----------------------------------------------------------------------------------------------------------------------
real invLerp(real _a, real _b, real _l);
//----------------------------------------------------------------------------------------------------------------------
/// @brief Returns a value at position xyz based on the values at the verticies
/// @param real * V : The values at each the verticies
/// @param real x : The x position of the point we are querying
/// @param real y : The y position of the point we are querying
/// @param real z : The z position of the point we are querying
//----------------------------------------------------------------------------------------------------------------------
real trilerp(real * V, real x, real y, real z);
//----------------------------------------------------------------------------------------------------------------------
/// @brief Checks if a point a is in the bounds of b and c
/// @param vec3 _a : The point we are querying
/// @param uvec3 _b : The bottom corner of the fluid cells (0,0,0)
/// @param uvec3 _c : The top corner of the fluid cells (1,1,1)
//----------------------------------------------------------------------------------------------------------------------
bool isInBounds(vec3 _a, vec3 _b, vec3 _c);
//----------------------------------------------------------------------------------------------------------------------
/// @brief Prints the vector inputted
/// @param uvec3 _x : The vector to print
//----------------------------------------------------------------------------------------------------------------------
void printvec(uvec3 _x);
//----------------------------------------------------------------------------------------------------------------------
/// @brief Prints the vector inputted
/// @param vec3 _x : The vector to print
//----------------------------------------------------------------------------------------------------------------------
void printvec(vec3 _x);
//----------------------------------------------------------------------------------------------------------------------
/// @brief Calculates the trilinear hat kernel
/// @param vec3 _dist : The distance used to calculate the hat funtion
/// @param real _dx : The width of a cell
//----------------------------------------------------------------------------------------------------------------------
real trilinearHatKernel(vec3 _dist, real _dx);
//----------------------------------------------------------------------------------------------------------------------
/// @brief Calculates the hat function of r
/// @param real _r : The input to the function
//----------------------------------------------------------------------------------------------------------------------
real hatFunction(real _r);
//----------------------------------------------------------------------------------------------------------------------
/// @brief Calculates the divergent velocity of v1 and v2
/// @param real _v1 : The first input velocity
/// @param real _v2 : The second input velocity
/// @param real _dx : The distance we are diverging across
//----------------------------------------------------------------------------------------------------------------------
real divergentVelocity(real _v1, real _v2, real _dx);
//----------------------------------------------------------------------------------------------------------------------
//Reference https://www.physicsforums.com/threads/runge-kutta-in-c.1448/
real runge(real x, real y, real _dt);
real f(real x, real y);
//End Reference
} //End namespace
//----------------------------------------------------------------------------------------------------------------------




//----------------------------------------------------------------------------------------------------------------------
struct Particle
{
    vec3 pos;

    vec3 vel;

    vec3 cellCol;

    uint idx;

    Particle()
    {
        pos = vec3(0,0,0);
        vel = vec3(0,0,0);
    }

    void updatePos(vec3 _newPos)
    {
        pos = _newPos;
    }
    void updateVel(vec3 _newVel, vec3 _cellCol)
    {
        vel = _newVel;
        cellCol = _cellCol;
    }
    void advect(real _dt)
    {
//        pos.x += utils::runge(vel.x, pos.x, _dt);
//        pos.y += utils::runge(vel.y, pos.y, _dt);
//        pos.z += utils::runge(vel.z, pos.z, _dt);
        pos.x += vel.x * _dt;
        pos.y += vel.y * _dt;
        pos.z += vel.z * _dt;
    }
};
//----------------------------------------------------------------------------------------------------------------------




//Empty doc tags for copy pasting!
//----------------------------------------------------------------------------------------------------------------------
/// @brief
//----------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------
/// @brief
/// @param
//----------------------------------------------------------------------------------------------------------------------

#endif // UTILS_H
