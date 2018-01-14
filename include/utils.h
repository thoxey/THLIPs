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
//----------------------------------------------------------------------------------------------------------------------

typedef glm::uvec3 uvec3;
typedef glm::ivec3 ivec3;
typedef unsigned int uint;

enum cellType {FLUID, AIR, SOLID};

enum neighbors {RIGHT, LEFT, UP, DOWN, FORWARD, BACKWARD};

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

namespace utils
{
uint getIndex(uint _length, uvec3 _pos);

real randRange(real _max);

real randRange(real _min, real _max);

real lerp(real _a, real _b, real _x);

real invLerp(real _a, real _b, real _l);

real trilerp(real * V, real x, real y, real z);

bool isInBounds(vec3 _a, vec3 _b, vec3 _c);

void printvec(uvec3 _x);

void printvec(vec3 _x);

real trilinearHatKernel(vec3 _dist, real _dx);

real hatFunction(real _r);

real divergentVelocity(real _v1, real _v2, real _dx);

//Reference https://www.physicsforums.com/threads/runge-kutta-in-c.1448/
real runge(real x, real y, real _dt);

real f(real x, real y);
//End Reference
}

struct Particle
{
    vec3 pos;

    vec3 vel;

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
    void updateVel(vec3 _newVel)
    {
        vel = _newVel;
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

//Empty doc tags for copy pasting!
//----------------------------------------------------------------------------------------------------------------------
/// @brief
//----------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------
/// @brief
/// @param
//----------------------------------------------------------------------------------------------------------------------

#endif // UTILS_H
