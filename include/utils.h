#ifndef UTILS_H
#define UTILS_H

#include<algorithm>
#include<vector>
#include<unordered_map>
#include<cmath>
#include<random>

#include<glm/vec3.hpp>

#include<eigen3/Eigen/Sparse>
#include<eigen3/Eigen/SparseCholesky>

//----------------------------------------------------------------------------------------------------------------------
/// @file utils.h
/// @brief The utility header, contains definitions for grid components and universal functions
/// @author Tom Hoxey
/// @version 1.0
/// @date 19/01/17 Initial version
//----------------------------------------------------------------------------------------------------------------------

//#define USE_DOUBLE_PRECISION
//-fsingle-precision-constant
#ifdef USE_DOUBLE_PRECISION
typedef double real;
typedef Eigen::VectorXd VectorX;
typedef glm::dvec3 vec3;
#else
typedef float real;
typedef Eigen::VectorXf VectorX;
typedef glm::vec3 vec3;
#endif

typedef glm::uvec3 uvec3;
typedef unsigned int uint;

enum cellType {FLUID, AIR, SOLID};

enum neighbors {RIGHT, LEFT, UP, DOWN, FORWARD, BACKWARD};

struct MG_Cell
{
    uvec3 gridPos;

    //Vel field should be half a cell offset from the pressure
    vec3 velField, oldVelField;

    //Pressure
    real p;

    cellType type;

    uint key;
    //Fig 5.3 Bridsons book
    real rhs;

    MG_Cell(uvec3 _pos, vec3 _velField, real _p)
    {
        gridPos = _pos;
        key = 541*_pos.x+79*_pos.y+31*_pos.z;

        velField = _velField;

        _p = p;
    }
    MG_Cell(){}
    real u()
    {
        return velField.x;
    }
    real w()
    {
        return velField.y;
    }
    real v()
    {
        return velField.z;
    }

};

struct MG_Particle
{
    vec3 pos;

    vec3 vel;

    uint cellidx;

};

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
const uvec3 ubackwardVec = uvec3(0,0,1);

namespace utility
{
uint getIndex(uint _length, MG_Cell _c);

real randRange(real _max);

real randRange(real _min, real _max);

real lerp(real _a, real _b, real _x);

real invLerp(real _a, real _b, real _l);

real trilerp(std::vector<real> V, real x, real y, real z);

bool isInBounds(uvec3 _a, uvec3 _b, uvec3 _c);

}



//Empty doc tags for copy pasting!
//----------------------------------------------------------------------------------------------------------------------
/// @brief
//----------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------
/// @brief
/// @param
//----------------------------------------------------------------------------------------------------------------------

#endif // UTILS_H
