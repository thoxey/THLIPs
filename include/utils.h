#ifndef UTILS_H
#define UTILS_H

#include<algorithm>
#include<vector>
#include<unordered_map>
#include<cmath>

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

struct MG_Cell
{
    uvec3 gridPos;

    vec3 velField;

    //Pressure
    real p;

    int layer;

    cellType type;

    uint key;

    MG_Cell(uvec3 _pos, vec3 _velField, real _p)
    {
        gridPos = _pos;
        key = 541*_pos.x+79*_pos.y+31*_pos.z;

        velField = _velField;

        _p = p;
    }
    MG_Cell()
    {
        layer = -1;
    }

};

struct MG_Particle
{
    vec3 pos;

    vec3 vel;

    uint cellidx;
};

uint getIndex(uint _length, MG_Cell _c);

#endif // UTILS_H
