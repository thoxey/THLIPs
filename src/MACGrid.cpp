#include "MACGrid.h"

std::vector<MG_Cell> MACGrid::getNeighbors()
{
    return std::vector<MG_Cell>();
}

bool MACGrid::checkInBounds(MG_Cell _c)
{
    return true;
}

float MACGrid::getMaxSpeed()
{
    float maxSpeed = 0.0f;
    for(MG_Particle p : m_particles)
    {
        maxSpeed = std::max(maxSpeed, p.vel.x);
        maxSpeed = std::max(maxSpeed, p.vel.y);
        maxSpeed = std::max(maxSpeed, p.vel.z);
    }
}
/*
// Trace a particle from point (x, y, z) for t time using RK2.
Point traceParticle(float x, float y, float z, float t)
    Vector V = getVelocity(x, y, z);
    V = getVelocity(x+0.5*t*V.x, y+0.5*t*V.y, z+0.5*t*V.z);
    return Point(x, y, z) + t*V;
 */
glm::vec3 MACGrid::tracePoint(glm::vec3 _p, float _t)
{
    glm::vec3 V = getVelocity(_p);

    V = getVelocity(glm::vec3(_p.x+0.5f*_t*V.x, _p.y+0.5f*_t*V.y, _p.z+0.5f*_t*V.z));

    return _p + _t*V;
}
/*
// Get the interpolated velocity at a point in space.
Vector getVelocity(float x, float y, float z)
    Vector V;
    V.x = getInterpolatedValue(x/h, y/h-0.5, z/h-0.5, 0);
    V.y = getInterpolatedValue(x/h-0.5, y/h, z/h-0.5, 1);
    V.z = getInterpolatedValue(x/h-0.5, y/h-0.5, z/h, 2);
    return V;
 */
glm::vec3 MACGrid::getVelocity(glm::vec3 _v)
{
    _v.x = getInterpolatedValue(_v, 0);
    _v.y = getInterpolatedValue(_v, 1);
    _v.z = getInterpolatedValue(_v, 2);
    return _v;
}
/*
// Get an interpolated data value from the grid.
float getInterpolatedValue(float x, float y, float z, int index)
    int i = floor(x);
    int j = floor(y);
    int k = floor(z);
    return
    (i+1-x) * (j+1-y) * (k+1-z) * cell(i, j, k).u[index] +
    (x-i) * (j+1-y) * (k+1-z) * cell(i+1, j, k).u[index] +
    (i+1-x) * (y-j) * (k+1-z) * cell(i, j+1, k).u[index] +
    (x-i) * (y-j) * (k+1-z) * cell(i+1, j+1, k).u[index] +
    (i+1-x) * (j+1-y) * (z-k) * cell(i, j, k+1).u[index] +
    (x-i) * (j+1-y) * (z-k) * cell(i+1, j, k+1).u[index] +

(i+1-x) * (y-j) * (z-k) * cell(i, j+1, k+1).u[index] +

(x-i) * (y-j) * (z-k) * cell(i+1, j+1, k+1).u[index];
 */
float MACGrid::getInterpolatedValue(glm::vec3 _v, unsigned int idx)
{
    int i = std::floor(_v.x);
    int j = std::floor(_v.y);
    int k = std::floor(_v.z);

    return (i+1-_v.x) * (j+1-_v.y) * (k+1-_v.z) * getCell(i, j, k).velField[idx] +
           (_v.x - i) * (j+1 - _v.y) * (k+1-_v.z) * getCell(i+1, j, k).velField[idx] +
           (i+1-_v.x) * (_v.y - j) * (k+1-_v.z) * getCell(i, j, k).velField[idx] +
           (_v.x - i) * (_v.y - j) * (k+1-_v.z) * getCell(i+1, j+1, k).velField[idx] +
           (i+1 - _v.x) * (j+1 - _v.y) * (_v.z - k) * getCell(i, j, k+1).velField[idx] +
           (_v.x - i) * (j+1 - _v.y) * (_v.z - k) * getCell(i+1, j, k+1).velField[idx] +
           (i+1 - _v.x) * (_v.y - j) * (_v.z - k) * getCell(i, j+1, k+1).velField[idx] +
           (_v.x - i) * (_v.y - j) * (_v.z - k) * getCell(i+1, j+1, k+1).velField[idx];
}



MG_Cell MACGrid::getCell(unsigned int _i, unsigned int _j, unsigned int _k)
{
    std::unordered_map<int, MG_Cell>::const_iterator ret = m_hashTable.find(generateKey(_i,_j,_k));
    return ret->second;
}

unsigned int MACGrid::generateKey(unsigned int _i, unsigned int _j, unsigned int _k)
{
    return 541*_i+79*_j+31*_k;
}

void MACGrid::insertCellInHashTable(MG_Cell _c)
{
    m_hashTable.emplace(_c.key, _c);
}
