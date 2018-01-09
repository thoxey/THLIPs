#include "MACGrid.h"

MACGrid::MACGrid(uvec3 _size)
{
    m_i_length = _size.x;
    m_j_length = _size.y;
    m_k_length = _size.y;
}

MACGrid::MACGrid(uint _i, uint _j, uint _k)
{
    m_i_length = _i;
    m_j_length = _j;
    m_k_length = _k;
}

vec3 MACGrid::getJitteredPos(MG_Cell _c)
{
    real x = randRange(h);
    real y = randRange(h);
    real z = randRange(h);

    return vec3((_c.gridPos.x*h) + x, (_c.gridPos.y*h) + y, (_c.gridPos.z*h) + z);
}

uint MACGrid::generateKey(uint _i, uint _j, uint _k)
{
    return 541*_i+79*_j+31*_k;
}

uint MACGrid::generateKey(uvec3 _pos)
{
    return 541*_pos.x+79*_pos.y+31*_pos.z;
}

void MACGrid::initialiseCell(MG_Cell _c, uvec3 _pos)
{
    _c.gridPos = _pos;
    _c.key = generateKey(_pos);
}

void MACGrid::initialiseCellWithFluid(MG_Cell _c, uvec3 _pos)
{
    _c.gridPos = _pos;
    _c.key = generateKey(_pos);
    for(uint i = 0; i < 8; i++)
    {
        MG_Particle p;
        p.cellidx = _c.key;
        p.pos = getJitteredPos(_c);
        m_particles.push_back(p);
    }
}

std::vector<MG_Cell> MACGrid::getNeighbors(MG_Cell _c)
{
    std::vector<MG_Cell> ret;
    ret.reserve(6);
    uvec3 centre = _c.gridPos;
    MG_Cell tmp;
    if(checkForCell(centre + uleftVec, tmp))
        ret[LEFT]=tmp;
    if(checkForCell(centre + urightVec, tmp))
        ret[RIGHT]=tmp;
    if(checkForCell(centre + uupVec, tmp))
        ret[UP]=tmp;
    if(checkForCell(centre + udownVec, tmp))
        ret[DOWN]=tmp;
    if(checkForCell(centre + uforwardVec, tmp))
        ret[FORWARD]=tmp;
    if(checkForCell(centre + ubackwardVec, tmp))
        ret[BACKWARD]=tmp;

    return ret;
}

bool MACGrid::checkInBounds(MG_Cell _c)
{
    return true;
}

real MACGrid::getMaxSpeed()
{
    real maxSpeed = 0.0;
    for(MG_Particle p : m_particles)
    {
        maxSpeed = std::max(maxSpeed, p.vel.x);
        maxSpeed = std::max(maxSpeed, p.vel.y);
        maxSpeed = std::max(maxSpeed, p.vel.z);
    }
}
/*
                // Trace a particle from point (x, y, z) for t time using RK2.
                Point traceParticle(real x, real y, real z, real t)
                    Vector V = getVelocity(x, y, z);
                    V = getVelocity(x+0.5*t*V.x, y+0.5*t*V.y, z+0.5*t*V.z);
                    return Point(x, y, z) + t*V;
                 */
vec3 MACGrid::tracePoint(vec3 _p, real _t)
{
    vec3 V = getVelocity(_p);

    V = getVelocity(vec3(_p.x+0.5f*_t*V.x, _p.y+0.5f*_t*V.y, _p.z+0.5f*_t*V.z));

    return _p + _t*V;
}
/*
                // Get the interpolated velocity at a point in space.
                Vector getVelocity(real x, real y, real z)
                    Vector V;
                    V.x = getInterpolatedValue(x/h, y/h-0.5, z/h-0.5, 0);
                    V.y = getInterpolatedValue(x/h-0.5, y/h, z/h-0.5, 1);
                    V.z = getInterpolatedValue(x/h-0.5, y/h-0.5, z/h, 2);
                    return V;
                 */
vec3 MACGrid::getVelocity(vec3 _v)
{
    _v.x = getInterpolatedValue(_v, 0);
    _v.y = getInterpolatedValue(_v, 1);
    _v.z = getInterpolatedValue(_v, 2);
    return _v;
}
/*
                // Get an interpolated data value from the grid.
                real getInterpolatedValue(real x, real y, real z, int index)
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
real MACGrid::getInterpolatedValue(vec3 _v, uint idx)
{
    int i = std::floor(_v.x);
    int j = std::floor(_v.y);
    int k = std::floor(_v.z);

    //I've tripled checked, but come back here if there are mistakes
    return (i+1-_v.x) * (j+1-_v.y) * (k+1-_v.z) * getCell(i, j, k).velField[idx] +
            (_v.x - i) * (j+1 - _v.y) * (k+1-_v.z) * getCell(i+1, j, k).velField[idx] +
            (i+1-_v.x) * (_v.y - j) * (k+1-_v.z) * getCell(i, j+1, k).velField[idx] +
            (_v.x - i) * (_v.y - j) * (k+1-_v.z) * getCell(i+1, j+1, k).velField[idx] +
            (i+1 - _v.x) * (j+1 - _v.y) * (_v.z - k) * getCell(i, j, k+1).velField[idx] +
            (_v.x - i) * (j+1 - _v.y) * (_v.z - k) * getCell(i+1, j, k+1).velField[idx] +
            (i+1 - _v.x) * (_v.y - j) * (_v.z - k) * getCell(i, j+1, k+1).velField[idx] +
            (_v.x - i) * (_v.y - j) * (_v.z - k) * getCell(i+1, j+1, k+1).velField[idx];
}

bool MACGrid::checkForCell(uint _i, uint _j, uint _k, MG_Cell& _c)
{
    std::unordered_map<int, MG_Cell>::const_iterator ret = m_hashTable.find(generateKey(_i,_j,_k));
    if(ret == m_hashTable.end())
        return false;
    else
    {
        _c = ret->second;
        return true;
    }

}

bool MACGrid::checkForCell(uvec3 _pos, MG_Cell& _c)
{
    std::unordered_map<int, MG_Cell>::const_iterator ret = m_hashTable.find(generateKey(_pos.x,_pos.y,_pos.z));
    if(ret == m_hashTable.end())
        return false;
    else
    {
        _c = ret->second;
        return true;
    }
}

bool MACGrid::checkForCell(MG_Particle _p)
{
    return true;
}

MG_Cell MACGrid::getCell(uint _i, uint _j, uint _k)
{
    std::unordered_map<int, MG_Cell>::const_iterator ret = m_hashTable.find(generateKey(_i,_j,_k));
    return ret->second;
}

MG_Cell MACGrid::getCell(uvec3 _pos)
{
    std::unordered_map<int, MG_Cell>::const_iterator ret = m_hashTable.find(generateKey(_pos.x,_pos.y,_pos.z));
    return ret->second;
}

void MACGrid::insertCellInHashTable(MG_Cell _c)
{
    m_hashTable.emplace(_c.key, _c);
}
