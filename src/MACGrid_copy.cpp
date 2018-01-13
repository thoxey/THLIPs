/*
#include "MACGrid.h"

MACGrid::MACGrid(uint _size, real _cellWidth)
{
    m_i_length = _size;
    m_j_length = _size;
    m_k_length = _size;

    h = _cellWidth;
}

vec3 MACGrid::getJitteredPos(MG_Cell _c, uint _count)
{
    real ho2 = h/2;
    real x, y, z;

    enum corners {OOO, IOO, OIO, OOI, IOI, OII, IIO, III};
    switch (_count) {
    case OOO:
        x = utility::randRange(ho2);
        y = utility::randRange(ho2);
        z = utility::randRange(ho2);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case IOO:
        x = utility::randRange(ho2, h);
        y = utility::randRange(ho2);
        z = utility::randRange(ho2);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case OIO:
        x = utility::randRange(ho2);
        y = utility::randRange(ho2, h);
        z = utility::randRange(ho2);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case OOI:
        x = utility::randRange(ho2);
        y = utility::randRange(ho2);
        z = utility::randRange(ho2, h);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case IOI:
        x = utility::randRange(ho2, h);
        y = utility::randRange(ho2);
        z = utility::randRange(ho2, h);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case OII:
        x = utility::randRange(ho2);
        y = utility::randRange(ho2, h);
        z = utility::randRange(ho2, h);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case IIO:
        x = utility::randRange(ho2, h);
        y = utility::randRange(ho2, h);
        z = utility::randRange(ho2);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case III:
        x = utility::randRange(ho2, h);
        y = utility::randRange(ho2, h);
        z = utility::randRange(ho2, h);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    default:
        return getCellPos(_c);
        break;
    }
}



uint MACGrid::generateKey(uint _i, uint _j, uint _k)
{
    return 541*_i+79*_j+31*_k;
}

uint MACGrid::generateKey(uvec3 _pos)
{
    return 541*_pos.x+79*_pos.y+31*_pos.z;
}


void MACGrid::initialiseCells(uvec3 _b, uvec3 _c)
{

}

void MACGrid::initialiseCellWithFluid(MG_Cell& _c, uvec3 _pos)
{
    _c.gridPos = _pos;
    _c.key = generateKey(_pos);
    _c.type = FLUID;
    for(uint i = 0; i < 8; i++)
    {
        MG_Particle p;
        p.cellidx = _c.key;
        p.pos = getJitteredPos(_c, i);
        m_particles.push_back(p);
    }
}

std::vector<MG_Cell> MACGrid::getNeighbors(MG_Cell _c)
{
    std::vector<MG_Cell> ret;
    ret.reserve(6);
    uvec3 centre = _c.gridPos;
    MG_Cell tmp;
    MG_Cell empty;
    if(checkForCell(centre + uleftVec, tmp))
        ret.push_back(tmp);
    else
        ret.push_back(empty);
    if(checkForCell(centre + urightVec, tmp))
        ret.push_back(tmp);
    else
        ret.push_back(empty);
    if(checkForCell(centre + uupVec, tmp))
        ret.push_back(tmp);
    else
        ret.push_back(empty);
    if(checkForCell(centre + udownVec, tmp))
        ret.push_back(tmp);
    else
        ret.push_back(empty);
    if(checkForCell(centre + uforwardVec, tmp))
        ret.push_back(tmp);
    else
        ret.push_back(empty);
    if(checkForCell(centre + ubackwardVec, tmp))
        ret.push_back(tmp);
    else
        ret.push_back(empty);

    return ret;
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
    return maxSpeed;
}

vec3 MACGrid::tracePoint(vec3 _p, real _t)
{
    vec3 V = getVelocity(_p);

    //utility::printvec(V);
    V = getVelocity(vec3(_p.x+0.5f*_t*V.x, _p.y+0.5f*_t*V.y, _p.z+0.5f*_t*V.z));
    //utility::printvec(V);

    return _p + _t*V;
}

vec3 MACGrid::getVelocity(vec3 _v)
{
    _v.x = getInterpolatedValue(_v, 0);
    _v.y = getInterpolatedValue(_v, 1);
    _v.z = getInterpolatedValue(_v, 2);
    return _v;
}

real MACGrid::getInterpolatedValue(vec3 _v, uint idx)
{
    int i = std::floor(_v.x);
    int j = std::floor(_v.y);
    int k = std::floor(_v.z);

    //I've tripled checked, but come back here if there are mistakes
    real ret = (i+1-_v.x) * (j+1-_v.y) * (k+1-_v.z) * getCell(i, j, k).velField[idx] +
            (_v.x - i) * (j+1 - _v.y) * (k+1-_v.z) * getCell(i+1, j, k).velField[idx] +
            (i+1-_v.x) * (_v.y - j) * (k+1-_v.z) * getCell(i, j+1, k).velField[idx] +
            (_v.x - i) * (_v.y - j) * (k+1-_v.z) * getCell(i+1, j+1, k).velField[idx] +
            (i+1 - _v.x) * (j+1 - _v.y) * (_v.z - k) * getCell(i, j, k+1).velField[idx] +
            (_v.x - i) * (j+1 - _v.y) * (_v.z - k) * getCell(i+1, j, k+1).velField[idx] +
            (i+1 - _v.x) * (_v.y - j) * (_v.z - k) * getCell(i, j+1, k+1).velField[idx] +
            (_v.x - i) * (_v.y - j) * (_v.z - k) * getCell(i+1, j+1, k+1).velField[idx];
    return ret;
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


MG_Cell MACGrid::getCell(uint _key)
{
    auto ret = m_hashTable.find(_key);
    if(ret != m_hashTable.end())
        return ret->second;
    else
    {
        MG_Cell c;
        return c;
    }
}

MG_Cell MACGrid::getCell(uint _i, uint _j, uint _k)
{
    auto ret = m_hashTable.find(generateKey(_i,_j,_k));
    if(ret != m_hashTable.end())
        return ret->second;
    else
    {
        MG_Cell c;
        return c;
    }
}

MG_Cell MACGrid::getCell(uvec3 _pos)
{
    auto ret = m_hashTable.find(generateKey(_pos.x,_pos.y,_pos.z));
    if(ret != m_hashTable.end())
        return ret->second;
    else
    {
        MG_Cell c;
        return c;
    }
}

void MACGrid::insertCellInHashTable(MG_Cell _c)
{
    auto test = m_hashTable.find(_c.key);
    if(test != m_hashTable.end())
        test->second = _c;
    else
        m_hashTable.emplace(_c.key, _c);
}

vec3 MACGrid::getCellPos(Cell _c)
{
    return vec3(_c.gridPos.x * h, _c.gridPos.y * h, _c.gridPos.z * h);
}
*/
