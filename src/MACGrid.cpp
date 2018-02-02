#include "MACGrid.h"

MACGrid::MACGrid(uint _size, real _cellWidth):
    m_gridLength(_size),
    m_h(_cellWidth)
{}
//----------------------------------------------------------------------------------------------------------------------
vec3 MACGrid::getJitteredPos(Cell _c, uint _count)
{
    real ho2 = m_h/2;
    real x, y, z;

    enum corners {OOO, IOO, OIO, OOI, IOI, OII, IIO, III};
    switch (_count) {
    case OOO:
        x = utils::randRange(ho2);
        y = utils::randRange(ho2);
        z = utils::randRange(ho2);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case IOO:
        x = utils::randRange(ho2, m_h);
        y = utils::randRange(ho2);
        z = utils::randRange(ho2);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case OIO:
        x = utils::randRange(ho2);
        y = utils::randRange(ho2, m_h);
        z = utils::randRange(ho2);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case OOI:
        x = utils::randRange(ho2);
        y = utils::randRange(ho2);
        z = utils::randRange(ho2, m_h);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case IOI:
        x = utils::randRange(ho2, m_h);
        y = utils::randRange(ho2);
        z = utils::randRange(ho2, m_h);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case OII:
        x = utils::randRange(ho2);
        y = utils::randRange(ho2, m_h);
        z = utils::randRange(ho2, m_h);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case IIO:
        x = utils::randRange(ho2, m_h);
        y = utils::randRange(ho2, m_h);
        z = utils::randRange(ho2);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    case III:
        x = utils::randRange(ho2, m_h);
        y = utils::randRange(ho2, m_h);
        z = utils::randRange(ho2, m_h);
        return getCellPos(_c) + vec3(x, y, z);
        break;
    default:
        return getCellPos(_c);
        break;
    }
}
//----------------------------------------------------------------------------------------------------------------------
vec3 MACGrid::getCellPos(Cell _c)
{
    return vec3(_c.gridPos.x * m_h, _c.gridPos.y * m_h, _c.gridPos.z * m_h);
}
//----------------------------------------------------------------------------------------------------------------------
void MACGrid::reclassifyCells()
{
    for(auto&& c : m_cells)
    {
        //Solid cells wont change so skip
        if(c.type == SOLID)
            continue;

        //Set to AIR and clear particle indicies
        c.m_paticleIDXs.clear();
        c.type = AIR;

        for(auto&& p : m_particles)
        {
            vec3 cell000 = getCellPos(c);
            vec3 cell111 = getCellPos(c)+vec3(m_h);

            if(utils::isInBounds(p.pos, cell000, cell111))
            {
                c.type = FLUID;
                c.m_paticleIDXs.push_back(p.idx);
            }
        }
    }
}
//----------------------------------------------------------------------------------------------------------------------
void MACGrid::initialiseCells(uvec3 _b, uvec3 _c)
{
    for(uint k = 0; k < m_gridLength; k++)
        {
            for(uint j = 0; j < m_gridLength; j++)
            {
                for(uint i = 0; i < m_gridLength; i++)
                {
                    Cell c = Cell(uvec3(i,j,k));
                    if(i == 0 || i == m_gridLength-1 || j == 0 || j == m_gridLength-1 || k == 0 || k == m_gridLength-1)
                        c.type = SOLID;
                    else if(utils::isInBounds(uvec3(i,j,k), _b, _c))
                    {
                        c.type = FLUID;
                        for(uint i = 0; i < 8; i++)
                        {
                            Particle p;
                            p.idx = m_particleCount++;
                            p.pos = getJitteredPos(c, i);
                            m_particles.push_back(p);
                            c.m_paticleIDXs.push_back(p.idx);
                        }
                    }
                    else
                        c.type = AIR;

                    m_cells.push_back(c);
                }
            }
        }

}
//----------------------------------------------------------------------------------------------------------------------
int * MACGrid::getNeighbors(Cell _c)
{
    static int ret[6];
    uvec3 refPos  = _c.gridPos;
    ret[LEFT]     = utils::getIndex(m_gridLength, getCell(refPos - uleftVec).gridPos);
    ret[RIGHT]    = utils::getIndex(m_gridLength, getCell(refPos - urightVec).gridPos);
    ret[UP]       = utils::getIndex(m_gridLength, getCell(refPos - uupVec).gridPos);
    ret[DOWN]     = utils::getIndex(m_gridLength, getCell(refPos - udownVec).gridPos);
    ret[FORWARD]  = utils::getIndex(m_gridLength, getCell(refPos - uforwardVec).gridPos);
    ret[BACKWARD] = utils::getIndex(m_gridLength, getCell(refPos - ubackwardVec).gridPos);
    return ret;
}
//----------------------------------------------------------------------------------------------------------------------
real MACGrid::getMaxSpeed()
{
    real maxSpeed = 1.0;//Stops divide by zero error later
    for(auto&& p : m_particles)
    {
        maxSpeed = std::max(maxSpeed, p.vel.x);
        maxSpeed = std::max(maxSpeed, p.vel.y);
        maxSpeed = std::max(maxSpeed, p.vel.z);
    }
    return maxSpeed;
}
//----------------------------------------------------------------------------------------------------------------------
Cell & MACGrid::getCell(uint _i, uint _j, uint _k)
{
    uint i = utils::getIndex(m_gridLength, uvec3(_i,_j,_k));
    if(i < m_cells.size())
        return m_cells[i];
    else
    {
        Cell c = Cell(uvec3(_i,_j,_k));
        return c;
    }
}
//----------------------------------------------------------------------------------------------------------------------
Cell MACGrid::getCell(uvec3 _pos)
{
    int i = utils::getIndex(m_h, _pos);
    if(m_cells.size() < i && i > 0)
        return m_cells[i];
    else
    {
        Cell c = Cell(_pos);
        return c;
    }
}

