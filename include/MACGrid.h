#ifndef __MACGRID__
#define __MACGRID__

#include "utils.h"

class MACGrid
{
public:

    MACGrid(uvec3 _size);
    MACGrid(uint _i, uint _j, uint _k);

    std::vector<MG_Cell> m_cells;
    std::vector<MG_Particle> m_particles;

    std::unordered_map<int, MG_Cell> m_hashTable;

    std::vector<MG_Cell> getNeighbors(MG_Cell _c);

    real getMaxSpeed();

    bool checkInBounds(MG_Cell _c);

    vec3 tracePoint(vec3 _p, real _t);

    vec3 getVelocity(vec3 _v);

    real getInterpolatedValue(vec3 _v, uint idx);

    MG_Cell getCell(uint _i, uint _j, uint _k);
    MG_Cell getCell(uvec3 _pos);

    uint generateKey(uint _i, uint _j, uint _k);
    uint generateKey(uvec3 _pos);

    bool checkForCell(uint _i, uint _j, uint _k);
    bool checkForCell(uvec3 _pos);
    bool checkForCell(MG_Particle _p);

    void insertCellInHashTable(MG_Cell _c);



    //Width, using notaion from H&W paper
    real h;

private:
    uint m_i_length, m_j_length, m_k_length;


};

#endif
