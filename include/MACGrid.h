#ifndef __MACGRID__
#define __MACGRID__

#include<vector>
#include<unordered_map>

#include<glm/vec3.hpp>

enum cellType {FLUID, AIR, SOLID};

struct MG_Cell
{
    glm::uvec3 gridPos;

    glm::vec3 velField;

    //Pressure
    float p;

    int layer;

    cellType type;

    unsigned int key;

    MG_Cell(glm::uvec3 _pos, glm::vec3 _velField, float _p)
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
    glm::vec3 pos;

    glm::vec3 vel;

    unsigned int cellidx;
};


class MACGrid
{
public:
    std::vector<MG_Cell> m_cells;
    std::vector<MG_Particle> m_particles;

    std::unordered_map<int, MG_Cell> m_hashTable;

    std::vector<MG_Cell> getNeighbors();

    float getMaxSpeed();

    bool checkInBounds(MG_Cell _c);

    glm::vec3 tracePoint(glm::vec3 _p, float _t);

    glm::vec3 getVelocity(glm::vec3 _v);

    float getInterpolatedValue(glm::vec3 _v, unsigned int idx);

    MG_Cell getCell(unsigned int _i, unsigned int _j, unsigned int _k);

    unsigned int generateKey(unsigned int _i, unsigned int _j, unsigned int _k);

    void insertCellInHashTable(MG_Cell _c);

    //Width, using notaion from H&W paper
    float h;

private:
    int m_i_length, m_j_length, m_k_length;


};

#endif
