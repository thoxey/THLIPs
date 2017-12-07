#ifndef __FLIPSIM__
#define __FLIPSIM__

#include<algorithm>
#include<cmath>

#include<glm/vec3.hpp>

#include "MACGrid.h"

class FlipSim
{
public:
    void initialise();
    void step(float dt);
private:
    void initGrid();

    void updateGrid();
    void advectVelocityField();
    void addBodyForce();
    void project();

    bool checkForCell(MG_Particle _p);

    float cfl();



    MACGrid m_MACGrid;

    //Grid Dimensions
    unsigned int m_iSize, m_jSize, m_kSize;

    float m_k_cfl = 1;

};

#endif
