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
    void calculatePressure();
    void advectVelocityField();
    void addBodyForce();
    void project();



    float cfl();



    MACGrid m_MACGrid;

    //Grid Dimensions
    unsigned int m_iSize, m_jSize, m_kSize;

    float m_k_cfl = 1;

};

#endif
