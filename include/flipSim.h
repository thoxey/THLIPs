#ifndef __FLIPSIM__
#define __FLIPSIM__

#include "MACGrid.h"



class FlipSim
{
public:
    void initialise();
    void step(real dt);
private:
    void initGrid();

    void updateGrid();
    void calculatePressure();
    void applyPressure();
    void advectVelocityField();
    void addBodyForce();
    void project();



    real cfl();



    MACGrid m_MACGrid;

    //Grid Dimensions
    uint m_iSize, m_jSize, m_kSize;

    real m_k_cfl = 1.0;

};

#endif
