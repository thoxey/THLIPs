#ifndef __FLIPSIM__
#define __FLIPSIM__

#include<glm/vec3.hpp>

struct Particle
{
    glm::vec3 position;
};
class FlipSim
{
public:
    void initialise();
    void step(float dt);
private:
    void initGrid();

    void advectVelocityField();
    void addBodyForce();
    void project();

    float cfl();


    //Grid Dimensions
    unsigned int m_iSize, m_jSize, m_kSize;

};

#endif
