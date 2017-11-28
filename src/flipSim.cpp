#include "flipSim.h"

void FlipSim::initialise()
{
    //Place 8 jittered particles per voxel into a voxel grid
    initGrid();

}

void FlipSim::initGrid()
{

}


void FlipSim::step(float dt)
{
    /*
    Basic algorithm from Bridson's Course Notes:
    (U is a velocity field)
    >for timeStep n = 0,1,2...
        >Determine Time Step dt to go from tn to tn+1

        Advect quantity q, through vector field Un for time interval dt
        >set Ua = advect(Un, dt, q)

        Add body forces
        >Ub = Ua +dtg

        To handle pressure/incompressibility we make a function called project
        >Un+1 = project(dt, Ub)
    */

    float t = 0;

    while(t < dt)
    {
        //Calculate our substep
        float subStep = cfl();

        //Per cell
//        sample();
//        {
//            updateGrid();
//        }


        //Per cell
        project();

        //Per particle
        //advect the velocity field
        advectVelocityField();

        //Add gravity and stuff
        //Per particle calculation
        addBodyForce();

        t += subStep;
    }
}

void FlipSim::advectVelocityField()
{
    for(unsigned int k = 0; k < m_kSize; k++)
        for(unsigned int j = 0; j < m_jSize; j++)
            for(unsigned int i = 0; i < m_iSize; i++)
            {
                //do advection
            }
}

void FlipSim::addBodyForce()
{

}

void FlipSim::project()
{

}

float FlipSim::cfl()
{
    return 0.0f;
}
