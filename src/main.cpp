#include <iostream>
#include"flipSim.h"
#include"exporter.h"

#include <fenv.h>

/*
  TO DO:
  -FINISH DOC TAGS (always)
  -Check the density in pressure solve
  -Check Advection Algorithm
  -Look into level sets
  -Add body forces
  -Transfer velocities from grid to particles and vice versa
  -Actually make a chain of functions and export some fngo
 */

int main()
{
    feenableexcept(FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW);
    std::cout<<"Starting THLIPs... \n";

    FlipSim flipSim = FlipSim(10, 0.1, uvec3(0, 5, 0), uvec3(10,10,10));

    for(uint i = 0; i < 50; i++)
    {
        exporter::exportToHoudini(i, flipSim.getParticles());
        flipSim.step(0.1);
        std::cout<<"--------------------------\n"<<"Finished Exporting Frame: "<<i+1<<"\n"<<"--------------------------"<<std::endl;
    }

    std::cout<<"THLIPs Finished! \n";
}
