#include <iostream>
#include"flipSim.h"

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
    std::cout<<"Starting THLIPs... \n";

    FlipSim flipSim = FlipSim(uvec3(100, 100, 100), 0.1, uvec3(90, 90, 90), uvec3(100,100,100));

    for(uint i = 0; i < 100; i++)
    {
        flipSim.step(i);
    }

    std::cout<<"THLIPs Finished! \n";
}
