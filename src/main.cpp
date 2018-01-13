#include <iostream>
#include"flipSim.h"
#include"exporter.h"
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

    FlipSim flipSim = FlipSim(10, 1.0, uvec3(0, 0, 0), uvec3(5,5,5));

//    std::vector<MG_Particle> test;
//    for(uint i = 0; i < 500; i++)
//    {
//        real x = utility::randRange(20.0);
//        real y = utility::randRange(20.0);
//        real z = utility::randRange(20.0);
//        MG_Particle p;
//        p.pos = vec3(x,y,z);
//        test.push_back(p);
//    }

    for(uint i = 0; i < 20; i++)
    {
        flipSim.step(0.0025);
        exporter::exportToHoudini(i, flipSim.getParticles());
//        for(uint i = 0; i < test.size(); i++)
//        {
//            real x = utility::randRange(20.0);
//            real y = utility::randRange(20.0);
//            real z = utility::randRange(20.0);
//            test[i].pos = vec3(x,y,z);
//        }
//        exporter::exportToHoudini(i, test);
        std::cout<<"Finished Exporting Frame: "<<i<<std::endl;
    }

    std::cout<<"THLIPs Finished! \n";
}
