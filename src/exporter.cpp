#include "exporter.h"
namespace exporter
{
/// This function was originally written by Jon Macey
void exportToHoudini(uint _frameNumber, std::vector<Particle> _particles)
{
        char fname[50];
        std::sprintf(fname,"geo/THLIPSparticle.%03d.geo",++_frameNumber);
        // we will use a stringstream as it may be more efficient
        std::stringstream ss;
        std::ofstream file;
        file.open(fname);
        if (!file.is_open())
        {
            std::cerr << "failed to Open file "<<fname<<'\n';
            exit(EXIT_FAILURE);
        }
        // write header see here http://www.sidefx.com/docs/houdini15.0/io/formats/geo
        ss << "PGEOMETRY V5\n";
        ss << "NPoints " << _particles.size() << " NPrims 1\n";
        ss << "NPointGroups 0 NPrimGroups 1\n";
        // this is hard coded but could be flexible we have 1 attrib which is Colour
        ss << "NPointAttrib 1  NVertexAttrib 0 NPrimAttrib 2 NAttrib 0\n";
        // now write out our point attrib this case Cd for diffuse colour
        ss <<"PointAttrib \n";
        // default the colour to white
        ss <<"Cd 3 float 1 1 1\n";
        // now we write out the particle data in the format
        // x y z 1 (attrib so in this case colour)
        for(unsigned int i=0; i<_particles.size(); ++i)
        {
            ss<<_particles[i].pos.x<<" "<<_particles[i].pos.y<<" "<<_particles[i].pos.z << " 1 ";
            //ss<<"("<<_particles[i].cellCol.x<<" "<<_particles[i].cellCol.y<<" "<< _particles[i].cellCol.z<<")\n";
            ss<<"("<<std::abs(_particles[i].vel.x)<<" "<<std::abs(_particles[i].vel.y)<<" "<<std::abs(_particles[i].vel.z)<<")\n";
        }

        // now write out the index values
        ss<<"PrimitiveAttrib\n";
        ss<<"generator 1 index 1 location1\n";
        ss<<"dopobject 1 index 1 /obj/AutoDopNetwork:1\n";
        ss<<"Part "<<_particles.size()<<" ";
        for(size_t i=0; i<_particles.size(); ++i)
        {
            ss<<i<<" ";
        }
        ss<<" [0	0]\n";
        ss<<"box_object1 unordered\n";
        ss<<"1 1\n";
        ss<<"beginExtra\n";
        ss<<"endExtra\n";
        // dump string stream to disk;
        file<<ss.rdbuf();
        file.close();
}
/// end of function

void exportToOBJ(uint _frameNumber, std::vector<Particle> _particles)
{
        char fname[50];
        std::sprintf(fname,"geo/THLIPSparticle.%03d.obj",_frameNumber++);
        // we will use a stringstream as it may be more efficient
        std::stringstream ss;
        std::ofstream file;
        file.open(fname);
        if (!file.is_open())
        {
            std::cerr << "failed to Open file "<<fname<<'\n';
            exit(EXIT_FAILURE);
        }
        for(Particle p : _particles)
        ss<<"v" <<" "<< p.pos.x <<" "<< p.pos.y <<" "<< p.pos.z<<"\n";
        // dump string stream to disk;
        file<<ss.rdbuf();
        file.close();
}
}
