#ifndef __EXPORTER__
#define __EXPORTER__

#include <sstream>
#include <fstream>
#include <iostream>

#include "utils.h"

class Exporter
{
public:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief ctor
    //----------------------------------------------------------------------------------------------------------------------
    Exporter(std::vector<MG_Particle> _particles);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Export particles to Houdini Geo format, code from ///https://github.com/NCCA/SimulationExports/blob/master/HoudiniGeo/src/Emitter.cpp
    //----------------------------------------------------------------------------------------------------------------------
    void exportToHoudini(uint _frameNumber);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Export particles to OBJ format
    //----------------------------------------------------------------------------------------------------------------------
    void exportToOBJ(uint _frameNumber);

private:
    std::vector<MG_Particle> m_particles;
};
#endif
