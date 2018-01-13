#ifndef __EXPORTER__
#define __EXPORTER__

#include <sstream>
#include <fstream>

#include "utils.h"

namespace exporter
{
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Export particles to Houdini Geo format, code from ///https://github.com/NCCA/SimulationExports/blob/master/HoudiniGeo/src/Emitter.cpp
    //----------------------------------------------------------------------------------------------------------------------
    void exportToHoudini(uint _frameNumber, std::vector<Particle> _particles);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Export particles to OBJ format
    //----------------------------------------------------------------------------------------------------------------------
    void exportToOBJ(uint _frameNumber, std::vector<Particle> _particles);

};
#endif
