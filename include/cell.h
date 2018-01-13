#ifndef __CELL__
#define __CELL__

#include "utils.h"

//----------------------------------------------------------------------------------------------------------------------
/// @file gridCell.h
/// @brief The cells in te MAC Grid
/// @author Tom Hoxey
/// @version 1.0
/// @date 19/01/17 Initial version
//----------------------------------------------------------------------------------------------------------------------

class Cell
{
public:
    //Constructors
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief ctor
    //----------------------------------------------------------------------------------------------------------------------
    Cell(ivec3 _pos);

    void updateVel(vec3 _newVel);

    void increaseVel(vec3 _addition);

    real U();

    void setU(real _newU);

    real V();

    void setV(real _newV);

    real W();

    void setW(real _newW);

    cellType type;

    ivec3 gridPos;

    std::vector<uint> m_paticleIDXs;

    uint getIDX(uint _length);

    vec3 getDeltaVel();

    vec3 getVelField() const;

private:

    real p;


    //Vel field should be half a cell offset from the pressure
    vec3 velField, oldVelField;

    real negativeDivergence;
};

#endif
