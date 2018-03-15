#include "cell.h"

Cell::Cell(ivec3 _pos):
    gridPos(_pos),
    velField(vec3(0.0,0.0,0.0)),
    oldVelField(velField),
    p(0.0),
    negativeDivergence(0.0),
    type(SOLID)
{
}

real Cell::U()
{
    return velField.x;
}

void Cell::setU(real _newU)
{
    velField.x = _newU;
}

real Cell::W()
{
    return velField.z;
}

void Cell::setV(real _newV)
{
    velField.y = _newV;
}

real Cell::V()
{
    return velField.y;
}

void Cell::setW(real _newW)
{
    velField.z = _newW;
}

void Cell::updateVel(vec3 _newVel)
{
    velField = _newVel;
}

void Cell::increaseVel(vec3 _addition)
{
    velField += _addition;
}

uint Cell::getIDX(uint _length)
{
    return utils::getIndex(_length, gridPos);
}

vec3 Cell::getDeltaVel()
{
    return velField - oldVelField;
}

vec3 Cell::getVelField() const
{
    return velField;
}

real Cell::getDensity() const
{
    return density;
}

void Cell::setDensity(const real &value)
{
    density = value;
}
