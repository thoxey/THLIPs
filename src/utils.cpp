#include "utils.h"



uint getIndex(uint _length, MG_Cell _c)
{
    return _c.gridPos.x*_length*_length+_c.gridPos.y*_length+_c.gridPos.z;
}

real lerp(real a, real b, real f)
{
    return (a * (1.0 - f)) + (b * f);
}

real randRange(real _max)
{
    std::random_device r;

    std::mt19937 e(r());

    std::uniform_real_distribution<> uniform_dist(0.0, _max);

    return uniform_dist(e);
}
