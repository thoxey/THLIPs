#include "utils.h"

namespace utility
{
uint getIndex(uint _length, MG_Cell _c)
{
    return _c.gridPos.x*_length*_length+_c.gridPos.y*_length+_c.gridPos.z;
}

real lerp(real _a, real _b, real _x)
{
    return (_a * (1.0 - _x)) + (_b * _x);
}

real invLerp(real _a, real _b, real _l)
{
    return -(_l - _a)/(_a+_b);
}

real trilerp(std::vector<real> V, real x, real y, real z)
{
    //Is are 1s and Os are 0s, based on notation here:
    //http://paulbourke.net/miscellaneous/interpolation/
    enum corners {OOO, IOO, OIO, OOI, IOI, OII, IIO, III};
    return (V[OOO] * (1-x) * (1-y) *1-z) +
           (V[IOO] * (1-y) * (1-z)) +
           (V[OIO] * (1-x) * y * (1-z)) +
           (V[OOI] * (1-x) * (1-y) * z) +
           (V[IOI] * x * (1-y) * z) +
           (V[OII] * (1-x) * y * z) +
           (V[IIO] * x * y * (1-z)) +
           (V[III] * x * y *z);
}

real randRange(real _max)
{
    std::random_device r;

    std::mt19937 e(r());

    std::uniform_real_distribution<> uniform_dist(0.0, _max);

    return uniform_dist(e);
}
real randRange(real _min, real _max)
{
    std::random_device r;

    std::mt19937 e(r());

    std::uniform_real_distribution<> uniform_dist(_min, _max);

    return uniform_dist(e);
}
}
