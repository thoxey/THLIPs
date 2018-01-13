#include "utils.h"

namespace utils
{
uint getIndex(uint _length, uvec3 _pos)
{
    return (_length*_length*_pos.x)+(_length*_pos.y)+_pos.z;
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

//....c
//.....
//b....
bool isInBounds(uvec3 _a, uvec3 _b, uvec3 _c)
{
    if(_a.x > _b.x && _a.x < _c.x &&
       _a.y > _b.y && _a.y < _c.y &&
       _a.z > _b.y && _a.z < _c.y)
        return true;
    else
        return false;
}

void printvec(uvec3 _x)
{
    std::cout<<"X:"<<_x.x<<" Y: "<<_x.y<<" Z: "<<_x.z<<"\n";
}

void printvec(vec3 _x)
{
    std::cout<<"X:"<<_x.x<<" Y: "<<_x.y<<" Z: "<<_x.z<<"\n";
}

real trilinearHatKernel(vec3 _dist, real _dx)
{
    return hatFunction(_dist.x / _dx) * hatFunction(_dist.y / _dx) * hatFunction(_dist.z / _dx);
}

real hatFunction(real _r) {
    real rAbs = std::abs(_r);
    if (rAbs <= 1) {
        return 1.0 - rAbs;
    } else {
        return 0.0;
    }
}

real divergentVelocity(real _v1, real _v2, real _dx)
{
    return (_v1-_v2)/_dx;
}

}
