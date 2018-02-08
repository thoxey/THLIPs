#include "utils.h"

namespace utils
{
uint getIndex(uint _length, uvec3 _pos)
{
    return _pos.x+(_length*_pos.y)+(_length*_length*_pos.z);
}
//----------------------------------------------------------------------------------------------------------------------
real lerp(real _a, real _b, real _x)
{
    return (_a * (1.0 - _x)) + (_b * _x);
}
//----------------------------------------------------------------------------------------------------------------------
real invLerp(real _a, real _b, real _l)
{
    return -(_l - _a)/(_a+_b);
}
//----------------------------------------------------------------------------------------------------------------------
real trilerp(real *V, real x, real y, real z)
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
            (V[III] * x * y * z);
}
//----------------------------------------------------------------------------------------------------------------------
real randRange(real _max)
{
    std::random_device r;

    std::mt19937 e(r());

    std::uniform_real_distribution<> uniform_dist(0.0, _max);

    return uniform_dist(e);
}
//----------------------------------------------------------------------------------------------------------------------
real randRange(real _min, real _max)
{
    std::random_device r;

    std::mt19937 e(r());

    std::uniform_real_distribution<> uniform_dist(_min, _max);

    return uniform_dist(e);
}
//----------------------------------------------------------------------------------------------------------------------
//....c
//.....
//b....
bool isInBounds(vec3 _a, vec3 _b, vec3 _c)
{
    if(_a.x > _b.x && _a.x < _c.x &&
            _a.y > _b.y && _a.y < _c.y &&
            _a.z > _b.z && _a.z < _c.z)
        return true;
    else
        return false;
}
//----------------------------------------------------------------------------------------------------------------------
void printvec(uvec3 _x, std::string _message)
{
    std::cout<<"X:"<<_x.x<<" Y: "<<_x.y<<" Z: "<<_x.z<<_message<<"\n";
}
//----------------------------------------------------------------------------------------------------------------------
void printvec(vec3 _x, std::string _message)
{
    std::cout<<"X:"<<_x.x<<" Y: "<<_x.y<<" Z: "<<_x.z<<_message<<"\n";
}
//----------------------------------------------------------------------------------------------------------------------
real trilinearHatKernel(vec3 _dist, real _dx)
{
    return hatFunction(_dist.x / _dx) * hatFunction(_dist.y / _dx) * hatFunction(_dist.z / _dx);
}
//----------------------------------------------------------------------------------------------------------------------
real hatFunction(real _r)
{
    real rAbs = std::abs(_r);
    if (rAbs <= 1)
        return 1.0 - rAbs;
    else
        return 0.0;
}
//----------------------------------------------------------------------------------------------------------------------
real divergentVelocity(real _v1, real _v2, real _dx)
{
    return (_v1-_v2)/_dx;
}
//----------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------
//A Runge-Kutta Method for solving Differential Equations
//of the form y'=f(x,y)
//--------------------------------------------------

real runge(real x, real y, real _dt)
{
    real K1    = (_dt * f(x,y));
    real K2    = (_dt * f((x + 1 / 2 * _dt), (y + 1 / 2 * K1)));
    real K3    = (_dt * f((x + 1 / 2 * _dt), (y + 1 / 2 * K2)));
    real K4    = (_dt * f((x + _dt), (y + K3)));
    real runge = (y + (1 / 6) * (K1 + 2 * K2 + 2 * K3 + K4));
    return runge;
}
//----------------------------------------------------------------------------------------------------------------------
real f(real x, real y)
{
    real f = x+y;
    return f;
}
//----------------------------------------------------------------------------------------------------------------------

}
