#include "utils.h"
uint getIndex(uint _length, MG_Cell _c)
{
    return _c.gridPos.x*_length*_length+_c.gridPos.y*_length+_c.gridPos.z;
}
