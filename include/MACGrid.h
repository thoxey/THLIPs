#ifndef __MACGRID__
#define __MACGRID__

#include "utils.h"
#include "cell.h"

//----------------------------------------------------------------------------------------------------------------------
/// @file MACGrid.h
/// @brief The Marker and Cell Grid class, based on the Robert Bridson model
/// @author Tom Hoxey
/// @version 1.0
/// @date 19/01/17 Initial version
//----------------------------------------------------------------------------------------------------------------------

class MACGrid
{
public:

    //Constructors
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief ctor
    //----------------------------------------------------------------------------------------------------------------------
    MACGrid(uint _size, real _cellWidth);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A container to store all the cells
    //----------------------------------------------------------------------------------------------------------------------
    std::vector<Cell> m_cells;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A container to store all the particles
    //----------------------------------------------------------------------------------------------------------------------
    std::vector<Particle> m_particles;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Returns the neighbors of a particular cell
    //----------------------------------------------------------------------------------------------------------------------
    std::vector<Cell> getNeighbors(Cell _c);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Returns the speed of the fastest particle
    //----------------------------------------------------------------------------------------------------------------------
    real getMaxSpeed();

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Returns a cell based on grid pos
    //----------------------------------------------------------------------------------------------------------------------
    Cell getCell(uint _i, uint _j, uint _k);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Returns a cell based on grid pos
    //----------------------------------------------------------------------------------------------------------------------
    Cell getCell(uvec3 _pos);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Sets the bounds of where the fluid is, where the solids are and where the air is
    //----------------------------------------------------------------------------------------------------------------------
    void initialiseCells(uvec3 _b, uvec3 _c);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    void reclassifyCells();

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    vec3 getCellPos(Cell _c);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Width of a cell in real units, using notaion from H&W paper
    //----------------------------------------------------------------------------------------------------------------------
    real m_h;

private:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Size of the grid in cells
    //----------------------------------------------------------------------------------------------------------------------
    uint m_gridLength;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Size of the grid in cells
    //----------------------------------------------------------------------------------------------------------------------
    uint m_cellCount;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    uint m_particleCount = 0;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    void initialiseCell(Cell _c, uvec3 _pos, cellType _t);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    void initialiseCellWithFluid(Cell _c);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    vec3 getJitteredPos(Cell _c, uint _count);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Gets the velocity of a particle, used in trace point
    //----------------------------------------------------------------------------------------------------------------------
    vec3 getVelocity(vec3 _v);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Used to get a 1D velocity using RK2
    //----------------------------------------------------------------------------------------------------------------------
    real getInterpolatedValue(real _a, real _b, real _c);
};

#endif
