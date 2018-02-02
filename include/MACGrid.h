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
    /// @param uint _size : Amount of cells per axis
    /// @param real _cellWidth : The width (and therefore height and depth) of each cell in meters
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
    /// @brief Returns the indicies neighbors of a particular cell
    /// @param Cell _c : The cell whose neighbors we want to return
    //----------------------------------------------------------------------------------------------------------------------
    int *getNeighbors(Cell _c);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Returns the speed of the fastest particle
    //----------------------------------------------------------------------------------------------------------------------
    real getMaxSpeed();

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Returns a cell based on grid pos
    /// @param uint _i : The i pos of the cell
    /// @param uint _j : The j pos of the cell
    /// @param uint _k : The k pos of the cell
    //----------------------------------------------------------------------------------------------------------------------
    Cell &getCell(uint _i, uint _j, uint _k);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Returns a cell based on grid pos
    /// @param The position of the cell in te grid
    //----------------------------------------------------------------------------------------------------------------------
    Cell getCell(uvec3 _pos);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Sets the bounds of where the fluid is, where the solids are and where the air is
    /// @param uvec3 _b : The bottom corner of the fluid cells (0,0,0)
    /// @param uvec3 _c : The top corner of the fluid cells (1,1,1)
    //----------------------------------------------------------------------------------------------------------------------
    void initialiseCells(uvec3 _b, uvec3 _c);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Updates the types of each cell depending on the marker particles
    //----------------------------------------------------------------------------------------------------------------------
    void reclassifyCells();

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Returns the position of the cell in the grid
    /// @param Cell _c : The cell whose position we are querying
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
    /// @brief The amount of particles in the simulation
    //----------------------------------------------------------------------------------------------------------------------
    uint m_particleCount = 0;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Returns a random position in the bounds of a cell octant
    /// @param Cell _c : The cell we want the position in
    /// @param uint _count : The count of the particle (0-7) used to determine which octant we are in
    //----------------------------------------------------------------------------------------------------------------------
    vec3 getJitteredPos(Cell _c, uint _count);
};

#endif
