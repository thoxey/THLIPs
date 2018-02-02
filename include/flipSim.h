#ifndef __FLIPSIM__
#define __FLIPSIM__

#include "MACGrid.h"

//----------------------------------------------------------------------------------------------------------------------
/// @file flipSim.h
/// @brief The simulation header, contains the main simulation code
/// @author Tom Hoxey
/// @version 1.0
/// @date 19/01/17 Initial version
//----------------------------------------------------------------------------------------------------------------------

class FlipSim
{
public:
    //Constructors
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief ctor
    //----------------------------------------------------------------------------------------------------------------------
    FlipSim(uint _size, real _cellWidth, uvec3 _b, uvec3 _c);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Steps the simulation forward once
    /// @param real _dt : The real time passed since the last step
    //----------------------------------------------------------------------------------------------------------------------
    void step(real _dt);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Returns the particles active in the simulation
    //----------------------------------------------------------------------------------------------------------------------
    std::vector<Particle> getParticles() const;

private:
    //The steps of the fluid simulation
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Updates the grid cells depending on the particles
    //----------------------------------------------------------------------------------------------------------------------
    void updateGrid();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Calculate the negative divergence
    //----------------------------------------------------------------------------------------------------------------------
    void calculateNegativeDivergence();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Enforces the dirichlet condition i.e fluid cant move through a solid
    //----------------------------------------------------------------------------------------------------------------------
    void enforceDirichlet();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Calculate the pressure in the cells
    //----------------------------------------------------------------------------------------------------------------------
    void calculatePressure(real _dt);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Applies the pressure calculations to each of the cells depending on thier contents
    //----------------------------------------------------------------------------------------------------------------------
    void applyPressure(real _dt);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief The advection portion of the semi-lagrangian algorithm
    //----------------------------------------------------------------------------------------------------------------------
    void advectVelocityField(real _dt);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Updates based on the body forces, such as gravity
    //----------------------------------------------------------------------------------------------------------------------
    void addBodyForce(real _dt);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief The projection step, handles the pressure and incrompressibility portion of the fluid algoritm
    //----------------------------------------------------------------------------------------------------------------------
    void project(real _dt);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    void updateParticles();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    void wrangleParticles();

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Calculates the time step of the fluid equation (how many steps per frame)
    //----------------------------------------------------------------------------------------------------------------------
    real cfl();

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    vec3 interpVelocity();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    vec3 getSampledVelocity(const vec3 &_p);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief The grid that we store the particles in, and calculate based on
    //----------------------------------------------------------------------------------------------------------------------
    MACGrid m_Grid;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief The dimensions of the grid, cells in the x,y,z directions
    //----------------------------------------------------------------------------------------------------------------------
    uint m_gridLength;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief The user defined scalar for the cfl function
    //----------------------------------------------------------------------------------------------------------------------
    real m_k_cfl = 1.0;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    real m_dx = 0.0;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Density of the fluid, 1000 by default as in water
    //----------------------------------------------------------------------------------------------------------------------
    real m_density = 1000.0;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Density of the air surrounding the fluid
    //----------------------------------------------------------------------------------------------------------------------
    real m_airDensity = 1.3;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    const vec3 m_g = vec3(0.0, -9.81, 0.0);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    bool m_firstStep = true;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief
    //----------------------------------------------------------------------------------------------------------------------
    real m_xMin, m_xMax, m_yMin, m_yMax, m_zMin, m_zMax;
};

#endif
