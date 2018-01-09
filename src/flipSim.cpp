#include "flipSim.h"

FlipSim::FlipSim(uvec3 _size):m_MACGrid(_size)
{
    m_dx = m_MACGrid.h;

}

FlipSim::FlipSim(uint _i, uint _j, uint _k): m_MACGrid(_i, _j, _k)
{
    ;
}


void FlipSim::updateGrid()
{
    for(MG_Particle p : m_MACGrid.m_particles)
    {

    }

}

void FlipSim::calculateNegativeDivergence()
{
    real scale = -1/m_dx;
    for(uint k = 0; k < m_kSize; k++)
    {
        for(uint j = 0; j < m_jSize; j++)
        {
            for(uint i = 0; i < m_iSize; i++)
            {
                MG_Cell c;
                if(m_MACGrid.checkForCell(i, j, k, c))
                    if(c.type == FLUID)
                    {
                        real uip1 = 0.0;
                        real vjp1 = 0.0;
                        real wkp1 = 0.0;
                        std::vector<MG_Cell> neighbors;

                        if(neighbors[RIGHT].type == FLUID)
                            uip1 = neighbors[RIGHT].u();
                        if(neighbors[UP].type == FLUID)
                            vjp1 = neighbors[UP].w();
                        if(neighbors[FORWARD].type == FLUID)
                            wkp1 = neighbors[FORWARD].v();

                        c.rhs = scale * (uip1 - c.u() + vjp1 - c.w() + wkp1 - c.v());

                        //figure 5.4 Bridon's book
                        if(neighbors[LEFT].type == SOLID)
                            c.rhs -= scale * c.u();
                        if(neighbors[RIGHT].type == SOLID)
                            c.rhs += scale * neighbors[RIGHT].u();

                        if(neighbors[DOWN].type == SOLID)
                            c.rhs -= scale * neighbors[DOWN].v();
                        if(neighbors[UP].type == SOLID)
                            c.rhs += scale * neighbors[UP].v();

                        if(neighbors[BACKWARD].type == SOLID)
                            c.rhs -= scale * neighbors[BACKWARD].w();
                        if(neighbors[FORWARD].type == SOLID)
                            c.rhs += scale * neighbors[FORWARD].w();

                    }
            }
        }
    }

}

void FlipSim::calculatePressure(real _dt)
{
    real scale = _dt / (m_density * m_dx * m_dx);
    //Sum up the fluid cells
    uint n_fluidCells = 0;
    std::vector<MG_Cell> fluidCells;
    for(MG_Cell c : m_MACGrid.m_cells)
        if(c.type == FLUID)
        {
            n_fluidCells++;
            fluidCells.push_back(c);
        }

    //Coeff Matrix
    Eigen::SparseMatrix<real> A(n_fluidCells, n_fluidCells);
    A.reserve(VectorX::Constant(7, n_fluidCells));

    //Fill with rhs values
    //Divergence
    VectorX b(n_fluidCells);
    for(uint i = 0; i < n_fluidCells; i++)
        b[i] = fluidCells[i].rhs;

    uint length = m_MACGrid.m_cells.size();

    //Bridsons book figure 5.5
    for(uint k = 0; k < m_kSize; k++)
    {
        for(uint j = 0; j < m_jSize; j++)
        {
            for(uint i = 0; i < m_iSize; i++)
            {
                MG_Cell c;
                uint aDiagIdx, aXidx, aYidx, aZidx;
                std::vector<MG_Cell> neighbors;
                if(m_MACGrid.checkForCell(i,j,k,c))
                {
                    aDiagIdx = getIndex(length, c);
                    neighbors = m_MACGrid.getNeighbors(c);
                }
                else
                    continue;

                aXidx = getIndex(length, neighbors[RIGHT]);

                aYidx = getIndex(length, neighbors[UP]);

                aZidx = getIndex(length, neighbors[FORWARD]);

                if(c.type == FLUID)
                {
                    if(neighbors[LEFT].type == FLUID)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                    }
                    if(neighbors[RIGHT].type == FLUID)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                        A.coeffRef(aDiagIdx, aXidx) -= scale;
                    }
                    else if(neighbors[RIGHT].type == AIR)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                    }

                    if(neighbors[DOWN].type == FLUID)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                    }
                    if(neighbors[UP].type == FLUID)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                        A.coeffRef(aDiagIdx, aYidx) -= scale;
                    }
                    else if(neighbors[UP].type == AIR)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                    }

                    if(neighbors[BACKWARD].type == FLUID)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                    }
                    if(neighbors[FORWARD].type == FLUID)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                        A.coeffRef(aDiagIdx, aZidx) -= scale;
                    }
                    else if(neighbors[FORWARD].type == AIR)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                    }
                }
            }
        }
    }

    //A sparse solver: time and/or sanity saver
    Eigen::SimplicialLLT<Eigen::SparseMatrix<real>> solver;
    //To store the result in
    VectorX p(n_fluidCells);
    solver.compute(A);
    p = solver.solve(b);

    for(uint i = 0; i < n_fluidCells; i++)
        fluidCells[i].p = p[i];
}


void FlipSim::applyPressure(real _dt)
{
    //fig 5.2
    real scale = _dt / (m_density*m_dx);

    for(uint k = 0; k < m_kSize; k++)
    {
        for(uint j = 0; j < m_jSize; j++)
        {
            for(uint i = 0; i < m_iSize; i++)
            {
                real u = 0.0;
                real v = 0.0;
                real w = 0.0;

                MG_Cell c;
                std::vector<MG_Cell> neighbors;
                if(m_MACGrid.checkForCell(i,j,k,c))
                    neighbors = m_MACGrid.getNeighbors(c);
                else
                    continue;

                //Update u
                if(neighbors[LEFT].type == FLUID || c.type == FLUID)
                {
                    if(!(neighbors[LEFT].type == SOLID || c.type == SOLID))
                    {
                        u -= scale * (c.p - neighbors[LEFT].p);
                    }
                }
                else
                {
                    //mark u(i,j,k) as unknown
                }

                //update v
                if(neighbors[DOWN].type == FLUID || c.type == FLUID)
                {
                    if(!(neighbors[DOWN].type == SOLID || c.type == SOLID))
                    {
                        v -= scale * (c.p - neighbors[DOWN].p);
                    }
                }
                else
                {
                    //mark v(i,j,k) as unknown
                }
                //update w
                if(neighbors[BACKWARD].type == FLUID || c.type == FLUID)
                {
                    if(!(neighbors[BACKWARD].type == SOLID || c.type == SOLID))
                    {
                        w -= scale * (c.p - neighbors[BACKWARD].p);
                    }
                }
                else
                {
                    //mark w(i,j,k) as unknown
                }
                c.velField = vec3(u,v,w);
            }
        }
    }
}

void FlipSim::step(real _dt)
{
    /*
    Basic algorithm from Bridson's Course Notes:
    (U is a velocity field)
    >for timeStep n = 0,1,2...
        >Determine Time Step dt to go from tn to tn+1

        Advect quantity q, through vector field Un for time interval dt
        >set Ua = advect(Un, dt, q)

        Add body forces
        >Ub = Ua +dtg

        To handle pressure/incompressibility we make a function called project
        >Un+1 = project(dt, Ub)
    */

    real t = 0;

    while(t < _dt)
    {
        //Calculate our substep
        real subStep = cfl();

        //Add gravity and stuff
        //Per particle calculation
        addBodyForce();

        updateGrid();

        //Per cell
        project(_dt);

        //Per particle
        //advect the velocity field
        advectVelocityField(_dt);


        t += subStep;
    }
}

void FlipSim::advectVelocityField(real _dt)
{
    for(MG_Particle p : m_MACGrid.m_particles)
        m_MACGrid.tracePoint(p.pos, _dt);
}

void FlipSim::addBodyForce()
{

}

void FlipSim::project(real _dt)
{
    calculateNegativeDivergence();
    calculatePressure(_dt);
    applyPressure(_dt);
}


real FlipSim::cfl()
{
    real ret = 0.0;
    ret = m_MACGrid.h/m_MACGrid.getMaxSpeed();
    return ret * m_k_cfl;
}
