#include "flipSim.h"

FlipSim::FlipSim(uvec3 _size):m_MACGrid(_size)
{
    ;

}

FlipSim::FlipSim(uint _i, uint _j, uint _k): m_MACGrid(_i, _j, _k)
{
    ;
}


void FlipSim::updateGrid()
{

    for(MG_Cell c : m_MACGrid.m_cells)
        c.layer = -1;
    /*
      // update cells that currently have fluid in them
      for each marker particle, P
        if the cell, C, containing the center of P does not exist
            if C is within the simulation bounds
                create C and put it in the hash table
                set the cell type of C to “fluid”
                C.layer = 0
        else if C is not part of a solid object
            Set the cell type for C to “fluid”
            C.layer = 0
    */
    for(MG_Particle p : m_MACGrid.m_particles)
    {
        MG_Cell c;
        if(p.cellidx < m_MACGrid.m_cells.size())
            c = m_MACGrid.m_cells[p.cellidx];
        if(!m_MACGrid.checkForCell(p))
        {
            if(m_MACGrid.checkInBounds(c))
            {
                //Create new cell c
                c.type = FLUID;
                c.layer = 0;
            }
        }
        else if(c.type != SOLID)
        {
            c.type = FLUID;
            c.layer = 0;
        }
    }

    /*
    // create a buffer zone around the fluid
    for i = 1 to max(2, ⌈kc f l ⌉)
        for each liquid or air cell, C, such that C.layer == i−1
            for each of the six neighbors of C, N if N already exists in the hash table
                if N.layer == −1 and N is not solid
                    set the cell type of N to “air”
                    N.layer = i
                else
                    create N and put it in the hash table
                    N.layer = i
                    if N is in the simulation bounds
                        set the cell type of N to “air”
                    else
                        set the cell type of N to “solid”

    delete any cells with layer == −1
     */
    for(int i = 0; std::max(2, (int)std::ceil(m_k_cfl)); i++)
    {
        for(MG_Cell c : m_MACGrid.m_cells)
        {
            if(c.layer != -1)
                continue;

            std::vector<MG_Cell> neighbors = m_MACGrid.getNeighbors(c);
            for(MG_Cell neighbor : neighbors)
            {
                if(true)//if N exists in the hash table
                {
                    if(neighbor.layer == -1 && neighbor.type != SOLID)
                    {
                        neighbor.type = AIR;
                        neighbor.layer = i;
                    }
                }
                else
                {
                    //Create neighbor
                    neighbor.layer = i;
                    if(m_MACGrid.checkInBounds(neighbor))
                    {
                        neighbor.type = AIR;
                    }
                    else
                    {
                        neighbor.type = SOLID;
                    }
                }
            }
        }
    }

    //delete any cells with layer == -1
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
                MG_Cell c = m_MACGrid.getCell(i, j, k);
                if(c.type == FLUID)
                {
                    real uip1 = 0.0;
                    real vjp1 = 0.0;
                    real wkp1 = 0.0;
                    if(m_MACGrid.getCell(i+1, j, k).type == FLUID)
                        uip1 = m_MACGrid.getCell(i+1, j, k).u();
                    if(m_MACGrid.getCell(i, j+1, k).type == FLUID)
                        vjp1 = m_MACGrid.getCell(i, j+1, k).w();
                    if(m_MACGrid.getCell(i, j, k+1).type == FLUID)
                        wkp1 = m_MACGrid.getCell(i, j, k+1).v();

                    c.rhs = scale * (uip1 - c.u() + vjp1 - c.w() + wkp1 - c.v());

                    //figure 5.4 Bridon's book
                    if(m_MACGrid.getCell(i-1, j, k).type == SOLID)
                        c.rhs -= scale * c.u();
                    if(m_MACGrid.getCell(i+1, j, k).type == SOLID)
                        c.rhs += scale * m_MACGrid.getCell(i+1, j, k).u();

                    if(m_MACGrid.getCell(i, j-1, k).type == SOLID)
                        c.rhs -= scale * m_MACGrid.getCell(i, j, k).v();
                    if(m_MACGrid.getCell(i, j+1, k).type == SOLID)
                        c.rhs += scale * m_MACGrid.getCell(i, j+1, k).v();

                    if(m_MACGrid.getCell(i,j,k-1).type == SOLID)
                        c.rhs -= scale * m_MACGrid.getCell(i,j,k).w();
                    if(m_MACGrid.getCell(i,j,k+1).type == SOLID)
                        c.rhs += scale * m_MACGrid.getCell(i, j, k+1).w();

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
    std::vector<real> rhsVec;
    for(MG_Cell c : m_MACGrid.m_cells)
        if(c.type == FLUID)
        {
            n_fluidCells++;
            rhsVec.push_back(c.rhs);
        }

    //Coeff Matrix
    Eigen::SparseMatrix<real> A(n_fluidCells, n_fluidCells);
    A.reserve(VectorX::Constant(7, n_fluidCells));

    //Fill with rhs values
    //Divergence
    VectorX b(n_fluidCells);
    for(uint i = 0; i < n_fluidCells; i++)
        b[i] = rhsVec[i];

    uint length = m_MACGrid.m_cells.size();

    //Bridsons book figure 5.5
    for(uint k = 0; k < m_kSize; k++)
    {
        for(uint j = 0; j < m_jSize; j++)
        {
            for(uint i = 0; i < m_iSize; i++)
            {
                uint aDiagIdx = getIndex(length, m_MACGrid.getCell(i,j,k));
                uint aXidx = getIndex(length, m_MACGrid.getCell(i+1,j,k));
                uint aYidx = getIndex(length, m_MACGrid.getCell(i,j+1,k));
                uint aZidx = getIndex(length, m_MACGrid.getCell(i,j,k+1));
                if(m_MACGrid.getCell(i,j,k).type == FLUID)
                {
                    if(m_MACGrid.getCell(i-1,j,k).type == FLUID)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                    }
                    if(m_MACGrid.getCell(i+1,j,k).type == FLUID)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                        A.coeffRef(aDiagIdx, aXidx) -= scale;
                    }
                    else if(m_MACGrid.getCell(i+1,j,k).type == AIR)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                    }

                    if(m_MACGrid.getCell(i,j-1,k).type == FLUID)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                    }
                    if(m_MACGrid.getCell(i,j+1,k).type == FLUID)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                        A.coeffRef(aDiagIdx, aYidx) -= scale;
                    }
                    else if(m_MACGrid.getCell(i,j+1,k).type == AIR)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                    }

                    if(m_MACGrid.getCell(i,j,k-1).type == FLUID)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                    }
                    if(m_MACGrid.getCell(i,j,k+1).type == FLUID)
                    {
                        A.coeffRef(aDiagIdx, aDiagIdx) += scale;
                        A.coeffRef(aDiagIdx, aZidx) -= scale;
                    }
                    else if(m_MACGrid.getCell(i,j,k+1).type == AIR)
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
                //Update u
                if(m_MACGrid.getCell(i-1, j, k).type == FLUID || m_MACGrid.getCell(i, j, k).type == FLUID)
                {
                    if(!(m_MACGrid.getCell(i-1, j, k).type == SOLID || m_MACGrid.getCell(i, j, k).type == SOLID))
                    {
                        u -= scale * (m_MACGrid.getCell(i,j,k).p - m_MACGrid.getCell(i-1,j,k).p);
                    }
                }
                else
                {
                    //mark u(i,j,k) as unknown
                }

                //update v
                if(m_MACGrid.getCell(i, j-1, k).type == FLUID || m_MACGrid.getCell(i, j, k).type == FLUID)
                {
                    if(!(m_MACGrid.getCell(i, j-1, k).type == SOLID || m_MACGrid.getCell(i, j, k).type == SOLID))
                    {
                        v -= scale * (m_MACGrid.getCell(i,j,k).p - m_MACGrid.getCell(i,j-1,k).p);
                    }
                }
                else
                {
                    //mark v(i,j,k) as unknown
                }
                //update w
                if(m_MACGrid.getCell(i, j, k-1).type == FLUID || m_MACGrid.getCell(i, j, k).type == FLUID)
                {
                    if(!(m_MACGrid.getCell(i, j, k-1).type == SOLID || m_MACGrid.getCell(i, j, k).type == SOLID))
                    {
                        w -= scale * (m_MACGrid.getCell(i,j,k).p - m_MACGrid.getCell(i,j,k-1).p);
                    }
                }
                else
                {
                    //mark w(i,j,k) as unknown
                }
                m_MACGrid.getCell(i, j, k).velField = vec3(u,v,w);
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

        updateGrid();



        //Per cell
        project(_dt);

        //Per particle
        //advect the velocity field
        advectVelocityField();

        //Add gravity and stuff
        //Per particle calculation
        addBodyForce();

        t += subStep;
    }
}

void FlipSim::advectVelocityField()
{
    for(uint k = 0; k < m_kSize; k++)
        for(uint j = 0; j < m_jSize; j++)
            for(uint i = 0; i < m_iSize; i++)
            {
                //do advection
            }
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
