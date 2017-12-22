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
void FlipSim::calculatePressure()
{
    using namespace Eigen;
    {
        //Sum up the fluid cells
        uint n_fluidCells = 0;
        for(MG_Cell c : m_MACGrid.m_cells)
            if(c.type == FLUID)
                n_fluidCells++;


        SparseMatrix<int> pressureMat(n_fluidCells, n_fluidCells);
        pressureMat.reserve(VectorXi::Constant(6, n_fluidCells));

        VectorX divergence(n_fluidCells);

        uint length = m_MACGrid.m_cells.size();
        for(MG_Cell c : m_MACGrid.m_cells)
        {
            int p = 6;
            std::vector<MG_Cell> neighbors = m_MACGrid.getNeighbors(c);
            for(MG_Cell neighbor : neighbors)
            {
                int nidx = getIndex(length, neighbor);
                if(neighbor.type != SOLID)
                    pressureMat.insert(nidx, nidx) = -1;
                else
                    p--;
            }
            uint i = getIndex(length, c);
            pressureMat.insert(i,i) = p;
        }

        //A sparse solver: time and/or sanity saver
        SimplicialLLT<SparseMatrix<int>> solver;
        //To store the result in
        VectorX x(n_fluidCells);

        solver.compute(pressureMat);
        //Commented out for now, causes error poss bc it is empty?
        //x = solver.solve(divergence);

    }

}

void FlipSim::applyPressure()
{
    for(uint k = 0; k < m_kSize; k++)
        for(uint j = 0; j < m_jSize; j++)
            for(uint i = 0; i < m_iSize; i++)
            {
                real x = 0.0;
                real y = 0.0;
                real z = 0.0;
                //Update u
                if(m_MACGrid.getCell(i-1, j, k).type == FLUID || m_MACGrid.getCell(i, j, k).type == FLUID)
                {
                    if(!(m_MACGrid.getCell(i-1, j, k).type == SOLID || m_MACGrid.getCell(i, j, k).type == SOLID))
                    {
                        //u(i,j,k) -= scale * (p(i,j,k) - p(i-1,j,k))
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
                        //v(i,j,k) -= scale * (p(i,j,k) - p(i,j-1,k))
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
                        //w(i,j,k) -= scale * (p(i,j,k) - p(i,j,k-1))
                    }
                }
                else
                {
                    //mark w(i,j,k) as unknown
                }
                m_MACGrid.getCell(i, j, k).velField = vec3(x,y,z);
            }
}

void FlipSim::step(real dt)
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

    while(t < dt)
    {
        //Calculate our substep
        real subStep = cfl();

        updateGrid();



        //Per cell
        project();

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

void FlipSim::project()
{

}


real FlipSim::cfl()
{
    real ret = 0.0;
    ret = m_MACGrid.h/m_MACGrid.getMaxSpeed();
    return ret * m_k_cfl;
}
