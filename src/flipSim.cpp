#include "flipSim.h"

FlipSim::FlipSim(uint _size, real _cellWidth, uvec3 _b, uvec3 _c):
    m_Grid(_size, _cellWidth),
    m_gridLength(_size)
{
    m_Grid.initialiseCells(_b, _c);
}

void FlipSim::step(real _dt)
{
    real subStep = cfl();
    m_Grid.reclassifyCells();
    updateGrid();
    addBodyForce(subStep);
    project(_dt);
    updateParticles();
    advectVelocityField(_dt);
}


real FlipSim::cfl()
{
    return 1.0;
}

void FlipSim::updateGrid()
{


    std::vector<real> uNum;
    uNum.reserve(m_gridLength*m_gridLength*m_gridLength);
    std::vector<real> uDen;
    uDen.reserve(m_gridLength*m_gridLength*m_gridLength);
    std::vector<real> vNum;
    vNum.reserve(m_gridLength*m_gridLength*m_gridLength);
    std::vector<real> vDen;
    vDen.reserve(m_gridLength*m_gridLength*m_gridLength);
    std::vector<real> wNum;
    wNum.reserve(m_gridLength*m_gridLength*m_gridLength);
    std::vector<real> wDen;
    wDen.reserve(m_gridLength*m_gridLength*m_gridLength);

    for(uint k = 0; k < m_gridLength; k++)
    {
        for(uint j = 0; j < m_gridLength; j++)
        {
            for(uint i = 0; i < m_gridLength; i++)
            {
                uNum[utils::getIndex(m_gridLength, uvec3(i,j,k))] = 0.0;
                uDen[utils::getIndex(m_gridLength, uvec3(i,j,k))] = 0.0;
                vNum[utils::getIndex(m_gridLength, uvec3(i,j,k))] = 0.0;
                vDen[utils::getIndex(m_gridLength, uvec3(i,j,k))] = 0.0;
                wNum[utils::getIndex(m_gridLength, uvec3(i,j,k))] = 0.0;
                wDen[utils::getIndex(m_gridLength, uvec3(i,j,k))] = 0.0;
            }
        }
    }

    for(Particle p : m_Grid.m_particles)
    {
        for(uint k = 0; k < m_gridLength; k++)
        {
            for(uint j = 0; j < m_gridLength; j++)
            {
                for(uint i = 0; i < m_gridLength; i++)
                {
                    vec3 cellPos = m_Grid.getCellPos(m_Grid.getCell(i,j,k));
                    real kernel = utils::trilinearHatKernel(p.pos - cellPos - vec3(0.5,0,0), m_Grid.m_h);

                    uNum[utils::getIndex(m_Grid.m_h, uvec3(i,j,k))] = p.vel.x * kernel;
                    uDen[utils::getIndex(m_Grid.m_h, uvec3(i,j,k))] = kernel;

                    vNum[utils::getIndex(m_Grid.m_h, uvec3(i,j,k))] = p.vel.y * kernel;
                    vDen[utils::getIndex(m_Grid.m_h, uvec3(i,j,k))] = kernel;

                    wNum[utils::getIndex(m_Grid.m_h, uvec3(i,j,k))] = p.vel.z * kernel;
                    wDen[utils::getIndex(m_Grid.m_h, uvec3(i,j,k))] = kernel;
                }
            }
        }
    }

    for(uint k = 0; k < m_gridLength; k++)
    {
        for(uint j = 0; j < m_gridLength; j++)
        {
            for(uint i = 0; i < m_gridLength; i++)
            {
                real newU = uNum[utils::getIndex(m_Grid.m_h, uvec3(i,j,k))] / uDen[utils::getIndex(m_Grid.m_h, uvec3(i,j,k))];
                real newV = vNum[utils::getIndex(m_Grid.m_h, uvec3(i,j,k))] / vDen[utils::getIndex(m_Grid.m_h, uvec3(i,j,k))];
                real newW = wNum[utils::getIndex(m_Grid.m_h, uvec3(i,j,k))] / wDen[utils::getIndex(m_Grid.m_h, uvec3(i,j,k))];
                m_Grid.getCell(i,j,k).updateVel(vec3(newU, newV, newW));
            }
        }
    }
}

void FlipSim::addBodyForce(real _dt)
{
    for(uint k = 0; k < m_gridLength; k++)
    {
        for(uint j = 0; j < m_gridLength; j++)
        {
            for(uint i = 0; i < m_gridLength; i++)
            {
                m_Grid.getCell(i,j,k).updateVel(m_g * _dt);
            }
        }
    }
}

void FlipSim::project(real _dt)
{
    enforceDirichlet();
    calculatePressure(_dt);
    enforceDirichlet();
}

void FlipSim::enforceDirichlet()
{
    for(uint k = 0; k < m_gridLength; k++)
    {
        for(uint j = 0; j < m_gridLength; j++)
        {
            for(uint i = 0; i < m_gridLength; i++)
            {
                Cell c = m_Grid.getCell(i,j,k);
                std::vector<Cell> neighbors = m_Grid.getNeighbors(c);
                //If we are neighboring/are a SOLID dont add velocity into the solid

                // X velocity
                if((c.type  == SOLID && c.U() > 0) || (neighbors[LEFT].type == SOLID && c.U() < 0))
                {
                    c.setU(0.0);
                }
                // Y velocity
                if((c.type  == SOLID && c.V() > 0) || (neighbors[DOWN].type == SOLID && c.V() < 0))
                {
                    c.setV(0.0);
                }
                // Z velocity
                if((c.type  == SOLID && c.W() > 0) || (neighbors[BACKWARD].type == SOLID && c.W() < 0))
                {
                    c.setW(0.0);
                }
            }
        }
    }
}

void FlipSim::calculatePressure(real _dt)
{
    uint n_fluidCells = 0;
    std::vector<int> fluidIDXs;
    for(Cell c : m_Grid.m_cells)
    {
        if(c.type == FLUID)
        {
            c.fluidIDX = n_fluidCells;
            fluidIDXs.push_back(n_fluidCells++);
        }
        else
        {
            fluidIDXs.push_back(-1);
        }
    }

    //Laplacian Matrix
    Eigen::SparseMatrix<real> A(n_fluidCells, n_fluidCells);
    A.reserve(Eigen::VectorXi::Constant(n_fluidCells, 5));

    //Divergence
    VectorX b(n_fluidCells);

    real scale = 1/(_dt*_dt);

    for(uint k = 0; k < m_gridLength; k++)
    {
        for(uint j = 0; j < m_gridLength; j++)
        {
            for(uint i = 0; i < m_gridLength; i++)
            {
                if(fluidIDXs[utils::getIndex(m_gridLength, uvec3(i,j,k))] != -1)
                {
                    uint idx = m_Grid.getCell(i,j,k).fluidIDX;
                    uint n_NonSolidNeighbors = 0;
                    std::vector<Cell> neighbors = m_Grid.getNeighbors(m_Grid.getCell(i, j, k));

                    if(neighbors[LEFT].type != SOLID)
                    {
                        if(neighbors[LEFT].type == FLUID)
                        {
                            A.insert(fluidIDXs[utils::getIndex(m_Grid.m_h, neighbors[LEFT].gridPos)], idx) = scale;
                        }
                        n_NonSolidNeighbors++;
                    }
                    if(neighbors[RIGHT].type != SOLID)
                    {
                        if(neighbors[RIGHT].type == FLUID)
                        {
                            A.insert(fluidIDXs[utils::getIndex(m_Grid.m_h, neighbors[RIGHT].gridPos)], idx) = scale;
                        }
                        n_NonSolidNeighbors++;
                    }

                    if(neighbors[DOWN].type != SOLID)
                    {
                        if(neighbors[DOWN].type == FLUID)
                        {
                            A.insert(fluidIDXs[utils::getIndex(m_Grid.m_h, neighbors[DOWN].gridPos)], idx) = scale;
                        }
                        n_NonSolidNeighbors++;
                    }
                    if(neighbors[UP].type != SOLID)
                    {
                        if(neighbors[UP].type == FLUID)
                        {
                            A.insert(fluidIDXs[utils::getIndex(m_Grid.m_h, neighbors[UP].gridPos)], idx) = scale;
                        }
                        n_NonSolidNeighbors++;
                    }

                    if(neighbors[BACKWARD].type != SOLID)
                    {
                        if(neighbors[BACKWARD].type == FLUID)
                        {
                            A.insert(fluidIDXs[utils::getIndex(m_Grid.m_h, neighbors[BACKWARD].gridPos)], idx) = scale;
                        }
                        n_NonSolidNeighbors++;
                    }
                    if(neighbors[FORWARD].type != SOLID)
                    {
                        if(neighbors[FORWARD].type == FLUID)
                        {
                            A.insert(fluidIDXs[utils::getIndex(m_Grid.m_h, neighbors[FORWARD].gridPos)], idx) = scale;
                        }
                        n_NonSolidNeighbors++;
                    }

                    A.insert(idx, idx) = -n_NonSolidNeighbors*scale;
                    real divergence = utils::divergentVelocity(neighbors[RIGHT].U(), m_Grid.getCell(i,j,k).U(), m_Grid.m_h);
                    divergence += utils::divergentVelocity(neighbors[UP].V(), m_Grid.getCell(i,j,k).V(), m_Grid.m_h);
                    divergence += utils::divergentVelocity(neighbors[FORWARD].W(), m_Grid.getCell(i,j,k).W(), m_Grid.m_h);
                    b[idx] = divergence;
                }
            }
        }
    }

    Eigen::SimplicialLLT<Eigen::SparseMatrix<real>> solver;
    //To store the result in
    VectorX p(n_fluidCells);
    solver.compute(A);
    p = solver.solve(b);

    for(uint k = 0; k < m_gridLength; k++)
    {
        for(uint j = 0; j < m_gridLength; j++)
        {
            for(uint i = 0; i < m_gridLength; i++)
            {
                uint idx = utils::getIndex(m_gridLength, uvec3(i,j,k));
                std::vector<Cell> neighbors = m_Grid.getNeighbors(m_Grid.getCell(i, j, k));

                real P = p(idx);
                real leftP = p(neighbors[LEFT].getIDX(m_gridLength));
                real downP = p(neighbors[DOWN].getIDX(m_gridLength));
                real backwardP = p(neighbors[BACKWARD].getIDX(m_gridLength));


                //vel_x - dt / density * pressure_diff_x / mac_grid.deltaX();
                real pressureDiffX = P-leftP;
                real pressureDiffY = P-downP;
                real pressureDiffZ = P-backwardP;

                real newU = _dt / m_density * pressureDiffX / m_Grid.m_h;
                real newV = _dt / m_density * pressureDiffY / m_Grid.m_h;
                real newW = _dt / m_density * pressureDiffZ / m_Grid.m_h;

                m_Grid.getCell(i,j,k).increaseVel(-vec3(newU, newV, newW));
            }
        }
    }
}

void FlipSim::updateParticles()
{
    for(Cell c : m_Grid.m_cells)
    {
        for(uint idx : c.m_paticleIDXs)
        {
            m_Grid.m_particles[idx].updateVel(c.getDeltaVel());
        }
    }
}

void FlipSim::advectVelocityField(real _dt)
{
    for(Particle p : m_Grid.m_particles)
        p.advect(_dt);
}

std::vector<Particle> FlipSim::getParticles()
{
    return m_Grid.m_particles;
}
