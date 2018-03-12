#include "flipSim.h"
//----------------------------------------------------------------------------------------------------------------------
FlipSim::FlipSim(uint _size, real _cellWidth, uvec3 _b, uvec3 _c):
    m_Grid(_size, _cellWidth),
    m_gridLength(_size)
{
    m_Grid.initialiseCells(_b, _c);

    m_xMin = _cellWidth;
    m_yMin = _cellWidth;
    m_zMin = _cellWidth;

    m_xMax = _cellWidth * _size;
    m_yMax = _cellWidth * _size;
    m_zMax = _cellWidth * _size;
}
//----------------------------------------------------------------------------------------------------------------------
void FlipSim::step(real _dt)
{
    real t = 0.0;
    while(t < _dt)
    {
        real subStep = cfl();
        //#define VERBOSE_OUTPUT
#ifdef VERBOSE_OUTPUT
        std::cout<<"Reclassifying cells... \n";
#endif
        m_Grid.reclassifyCells();
#ifdef VERBOSE_OUTPUT
        std::cout<<"Updating Grid... \n";
#endif
        updateGrid();
#ifdef VERBOSE_OUTPUT
        std::cout<<"Adding Body Force... \n";
#endif
        addBodyForce(t);
#ifdef VERBOSE_OUTPUT
        std::cout<<"Projecting... \n";
#endif
        project(t);
#ifdef VERBOSE_OUTPUT
        std::cout<<"Updating Particles... \n";
#endif
        updateParticles();
#ifdef VERBOSE_OUTPUT
        std::cout<<"Advecting... \n";
#endif
        advectVelocityField(t);
#ifdef VERBOSE_OUTPUT
        std::cout<<"Wrangling... \n";
#endif
        wrangleParticles();
        t += _dt;//subStep;
    }
}
//----------------------------------------------------------------------------------------------------------------------
real FlipSim::cfl()
{
    real umax = m_Grid.getMaxSpeed() + (std::sqrt(5*m_Grid.m_h*(-1*m_g.y)));
    real ret = (5*m_Grid.m_h)/umax;
    std::cout<<"CFL: "<<ret * m_k_cfl<<"\n";
    return ret * m_k_cfl;
}
//----------------------------------------------------------------------------------------------------------------------
void FlipSim::updateGrid()
{
    uint gridCubed = m_gridLength*m_gridLength*m_gridLength;
    real uNum[gridCubed];
    real uDen[gridCubed];
    real vNum[gridCubed];
    real vDen[gridCubed];
    real wNum[gridCubed];
    real wDen[gridCubed];

    for(uint k = 0; k < m_gridLength; k++)
    {
        for(uint j = 0; j < m_gridLength; j++)
        {
            for(uint i = 0; i < m_gridLength; i++)
            {
                uint idx = utils::getIndex(m_gridLength, uvec3(i,j,k));
                uNum[idx] = 0.0;
                uDen[idx] = 0.0;
                vNum[idx] = 0.0;
                vDen[idx] = 0.0;
                wNum[idx] = 0.0;
                wDen[idx] = 0.0;
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
                    //FIX THIS===========================================================<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    real kernel = 1.0f; //utils::trilinearHatKernel(p.pos - cellPos - vec3(0.5,0,0), m_Grid.m_h);
                    //std::cout<<kernel<<" - Kernel \n";
                    uint idx = utils::getIndex(m_gridLength, uvec3(i,j,k));

                    uNum[idx] = p.vel.x * kernel;
                    uDen[idx] = kernel;

                    vNum[idx] = p.vel.y * kernel;
                    vDen[idx] = kernel;

                    wNum[idx] = p.vel.z * kernel;
                    wDen[idx] = kernel;
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
                uint idx = utils::getIndex(m_gridLength, uvec3(i,j,k));
                real newU = uNum[idx] / uDen[idx];
                real newV = vNum[idx] / vDen[idx];
                real newW = wNum[idx] / wDen[idx];
                m_Grid.getCell(i,j,k).updateVel(vec3(newU, newV, newW));
            }
        }
    }
}
//----------------------------------------------------------------------------------------------------------------------
void FlipSim::addBodyForce(real _dt)
{
    for(auto&& c : m_Grid.m_cells)
        c.increaseVel(m_g * _dt);
}
//----------------------------------------------------------------------------------------------------------------------
void FlipSim::project(real _dt)
{
    enforceDirichlet();
    calculatePressure(_dt);
    enforceDirichlet();
}
//----------------------------------------------------------------------------------------------------------------------
void FlipSim::enforceDirichlet()
{
    for(uint k = 0; k < m_gridLength; k++)
    {
        for(uint j = 0; j < m_gridLength; j++)
        {
            for(uint i = 0; i < m_gridLength; i++)
            {
                uint idx = utils::getIndex(m_gridLength, uvec3(i,j,k));
                int * neighbors = m_Grid.getNeighbors(m_Grid.m_cells[idx]);
                //If we are neighboring/are a SOLID dont add velocity into the solid

                if(m_Grid.m_cells[idx].type  == SOLID)
                {
                    m_Grid.m_cells[idx].setU(0.0);
                    m_Grid.m_cells[idx].setV(0.0);
                    m_Grid.m_cells[idx].setW(0.0);
                    continue;
                }

                // X velocity
                if(m_Grid.m_cells[neighbors[LEFT]].type == SOLID && m_Grid.m_cells[idx].U() < 0.0)
                {
                    //utils::printvec(m_Grid.m_cells[idx].getVelField(), "Dirichlet enforced bc LEFT");
                    m_Grid.m_cells[idx].setU(0.0);
                }
                if(m_Grid.m_cells[neighbors[RIGHT]].type == SOLID && m_Grid.m_cells[idx].U() > 0.0)
                {
                    //utils::printvec(m_Grid.m_cells[idx].getVelField(), "Dirichlet enforced bc RIGHT");
                    m_Grid.m_cells[idx].setU(0.0);
                }
                // Y velocity
                if(m_Grid.m_cells[neighbors[DOWN]].type == SOLID && m_Grid.m_cells[idx].V() > 0.0)
                {
                    //utils::printvec(m_Grid.m_cells[idx].getVelField(), "Dirichlet enforced bc DOWN");
                    m_Grid.m_cells[idx].setV(0.0);
                }
                if(m_Grid.m_cells[neighbors[UP]].type == SOLID && m_Grid.m_cells[idx].V() < 0.0)
                {
                    //utils::printvec(m_Grid.m_cells[idx].getVelField(), "Dirichlet enforced bc UP");
                    m_Grid.m_cells[idx].setV(0.0);
                }
                // Z velocity
                if(m_Grid.m_cells[neighbors[BACKWARD]].type == SOLID && m_Grid.m_cells[idx].W() < 0.0)
                {
                    //utils::printvec(m_Grid.m_cells[idx].getVelField(), "Dirichlet enforced bc BACKWARD");
                    m_Grid.m_cells[idx].setW(0.0);
                }
                if(m_Grid.m_cells[neighbors[FORWARD]].type == SOLID && m_Grid.m_cells[idx].W() > 0.0)
                {
                    //utils::printvec(m_Grid.m_cells[idx].getVelField(), "Dirichlet enforced bc FORWARD");
                    m_Grid.m_cells[idx].setW(0.0);
                }
            }
        }
    }
}
//----------------------------------------------------------------------------------------------------------------------
void FlipSim::calculatePressure(real _dt)
{
    uint n_fluidCells = 0;
    std::vector<int> fluidIDXs;
    for(auto&& c : m_Grid.m_cells)
    {
        if(c.type == FLUID)
        {
            fluidIDXs.push_back(n_fluidCells++);
        }
        else
        {
            fluidIDXs.push_back(-1);
        }
    }

    assert(n_fluidCells > 0);

    //Laplacian Matrix
    Eigen::SparseMatrix<real> A(n_fluidCells, n_fluidCells);
    A.reserve(Eigen::VectorXi::Constant(n_fluidCells, 5));

    //Divergence
    VectorX b(n_fluidCells);

    real scale = _dt/(m_Grid.m_h*m_Grid.m_h);

    for(uint k = 0; k < m_gridLength; k++)
    {
        for(uint j = 0; j < m_gridLength; j++)
        {
            for(uint i = 0; i < m_gridLength; i++)
            {
                if(fluidIDXs[utils::getIndex(m_gridLength, uvec3(i,j,k))] != -1)
                {
                    uint idx = fluidIDXs[utils::getIndex(m_gridLength, uvec3(i,j,k))];
                    uint n_NonSolidNeighbors = 0;
                    int * neighbor = m_Grid.getNeighbors(m_Grid.getCell(i, j, k));

                    if(m_Grid.m_cells[neighbor[LEFT]].type != SOLID)
                    {
                        if(m_Grid.m_cells[neighbor[LEFT]].type == FLUID)
                        {
                            A.insert(fluidIDXs[neighbor[LEFT]], idx) = scale;
                        }
                        n_NonSolidNeighbors++;
                    }

                    if(m_Grid.m_cells[neighbor[RIGHT]].type != SOLID)
                    {
                        if(m_Grid.m_cells[neighbor[RIGHT]].type == FLUID)
                        {
                            A.insert(fluidIDXs[neighbor[RIGHT]], idx) = scale;
                        }
                        n_NonSolidNeighbors++;
                    }

                    if(m_Grid.m_cells[neighbor[DOWN]].type != SOLID)
                    {
                        if(m_Grid.m_cells[neighbor[DOWN]].type == FLUID)
                        {
                            A.insert(fluidIDXs[neighbor[DOWN]], idx) = scale;
                        }
                        n_NonSolidNeighbors++;
                    }

                    if(m_Grid.m_cells[neighbor[UP]].type != SOLID)
                    {
                        if(m_Grid.m_cells[neighbor[UP]].type == FLUID)
                        {
                            A.insert(fluidIDXs[neighbor[UP]], idx) = scale;
                        }
                        n_NonSolidNeighbors++;
                    }

                    if(m_Grid.m_cells[neighbor[BACKWARD]].type != SOLID)
                    {
                        if(m_Grid.m_cells[neighbor[BACKWARD]].type == FLUID)
                        {
                            A.insert(fluidIDXs[neighbor[BACKWARD]], idx) = scale;
                        }
                        n_NonSolidNeighbors++;
                    }

                    if(m_Grid.m_cells[neighbor[FORWARD]].type != SOLID)
                    {
                        if(m_Grid.m_cells[neighbor[FORWARD]].type == FLUID)
                        {
                            A.insert(fluidIDXs[neighbor[FORWARD]], idx) = scale;
                        }
                        n_NonSolidNeighbors++;
                    }

                    A.coeffRef(idx, idx) = -n_NonSolidNeighbors*scale;

                    real uRight = m_Grid.m_cells[neighbor[RIGHT]].U();
                    real u = m_Grid.m_cells[utils::getIndex(m_gridLength, uvec3(i,j,k))].U();
                    real vUp = m_Grid.m_cells[neighbor[UP]].V();
                    real v = m_Grid.m_cells[utils::getIndex(m_gridLength, uvec3(i,j,k))].V();
                    real wForward = m_Grid.m_cells[neighbor[FORWARD]].W();
                    real w = m_Grid.m_cells[utils::getIndex(m_gridLength, uvec3(i,j,k))].W();

                    real divergence = scale * (uRight - u + vUp - v + wForward - w);

                    //figure 5.4 Bridon's book
                    if(m_Grid.m_cells[neighbor[LEFT]].type == SOLID)
                        divergence -= scale * u;
                    if(m_Grid.m_cells[neighbor[RIGHT]].type == SOLID)
                        divergence += scale * uRight;

                    if(m_Grid.m_cells[neighbor[DOWN]].type == SOLID)
                        divergence -= scale * v;
                    if(m_Grid.m_cells[neighbor[UP]].type == SOLID)
                        divergence += scale * vUp;

                    if(m_Grid.m_cells[neighbor[BACKWARD]].type == SOLID)
                        divergence -= scale * w;
                    if(m_Grid.m_cells[neighbor[FORWARD]].type == SOLID)
                        divergence += scale * wForward;

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

//    for(uint i = 0; i<p.rows(); i++)
//    {
//        if(isnan(p(i)))
//            p(i) = 0;
//    }

    for(uint k = 0; k < m_gridLength; k++)
    {
        for(uint j = 0; j < m_gridLength; j++)
        {
            for(uint i = 0; i < m_gridLength; i++)
            {
                uint idx = fluidIDXs[utils::getIndex(m_gridLength, uvec3(i,j,k))];
                int * neighbors = m_Grid.getNeighbors(m_Grid.getCell(i, j, k));

                if(idx < p.rows())
                {
                    real P = p(idx);

                    real leftP = 0.0;
                    uint leftIDX = m_Grid.m_cells[neighbors[LEFT]].getIDX(m_gridLength);
                    if(leftIDX < n_fluidCells)
                        leftP = p(leftIDX);

                    real downP = 0.0;
                    uint downIDX = m_Grid.m_cells[neighbors[DOWN]].getIDX(m_gridLength);
                    if(downIDX < n_fluidCells)
                        downP = p(downIDX);

                    real backwardP = 0.0;
                    uint backwardIDX = m_Grid.m_cells[neighbors[BACKWARD]].getIDX(m_gridLength);
                    if(backwardIDX < n_fluidCells)
                        backwardP = p(backwardIDX);

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
}
//----------------------------------------------------------------------------------------------------------------------
vec3 FlipSim::getSampledVelocity(const vec3 &_p)
{
    vec3 p = _p;
    real * Us = (real *) malloc(8*sizeof(real));
    real * Vs = (real *) malloc(8*sizeof(real));
    real * Ws = (real *) malloc(8*sizeof(real));
    for(int k = 0; k <2; k++)
    {
        for(int j = 0; j <2; j++)
        {
            for(int i = 0; i <2; i++)
            {
                p.x += (i+i-1) * m_Grid.m_h;
                p.y += (j+j-1) * m_Grid.m_h;
                p.z += (k+k-1) * m_Grid.m_h;

                int x = std::floor(p.x);
                x = std::max(x, 0);
                x = std::min(x, (int)m_gridLength-1);
                int y = std::floor(p.y);
                y = std::max(y, 0);
                y = std::min(y, (int)m_gridLength-1);
                int z = std::floor(p.z);
                z = std::max(z, 0);
                z = std::min(z, (int)m_gridLength-1);

                uint idx = utils::getIndex(m_gridLength, uvec3(x,y,z));

                vec3 tmpVel = m_Grid.m_cells[idx].getDeltaVel();

                real tmpU = tmpVel.x;
                real tmpV = tmpVel.y;
                real tmpW = tmpVel.z;

                Us[i+(j*2)+(k*2*2)] = tmpU;
                Vs[i+(j*2)+(k*2*2)] = tmpV;
                Ws[i+(j*2)+(k*2*2)] = tmpW;
            }
        }
    }
    real newU = utils::trilerp(Us, p.x, p.y, p.z);
    real newV = utils::trilerp(Vs, p.x, p.y, p.z);
    real newW = utils::trilerp(Ws, p.x, p.y, p.z);

    free(Us);
    free(Vs);
    free(Ws);

    return vec3(newU, newV, newW);
}
//----------------------------------------------------------------------------------------------------------------------
std::vector<Particle> FlipSim::getParticles() const
{
    return m_Grid.m_particles;
}
//----------------------------------------------------------------------------------------------------------------------
void FlipSim::updateParticles()
{
    for(auto&& c : m_Grid.m_cells)
    {
        for(uint idx : c.m_paticleIDXs)
        {
            //            vec3 sampledVel = getSampledVelocity(m_Grid.m_particles[idx].pos);
            //            utils::printvec(sampledVel);
            //            vec3 cellCol = c.getVelField();
            //            m_Grid.m_particles[idx].updateVel(sampledVel, cellCol);

            vec3 cellCol = c.getVelField();
            m_Grid.m_particles[idx].updateVel(cellCol, cellCol);
        }
    }
}
//----------------------------------------------------------------------------------------------------------------------
void FlipSim::advectVelocityField(real _dt)
{
    for(auto&& p : m_Grid.m_particles)
        p.advect(_dt);
}
//----------------------------------------------------------------------------------------------------------------------
void FlipSim::wrangleParticles()
{
    for(auto&& p : m_Grid.m_particles)
    {
        if(p.pos.x < m_xMin)
        {
            p.pos.x = m_xMin;
        }
        if(p.pos.x > m_xMax)
        {
            p.pos.x = m_xMax;
        }

        if(p.pos.y < m_yMin)
        {
            p.pos.y = m_yMin;
        }
        if(p.pos.y > m_yMax)
        {
            p.pos.y = m_yMax;
        }

        if(p.pos.z < m_zMin)
        {
            p.pos.z = m_zMin;
        }
        if(p.pos.z > m_zMax)
        {
            p.pos.z = m_zMax;
        }

    }
    //enforceDirichlet();
}
//----------------------------------------------------------------------------------------------------------------------
