#include "flipSim.h"

FlipSim::FlipSim(uvec3 _size, real _cellWidth, uvec3 _b, uvec3 _c):m_MACGrid(_size, _cellWidth)
{
    m_dx = _cellWidth;

    m_iSize = _size.x;
    m_jSize = _size.y;
    m_kSize = _size.z;

    m_MACGrid.initialiseCells(_b, _c);
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

        //Add gravity and stuff
        //Per particle calculation
        addBodyForce(subStep);

        //Per cell
        project(subStep);

        updateParticles();

        //Per particle
        //advect the velocity field
        advectVelocityField(_dt);

        for(uint i = 0; i < m_MACGrid.m_cells.size(); i++)
        {
            MG_Cell c = m_MACGrid.getCell(m_MACGrid.m_cells[i]);
            //utility::printvec(c.oldVelField);
            c.oldVelField = c.velField;
            m_MACGrid.insertCellInHashTable(c);
            //utility::printvec(c.oldVelField);
        }

        t += subStep;
    }
}


real FlipSim::cfl()
{
    real umax = m_MACGrid.getMaxSpeed() + (std::sqrt(5*m_dx*(-1*m_g.y)));

    real ret = (5*m_dx)/std::max(umax, 0.000001);
    return ret * m_k_cfl;
}

void FlipSim::addBodyForce(real _dt)
{
    for(uint k = 0; k < m_kSize; k++)
    {
        for(uint j = 0; j < m_jSize; j++)
        {
            for(uint i = 0; i < m_iSize; i++)
            {
                //utility::printvec(m_MACGrid.getCell(i,j,k).velField);
                MG_Cell c = m_MACGrid.getCell(i,j,k);
                c.velField += _dt * m_g;
                m_MACGrid.insertCellInHashTable(c);
                //utility::printvec(m_MACGrid.getCell(i,j,k).velField);
            }
        }
    }
}

void FlipSim::updateGrid()
{
    for(uint k = 0; k < m_kSize; k++)
    {
        for(uint j = 0; j < m_jSize; j++)
        {
            for(uint i = 0; i < m_iSize; i++)
            {
                MG_Cell c;
                if(m_MACGrid.checkForCell(i, j, k, c))
                {
                    //utility::printvec(c.velField);
                    std::vector<MG_Cell> neighbors = m_MACGrid.getNeighbors(c);
                    std::vector<MG_Particle> tmpP;
                    vec3 vel;
                    for(MG_Particle p : m_MACGrid.m_particles)
                    {
                        if(c.key == p.cellidx)
                            tmpP.push_back(p);
                    }
                    if(tmpP.size() == 0 && c.type != SOLID)
                    {
                        c.type = AIR;
                    }
                    else if(!tmpP.size() == 0)
                    {
                        c.type = FLUID;

                        for(MG_Particle p : tmpP)
                        {
                            vec3 cellpos = m_MACGrid.getCellPos(c);
                            real newU = utility::lerp(neighbors[LEFT].u(), c.u(), p.pos.x-cellpos.x);
                            real newV = utility::lerp(neighbors[DOWN].v(), c.v(), p.pos.y-cellpos.y);
                            real newW = utility::lerp(neighbors[BACKWARD].w(), c.w(), p.pos.z-cellpos.z);
                            vel += vec3(newU, newV, newW);
                        }
                        vel /= tmpP.size();
                        c.velField = vel;
                        m_MACGrid.insertCellInHashTable(c);
                    }
                    //utility::printvec(c.velField);
                }
            }
        }
    }
}

void FlipSim::project(real _dt)
{
    calculateNegativeDivergence();
    calculatePressure(_dt);
    applyPressure(_dt);
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
                        std::vector<MG_Cell> neighbors = m_MACGrid.getNeighbors(c);

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
    std::vector<int> fluidCellKeys;
    uint fluidIDX = 0;
    for(uint c : m_MACGrid.m_cells)
        if(m_MACGrid.getCell(c).type == FLUID)
        {
            n_fluidCells++;
            fluidCellKeys.push_back(c);
//            MG_Cell tmp = m_MACGrid.getCell(c);
//            tmp.fluidIDX = fluidIDX++;
//            m_MACGrid.insertCellInHashTable(tmp);
            m_MACGrid.getCell(c).updateFluidIDX(fluidIDX++);
        }

    //Coeff Matrix
    Eigen::SparseMatrix<real> A(n_fluidCells, n_fluidCells);
    A.reserve(Eigen::VectorXi::Constant(n_fluidCells, 7));

    //Fill with rhs values
    //Divergence
    VectorX b(n_fluidCells);
    for(uint i = 0; i < n_fluidCells; i++)
        b[i] = m_MACGrid.getCell(fluidCellKeys[i]).rhs;

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
                    aDiagIdx = c.fluidIDX;
                    neighbors = m_MACGrid.getNeighbors(c);
                }
                else
                    continue;

                aXidx = neighbors[RIGHT].fluidIDX;

                aYidx = neighbors[UP].fluidIDX;

                aZidx = neighbors[FORWARD].fluidIDX;

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
    {
        MG_Cell c = m_MACGrid.getCell(fluidCellKeys[i]);
        //std::cout<<c.p<<"\n";
        c.p = p[i];
        m_MACGrid.insertCellInHashTable(c);
        //std::cout<<c.p<<"\n";
    }
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
                {
                    neighbors = m_MACGrid.getNeighbors(c);
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
                    //utility::printvec(c.velField);
                    c.velField = vec3(u,v,w);
                    m_MACGrid.insertCellInHashTable(c);
                    if(c.type == FLUID)
                        real t = c.p;
                    //utility::printvec(c.velField);
                }
            }
        }
    }
}

void FlipSim::updateParticles()
{
    std::vector<MG_Particle> newP;
    for(MG_Particle p : m_MACGrid.m_particles)
    {
        p.vel = m_MACGrid.getCell(p.cellidx).oldVelField - m_MACGrid.getCell(p.cellidx).velField;
        newP.push_back(p);
//        utility::printvec(m_MACGrid.getCell(p.cellidx).velField);
    }
    m_MACGrid.m_particles = newP;

}

void FlipSim::advectVelocityField(real _dt)
{
//    for(MG_Particle p : m_MACGrid.m_particles)
//    {
//        p.pos = m_MACGrid.tracePoint(p.pos, _dt);
//    }
    for(uint i = 0; i < m_MACGrid.m_particles.size(); i++)
    {
//        utility::printvec(m_MACGrid.m_particles[i].pos);
        m_MACGrid.m_particles[i].pos = m_MACGrid.tracePoint(m_MACGrid.m_particles[i].pos, _dt);
//        utility::printvec(m_MACGrid.m_particles[i].pos);
    }
}

std::vector<MG_Particle> FlipSim::getParticles()
{
    return m_MACGrid.m_particles;
}
