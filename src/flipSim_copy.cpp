/*
#include "flipSim.h"

FlipSim::FlipSim(uint _size, real _cellWidth, uvec3 _b, uvec3 _c):m_MACGrid(_size, _cellWidth)
{
    m_dx = _cellWidth;

    m_iSize = _size;
    m_jSize = _size;
    m_kSize = _size;

    m_MACGrid.initialiseCells(_b, _c);
}

void FlipSim::step(real _dt)
{
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
        advectVelocityField(subStep);

        for(uint i = 0; i < m_MACGrid.m_cells.size(); i++)
        {
            MG_Cell c = m_MACGrid.getCell(m_MACGrid.m_cells[i]);
            //utility::printvec(c.oldVelField);
            c.oldVelField = c.velField;
            m_MACGrid.insertCellInHashTable(c);
            //utility::printvec(c.oldVelField);
        }

        t += subStep;
        if(t>_dt)
            t = _dt;
    }
}


real FlipSim::cfl()
{
    real umax = m_MACGrid.getMaxSpeed() + (std::sqrt(5*m_dx*(-1*m_g.y)));

    real ret = (5*m_dx)/std::max(umax, 1.0);
    ret = 1.0;
    return ret * m_k_cfl;
}

void FlipSim::classifyParticles()
{

}

void FlipSim::updateGrid()
{
    for(uint k = 0; k < m_kSize; k++)
    {
        for(uint j = 0; j < m_jSize; j++)
        {
            for(uint i = 0; i < m_iSize; i++)
            {
                MG_Cell c = m_MACGrid.getCell(i,j,k);
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
                else if(tmpP.size() != 0)
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

void FlipSim::addBodyForce(real _dt)
{
    for(uint k = 0; k < m_kSize; k++)
    {
        for(uint j = 0; j < m_jSize; j++)
        {
            for(uint i = 0; i < m_iSize; i++)
            {
                //                utility::printvec(m_MACGrid.getCell(i,j,k).velField);
                MG_Cell c = m_MACGrid.getCell(i,j,k);
                c.velField += _dt * m_g;
                m_MACGrid.insertCellInHashTable(c);
                //                utility::printvec(m_MACGrid.getCell(i,j,k).velField);
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
                MG_Cell c = m_MACGrid.getCell(uvec3(i,j,k));
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
    uint n_cells = m_iSize, m_jSize, m_kSize;

    //Coeff Matrix
    Eigen::SparseMatrix<real> A(n_cells, n_cells);
    A.reserve(Eigen::VectorXi::Constant(n_cells, 7));

    //Fill with rhs values
    //Divergence
    VectorX b(n_cells);
    for(uint i = 0; i < n_cells; i++)
        b[i] = m_MACGrid.getCell(m_MACGrid.m_cells[i]).rhs;

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
                    aDiagIdx = utility::getIndex(m_iSize, c);
                    neighbors = m_MACGrid.getNeighbors(c);
                }
                else
                    continue;

                aXidx = utility::getIndex(m_iSize, neighbors[RIGHT]);

                aYidx = utility::getIndex(m_iSize, neighbors[UP]);

                aZidx = utility::getIndex(m_iSize, neighbors[FORWARD]);

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
    VectorX p(n_cells);
    solver.compute(A);
    p = solver.solve(b);

    for(uint i = 0; i < n_cells; i++)
    {
        MG_Cell c = m_MACGrid.getCell(m_MACGrid.m_cells[i]);
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

                MG_Cell c = m_MACGrid.getCell(i,j,k);
                std::vector<MG_Cell> neighbors;
                neighbors = m_MACGrid.getNeighbors(c);
                //Update u
                if(neighbors[LEFT].type == FLUID || c.type == FLUID)
                {
                    if(!(neighbors[LEFT].type == SOLID || c.type == SOLID))
                    {
                        u = 0.0;
                    }
                    else
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
                        v = 0.0;
                    }
                    else
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
                        w = 0.0;
                    }
                    else
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
                //utility::printvec(c.velField);
            }
        }
    }
}

void FlipSim::updateParticles()
{
    for(MG_Particle p : m_MACGrid.m_particles)
    {
        p.updateVel(p.vel + m_MACGrid.getCell(p.cellidx).oldVelField - m_MACGrid.getCell(p.cellidx).velField);
        //        p.updateVel(m_MACGrid.getCell(p.cellidx).velField);
        utility::printvec(m_MACGrid.getCell(p.cellidx).velField);
    }

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
        m_MACGrid.m_particles[i].pos = m_MACGrid.tracePoint(m_MACGrid.m_particles[i].vel, _dt);
        //        utility::printvec(m_MACGrid.m_particles[i].pos);
    }
}

std::vector<MG_Particle> FlipSim::getParticles()
{
    return m_MACGrid.m_particles;
}
*/
