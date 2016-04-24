/*************************************************************************
    > File Name: simulator.hpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Thu 23 Jul 2015 11:59:50 PM EDT
 ************************************************************************/

#ifndef SIMULATOR
#define SIMULATOR

#include<functional>

#include<Eigen/Dense>
#include<vector>

#include<mat.h>

// System: system paramter
// State: state of the system
// N: dimension of control inputs

template<typename Robot> class Sim
{
	public:
		typedef typename Robot::System System;

		typedef typename Robot::Ref Ref;

		typedef typename Robot::State State;
		typedef typename Robot::DState DState;
		
		typedef typename Robot::V V;
		typedef typename Robot::U U;

	public:
		const System system;
		const double dt; // simulation time step, dt>0

	protected:
		std::vector<State> states[4];
		std::vector<U> inputs;

	public:
		Sim(System system_, double dt_):system(system_), dt(dt_)
		{
		}

		State get_state() const
		{
			return states[0].back();
		}

		std::vector<State> const & get_states() const
		{
			return states[0];
		}

		std::vector<U> const & get_inputs() const
		{
			return inputs;
		}

		void init(State const & state0, size_t const & n)
		{
			states[0].clear();
			states[1].clear();
			states[2].clear();
			states[3].clear();

			states[0].reserve(n);
			states[1].reserve(n);
			states[2].reserve(n);
			states[3].reserve(n);

			inputs.clear();
			inputs.reserve(n);
			states[0].push_back(state0);
		}

		State update(U u) // dt>0, simulate states forwards by RK4
		{
			State state[4];

			state[0]=states[0].back();
			DState k1(system,state[0],u);

			state[1]=state[0].update(k1,dt/2);
			DState k2(system,state[1],u);

			state[2]=state[0].update(k2,dt/2);
			DState k3(system,state[2],u);

			state[3]=state[0].update(k3,dt);
			DState k4(system,state[3],u);

			DState k(k1,k2,k3,k4);
			State state_curr=state[0].update(k,dt);

			states[0].push_back(state_curr);
			states[1].push_back(state[1]);
			states[2].push_back(state[2]);
			states[3].push_back(state[3]);

			inputs.push_back(u);

			return state_curr;
		}

		void clear()
		{
			states[0].clear();
			states[1].clear();
			states[2].clear();
			states[3].clear();

			inputs.clear();
		}

		void save(std::string path) const
		{
			State::save(states[0],path);
			
			std::cout<<path<<std::endl;
			MATFile *file;
			void *p;

			mwSize *dims=new mwSize[3];
			dims[0]=Robot::N;
			dims[1]=1;
			dims[2]=inputs.size();

			file=matOpen(path.c_str(),"u");

			typename std::vector<Eigen::Matrix<double,Robot::N,1> >::const_iterator itr_u=inputs.begin();

			mxArray *pU=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);

			for(p=mxGetPr(pU);itr_u!=inputs.end();itr_u++)
			{
				memcpy(p,itr_u->data(),sizeof(double)*Robot::N);
				p+=sizeof(double)*Robot::N;
			}

			matPutVariable(file,"u",pU);
		}
};

#endif
