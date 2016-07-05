/*Raul P. Pelaez 2016. Two step velocity VerletNVT Integrator derived class implementation

  An Integrator is intended to be a separated module that handles the update of positions given the forces

  It takes care of keeping the positions updated.
  The positions must be provided, they are not created by the module.
  Also takes care of writing to disk
 
  Currently uses a BBK thermostat to maintain the temperature.
  Solves the following differential equation using a two step velocity verlet algorithm, see GPU code:
      X[t+dt] = X[t] + v[t]·dt
      v[t+dt]/dt = -gamma·v[t] - F(X) + sigma·G
  gamma is a damping constant, sigma = sqrt(2·gamma·T) and G are normal distributed random numbers with var 1.0 and mean 0.0.

TODO:
100- Allow velocities from outside
90-  Implement more thermostats
80-  Gamma should be chosen by the user
70- Improve the intial conditions for the velocity, should be normal, not uniform.
*/
#include "VerletNVT.h"


VerletNVT::VerletNVT():
  Integrator(), /*After initializing the base class, you have access to things like N, L...*/
  noise(N +((3*N)%2)), vel(N)
{
  cerr<<"Initializing Verlet NVT Integrator..."<<endl;

  cerr<<"\tSet T="<<gcnf.T<<endl;
  
  gamma = 0.1f;

  /*Set params and init GPU parameters*/
  params.gamma = gamma;
  params.dt = dt;
  params.T = gcnf.T;
  params.noiseAmp = sqrt(dt*0.5f)*sqrt(2.0f*gamma*params.T);
  initVerletNVTGPU(params);

  
  /*Init rng*/
  curandCreateGenerator(&rng, CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(rng, 1234ULL);  //WARNING, THE SEED!!!
  /*This shit is obscure, curand will only work with an even number of elements*/
  curandGenerateNormal(rng, (float*) noise.d_m, 3*N + ((3*N)%2), 0.0f, 1.0f);
  
  /*Distribute the velocities according to the temperature*/
  float vamp = sqrt(3.0f*params.T);
  /*Create velocities*/
  vel.fill_with(make_float3(0.0f));
  fori(0,N){
    vel[i].x = vamp*(2.0f*(rand()/(float)RAND_MAX)-1.0f);
    vel[i].y = vamp*(2.0f*(rand()/(float)RAND_MAX)-1.0f);
    vel[i].z = vamp*(2.0f*(rand()/(float)RAND_MAX)-1.0f);
  }
  vel.upload();
  
  cerr<<"Verlet NVT Integrator\t\tDONE!!\n"<<endl;
}

VerletNVT::~VerletNVT(){}

//The integration process can have two steps
void VerletNVT::update(){
  if(steps==0)
    for(auto forceComp: interactors) forceComp->sumForce();
  
  steps++;
  if(steps%1000==0) cerr<<"\rComputing step: "<<steps<<"   ";

  /**First integration step**/
  /*Gen noise*/
  curandGenerateNormal(rng, (float*) noise.d_m, 3*N + ((3*N)%2), 0.0f, 1.0f);
  integrateVerletNVTGPU(pos->d_m, vel, force->d_m, noise, N, 1);
  /**Compute all the forces**/
  
  /*Reset forces*/
  cudaMemset((float *)force->d_m, 0.0f, 4*N*sizeof(float));
  for(auto forceComp: interactors) forceComp->sumForce();
  
  /**Second integration step**/
  /*Gen noise*/
  curandGenerateNormal(rng, (float*) noise.d_m, 3*N + ((3*N)%2), 0.0f, 1.0f);
  integrateVerletNVTGPU(pos->d_m, vel, force->d_m, noise, N, 2);
}


float VerletNVT::sumEnergy(){
  return computeKineticEnergyVerletNVT(vel, N);
}
