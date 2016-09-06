/*Raul P. Pelaez 2016. Brownian Euler Maruyama with hydrodynamics integrator GPU kernels and callers

  Functions to integrate movement. The integration is done via a functor wich accesor ()
     takes a thrust::Tuple containing positions, velocities and forces on each particle. 
*/

#ifndef INTEGRATORBROWNIANHYDRODYNAMICSEULERMARUYAMAGPU_CUH
#define INTEGRATORBROWNIANHYDRODYNAMICSEULERMARUYAMAGPU_CUH

#include<thrust/device_vector.h>
struct BrownianHydrodynamicsEulerMaruyamaParameters{
  float sqrtdt;
  float dt;
  float sqrt2dt;
  float4* K;
  float D0, rh;
};


void copy_pos(float4 *force, float3* force3, uint N);

void rodne_call(float *d_D, float4 *d_R, cudaStream_t stream, uint N);

//void fix_B_call(float *B, uint N, cudaStream_t stream);


void initBrownianHydrodynamicsEulerMaruyamaGPU(BrownianHydrodynamicsEulerMaruyamaParameters m_params);

void integrateBrownianHydrodynamicsEulerMaruyamaGPU(float4 *pos, float3* DF, float3* BdW, uint N);

#endif
