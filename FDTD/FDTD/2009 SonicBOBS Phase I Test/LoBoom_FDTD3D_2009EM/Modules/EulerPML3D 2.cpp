/*******************************************************************************
 *
 *  EulerPML3D.cpp
 *
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  7 / 22 / 2011
 *  
 *    This module performs a second-order centered-difference finite-
 *  differencing on spatial derivatives of the field variables for 
 *  LoBoomFDTD3D_SingleBldg program.  The grid points inside the PML and points 
 *  adjacent to rigid boundaries at the ground and the building structure are 
 *  separately treated.  The rigid boundaries are modeled to be at the mid-point
 *  between two points.  Spatial finite-differencing operated on the field 
 *  prepared with this module will produce reflections at these boundaries.
 *    
 *******************************************************************************/     

#include "Headers.h"

void EulerPML3D(Grid3D &U, Grid3D &dU, Array &sigma)
{
  int i,j,k;
  // Precalculation of equation constants for reducing computation time
  float CFLz0o2 = -CFL/2*z0;
  float CFLo2z0 = -CFL/2/z0;
  float CFLo2 = -CFL/2;
  
  #pragma omp for
  for (i=1; i<I_MAX+1; i++) {
    for (j=1; j<J_MAX+1; j++) {
      for (k=1; k<K_MAX+1; k++) {
        // Perform FDTD calculation only on grid points exterior to the building
        if (i<wall_loc_i1 || i>wall_loc_i2 || j<wall_loc_j1 || j>wall_loc_j2 || k<wall_loc_k1 || k>wall_loc_k2) {
          // For x-direction
          if (i==1) {           // PML region @ domain boundary
            dU[1][j][k].px = CFLo2*(z0*(U[2][j][k].u+U[1][j][k].u)+sigma[0]*U[1][j][k].px);
            dU[1][j][k].u = CFLo2*((U[2][j][k].px-U[1][j][k].px+U[2][j][k].py-U[1][j][k].py+U[2][j][k].pz-U[1][j][k].pz)/z0+sigma[0]*U[1][j][k].u);
          }
          else if (i<=D) {      // PML region @ x < 0
            dU[i][j][k].px = CFLo2*(z0*(U[i+1][j][k].u-U[i-1][j][k].u)+sigma[i-1]*U[i][j][k].px);
            dU[i][j][k].u = CFLo2*((U[i+1][j][k].px-U[i-1][j][k].px+U[i+1][j][k].py-U[i-1][j][k].py+U[i+1][j][k].pz-U[i-1][j][k].pz)/z0+sigma[i-1]*U[i][j][k].u);
          }
          else if (i==I_MAX) {  // PML region @ domain boundary
            dU[I_MAX][j][k].px = CFLo2*(z0*(-U[I_MAX][j][k].u-U[I_MAX-1][j][k].u)+sigma[0]*U[I_MAX][j][k].px);
            dU[I_MAX][j][k].u = CFLo2*((U[I_MAX][j][k].px-U[I_MAX-1][j][k].px+U[I_MAX][j][k].py-U[I_MAX-1][j][k].py+U[I_MAX][j][k].pz-U[I_MAX-1][j][k].pz)/z0+sigma[0]*U[I_MAX][j][k].u);
          }          
          else if (i>I_MAX-D) { // PML region @ x > x_max
            dU[i][j][k].px = CFLo2*(z0*(U[i+1][j][k].u-U[i-1][j][k].u)+sigma[I_MAX-i]*U[i][j][k].px);
            dU[i][j][k].u = CFLo2*((U[i+1][j][k].px-U[i-1][j][k].px+U[i+1][j][k].py-U[i-1][j][k].py+U[i+1][j][k].pz-U[i-1][j][k].pz)/z0+sigma[I_MAX-i]*U[i][j][k].u);
          }
          else {                // Euler region in x
            // West wall
            if (i==wall_loc_i1-1 && j>=wall_loc_j1 && j<=wall_loc_j2 && k>=wall_loc_k1 && k<=wall_loc_k2) {
              dU[i][j][k].px = CFLz0o2*(-U[i][j][k].u-U[i-1][j][k].u);
              dU[i][j][k].u = CFLo2z0*(U[i][j][k].px-U[i-1][j][k].px+U[i][j][k].py-U[i-1][j][k].py+U[i][j][k].pz-U[i-1][j][k].pz);
            }
            // East wall
            else if (i==wall_loc_i2+1 && j>=wall_loc_j1 && j<=wall_loc_j2 && k>=wall_loc_k1 && k<=wall_loc_k2) {
                dU[i][j][k].px = CFLz0o2*(U[i+1][j][k].u+U[i][j][k].u);
                dU[i][j][k].u = CFLo2z0*(U[i+1][j][k].px-U[i][j][k].px+U[i+1][j][k].py-U[i][j][k].py+U[i+1][j][k].pz-U[i][j][k].pz);
            }
            // Euler region
            else {
              dU[i][j][k].px = CFLz0o2*(U[i+1][j][k].u-U[i-1][j][k].u);
              dU[i][j][k].u = CFLo2z0*(U[i+1][j][k].px-U[i-1][j][k].px+U[i+1][j][k].py-U[i-1][j][k].py+U[i+1][j][k].pz-U[i-1][j][k].pz);
            }
          }
          
          // For y-direction
          if (j==1) {           // PML region @ domain boundary
            dU[i][1][k].py = CFLo2*(z0*(U[i][2][k].v+U[i][1][k].v)+sigma[0]*U[i][1][k].py);
            dU[i][1][k].v = CFLo2*((U[i][2][k].px-U[i][1][k].px+U[i][2][k].py-U[i][1][k].py+U[i][2][k].pz-U[i][1][k].pz)/z0+sigma[0]*U[i][1][k].v);
          }
          else if (j<=D) {      // PML region @ y < 0
            dU[i][j][k].py = CFLo2*(z0*(U[i][j+1][k].v-U[i][j-1][k].v)+sigma[j-1]*U[i][j][k].py);
            dU[i][j][k].v = CFLo2*((U[i][j+1][k].px-U[i][j-1][k].px+U[i][j+1][k].py-U[i][j-1][k].py+U[i][j+1][k].pz-U[i][j-1][k].pz)/z0+sigma[j-1]*U[i][j][k].v);
          }
          else if (j==J_MAX) {  // PML region @ domain boundary
            dU[i][J_MAX][k].py = CFLo2*(z0*(-U[i][J_MAX][k].v-U[i][J_MAX-1][k].v)+sigma[0]*U[i][J_MAX][k].py);
            dU[i][J_MAX][k].v = CFLo2*((U[i][J_MAX][k].px-U[i][J_MAX-1][k].px+U[i][J_MAX][k].py-U[i][J_MAX-1][k].py+U[i][J_MAX][k].pz-U[i][J_MAX-1][k].pz)/z0+sigma[0]*U[i][J_MAX][k].v);
          }
          else if (j>J_MAX-D) { // PML region @ y > y_max
            dU[i][j][k].py = CFLo2*(z0*(U[i][j+1][k].v-U[i][j-1][k].v)+sigma[J_MAX-j]*U[i][j][k].py);
            dU[i][j][k].v = CFLo2*((U[i][j+1][k].px-U[i][j-1][k].px+U[i][j+1][k].py-U[i][j-1][k].py+U[i][j+1][k].pz-U[i][j-1][k].pz)/z0+sigma[J_MAX-j]*U[i][j][k].v);
          }
          else {                // Euler region in y
            // South wall
            if (j==wall_loc_j1-1 && i>=wall_loc_i1 && i<=wall_loc_i2 && k>=wall_loc_k1 && k<=wall_loc_k2) {
              dU[i][j][k].py = CFLz0o2*(-U[i][j][k].v-U[i][j-1][k].v);
              dU[i][j][k].v = CFLo2z0*(U[i][j][k].px-U[i][j-1][k].px+U[i][j][k].py-U[i][j-1][k].py+U[i][j][k].pz-U[i][j-1][k].pz);
            }
            // North wall
            else if (j==wall_loc_j2+1 && i>=wall_loc_i1 && i<=wall_loc_i2 && k>=wall_loc_k1 && k<=wall_loc_k2) {
              dU[i][j][k].py = CFLz0o2*(U[i][j+1][k].v+U[i][j][k].v);
              dU[i][j][k].v = CFLo2z0*(U[i][j+1][k].px-U[i][j][k].px+U[i][j+1][k].py-U[i][j][k].py+U[i][j+1][k].pz-U[i][j][k].pz);
            }
            // Euler region
            else {
              dU[i][j][k].py = CFLz0o2*(U[i][j+1][k].v-U[i][j-1][k].v);
              dU[i][j][k].v = CFLo2z0*(U[i][j+1][k].px-U[i][j-1][k].px+U[i][j+1][k].py-U[i][j-1][k].py+U[i][j+1][k].pz-U[i][j-1][k].pz);
            }
          }

          // For z-direction
          if (k==1) {           // Euler region @ domain boundary
            dU[i][j][1].pz = CFLz0o2*(U[i][j][2].w+U[i][j][1].w);
            dU[i][j][1].w = CFLo2z0*(U[i][j][2].px-U[i][j][1].px+U[i][j][2].py-U[i][j][1].py+U[i][j][2].pz-U[i][j][1].pz);
          }
          else if (k==K_MAX) {  // PML region @ domain boundary
            dU[i][j][K_MAX].pz = CFLo2*(z0*(-U[i][j][K_MAX].w-U[i][j][K_MAX-1].w)+sigma[0]*U[i][j][K_MAX].pz);
            dU[i][j][K_MAX].w = CFLo2*((U[i][j][K_MAX].px-U[i][j][K_MAX-1].px+U[i][j][K_MAX].py-U[i][j][K_MAX-1].py+U[i][j][K_MAX].pz-U[i][j][K_MAX-1].pz)/z0+sigma[0]*U[i][j][K_MAX].w);
          }
          else if (k>K_MAX-D) { // PML region @ z > z_max
            dU[i][j][k].pz = CFLo2*(z0*(U[i][j][k+1].w-U[i][j][k-1].w)+sigma[K_MAX-k]*U[i][j][k].pz);
            dU[i][j][k].w = CFLo2*((U[i][j][k+1].px-U[i][j][k-1].px+U[i][j][k+1].py-U[i][j][k-1].py+U[i][j][k+1].pz-U[i][j][k-1].pz)/z0+sigma[K_MAX-k]*U[i][j][k].w);
          }
          else {                // Euler region in z
            // Roof
            if (k==wall_loc_k2+1 && i>=wall_loc_i1 && i<=wall_loc_i2 && j>=wall_loc_j1 && j<=wall_loc_j2) {
              dU[i][j][k].pz = CFLz0o2*(U[i][j][k+1].w+U[i][j][k].w);
              dU[i][j][k].w = CFLo2z0*(U[i][j][k+1].px-U[i][j][k].px+U[i][j][k+1].py-U[i][j][k].py+U[i][j][k+1].pz-U[i][j][k].pz);
            }
            // Euler region
            else {
              dU[i][j][k].pz = CFLz0o2*(U[i][j][k+1].w-U[i][j][k-1].w);
              dU[i][j][k].w = CFLo2z0*(U[i][j][k+1].px-U[i][j][k-1].px+U[i][j][k+1].py-U[i][j][k-1].py+U[i][j][k+1].pz-U[i][j][k-1].pz);
            }
          }
        }
      }
    }
  }
  
  return;
}
