/******************************************************************************
 *
 *  Headers.h
 *
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  7 / 14 / 2011
 *     
 *    This headear module contains all of the library headers, variable type 
 *  definitions, and declaration of global variables, file input/output streams, 
 *  and function prototypes for LoBoomFDTD2D_2006 program.  This header must be 
 *  included in the main function and every member modules for correct 
 *  compilation of the program.
 *    
 ******************************************************************************/  
#include <iostream>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <vector>
#include <string>

using namespace std;

#define PI 3.14159265

// Definition of Custom Variable Types
struct Mic {
  int i,j;
};

// Variables
extern string input_name;
extern int input_length_i, num_floors, x_indmax, z_indmax, t_indmax;
extern double input_fs, x_max, z_max, x_intersect, x0, dx, tmax, c0, theta_inc, dt, cXcXdtXdt, dxXdx, theta_ref;
extern int num_frame, mov_skip, x_range_ind_lo, x_range_ind_hi, z_range_ind_hi;
extern double x_range_lo, x_range_hi, z_range_hi;
extern double sin_inc, sin_ref, cos_inc, cos_ref, tan14, z1, z2, z3, z4;
extern int prop_time_inc, prop_time_ref, input_count;
extern int num_mics;
extern vector<Mic> mics;

// File I/O stream
extern ofstream p_2D_initial, p_domain, p_mic, p_2D_final;

// Function Prototypes
extern void ReadParameters();
extern void InitializeDataFiles();