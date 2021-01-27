/*******************************************************************************
 *
 *  PrepareDomain_2009_Museum.cpp
 *  
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  9 / 28 / 2011
 * 
 *    This program generates the simulation parameter and geometry input files 
 *  necessary for LoBoomFDTD3D_MultiBldg program simulating the museum buildings
 *  from 2009 SonicBOBS experiment.  The computational domain grid including the
 *  building structures is  defined, where the location coordinates of the 
 *  building corners and mics, calculated with the GPS survey information stored
 *  in In_coordinates.txt, are used to identify the interior and exterior points
 *  and also the boundary points.
 *           
 *******************************************************************************/

#include "Headers.h"

// Definition of Custom Variable Types
struct Point
{
  double x,y,z;
};
struct Plane
{
  double a,b,c,d,D;
};
typedef vector<vector<double> > Matrix2D;

// Function Prototypes
double Determinant(double Mat[3][3]);

void PrepareDomain(Grid3D &U)
{
  int i, j, k;
  
  /*****************************************************************************
   *  Coordinate Rotating & Shifting
   *****************************************************************************/
  cout << "Rotating & Shifting coordinates... ";  
  
  // Original coordinates of buildings and microphones
  int tot_num_pt = 28;    // Total number of coordinate points stored in In_coordinates_polar.txt
  ifstream f_coord("In_coordinates_polar.txt");
  if (!f_coord) {
    cout << "ERROR: Could not open In_coordinates_polar.txt" <<endl;
    exit(1);
  }
  
  Matrix2D Coord_shift(tot_num_pt,vector<double>(2));
  vector<double> Distance(tot_num_pt), Azimuth(tot_num_pt), Azimuth_rot(tot_num_pt,0);
  vector<double> Coord_x(tot_num_pt), Coord_y(tot_num_pt);  

  for (i=0; i<tot_num_pt; i++) {
    f_coord >> Distance[i];
    f_coord >> Azimuth[i];
  }

  // Rotate to align the coordinate system with the incoming boom azimuth
  for (i=1; i<tot_num_pt; i++)  Azimuth_rot[i] = Azimuth[i]+(270-angle_azm);

  for (i=0; i<tot_num_pt; i++) {
    Coord_x[i] = Distance[i]*sin(Azimuth_rot[i]*PI/180);
    Coord_y[i] = Distance[i]*cos(Azimuth_rot[i]*PI/180);
  }

  // Shift to accommodate extra_dim_x,y around the building geometry
  double x_shift = -*min_element(Coord_x.begin(),Coord_x.begin()+8)+extra_dim_x;
  double y_shift = -*min_element(Coord_y.begin(),Coord_y.begin()+8)+extra_dim_y;
  
  for (i=0; i<tot_num_pt; i++) {
    Coord_shift[i][0] = Coord_x[i]+x_shift;  // The smallest x-coord of the building corners is set at extra_dim_x
    Coord_shift[i][1] = Coord_y[i]+y_shift;  // The smallest y-coord of the building corners is set at extra_dim_y
  }
  
  Matrix2D Bldg1(Coord_shift.begin(),Coord_shift.begin()+4);      // Coordinates for building #1
  Matrix2D Bldg2(Coord_shift.begin()+4,Coord_shift.begin()+8);    // Coordinates for building #2
  Matrix2D CMics(Coord_shift.begin()+8,Coord_shift.begin()+18);   // Coordinates for C mics
  Matrix2D WMics(Coord_shift.begin()+18,Coord_shift.begin()+28);  // Coordinates for W mics
  
  cout << "DONE" << endl;
  
  
  /*****************************************************************************
   *  Defining Planes
   *****************************************************************************/  
  cout << "Defining planes... ";
    
  Point pt[20];
  Plane P[12];
  // Define the location of the building corner points
  pt[0].x = Bldg1[2][0];
  pt[0].y = Bldg1[2][1];
  pt[0].z = 0;
  pt[1].x = Bldg1[3][0];
  pt[1].y = Bldg1[3][1];
  pt[1].z = 0;
  pt[2].x = Bldg1[1][0];
  pt[2].y = Bldg1[1][1];
  pt[2].z = 0;
  pt[3].x = Bldg1[0][0];
  pt[3].y = Bldg1[0][1];
  pt[3].z = 0;
  pt[4].x = Bldg1[2][0];
  pt[4].y = Bldg1[2][1];
  pt[4].z = 6.54;
  pt[5].x = Bldg1[3][0];
  pt[5].y = Bldg1[3][1];
  pt[5].z = 6.54;
  pt[6].x = (Bldg1[1][0]+Bldg1[3][0])/2;
  pt[6].y = (Bldg1[1][1]+Bldg1[3][1])/2;
  pt[6].z = 7.479;
  pt[7].x = Bldg1[1][0];
  pt[7].y = Bldg1[1][1];
  pt[7].z = 6.54;
  pt[8].x = Bldg1[0][0];
  pt[8].y = Bldg1[0][1];
  pt[8].z = 6.54;
  pt[9].x = (Bldg1[0][0]+Bldg1[2][0])/2;
  pt[9].y = (Bldg1[0][1]+Bldg1[2][1])/2;
  pt[9].z = 7.479;
  pt[10].x = Bldg2[2][0];
  pt[10].y = Bldg2[2][1];
  pt[10].z = 0;
  pt[11].x = Bldg2[3][0];
  pt[11].y = Bldg2[3][1];
  pt[11].z = 0;
  pt[12].x = Bldg2[1][0];
  pt[12].y = Bldg2[1][1];
  pt[12].z = 0;
  pt[13].x = Bldg2[0][0];
  pt[13].y = Bldg2[0][1];
  pt[13].z = 0;
  pt[14].x = Bldg2[2][0];
  pt[14].y = Bldg2[2][1];
  pt[14].z = 4.519;
  pt[15].x = Bldg2[3][0];
  pt[15].y = Bldg2[3][1];
  pt[15].z = 4.519;
  pt[16].x = (Bldg2[3][0]+Bldg2[1][0])/2;
  pt[16].y = (Bldg2[3][1]+Bldg2[1][1])/2;
  pt[16].z = 6.66;
  pt[17].x = Bldg2[1][0];
  pt[17].y = Bldg2[1][1];
  pt[17].z = 4.519;
  pt[18].x = Bldg2[0][0];
  pt[18].y = Bldg2[0][0];
  pt[18].z = 4.519;
  pt[19].x = (Bldg2[0][0]+Bldg2[2][0])/2;
  pt[19].y = (Bldg2[0][1]+Bldg2[2][1])/2;
  pt[19].z = 6.66;
  
  // Define the plane equation coefficient using the location of three corner points that lie on that respective plane
  double Mat0[3][3] = {{pt[0].x, pt[0].y, pt[0].z}, {pt[3].x, pt[3].y, pt[3].z}, {pt[4].x, pt[4].y, pt[4].z}};
  double Mat0a[3][3] = {{1, pt[0].y, pt[0].z}, {1, pt[3].y, pt[3].z}, {1, pt[4].y, pt[4].z}};
  double Mat0b[3][3] = {{pt[0].x, 1, pt[0].z}, {pt[3].x, 1, pt[3].z}, {pt[4].x, 1, pt[4].z}};
  double Mat0c[3][3] = {{pt[0].x, pt[0].y, 1}, {pt[3].x, pt[3].y, 1}, {pt[4].x, pt[4].y, 1}};
  P[0].D = Determinant(Mat0);
  P[0].d = 1;
  P[0].a = -P[0].d/P[0].D*Determinant(Mat0a);
  P[0].b = -P[0].d/P[0].D*Determinant(Mat0b);
  P[0].c = -P[0].d/P[0].D*Determinant(Mat0c);
  
  double Mat1[3][3] = {{pt[0].x, pt[0].y, pt[0].z}, {pt[1].x, pt[1].y, pt[1].z}, {pt[4].x, pt[4].y, pt[4].z}};
  double Mat1a[3][3] = {{1, pt[0].y, pt[0].z}, {1, pt[1].y, pt[1].z}, {1, pt[4].y, pt[4].z}};
  double Mat1b[3][3] = {{pt[0].x, 1, pt[0].z}, {pt[1].x, 1, pt[1].z}, {pt[4].x, 1, pt[4].z}};
  double Mat1c[3][3] = {{pt[0].x, pt[0].y, 1}, {pt[1].x, pt[1].y, 1}, {pt[4].x, pt[4].y, 1}};
  P[1].D = Determinant(Mat1);
  P[1].d = 1;
  P[1].a = -P[1].d/P[1].D*Determinant(Mat1a);
  P[1].b = -P[1].d/P[1].D*Determinant(Mat1b);
  P[1].c = -P[1].d/P[1].D*Determinant(Mat1c);
  
  double Mat2[3][3] = {{pt[1].x, pt[1].y, pt[1].z}, {pt[2].x, pt[2].y, pt[2].z}, {pt[5].x, pt[5].y, pt[5].z}};
  double Mat2a[3][3] = {{1, pt[1].y, pt[1].z}, {1, pt[2].y, pt[2].z}, {1, pt[5].y, pt[5].z}};
  double Mat2b[3][3] = {{pt[1].x, 1, pt[1].z}, {pt[2].x, 1, pt[2].z}, {pt[5].x, 1, pt[5].z}};
  double Mat2c[3][3] = {{pt[1].x, pt[1].y, 1}, {pt[2].x, pt[2].y, 1}, {pt[5].x, pt[5].y, 1}};
  P[2].D = Determinant(Mat2);
  P[2].d = 1;
  P[2].a = -P[2].d/P[2].D*Determinant(Mat2a);
  P[2].b = -P[2].d/P[2].D*Determinant(Mat2b);
  P[2].c = -P[2].d/P[2].D*Determinant(Mat2c);
  
  double Mat3[3][3] = {{pt[7].x, pt[7].y, pt[7].z}, {pt[2].x, pt[2].y, pt[2].z}, {pt[3].x, pt[3].y, pt[3].z}};
  double Mat3a[3][3] = {{1, pt[7].y, pt[7].z}, {1, pt[2].y, pt[2].z}, {1, pt[3].y, pt[3].z}};
  double Mat3b[3][3] = {{pt[7].x, 1, pt[7].z}, {pt[2].x, 1, pt[2].z}, {pt[3].x, 1, pt[3].z}};
  double Mat3c[3][3] = {{pt[7].x, pt[7].y, 1}, {pt[2].x, pt[2].y, 1}, {pt[3].x, pt[3].y, 1}};
  P[3].D = Determinant(Mat3);
  P[3].d = 1;
  P[3].a = -P[3].d/P[3].D*Determinant(Mat3a);
  P[3].b = -P[3].d/P[3].D*Determinant(Mat3b);
  P[3].c = -P[3].d/P[3].D*Determinant(Mat3c);
  
  double Mat4[3][3] = {{pt[5].x, pt[5].y, pt[5].z}, {pt[6].x, pt[6].y, pt[6].z}, {pt[9].x, pt[9].y, pt[9].z}};
  double Mat4a[3][3] = {{1, pt[5].y, pt[5].z}, {1, pt[6].y, pt[6].z}, {1, pt[9].y, pt[9].z}};
  double Mat4b[3][3] = {{pt[5].x, 1, pt[5].z}, {pt[6].x, 1, pt[6].z}, {pt[9].x, 1, pt[9].z}};
  double Mat4c[3][3] = {{pt[5].x, pt[5].y, 1}, {pt[6].x, pt[6].y, 1}, {pt[9].x, pt[9].y, 1}};
  P[4].D = Determinant(Mat4);
  P[4].d = 1;
  P[4].a = -P[4].d/P[4].D*Determinant(Mat4a);
  P[4].b = -P[4].d/P[4].D*Determinant(Mat4b);
  P[4].c = -P[4].d/P[4].D*Determinant(Mat4c);
  
  double Mat5[3][3] = {{pt[7].x, pt[7].y, pt[7].z}, {pt[6].x, pt[6].y, pt[6].z}, {pt[9].x, pt[9].y, pt[9].z}};
  double Mat5a[3][3] = {{1, pt[7].y, pt[7].z}, {1, pt[6].y, pt[6].z}, {1, pt[9].y, pt[9].z}};
  double Mat5b[3][3] = {{pt[7].x, 1, pt[7].z}, {pt[6].x, 1, pt[6].z}, {pt[9].x, 1, pt[9].z}};
  double Mat5c[3][3] = {{pt[7].x, pt[7].y, 1}, {pt[6].x, pt[6].y, 1}, {pt[9].x, pt[9].y, 1}};
  P[5].D = Determinant(Mat5);
  P[5].d = 1;
  P[5].a = -P[5].d/P[5].D*Determinant(Mat5a);
  P[5].b = -P[5].d/P[5].D*Determinant(Mat5b);
  P[5].c = -P[5].d/P[5].D*Determinant(Mat5c);
  
  double Mat6[3][3] = {{pt[10].x, pt[10].y, pt[10].z}, {pt[14].x, pt[14].y, pt[14].z}, {pt[13].x, pt[13].y, pt[13].z}};
  double Mat6a[3][3] = {{1, pt[10].y, pt[10].z}, {1, pt[14].y, pt[14].z}, {1, pt[13].y, pt[13].z}};
  double Mat6b[3][3] = {{pt[10].x, 1, pt[10].z}, {pt[14].x, 1, pt[14].z}, {pt[13].x, 1, pt[13].z}};
  double Mat6c[3][3] = {{pt[10].x, pt[10].y, 1}, {pt[14].x, pt[14].y, 1}, {pt[13].x, pt[13].y, 1}};
  P[6].D = Determinant(Mat6);
  P[6].d = 1;
  P[6].a = -P[6].d/P[6].D*Determinant(Mat6a);
  P[6].b = -P[6].d/P[6].D*Determinant(Mat6b);
  P[6].c = -P[6].d/P[6].D*Determinant(Mat6c);
  
  double Mat7[3][3] = {{pt[10].x, pt[10].y, pt[10].z}, {pt[14].x, pt[14].y, pt[14].z}, {pt[11].x, pt[11].y, pt[11].z}};
  double Mat7a[3][3] = {{1, pt[10].y, pt[10].z}, {1, pt[14].y, pt[14].z}, {1, pt[11].y, pt[11].z}};
  double Mat7b[3][3] = {{pt[10].x, 1, pt[10].z}, {pt[14].x, 1, pt[14].z}, {pt[11].x, 1, pt[11].z}};
  double Mat7c[3][3] = {{pt[10].x, pt[10].y, 1}, {pt[14].x, pt[14].y, 1}, {pt[11].x, pt[11].y, 1}};
  P[7].D = Determinant(Mat7);
  P[7].d = 1;
  P[7].a = -P[7].d/P[7].D*Determinant(Mat7a);
  P[7].b = -P[7].d/P[7].D*Determinant(Mat7b);
  P[7].c = -P[7].d/P[7].D*Determinant(Mat7c);
  
  double Mat8[3][3] = {{pt[12].x, pt[12].y, pt[12].z}, {pt[17].x, pt[17].y, pt[17].z}, {pt[11].x, pt[11].y, pt[11].z}};
  double Mat8a[3][3] = {{1, pt[12].y, pt[12].z}, {1, pt[17].y, pt[17].z}, {1, pt[11].y, pt[11].z}};
  double Mat8b[3][3] = {{pt[12].x, 1, pt[12].z}, {pt[17].x, 1, pt[17].z}, {pt[11].x, 1, pt[11].z}};
  double Mat8c[3][3] = {{pt[12].x, pt[12].y, 1}, {pt[17].x, pt[17].y, 1}, {pt[11].x, pt[11].y, 1}};
  P[8].D = Determinant(Mat8);
  P[8].d = 1;
  P[8].a = -P[8].d/P[8].D*Determinant(Mat8a);
  P[8].b = -P[8].d/P[8].D*Determinant(Mat8b);
  P[8].c = -P[8].d/P[8].D*Determinant(Mat8c);
  
  double Mat9[3][3] = {{pt[12].x, pt[12].y, pt[12].z}, {pt[17].x, pt[17].y, pt[17].z}, {pt[13].x, pt[13].y, pt[13].z}};
  double Mat9a[3][3] = {{1, pt[12].y, pt[12].z}, {1, pt[17].y, pt[17].z}, {1, pt[13].y, pt[13].z}};
  double Mat9b[3][3] = {{pt[12].x, 1, pt[12].z}, {pt[17].x, 1, pt[17].z}, {pt[13].x, 1, pt[13].z}};
  double Mat9c[3][3] = {{pt[12].x, pt[12].y, 1}, {pt[17].x, pt[17].y, 1}, {pt[13].x, pt[13].y, 1}};
  P[9].D = Determinant(Mat9);
  P[9].d = 1;
  P[9].a = -P[9].d/P[9].D*Determinant(Mat9a);
  P[9].b = -P[9].d/P[9].D*Determinant(Mat9b);
  P[9].c = -P[9].d/P[9].D*Determinant(Mat9c);
  
  double Mat10[3][3] = {{pt[15].x, pt[15].y, pt[15].z}, {pt[16].x, pt[16].y, pt[16].z}, {pt[19].x, pt[19].y, pt[19].z}};
  double Mat10a[3][3] = {{1, pt[15].y, pt[15].z}, {1, pt[16].y, pt[16].z}, {1, pt[19].y, pt[19].z}};
  double Mat10b[3][3] = {{pt[15].x, 1, pt[15].z}, {pt[16].x, 1, pt[16].z}, {pt[19].x, 1, pt[19].z}};
  double Mat10c[3][3] = {{pt[15].x, pt[15].y, 1}, {pt[16].x, pt[16].y, 1}, {pt[19].x, pt[19].y, 1}};
  P[10].D = Determinant(Mat10);
  P[10].d = 1;
  P[10].a = -P[10].d/P[10].D*Determinant(Mat10a);
  P[10].b = -P[10].d/P[10].D*Determinant(Mat10b);
  P[10].c = -P[10].d/P[10].D*Determinant(Mat10c);
  
  double Mat11[3][3] = {{pt[17].x, pt[17].y, pt[17].z}, {pt[16].x, pt[16].y, pt[16].z}, {pt[19].x, pt[19].y, pt[19].z}};
  double Mat11a[3][3] = {{1, pt[17].y, pt[17].z}, {1, pt[16].y, pt[16].z}, {1, pt[19].y, pt[19].z}};
  double Mat11b[3][3] = {{pt[17].x, 1, pt[17].z}, {pt[16].x, 1, pt[16].z}, {pt[19].x, 1, pt[19].z}};
  double Mat11c[3][3] = {{pt[17].x, pt[17].y, 1}, {pt[16].x, pt[16].y, 1}, {pt[19].x, pt[19].y, 1}};
  P[11].D = Determinant(Mat11);
  P[11].d = 1;
  P[11].a = -P[11].d/P[11].D*Determinant(Mat11a);
  P[11].b = -P[11].d/P[11].D*Determinant(Mat11b);
  P[11].c = -P[11].d/P[11].D*Determinant(Mat11c);

  cout << "DONE" << endl;
  

  /*****************************************************************************
   *  Calculate Computational Domain
   *****************************************************************************/  
  cout << "Calculating computational domain... ";
    
  // Domain Size
  double bldg_max_x, bldg_max_y, x_intersect, boom_length_x, X_MAX, Y_MAX, Z_MAX;

  vector<double> Coord_shift_col_x(8), Coord_shift_col_y(8);
  for (i=0; i<8; i++) {
    Coord_shift_col_x[i] = Coord_shift[i][0];
    Coord_shift_col_y[i] = Coord_shift[i][1];
  }  
  bldg_max_x = *max_element(Coord_shift_col_x.begin(),Coord_shift_col_x.end());  // Building corner with the largest x-coordinate value
  bldg_max_y = *max_element(Coord_shift_col_y.begin(),Coord_shift_col_y.end());  // Building corner with the largest y-coordinate value
  x_intersect = bldg_max_x+pt[8].z*tan(theta);    // Initial position of the boom at the ground [m], based on boom wavefront being initially aligned with the coordinate point at bldg_max_x
  i_intersect = int(ceil(x_intersect/dx))+D;      // Initial position of the boom at the ground [index]
  boom_length_x = boom_length_i/fs*c0;            // Length of the boom portion of the input boom file [m], determining how far the domain should stretch in x & z directions.
  
  X_MAX = x_intersect+boom_length_x/cos(theta)+buffer_dim_x; // Domain size in x-direction [m] (excluding the PML region)
  Y_MAX = bldg_max_y+extra_dim_y;                            // Domain size in y-direction [m] (excluding the PML region)
  Z_MAX = ((pt[16].z/tan(theta)-pt[16].x)+x_intersect+boom_length_x/cos(theta))*sin(theta)*cos(theta)+buffer_dim_z;    // Domain size in z-direction [m] (excluding the PML region)

  I_MAX = int(round(X_MAX/dx))+2*D;             // Domain size in x-direction [index] (including the PML region)
  J_MAX = int(round(Y_MAX/dx))+2*D;             // Domain size in y-direction [index] (including the PML region)
  K_MAX = int(round(Z_MAX/dx))+D;               // Domain size in z-direction [index] (including the PML region)

  cout << "DONE" << endl;
  
  
  /*****************************************************************************
  * Create Building Geometry
  ******************************************************************************/  
  cout << "Creating building geometry file... ";

  // Resize U according to the domain size calculated above
  U.resize(I_MAX+2);
  for (i=0; i<I_MAX+2; i++) {
    U[i].resize(J_MAX+2);
    for (j=0; j<J_MAX+2; j++) {
      U[i][j].resize(K_MAX+2);
    }
  }

  float wall_thickness_normalized = wall_thickness*sqrt(P[0].a*P[0].a+P[0].b*P[0].b+P[0].c*P[0].c);  // Wall thickness normalized to calculate distance to P[0]
  float X[I_MAX+2], Y[J_MAX+2], Z[K_MAX+2];

  for (i=1; i<I_MAX+1; i++) X[i] = (i-D)*dx;
  for (j=1; j<J_MAX+1; j++) Y[j] = (j-D)*dx; 
  for (k=1; k<K_MAX+1; k++) Z[k] = k*dx;

  for (k=0; k<K_MAX+1; k++) {
    for (j=0; j<J_MAX+1; j++) {
      for (i=0; i<I_MAX+1; i++) {
        // Interior points of Building 1
        if ( (P[0].a*X[i]+P[0].b*Y[j]+P[0].c*Z[k]+P[0].d)>0 && (P[1].a*X[i]+P[1].b*Y[j]+P[1].c*Z[k]+P[1].d)<0 && (P[2].a*X[i]+P[2].b*Y[j]+P[2].c*Z[k]+P[2].d)>0 && (P[3].a*X[i]+P[3].b*Y[j]+P[3].c*Z[k]+P[3].d)>0 && (P[4].a*X[i]+P[4].b*Y[j]+P[4].c*Z[k]+P[4].d)>0 && (P[5].a*X[i]+P[5].b*Y[j]+P[5].c*Z[k]+P[5].d)>0 )
          U[i][j][k].InOrOut = 1;
        // Interior points of Building 2
        else if ( (P[6].a*X[i]+P[6].b*Y[j]+P[6].c*Z[k]+P[6].d)>0 && (P[7].a*X[i]+P[7].b*Y[j]+P[7].c*Z[k]+P[7].d)<0 && (P[8].a*X[i]+P[8].b*Y[j]+P[8].c*Z[k]+P[8].d)<0 && (P[9].a*X[i]+P[9].b*Y[j]+P[9].c*Z[k]+P[9].d)>0 && (P[10].a*X[i]+P[10].b*Y[j]+P[10].c*Z[k]+P[10].d)<0 && (P[11].a*X[i]+P[11].b*Y[j]+P[11].c*Z[k]+P[11].d)>0 )
          U[i][j][k].InOrOut = 1;
        // Interior points of Wall
        else if ( (k*dx<=wall_height) && (P[0].a*X[i]+P[0].b*Y[j]+P[0].c*Z[k]+P[0].d)<=wall_thickness_normalized && (P[0].a*X[i]+P[0].b*Y[j]+P[0].c*Z[k]+P[0].d)>0 && (P[1].a*X[i]+P[1].b*Y[j]+P[1].c*Z[k]+P[1].d)>=0 && (P[6].a*X[i]+P[6].b*Y[j]+P[6].c*Z[k]+P[6].d)<=0)
          U[i][j][k].InOrOut = 1;
        // Exterior points
        else
          U[i][j][k].InOrOut = 0;
      }
    }
  }
  
  cout << "DONE" << endl;


  /*****************************************************************************
   *  Microphone selection
   *****************************************************************************/
  cout << "Reading microphone locations where pressure is to be recorded... ";
  
  ifstream f_par_mic("In_parameters_mic.txt");
  if (!f_par_mic) {
    cout << "ERROR: Could not open In_parameters_mic.txt" <<endl;
    exit(1);
  }
  Mic mic_temp;
  string line;
  while (f_par_mic.good()) {
    getline (f_par_mic,line);
    if (line=="CA") {
      mic_temp.i = int(round(CMics[0][0]/dx))+D;
      mic_temp.j = int(round(CMics[0][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="LB") {
      mic_temp.i = int(round(CMics[1][0]/dx))+D;
      mic_temp.j = int(round(CMics[1][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="CC") {
      mic_temp.i = int(round(CMics[2][0]/dx))+D;
      mic_temp.j = int(round(CMics[2][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="CD") {
      mic_temp.i = int(round(CMics[3][0]/dx))+D;
      mic_temp.j = int(round(CMics[3][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="CE") {
      mic_temp.i = int(round(CMics[4][0]/dx))+D;
      mic_temp.j = int(round(CMics[4][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="CF") {
      mic_temp.i = int(round(CMics[5][0]/dx))+D;
      mic_temp.j = int(round(CMics[5][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="CG") {
      mic_temp.i = int(round(CMics[6][0]/dx))+D;
      mic_temp.j = int(round(CMics[6][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="CH") {
      mic_temp.i = int(round(CMics[7][0]/dx))+D;
      mic_temp.j = int(round(CMics[7][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="CI") {
      mic_temp.i = int(round(CMics[8][0]/dx))+D;
      mic_temp.j = int(round(CMics[8][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="CJ") {
      mic_temp.i = int(round(CMics[9][0]/dx))+D;
      mic_temp.j = int(round(CMics[9][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="WA") {
      mic_temp.i = int(round(WMics[0][0]/dx))+D;
      mic_temp.j = int(round(WMics[0][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="WB") {
      mic_temp.i = int(round(WMics[1][0]/dx))+D;
      mic_temp.j = int(round(WMics[1][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="WC") {
      mic_temp.i = int(round(WMics[2][0]/dx))+D;
      mic_temp.j = int(round(WMics[2][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="WD") {
      mic_temp.i = int(round(WMics[3][0]/dx))+D;
      mic_temp.j = int(round(WMics[3][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="WE") {
      mic_temp.i = int(round(WMics[4][0]/dx))+D;
      mic_temp.j = int(round(WMics[4][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="WF") {
      mic_temp.i = int(round(WMics[5][0]/dx))+D;
      mic_temp.j = int(round(WMics[5][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="WG") {
      mic_temp.i = int(round(WMics[6][0]/dx))+D;
      mic_temp.j = int(round(WMics[6][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="WH") {
      mic_temp.i = int(round(WMics[7][0]/dx))+D;
      mic_temp.j = int(round(WMics[7][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="WY") {
      mic_temp.i = int(round(WMics[8][0]/dx))+D;
      mic_temp.j = int(round(WMics[8][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="WZ") {
      mic_temp.i = int(round(WMics[9][0]/dx))+D;
      mic_temp.j = int(round(WMics[9][1]/dx))+D;
      mic_temp.k = 1;
      mics.push_back(mic_temp);
    }
  }
  f_par_mic.close();
  num_mics = mics.size();   // Total number of microphones recording
  
  cout << "DONE" << endl;
  
  
  /*****************************************************************************
   *  Domain & Microphone Parameters
   *****************************************************************************/
  // Display the domain size 
  long long TotalNumPt = (long long)I_MAX*(long long)J_MAX*(long long)K_MAX;
  cout << endl;
  cout << "\t Size of the domain (X_MAX, Y_MAX, Z_MAX) = (" << (I_MAX-2*D)*dx << ", " << (J_MAX-2*D)*dx << ", " << (K_MAX-D)*dx << ")" << endl;
  cout << "\t Size of the domain (I_MAX, J_MAX, K_MAX) = (" << I_MAX << ", " << J_MAX << ", " << K_MAX << ")" << endl;
  cout << "\t Total number of grid points = " << TotalNumPt << endl;
  if (boom_length_i > input_length_i)
    cout << "WARNING: boom_length > input_length" << endl; 
    
      
  return;
}

double Determinant(double Mat[3][3])
{
  return Mat[0][0]*Mat[1][1]*Mat[2][2]+Mat[0][1]*Mat[1][2]*Mat[2][0]+Mat[0][2]*Mat[1][0]*Mat[2][1]-Mat[0][0]*Mat[1][2]*Mat[2][1]-Mat[0][1]*Mat[1][0]*Mat[2][2]-Mat[0][2]*Mat[1][1]*Mat[2][0];
}