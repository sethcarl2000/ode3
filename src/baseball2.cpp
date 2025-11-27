///
/// Starter template for second baseball problem
///
///  Do not change the interface for running the program
///  Fill in the value of vPitch in the print statement with your solution
///  at the end of main()
///

#include "RungeKutta.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <functional> 
#include <cmath> 

using namespace std;


int main(int argc, char **argv){

  // we have 6 initial conditions for this problem
  // y[0] = y[2] = y[4] = 0;  // init x,y,z
  // y[1] = v0*cos(theta0);   // vx  "x is line towards the plate
  // y[3] = 0;                // vy  "y" is measured as left/right divergence from line to plate
  // y[5] = v0*sin(theta0);   // vz  "z" is vertival measure
  vector<double> y0(6);

  bool showPlot=false;
  // pitches
  // slider ip=0
  // curve ip=1
  // screwball ip=2
  // fast ip=3
  int ip=1;    // default pitch
  int c;
  while ((c = getopt (argc, argv, "p:n")) != -1)
    switch (c) {
    case 'p':
      ip = atoi(optarg);
      break;
    case 'n':
      showPlot=false;
      break;
    }

  //angle the baseball's axis of rotation makes with the horizontal (+y) axis.
  // its assumed that the baseball's axis of rotation is perp. to the baseball's direction of travel. 
  double phi = 0.; // rad

  constexpr double pi = 3.14159265359; 
  TString title;
  if (ip==0){
    cout << "Setting up initial conditions for slider" << endl;
    phi = 0; 
  }
  else if (ip==1){
    cout << "Setting up initial conditions for curveball" << endl;
    phi = pi / 4.; 
  }
  else if (ip==2){
    cout << "Setting up initial conditions for screwball" << endl;
    phi = -pi / 4.; 
  }
  else {
    cout << "Setting up initial conditions for fastball" << endl;
    phi = 0.45 * pi; 
  }

  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  double xend=60;   // feet
  double yend=0;    // tbd
  double zend=0;    // tbd
  double vxend=0;
  double vyend=0;
  double vzend=0;

  // write code here
  
  //this is just a convienient renaming for the template-structs of PhasePoint_t<> and FirstDerivativeSystem<>
  constexpr unsigned int DOF = 6; //the number of dependent variables in our phase-space
  using XVpt = PhasePoint_t<DOF>; 
  using XVsystem = FirstDerivativeSystem<6>;


  const double g  =9.81;      // acceleration of gravity (m/s^2)
  const double m  =0.145;     // mass of baseball (kg)
  const double d  =0.075;     // diameter of baseball (m)
  const double Cd =0.25;      // drag / diameter^2 v^2
  const double Cb =1.6e-4;    // drag / diameter v ( N / v * m )
  const double B  =4.1e-4;    // magnitude of magnus force ( N / Hz )
  const double omega = 25. * 6.28318530718; // angular velocity of baseball (rad/s) 

  

  //a system of our 6 dependent variables, wich returns a 6-dimensional array of the fist time-derivatives of each dependent variable.
  //in this system of depedent variables, the velocites depend change from drag and the magnus force.  
  XVsystem sys_drag_and_magnus = [g,m,d,Cd,Cb,B,omega,phi](const XVpt& pt){
      
    double v = sqrt( pt.X[3]*pt.X[3] + pt.X[4]*pt.X[4] + pt.X[5]*pt.X[5] );

    double b = Cb * d; 
    double c = Cd * d*d;

    return array<double,6>{
        pt.X[3],    // dx/dt = vx
        pt.X[4],    // dy/dt = vy
        pt.X[5],    // dz/dt = vz
        -pt.X[3]*( b  +  c*v )/m  -  B*omega*( pt.X[5]*sin(phi) - pt.X[4]*cos(phi) ),            
        -pt.X[4]*( b  +  c*v )/m  -  B*omega*pt.X[3]*cos(phi),           
        -pt.X[5]*( b  +  c*v )/m  -  B*omega*pt.X[4]*sin(phi)   - g   
    };
  };

  const double z_strike = 0.9; //minimum height of ball once it reaches home-plate, in meters

  const double ymax{10.}, zmax{10.}; 
  //our stopping condition; we want to return if our baseball is out-of-bounds
  StoppingCondition<6> end_pitch_test = [xend,ymax,zmax](const XVpt& pt) 
  {
    //return 'true' if the stopping condition is met (so we should stop iterations), and 'false' if iterations should continue. 
    if (pt.X[0] > xend*0.3048 ||
        fabs(pt.X[1]) > ymax  ||
        fabs(pt.X[2]) > zmax) return true; 
    return false;  
  };

  //initial conditions (in m/s)
  XVpt starting_point{ 0., {0.,0.,0., 0.,1.,40.}}; 

  auto points = RungeKutta4N<6>(sys_drag_and_magnus, starting_point, 0.05, 1000, &end_pitch_test); 

  auto gs = MakeTGraphsFromPts<6>(points); 

  XVpt last_point = points.back(); 
  xend = last_point.X[0]; 
  yend = last_point.X[1]; 
  zend = last_point.X[2]; 

  vxend = last_point.X[3]; 
  vyend = last_point.X[4]; 
  vzend = last_point.X[5]; 

  // to compare to the plots in Fitzpatrick, output your results in **feet**
  // do not change these lines
  printf("********************************\n");
  printf("Coordinates when x=60 feet\n");
  printf("(x,y,z) = (%lf,%lf,%lf)\n",xend,yend,zend);
  printf("(vx,vy,vz) = (%lf,%lf,%lf)\n",vxend,vyend,vzend);
  printf("********************************\n");

  // plot the trajectory.  See Fitzpatrick for plot details
  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}

