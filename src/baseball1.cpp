///
/// Starter template for first baseball problem
/// Solve for the initial speed of the pitch given the initial parameters
/// xend : distance to home plate [18.5] m
/// z0 : height of release of ball [1.4] m
/// theta0 : angle of release above horizontal [1] degree
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
#include <cmath> 
#include <vector> 

using namespace std;

struct Params {
  double g;   // acceleration [m/s^2]
  double m;   // mass of object [kg], nb proj. In vacuum funcs do not depend on the mass
  double d;   // m diameter of ball
  double b;   // b,c params for air resistance
  double c;
};

int main(int argc, char **argv){

  // examples of parameters
  Params pars;
  pars.g=9.81;      // acceleration of gravity (m/s^2)
  pars.m=0.145;     // mass of baseball (kg)
  pars.d=0.0075;    // diameter of baseball (m)
  pars.b=1.6e-4;    // drag / diameter ( N / v * m )
  pars.c=0.25;      // drag / diameter^2 ( N / v^2 * m^2 )
  
  double xend=18.5;        // meters to plate
  double z0=1.4;           // height of release [m]
  double theta0=1;         // angle of velocity at release (degrees)
                                      // convert to radians before using!
  bool showPlot=false;    // keep this flag false by default
  
  // allow changing the parameters from the command line
  int c;
  while ((c = getopt (argc, argv, "x:z:t:p")) != -1)
    switch (c) {
    case 'x':
      xend = atof(optarg);
      break;
    case 'z':
      z0 = atof(optarg);
      break;
    case 't':
      theta0 = atof(optarg);
      break;
    case 'p':
      showPlot=true;
      break;
    case '?':
      fprintf (stderr, "Unknown option `%c'.\n", optopt);
    }
  TApplication theApp("App", &argc, argv); // init ROOT App for displays


  double vPitch = 100.;   // m/s of pitch needed to land in strike zone at 0.9 meters
  // write code to solve for vPitch here

  //so, the goal of this code is to determine the minimum velocity (vPitch) at which the ball will arrive over home-plate 
  // within the 'strike-zone' (z_final >= 0.9 m).
  
  // So, let's establish the system of equations: 
  //now, let's turn air-resistance back on. 

  //this is just a convienient renaming for the template-structs of PhasePoint_t<> and FirstDerivativeSystem<>
  constexpr unsigned int DOF = 6; //the number of dependent variables in our phase-space
  using XVpt = PhasePoint_t<DOF>; 
  using XVsystem = FirstDerivativeSystem<6>;

  //a system of our 6 dependent variables, wich returns a 6-dimensional array of the fist time-derivatives of each dependent variable. 
  XVsystem sys_with_drag = [&pars](const XVpt& pt){
      
      double v = sqrt( pt.X[3]*pt.X[3] + pt.X[4]*pt.X[4] + pt.X[5]*pt.X[5] );
      
      double b = pars.b * pars.d; 
      double c = pars.c * pars.d * pars.d; 

      return array<double,6>{
          pt.X[3],    // dx/dt = vx
          pt.X[4],    // dy/dt = vy
          pt.X[5],    // dz/dt = vz
          -pt.X[3]*( b  +  c*v )/pars.m,            // dvx/dt = -F(v)*vx/m
          -pt.X[4]*( b  +  c*v )/pars.m,            // dvy/dt = -F(v)*vy/m
          -pt.X[5]*( b  +  c*v )/pars.m  - pars.g   // dvz/dt = -F(v)*vz/m - g  
      };
  };

  const double z_strike = 0.9; //minimum height of ball once it reaches home-plate, in meters

  //our stopping condition; we want to return if our baseball is out-of-bounds
  StoppingCondition<6> end_pitch_test = [&pars,xend,z_strike](const XVpt& pt) 
  {
    //return 'true' if the stopping condition is met (so we should stop iterations), and 'false' if iterations should continue. 
    if (pt.X[2] < z_strike || pt.X[0] > xend) return true;
    return false; 
  };

  //now, let's attempt to change 'v', using the bisection method. 
  //returns 'true' if pitch is in the strike-zone, and false otherwise. 
  auto is_pitch_strike = [&end_pitch_test, &sys_with_drag,theta0,z0,xend,z_strike](const double v0)
  {
    //convert our angle to radians
    double thet = theta0 * 0.0174532925199; 

    //supply the initial position & velocities 
    XVpt initial_conditions{ 0., {0.,0.,z0,  v0*cos(thet),0.,v0*sin(thet)}}; 

    auto pts = RungeKutta4N<6>(sys_with_drag, initial_conditions, 0.01, 250, &end_pitch_test); 
    
    auto ending_pt = pts.back(); 
    if (ending_pt.X[0] < xend || ending_pt.X[2] < z_strike) return false; 
    return true; 
  };  

  //now, let's use the bisection method to determine if the baseball was a strke or not. 
  const double tolerance = 1e-3; //in m/s
  double v1 = vPitch;

  //check to make sure our first guess is fast enough
  while (!is_pitch_strike(v1)) v1 *= 2.; 

  double v0 = 0.;  
  while (v1-v0 > tolerance) {

    //check if the midpoint velocity is fast enough.
    double v_midpoint = 0.5*(v1 + v0); 
    if (is_pitch_strike(v_midpoint)) { 
      //if the v_midpoint IS fast enough, then the minimum velocity is LESS than v_midpoint. 
      v1 = v_midpoint; continue; 
    } else {
      //if the v_midpoint IS NOT fast enough, then the minimum velocity is MORE than v_midpoint. 
      v0 = v_midpoint; continue; 
    }
  }

  //this is our best guess, the 'true' answer should be: tolerance < | v_min - vPitch |, if all went well. 
  vPitch = 0.5 * (v1+v0); 

  // do not change these lines
  printf("********************************\n");
  printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n",xend,z0,theta0);
  printf("v_pitch = %lf m/s\n",vPitch);
  printf("********************************\n");

  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}

