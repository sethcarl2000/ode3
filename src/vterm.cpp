//Local headers
#include "RungeKutta.hpp"
//ROOT headers
#include <TCanvas.h>
#include <TGraph.h> 
//stdlib headers
#include <functional>
#include <vector>
#include <array>
#include <cstdio>
#include <iostream>  
#include <cmath> 

constexpr int DOF = 6; 

int main(int argc, char* argv[])
{
    using namespace std; 

    //compute the velocity of the ball in vacuum. 

    //our first system has a phase-space of 6 dependent variables: 
    // x,y,z vx,vy,vz   
    // therefore, we want a 6-dimensional phase-space point (and correspoding firt-deriv system): 
    using XVpt = PhasePoint_t<DOF>; 
    using XVsystem = FirstDerivativeSystem<6>; 


    //First, we're going to integrate the equaitons of motion for a projectile with NO air resistance, 
    // and an acceleration 'g' in the z-direction (and a mass of unity)
    const double m = 1.;    //kg 
    const double g = 9.81;  //m/s

    //now, let's define a stopping condition. this stops iterations if our particle 'hits the ground' (z <= 0). 
    StoppingCondition<6> stop_at_ground = [](const XVpt& pt){ return pt.X[2] <= 0.; };

    //this is our system of equations for the free particle (in vacuum). 
    //the 'XVpt' represents the input time (t) and phase-space coords (x,y,z,vx,vy,vz: our dependent vars.).
    //each element of the array it returns represents the first-deriviative of a dependent variable, computed
    //at input time 't'. Of course, these first derivatives for the free-particle don't have explicit 't'-depedance. 
    XVsystem sys_vacuum_with_grav = [m,g](const XVpt& pt){
        //double v = sqrt( pt.X[3]*pt.X[3] + pt.X[4]*pt.X[4] + pt.X[5]*pt.X[5] );
        return array<double,6>{
            pt.X[3],    // dx = vx
            pt.X[4],    // dy = vy
            pt.X[5],    // dz = vz
            0.,             // dvx = 0.
            0.,             // dvy = 0.
            0. - (m*g)/m    // dvz = 0. - g
        };
    };

    
    //now, let's try some iterations
    XVpt starting_point{ .t=0., .X={0.,0.,3., 0.,0.,2.} };   

    cout << "attempting integration..." << flush; 
    auto points = RungeKutta4N<6>( sys_vacuum_with_grav, starting_point, 0.01, 200, &stop_at_ground ); 
    cout << "done" << endl; 

    cout << "making graphs from points..." << flush; 
    auto gs = MakeTGraphsFromPts<6>(points); 
    cout << "done." << endl; 

    cout << "Drawing graph..." << flush; 
    auto c = new TCanvas; 
    gs[2]->Draw(); 
    gs[0]->SetLineColor(kRed);  gs[0]->Draw("SAME");
    gs[1]->SetLineColor(kBlue); gs[1]->Draw("SAME");
    cout << "done." << endl; 

    c->SaveAs("test1.png"); 
    
    return 0; 
}
