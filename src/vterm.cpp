//Local headers
#include "RungeKutta.hpp"
//ROOT headers
#include <TCanvas.h>
#include <TGraph.h>
#include <TPad.h>  
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
    double m = 0.145; //kg 
    double g = 9.81;  //m/s

    //now, let's define a stopping condition. this stops iterations if our particle 'hits the ground' (z <= 0). 
    StoppingCondition<6> stop_at_ground = [](const XVpt& pt){ return pt.X[2] <= 0.; };

    //this is our system of equations for the free particle (in vacuum). 
    //the 'XVpt' represents the input time (t) and phase-space coords (x,y,z,vx,vy,vz: our dependent vars.).
    //each element of the array it returns represents the first-deriviative of a dependent variable, computed
    //at input time 't'. Of course, these first derivatives for the free-particle don't have explicit 't'-depedance. 
    XVsystem sys_vacuum_with_grav = [m,g](const XVpt& pt){
        return array<double,6>{
            pt.X[3],    // dx = vx
            pt.X[4],    // dy = vy
            pt.X[5],    // dz = vz
            0.,             // dvx = 0.
            0.,             // dvy = 0.
            0. - (m*g)/m    // dvz = 0. - g
        };
    };

    //Let's integrate this with a bunch of different step-sizes, and see how energy conservation is (or is not) conserved.
    //The way we'll do that is this: 
    // Let a projectile go at z0 = 5m, and then record its kinetic energy when it hits the ground. 
    // Then, compare the sum of potential (U) and kinetic (KE) before and after the integration of motion. 
    const double h_min = 1e-2;
    const double h_max = 1.;
    int n_points = 10; 
    double h=h_min; 
    for (int i=0; i<n_points; i++) {
        
        //now, let's try some iterations
        XVpt starting_point{ .t=0., .X={0.,0.,5., 0.,0.,0.} };   

        auto points = RungeKutta4N<6>( sys_vacuum_with_grav, starting_point, h, 200, &stop_at_ground ); 

        XVpt ending_point = points.back(); 

        //check the last point. compute the kinetic energy, and also the potential energy. 
        
        //given an 'XVpt', computes the kinetic energy
        auto kinetic_energy = [m](const XVpt& pt){ 
            double v2 = pt.X[3]*pt.X[3] + pt.X[4]*pt.X[4] + pt.X[5]*pt.X[5]; 
            return 0.5*m*v2;
        }; 

        //compute potential energy of a point (U=0 at z=0)
        auto potential_energy = [m,g](const XVpt& pt){
            return m * g * pt.X[2];
        };

        //compute the difference in kinetic & potential energy
        double delta_KE = kinetic_energy(ending_point) - kinetic_energy(starting_point); 
        
        double delta_U  = potential_energy(ending_point) - potential_energy(starting_point); 

        printf(
            "Step size: %.3e s\n"
            "Final coordinate (x,y,z): (% .2f m, % .2f m, % .2f m)\n"
            "   Change in kinetic:     %+.6f J\n"
            "   Change in potential:   %+.6f J\n"
            "   Energy discrepancy:    %+.3e J\n",
            h, 
            ending_point.X[0],ending_point.X[1],ending_point.X[2], 
            delta_KE, 
            delta_U,
            delta_KE + delta_U
        );

        h *= pow(h_max/h_min, 1./((double)n_points-1)); 
    }

    const double diameter = 0.075; // (m) diameter of the baseball
    
    //coefficients of drag (for v^2 and v force, respectivley)
    const double Cd = 0.25;     // ( N t^2 / m^2 ) v^2, drag force coefficient
    const double Bd = 1.6e-4;   // ( N t / m ) drag-force coefficient

    
    //now, let's turn air-resistance back on. 
    XVsystem sys_with_drag = [&m,g,diameter,Bd,Cd](const XVpt& pt){
        
        double v = sqrt( pt.X[3]*pt.X[3] + pt.X[4]*pt.X[4] + pt.X[5]*pt.X[5] );
        
        double b = Bd * diameter; 
        double c = Cd * diameter * diameter; 

        return array<double,6>{
            pt.X[3],    // dx = vx
            pt.X[4],    // dy = vy
            pt.X[5],    // dz = vz
            -pt.X[3]*( b  +  c*v )/m,       // dvx = vx
            -pt.X[4]*( b  +  c*v )/m,       // dvy = vy
            -pt.X[5]*( b  +  c*v )/m  - g   // dvz = vz  
        };
    };

    //exit iterations when the magnitude of acceleration falls below this value
    const double min_acceleration = 1e-4; 
    StoppingCondition<6> stop_min_acceleration = [&sys_with_drag, min_acceleration](const XVpt& pt)
    {   
        //compute the magnitude of the acceleration
        auto&& dX_dt = sys_with_drag(pt); 
        double mag_acceleration = sqrt( dX_dt[3]*dX_dt[3] + dX_dt[4]*dX_dt[4] + dX_dt[5]*dX_dt[5] );
        return mag_acceleration < min_acceleration; 
    };  

    
    //return the magnitude of the velocity at a phase-space point    
    auto velocity_mag = [](const XVpt& pt) { return sqrt( pt.X[3]*pt.X[3] + pt.X[4]*pt.X[4] + pt.X[5]*pt.X[5] ); };


    XVpt starting_point{ 0., {0.,0.,0., 0.,0.,0.}};
    auto points = RungeKutta4N<6>(sys_with_drag, starting_point, 0.10, 10000, &stop_min_acceleration);
    printf(
        "Number of points to reach terminal velocity: %zi (stopping condition: |a| < %.1e m/s^2)\n"
        "   terminal velocity (default parameters): %.4f m/s \n",
        points.size(),
        min_acceleration,
        velocity_mag(points.back())
    );  

    TCanvas *c; 
    c = new TCanvas; 
    auto gs = MakeTGraphsFromPts<6>(points); 
    gs[5]->SetTitle("z-velocity vs time;time (s);z-velocity (m/s)"); 
    gs[5]->Draw(); 
    c->SaveAs("v_vs_time.png"); 

    //now, let's study the terminal velocity as a function of mass: 
    const double mass_min = 0.001; // kg 
    const double mass_max = 10.;   // kg 

    m = mass_min; 
    n_points = 40; 

    vector<double> pts_mass, pts_vt; 
    for (int i=0; i<n_points; i++) {

        starting_point = XVpt{ 0., {0.,0.,0., 0.,0.,0.}};
        points = RungeKutta4N<6>(sys_with_drag, starting_point, 0.10, 10000, &stop_min_acceleration);
        
        //record the mass, and the reported velocity at the last point. 
        pts_mass.push_back(m); 
        pts_vt  .push_back(velocity_mag(points.back()));
        //printf(" n.points = %-5zi, v-term: %.4f m/s\n", points.size(), pts_vt.back()); 

        m *= pow( mass_max/mass_min, 1./((double)n_points-1) );
    }
    
    c = new TCanvas; 
    auto graph = new TGraph(pts_mass.size(), pts_mass.data(), pts_vt.data()); 
    
    gPad->SetLogx(1);
    //gPad->SetLogy(1); 
    graph->SetTitle("Baseball mass vs. Terminal Velocity;mass (kg);Terminal velocity (m/s)");
    graph->SetMarkerStyle(kOpenCircle);
    graph->SetMarkerSize(0.7);  
    graph->Draw("ALP"); 

    c->SaveAs("vterm.png"); 

    return 0; 
}
