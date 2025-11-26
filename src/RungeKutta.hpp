#ifndef RungeKutta_HPP
#define RungeKutta_HPP

#include <functional> 
#include <vector> 
#include <array> 
#include <TGraph.h> 
#include <iostream> 

//A point in phase-space, with a phase-space of dependent variables of size 'DOF' 
template<unsigned int DOF> struct PhasePoint_t {
    
    //independent variable
    double t; 
    //dependent variables
    std::array<double,DOF> X; 

    //DoF of this phase-space point
    static unsigned int size() { return DOF; }
};

//system of funcitons with 'DOF' independent variables. 
// input is a single phase-space point (PhasePoint_t<DOF>), 
// output is an array of first time-derivatives for each dependent coordiante of the phase-space 
template<int DOF> using FirstDerivativeSystem =
    std::function<std::array<double,DOF>(const PhasePoint_t<DOF>&)>; 

// A function which takes as input a point in phase-space (a t-value, and a list of dependent-variable values)
// and returns 'true' if iterations should be stopped, and 'false' if iterations should continue. 
// To reiterate, design this function such that: 
// return true;     -   Stop iterations at this point
// return false;    -   Continue iterations
template<int DOF> using StoppingCondition = 
    std::function<bool(const PhasePoint_t<DOF>&)>; 

//Takes an initial phase-space point, and a system of first-derivative rules for each dependent variable, and 
// uses the 4th-order Runge-Kutta method for a N-system of coupled ODE's to integrate a solution. 
template<unsigned int DOF> std::vector<PhasePoint_t<DOF>> RungeKutta4N(
    const FirstDerivativeSystem<DOF>& first_derivs,     //system of first t-derivatives for each dependent variable
    const PhasePoint_t<DOF>& starting_point,            //starting-point in phase-space
    const double h=1e-2,                                //fixed time step-size
    const size_t max_iterations=100,                    //maximum iterations, after which iterations will stop
    StoppingCondition<DOF> *stop_iterations=nullptr     //stopping condition, returns 'true' if iterations should stop. 
                                                        // if this arg is 'nullptr', then no stopping condition is evaluated (iterations continue
                                                        // until 'max_iterations' is reached). Otherwise, if a ptr to a valid funciton is supplied, 
                                                        // iterations will stop after the function evaluates to 'true' for the first time. 
)
{
    using namespace std; 
    using Point = PhasePoint_t<DOF>;
    using uint = unsigned int; 

    vector<Point> points{starting_point}; 
    
    const bool use_stopping_condition = (stop_iterations!=nullptr);

    //here are the 4-th order, runge-kutta methods:
    while(points.size() < max_iterations) {

        //start from the most recently added point
        Point pt = points.back(); 

        //We must compute the runge-kutta method for each point. 
        array<double,DOF> K1, K2, K3, K4; 

        K1 = first_derivs(pt);

        Point pt_K2 = pt; 
        pt_K2.t += h/2.;
        for (uint i=0; i<pt.X.size(); i++) pt_K2.X[i] += K1[i] * h/2; 

        K2 = first_derivs(pt_K2); 

        Point pt_K3 = pt; 
        pt_K3.t += h/2.; 
        for (uint i=0; i<pt.X.size(); i++) pt_K3.X[i] += K2[i] * h/2;
        
        K3 = first_derivs(pt_K3); 

        Point pt_K4 = pt; 
        pt_K4.t += h; 
        for (uint i=0; i<pt.X.size(); i++) pt_K4.X[i] += K3[i] * h;

        K4 = first_derivs(pt_K4); 

        pt.t += h;
        for (uint i=0; i<pt.X.size(); i++) {
            pt.X[i] += (K1[i] + 2.*K2[i] + 2.*K3[i] + K4[i]) * h/6.; 
        }

        points.push_back(pt); 

        //if the 'stopping condition' function is applied, then check if this condition is met, 
        // and, if so, then quit iterations. 
        if (use_stopping_condition && (*stop_iterations)(pt)) return points; 
    }

    return points; 
}

//make an array of TGraphs from the arrays of points
template<int DOF> std::vector<TGraph*> MakeTGraphsFromPts(const std::vector<PhasePoint_t<DOF>>& points)
{
    using namespace std; 
    using uint = unsigned int; 

    vector<double> t{}; t.reserve(points.size()); 
    vector<vector<double>> X{};
    
    //put all time-points in a vector
    for (const auto& pt : points) t.push_back(pt.t); 

    for (uint i=0; i<DOF; i++) {
        X.push_back({});
        X.back().reserve(points.size());
    }

    //put all space-points in their own vectors
    for (const auto& pt : points) {
        for (uint i=0; i<DOF; i++) X[i].push_back(pt.X[i]);
    }

    //put all this data into different graphs
    vector<TGraph*> graphs; 
    for (uint i=0; i<DOF; i++) {
        graphs.push_back(new TGraph(points.size(), t.data(), X[i].data()));
    }
    return graphs; 
}

#endif 