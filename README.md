# ODE-3 

## Compiling 
Everything should compile the same way as in the original fork, using ```make``` at the top-level dir of the repo. All executables (```vterm```, ```baseball1```, ```baseball2```) are in the ```src``` directory. 

## custom version of ode-lib

I wrote my own 4th-order Runge-Kutta solver, you can see its imlpementation in ```src/RungeKutta.hpp```. I did this for a few reasons:  
- I wanted to understand the R-K method better
- I wanted practice using C++ templates, which my solver uses.
- Specifically, the Solver method ```RungeKutta4N``` has a template-argument which is an unsigned integer (as you can see written my my source-code files as ```RungeKutta4N<6>```). This number 6 specifies the number of dependent variables in our system of coupled ODEs. The reason I bothered using a template for this, rather than just using a ```std::vector``` (where I would not have to specify the DOF of my ODE system beforehand), is that when you use templates like this, you _must_ specify the DoF of your ODE system at compile-time. in other words, If I used vectors, I could write a lot of code which would compile just fine, but would break in weird and unpredictable ways if I mismatch the number of independent vars in different places. If you specify exactly how many DOF you have with a template, the code will refuse to compile if you accidentally use the wrong DOF somehwere in your code.

## RungeKutta solver 
The way the RungeKutta solver works is this: 
#### PhasePoint_t 
This is a small struct which contains a double for the independent-var (t), and an array of points, each of which is the value of different independent variables in phase-space. the template is the number of dependent vars a system has. for example, for a baseball with three positions and three velocities, the phase-space has 6 dependent vars (x,y,z + vx,vy,vz), and thus a phase-space point for this problem could be: 
```
PhasePoint_t<6> point;
```
#### FirstDerivativeSystem
This is just a re-name of a ```std::function<>``` template, which takes in a single PhasePoint_t sturct (a single point in phase-space), and computes the first derivate of each dependent variable at the given phase-space point. See the executables in the ```src``` directory for concrete examples of ```FirstDerivativeSystem``` objects. 

#### StoppingCondition
This is another rename of a ```std::function<>``` template, which takes a PhaePoint_t point, but it returns a bool - if this function returns 'true' then iterations are terminated at this phase-point, and if it returns 'false', then iterations continue. See the executables in the ```src``` directory for concrete examples of ```StoppingCondition``` objects. 

#### RungeKutta4N 
This is what actually applies the 4th-order Runge-Kutta integration method. The arguments are: 
+ ```FirstDerivativeSystem``` 'first_derivs': this is the systems of coupled first-derivatives; one for each dependent variable.
+ ```PhasePoint_t``` 'starting_point': this is the starting point in phase space, from which each dependent variable will be evolved forward in time
+ ```double h```: the fixed-size time-step
+ ```max_iterations``` the max number of time-steps to be integrated.
+ ```StoppingCondition*``` 'stop_iterations' A ptr to a ```StoppingCondition``` fcn object. If a ```nullptr``` is provided for this argument, then no stopping condition is evaluated; iterations will continue until ```max_iterations``` is reached. If it isn't null, then this function is evaluated at each integration step; if it ever returns ```true```, then iterations are ceased.

This funciton returns a vector of ```PhasePoint_t``` objects, which are the results of integrating from the starting point supplied; they are spaced in time by ```h```, and include the initial phase-point supplied. 

