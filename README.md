# Robotics2
Tools and exercises for the "robotics2" course

## Parametrization 
Parametrization.m is a script which estimates uknown parameters (mass, inertia, arm lengths..) of an RPR manipulator. Just launch it on MATLAB and have fun <br/>
**The Process**
I follow the Lagrangian formulation of the robot dynamics<img src="https://render.githubusercontent.com/render/math?math=M(q)+c(q,\dot{q})+g(q)=\tau"> <br/>
M,c and g contains unkown and known parameters, so the dynamics must be rewritten as 
