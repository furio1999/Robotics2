# Robotics2
Tools and exercises for the "robotics2" course

## Parametrization 
Parametrization.m is a script which estimates uknown parameters (mass, inertia, arm lengths..) of an RPR manipulator. Just launch it on MATLAB and have fun <br/>
## The process
I follow the Lagrangian formulation of the robot dynamics, with inertia matrix M and Coriolis terms c. 

I have 2 remarks:
1) For this particular case I assume g(q)=0, given that our RPR Robot lies on an horizontal plane <br/>
2) Known parameters: mass, length and CoM position of link3  <br/>
