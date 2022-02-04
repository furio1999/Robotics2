# Robotics2
Tools and exercises for the "robotics2" course

## Parametrization 
Parametrization.m is a script which estimates uknown parameters (mass, inertia, arm lengths..) of an RPR manipulator. Just launch it on MATLAB and have fun <br/>
## The process
I follow the Lagrangian formulation of the robot dynamics <img src="https://render.githubusercontent.com/render/math?math=M(q) \, c(q,\dot{q}) \, g(q)  \, \tau"> <br/>
M,c and g contains unkown and known parameters, so the dynamics must be rewritten as <img src="https://render.githubusercontent.com/render/math?math=M(q) \, c(q,\dot{q}) \,Y_ka_k \, Y_ua_u"> <br/>

I have 2 remarks:
1) For this particular case I assume g(q)=0, given that our RPR Robot lies on an horizontal plane <br/>
2) Known parameters: <img src="https://render.githubusercontent.com/render/math?math=m_3 \, dc3 \, l_3>  <br/>
