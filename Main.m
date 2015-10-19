clear; clc

% just a line of code to see if I can send stuff to the github rep
timestep = 1e-3;
duration = 300;
initial_x = 1;
initial_y = 2;
initial_z = 1;

% prealloc
time = 0:timestep:duration;
x = zeros(1,length(time)); y = x; z= x;

% parametrization
a = 1; b = 1; c = 0.5; d=0.9; e=0.7; f=0.9; g = 0.8; 

x(1) = initial_x; y(1) = initial_y; z(1) = initial_z;
for i = 2:duration/timestep+1
    i_ = i-1;
    %
    dxdt = a * x(i_) - b * x(i_) * y(i_);
    dydt = c * x(i_) * y(i_) - g * y(i_) - e * y(i_) * z(i_);
    dzdt = -f * z(i_) + g * y(i_) * z(i_);
    %
    x(i) = x(i_) + dxdt * timestep;
    y(i) = y(i_) + dydt * timestep;
    z(i) = z(i_) + dzdt * timestep;
end
figure
plot3(x,y,z)