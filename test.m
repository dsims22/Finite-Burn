clear all;
clc;

% Define Parameters
params.mu = 398600; % gravitational parameter, km^3/s^2 
params.mdot = 0.025; % kg/s (made up a number that sounds reasonable for electric propulsion. Will add a real number in the future)  
params.g = 9.81; % gravity m/s^2 (probably need to use gravity at altitude but this is fine for a first look)
mu = params.mu;
g = params.g;
mdot = params.mdot;

% Call functions
addpath("orbits")

% Initial Orbit
r1 = [-10063; -472; -12487]; % km
v1 = [-0.359; -4.950; 0.475]; % km/s
y01 = [r1; v1];

% Final Orbit
r2 = [-7704; -7223; -12077]; % km
v2 = [-2.4078; -4.249; -1.0066]; % km/s
y02 = [r2; v2];

% Spacecraft mass and ISP
Isp = 1960; % Sec using Busek USA	BIT-3	RF Ion	Iodine	https://www.nasa.gov/smallsat-institute/sst-soa/in-space_propulsion/
m0 = 50; % kg (made up, also very unrealistic for a satellite)
params.m0=m0;

% Ode 45 options
options = odeset(RelTol=1e-6, AbsTol=1e-8);

% Lambert Problem Velocities
t_goal= 0.45 * 60 * 60; % 12 hours for the transfer
N = 0;
D = 0;
[vi,vf,vH1,vH2,theta] = LambertBattin(r1,r2,t_goal,N,D,v1,mu);

% Calculate times of the burns
Ve = (Isp * g) / 1000; % km/s         

% Burn 1
delta_v1 = norm(vi - v1);       
mf1 = m0 * exp(-delta_v1 / Ve); 
delta_m1  = m0 - mf1;                       
t1   = delta_m1 / mdot;                      

% Burn 2
delta_v3 = norm(vf - v2);        
mf3 = delta_m1 * exp(-delta_v3 / Ve);
delta_m3  = delta_m1 - mf3;
t3   = delta_m3 / mdot;

% Coast
t2 = t_goal - t1 - t3;

fprintf('Time of first burn in hours %f\n', t1/3600)
fprintf('Time of coast in hours %f\n', t2/3600)
fprintf('Time of third burn in hours %f\n', t3/3600)

 
% Defining the P vector (also making up a number for thrust and times,
% thrust = 1 mN
T =   mdot*Ve;% 1e-6; % Changed to try to get it to converge
P = [t1 T T T t2 t3 T T T]';

% Simulate First Burn
tspan1 = [0 P(1)];
x01 = [r1; v1];
x = [r1; v1];
booleanT = "True";
[t_burn1, xout1] = ode45(@(t,x)sys_dynamics(x,P(2:4),m0,booleanT,params), tspan1, x01, options);

% Simulate Coast
tspan2 = [P(1) P(5)+P(1)];
x02 = xout1(end,:)';
booleanF = "False";
[t_burn2, xout2] = ode45(@(t,x)sys_dynamics(x,P(2:4),delta_m1,booleanF,params), tspan2, x02, options);

% Simulate Second Burn
tspan3 = [P(5)+P(1) P(5)+P(1)+P(6)];
x03 = xout2(end,:)';
[t_burn3, xout3] = ode45(@(t,x)sys_dynamics(x,P(7:9),delta_m3,booleanT,params), tspan3, x03, options);

% Simulate orbits
tspan = [0 3600];
tspan2 = [0 100000];
OE1 = orbits.RV2COE(r1,v1,params);
OE2 = orbits.RV2COE(r2,v2,params);
n = 1000; % Number of points for the plot
Data_Points = linspace(0, 2*pi, n);
[i1,j1,k1] = generate_orbit_points(OE1, Data_Points);
[i2, j2, k2] = generate_orbit_points(OE2, Data_Points);
