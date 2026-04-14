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
%addpath('/Finite burn')
load_subroutines


% Circular orbits from Vallado pg 326 
rA = 191.34411; % km Initial orbit 
rB = 35781.34857; % km Final Orbit
ER = 6378.137;

r_initial = (rA+ER);
r_final = (rB+ER);

v1 = sqrt(mu/r_initial);
v1 = [v1 0 0]';

y_A = [r_initial, 0, 0]';
y_B = [r_final, 0, 0]';

t = 315.402974 *60; % time for the Hohmann transfer
%N = 0;
%D = 0;
%[vi,vf,vH1,vH2,theta] = LambertBattin(y_A,y_B,t,N,D,v1,mu);
tol = 1e-5;
max_iter = 1000;
[vi, vf] = P_method_Lambert(y_A,y_B,t, max_iter, tol);


fprintf("The Initial Lambert velocity in km/s is %f\n", vi);
fprintf("The Final Lambert velocity in km/s  is %f\n", vf);

TwoD_Hohmann_transfer(rA, rB, norm(vi), norm(vf), params)