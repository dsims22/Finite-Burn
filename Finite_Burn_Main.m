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
addpath('/home/arsl/Desktop/Research/Finite burn')
load_subroutines

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
g = 9.81;
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
[t_burn1,t_burn2,t_burn3, xout1, xout2, xout3 ] = bcb_sim(r1, v1,m0, delta_m1, delta_m3,P,params);

% Simulate orbits
tspan = [0 3600];
tspan2 = [0 100000];
OE1 = RV2COE(r1,v1,params);
OE2 = RV2COE(r2,v2,params);
n = 1000; % Number of points for the plot
Data_Points = linspace(0, 2*pi, n);
[i1,j1,k1] = generate_orbit_points(OE1, Data_Points);
[i2, j2, k2] = generate_orbit_points(OE2, Data_Points);


figure();
hold on;
plot3(i1, j1, k1, 'DisplayName', 'Initial Orbit');
plot3(i2, j2, k2, 'DisplayName', 'Final Orbit');

plot3(xout1(:,1), xout1(:,2), xout1(:,3), 'DisplayName', 'Burn 1');
plot3(xout2(:,1), xout2(:,2), xout2(:,3), 'DisplayName', 'Coast');
plot3(xout3(:,1), xout3(:,2), xout3(:,3), 'DisplayName', 'Burn 2');

scatter3(xout1(1,1),   xout1(1,2),   xout1(1,3),   30, 'k', 'filled', 'DisplayName', 'Start');
scatter3(xout3(end,1), xout3(end,2), xout3(end,3), 30, 'r', 'filled', 'DisplayName', 'End');
scatter3(0, 0, 0,   80, 'b', 'filled', 'DisplayName', 'Earth');

xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
legend; grid on; axis equal;
view(3)
hold off;

figure();
subplot(3,1,1);
hold on;
plot(t_burn1, xout1(:,4)); 
plot(t_burn2, xout2(:,4)); 
plot(t_burn3, xout3(:,4));
ylabel('v_x (m/s)');
grid on;

subplot(3,1,2);
hold on;
plot(t_burn1, xout1(:,5)); 
plot(t_burn2, xout2(:,5)); 
plot(t_burn3, xout3(:,5));
ylabel('v_y (m/s)'); 
grid on;

subplot(3,1,3);
hold on;
plot(t_burn1, xout1(:,6)); 
plot(t_burn2, xout2(:,6)); 
plot(t_burn3, xout3(:,6));
xlabel('Time (s)');
ylabel('v_z (m/s)');  
grid on;

legend('Burn 1', 'Coast', 'Burn 2');
hold off;


%%%%%%%%%%%% Fmincon %%%%%%%%%%%

% Initial and final conditions
%x0 = [r1; v1];
des_val = [r2;v2];

% Linear inequality constraints
A = [-1 0 0 0 0 0 0 0 0 ; % P_1, P_5, P_6 A matrix leq
     0 0 0 0 -1 0 0 0 0;
     0 0 0 0 0 -1 0 0 0];

b = [0;0;0];

% Linear equality Constraints
Aeq = [];
beq = [];

% Upper and lower bounds 
ub = [];
lb = [];
%x_feas = linprog(zeros(size(P)), A, b, Aeq, beq, lb, ub);

options1 = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', ...
    'interior-point', 'MaxFunctionEvaluations', 500000, ...
    'MaxIterations', 50000,PlotFcn={@optimplotx,@optimplotfval,@optimplotfirstorderopt});%,    Algorithm="interior-point",...
    %EnableFeasibilityMode=true,...
    %SubproblemAlgorithm="cg");%'SpecifyObjectiveGradient',true); 
x01 = [r1; v1];

x_fmincon1 = fmincon(@(P)J(P, m0, delta_m3), P, A, b, Aeq, beq, lb, ub, @(x)mycon(x,P,x01,m0, options, des_val,delta_m1,delta_m3, params), options1);

% First Burn
%x_fmincon1 = fmincon(@(P)J(P, Ve, m0, mdot), P, A, b, Aeq, beq, lb, ub, @(x)mycon(x,P(2:4),x01,m0, booleanT, tspan1, options, des_val, params), options1);

% x_fmincon1 = fmincon(@(x)norm(sys_dynamics(x,P(2:4),m0,booleanT,params)), x01, A, b, Aeq, beq, lb, ub, @(x)mycon(x,P(2:4),x01,m0, booleanT, tspan1, options, des_val, params));

% Coast
%x_fmincon2 = fmincon(@(P)J(P, Ve, delta_m1, mdot), P, A, b, Aeq, beq, lb, ub, @(x)mycon(x,P(2:4),x02,delta_m1, booleanF, tspan2, options, des_val,params), options1);

% x_fmincon2 = fmincon(@(x)norm(sys_dynamics(x,P(2:4),delta_m1,booleanF,params)), x02, A, b, Aeq, beq, lb, ub, @(x)mycon(x,P(2:4),x02,delta_m1, booleanF, tspan2, options, des_val,params));

% Second Burn
%[t_f3,x_fmincon3] = fmincon(@(P)J(P, Ve, delta_m3, mdot), P, A, b, Aeq, beq, lb, ub, @(x)mycon(x,P(7:9),x03,delta_m3, booleanT, tspan3, options, des_val,params), options1);

%[t_f3,x_fmincon3] = fmincon(@(x)norm(sys_dynamics(x,P(7:9),delta_m3,booleanT,params)), x03, A, b, Aeq, beq, lb, ub, @(x)mycon(x,P(7:9),x03,delta_m3, booleanT, tspan3, options, des_val,params));
%%%%%%%% Function to be minimized %%%%%%%%

function [costf] = J(P, m0, delta_m3)
    
    costf = P(1)*norm(P(2:4))/m0 + P(6)*norm(P(7:9))/delta_m3;
    %if nargout>1
    %    costg = zeros(size(P));
    %end

end

%function costf = J(P, Ve, m0, mdot)
%    t1 = P(1);
%    t2 = P(6);
%    
%    mf1  = m0   - mdot * t1;
%    mf2  = mf1  - mdot * t2;
%    
%    dv1 = Ve * log(m0  / mf1);   
%    dv2 = Ve * log(mf1 / mf2);   
%    
%    costf = dv1 + dv2;            
%end

%%%%%%%%% Non-linear constraints function for fmincon%%%%%%%%%%%
%function [c, ceq] = mycon(x,T,x0,m0, boolean, tspan, options, des_val,params)
%       % Extract the values
%       x_des = des_val(1:3);
%       v_des = des_val(4:6);
%                
%       [t, xout] = ode45(@(t,x)sys_dynamics(x,T,m0,boolean,params), tspan, x0, options);
%
%       xf = xout(end, 1:3);
%       vf = xout(end, 4:6);
%
%       % Errors in position and velocity
%       e_x = xf - x_des;
%       e_v = vf - v_des;
%        
%       % Nonlinear inequalities at x
%       c = [e_x; e_v];
%       ceq = [];
%end
function [c, ceq] = mycon(x,P,x0,m0, options, des_val,delta_m1,delta_m3, params)
       % Extract the values
       x_des = des_val(1:3);
       v_des = des_val(4:6);
                
        [t_burn1,t_burn2,t_burn3, xout1, xout2, xout3 ] = bcb_sim(r1, v1,m0, delta_m1, delta_m3,P,params);

       xf = xout3(end, 1:3);
       vf = xout3(end, 4:6);

       % Errors in position and velocity
       e_x = xf - x_des;
       e_v = vf - v_des;
        
       % Nonlinear inequalities at x
       c = [e_x; e_v];
       ceq = [];
end










