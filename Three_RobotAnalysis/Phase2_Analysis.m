% Analysis of Phase 2 of Three Robot Deadlock. In this Phase, each robot
% has both of its collision avoidance constraints active

clear all
close all
clc

syms x y Dt Di Ds Dg kp gamma x1 y1 alpha t eta 'real'


assume(Di,'positive')                                       % Initial distance between robots at t=0
assume(Dg,'positive')                                       % Distance between robot 1,2,3 and its goal at t=0
assume(Ds,'positive')                                       % Desired safety margin
assume(kp,'positive')                                       % Proportional gain of all robots
assume(gamma,'positive')                                    % Barrier function parameter

p10 = [x;y];                                                % Equation 54, Initial position of robot 1
p20 = p10 + Di*[cos(alpha);sin(alpha)];                     % Equation 54, Initial position of robot 2
p30 = p20 + Di*[cos(alpha +(2*pi/3));sin(alpha+(2*pi/3))];  % Equation 54, Initial position of robot 3
c   = (p10+p20+p30)/3 ; 


pd1 = simplify(p10 + Dg*((c-p10)/norm(c-p10)));             % Equation 55, Initial position of robot 1
pd2 = simplify(p20 + Dg*((c-p20)/norm(c-p20)));             % Equation 55, Initial position of robot 2
pd3 = simplify(p30 + Dg*((c-p30)/norm(c-p30)));             % Equation 55, Initial position of robot 3


p1t = eta*(p10) + (1-eta)*pd1;                              % Equation 67, positions of robot 1 at time t
p2t = p1t + Dt*[cos(alpha) ; sin(alpha)];                   % Equation 67, positions of robot 2 at time t
p3t = p2t + Dt*[cos(alpha+(2*pi/3));sin(alpha + (2*pi/3))]; % Equation 67, positions of robot 3 at time t


u1cap_t = -kp*(p1t-pd1);                
u2cap_t = -kp*(p2t-pd2);
u3cap_t = -kp*(p3t-pd3);


h_12   = simplify(norm(p1t-p2t)^2) - Ds^2;
h_13   = simplify(norm(p1t-p3t)^2) - Ds^2;
h_23   = simplify(norm(p2t-p3t)^2) - Ds^2;


a12    = -(p1t-p2t);  b12 = (gamma*h_12/4);
a13    = -(p1t-p3t);  b13 = (gamma*h_13/4);
a23    = -(p2t-p3t);  b23 = (gamma*h_23/4);
a21    = -a12            ;  b21 = b12 ;
a31    = -a13            ;  b31 = b13 ; 
a32    = -a23            ;  b32 = b23 ;


% Robot 1 
G      = [a12'*a12 a12'*a13 ; a13'*a12 a13'*a13];           % Grammian for computing Lagrange multipliers for robot 1's constraints with 2 and 3                
f      = 2*[(a12'*u1cap_t) - b12 ; (a13'*u1cap_t) - b13];
mu     = simplify(G\f);                                     % Lagrange multipliers for robot 1's constraints with 2 and 3  
u1star = simplify(u1cap_t - 0.5*([a12 a13]*mu));            % Equation 66, both collision avoidance constraints of robot 1 are active


% Robot 2 
G      = [a21'*a21 a21'*a23 ; a23'*a21 a23'*a23];           % Grammian for computing Lagrange multipliers for robot 2's constraints with 1 and 3  
f      = 2*[(a21'*u2cap_t) - b21 ; (a23'*u2cap_t) - b23];
mu     = simplify(G\f);                                     % Lagrange multipliers for robot 2's constraints with 1 and 3       
u2star = simplify(u2cap_t - 0.5*([a21 a23]*mu));            % Equation 66, both collision avoidance constraints of robot 2 are active


% Robot 3 
G      = [a31'*a31 a31'*a32 ; a32'*a31 a32'*a32];           % Grammian for computing Lagrange multipliers for robot 3's constraints with 1 and 2  
f      = 2*[(a31'*u3cap_t) - b31 ; (a32'*u3cap_t) - b32];
mu     = simplify(G\f);                                     % Lagrange multipliers for robot 3's constraints with 1 and 2  
u3star = simplify(u3cap_t - 0.5*([a31 a32]*mu));            % Equation 66, both collision avoidance constraints of robot 3 are active



% Relative Dynamics
delp21 = p2t-p1t;
delv21 = simplify(u2star-u1star);
Dtdot  = simplify([cos(alpha) sin(alpha)]*delv21)           % Equation 70, Dt dot    

 