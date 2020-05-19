% This file does symbolic computation of Phase 2 of two robot deadlock. In
% this phase, collision avoidance constraints of both robot 1 and 2 are
% active
clear all
close all
clc

syms x1t y1t Ds t kp1 kp2 gamma x10 y10 x1t y1t Dt alpha Dinit DG1 DG2 kp DG 'real'
global Ds
assume(Dinit,'positive')                           % Initial distance between robots at t=0
assume(DG1,'positive')                             % Distance between robot 1 and its goal at t=0
assume(t,'positive')                               % Time variable
assume(DG2,'positive')                             % Distance between robot 1 and its goal at t=0
assume(Ds,'positive')                              % Desired safety margin
assume(Dt,'positive')                              % Distance between robot 1 and robot 2 at time t
assume(kp1,'positive')                             % Proportional gain of robot 1
assume(kp2,'positive')                             % Proportional gain of robot 2
assume(gamma,'positive')                           % Barrier function parameter

ealpha      = [cos(alpha) ; sin(alpha)];           % Unit vector along the line segment joining the robots, refer to Figure 1 in the paper
eperp       = [cos(alpha+pi/2) ; sin(alpha+pi/2)]; % Needed to verify u1perp=0 in Equation 40
p10         = [x10 ; y10];                         % Equation 17, Defining initial positions
p20         = p10 + Dinit*ealpha;                  % Equation 17, Defining initial positions
pd1         = p10 + DG1*ealpha ;                   % Equation 18, Defining goals, 
pd2         = p20 - DG2*ealpha;                    % Equation 18, Defining goals
p1t         = [x1t;y1t];                           % p1t = (x1t,y1t) = position of robot 1 at time t, it's a floating variable.


p2t         = p1t + Dt*ealpha;                     % Robot 1 position is written in terms of robot 2
u1cap       = -kp1*(p1t-pd1);                      % Nominal control of robot 1
u2cap       = -kp2*(p2t-pd2);                      % Nominal control of robot 2


a12         = -(p1t-p2t);
b12         = gamma*h(p1t,p2t)/4;
mu_12       = 2*((a12'*u1cap) - b12)/(a12'*a12);   % Equation 49, Lagrange Multiplier
u1star      = simplify(u1cap - (0.5*mu_12*a12));   % Equation 48, Robot 1 control
v1t         = u1star;                              % Equation 50, Robot 1 velocity = control

a21         = -(p2t-p1t);
b21         = gamma*h(p2t,p1t)/4;
mu_21       = 2*((a21'*u2cap) - b21)/(a21'*a21);   % Equation 49, Lagrange Multiplier
u2star      = simplify(u2cap - (0.5*mu_21*a21));   % Equation 48, Robot 2 control
v2t         = u2star;                              % Equation 50, Robot 1 velocity = control

deltap21    = simplify(p2t-p1t);
deltav21    = simplify(v2t-v1t);

Dtdot       = simplify(ealpha'*deltav21)          ;% Equation 52, rate at which distance between robots is decreasing 
dist        = simplify(sqrt(deltap21'*deltap21))  ;% This should simply be Dt


function hval = h(p1,p2)

global Ds
hval = ((p1-p2)'*(p1-p2)) - Ds^2;

end