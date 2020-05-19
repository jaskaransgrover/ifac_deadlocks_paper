% This file does symbolic computation of Phase 1 of two robot deadlock. In
% this phase, collision avoidance constraints of both robots stay inactive. 
clear all
close all
clc
syms Ds t kp1 kp2 gamma x10 y10 x1t y1t Dt alpha Dinit DG1 DG2 kp DG 'real'
global Ds
assume(Dinit,'positive')                        % Initial distance between robots at t=0
assume(DG1,'positive')                          % Distance between robot 1 and its goal at t=0
assume(t,'positive')                            % Time variable
assume(DG2,'positive')                          % Distance between robot 1 and its goal at t=0
assume(Ds,'positive')                           % Desired safety margin
assume(Dt,'positive')                           % Distance between robot 1 and robot 2 at time t
assume(kp1,'positive')                          % Proportional gain of robot 1
assume(kp2,'positive')                          % Proportional gain of robot 2
assume(gamma,'positive')                        % Barrier function parameter

ealpha   = [cos(alpha) ; sin(alpha)];           % Unit vector along the line segment joining the robots, refer to Figure 1 in the paper
p10      = [x10 ; y10];                         % Equation 17, Definiting initial positions
p20      = p10 + Dinit*ealpha;                  % Equation 17, Definiting initial positions
pd1      = p10 + DG1*ealpha ;                   % Equation 18, Defining goals, 
pd2      = p20 - DG2*ealpha;                    % Equation 18, Defining goals

eta_t1   = exp(-kp1*t);                         % Decaying exponential, used in p_1(t)
eta_t2   = exp(-kp2*t);                         % Decaying exponential, used in p_2(t)
p1t      = pd1 - (DG1*eta_t1*ealpha);           % Equation 29, Position of robot 1 at time t where robot 1 dynamics are p1dot = -kp1(p1-pd1)
p2t      = pd2 + (DG2*eta_t2*ealpha);           % Equation 29, Position of robot 2 at time t where robot 2 dynamics are p2dot = -kp2(p2-pd2)
v1t      = simplify(diff(p1t,t));               % velocity of robot 1 at time t
v2t      = simplify(diff(p2t,t));               % velocity of robot 2 at time t

deltap21 = simplify(p2t-p1t);                   % displacement vector along ealpha, this is just a sanity check, displacement along e_{alpha+pi/2} should be zero
deltav21 = simplify(v2t-v1t);                   % relative velocity vector along ealpha
dist_t   = simplify(ealpha'*deltap21);          % Equation 30, Scalar distance between robots

u1cap    = simplify(-kp1*(p1t-pd1));            % Equation 32, Control input  for robot 1
u2cap    = simplify(-kp2*(p1t-pd1));            % Equation 33, Control input  for robot 2

p2t      = p1t + Dt*ealpha;                     % Equation 25, Writing the positions another way to get symbolic equation for roots of flag f_12
a12      = -(p1t-p2t);                          % Equation 25, Writing the positions another way to get symbolic equation for roots of flag f_12
b12      = gamma*h(p1t,p2t)/4;                  % Equation 25, Writing the positions another way to get symbolic equation for roots of flag f_12
f12      = simplify((a12'*u1cap) - b12);        % Equation 25, Writing the positions another way to get symbolic equation for roots of flag f_12
c        = subs(f12,'Dt',0);                    % a,b,c are coefficients of the quadratic polynomial, ax^2 + bx + c, here x = Dt. 
b        = subs(diff(f12,Dt),'Dt',0);           % a,b,c are coefficients of the quadratic polynomial, ax^2 + bx + c, here x = Dt. 
a        = diff(diff(f12,Dt),Dt)/2 ;            % a,b,c are coefficients of the quadratic polynomial, ax^2 + bx + c, here x = Dt. 
p        = [a b c];
polyroots= simplify(roots(p))   ;               
beta1plus= polyroots(1)                         % Equation 34
dist_t                                          % Equation 30, Scalar distance between robots
    

function hval = h(p1,p2)
global Ds
hval = ((p1-p2)'*(p1-p2)) - Ds^2;
end