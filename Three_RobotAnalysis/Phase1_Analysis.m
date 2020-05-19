% Analysis of Phase 1 of Three Robot Deadlock. In this Phase, each robot
% has both of its collision avoidance constraints inactive

clear all
close all
clc

global Dt
syms x1t y1t Di Ds Dg kp t gamma Dt x10 y10 alpha 'real'

assume(Di,'positive')                                       % Initial distance between robots at t=0
assume(Dg,'positive')                                       % Distance between robot 1,2,3 and its goal at t=0
assume(Ds,'positive')                                       % Desired safety margin
assume(kp,'positive')                                       % Proportional gain of all robots
assume(gamma,'positive')                                    % Barrier function parameter



ealpha     = [cos(alpha) ; sin(alpha)];
ealpha2pi3 = [cos(alpha + (2*pi/3)) ; sin(alpha + (2*pi/3))];
p10        = [x10 ; y10];
p20        = p10 + Di*ealpha;
p30        = p20 + Di*ealpha2pi3;
c          = (p10+p20+p30)/3 ; 

pd1        = simplify(p10 + Dg*((c-p10)/norm(c-p10)));
pd2        = simplify(p20 + Dg*((c-p20)/norm(c-p20)));
pd3        = simplify(p30 + Dg*((c-p30)/norm(c-p30)));

p1t        = (exp(-kp*t)*p10) + (1-exp(-kp*t))*pd1;
p2t        = (exp(-kp*t)*p20) + (1-exp(-kp*t))*pd2;
p3t        = (exp(-kp*t)*p30) + (1-exp(-kp*t))*pd3;
ct         = (p1t+p2t+p3t)/3 ; 

simplify(c-ct)                                              % Equation 60, centroid remains invariant.


u1cap      = -kp*(p1t-pd1);
u2cap      = -kp*(p2t-pd2);
u3cap      = -kp*(p3t-pd3);

deltap21   = p2t  - p1t ; 
dist_t     = simplify(ealpha'*deltap21);                   % Equation 64
p2t        = p1t  + (Dt*ealpha);

h_12       = simplify(norm(p1t-p2t)^2) - Ds^2;
a12        = -(p1t-p2t);  b12 = (gamma*h_12/4);             % Equation 65
fa         = simplify((a12'*u1cap)-b12);
beta_a     = returnroots(fa);
beta1plus  = beta_a(1)

function beta = returnroots(expr)

global Dt
c      = subs(expr,'Dt',0);
b      = subs(diff(expr,Dt),'Dt',0);
a      = simplify(diff(diff(expr,Dt),Dt)/2) ;
p      = [a b c];
beta   = simplify(roots(p));

end



