% This file does symbolic computation of results of Lemma 2, i.e. at t = 0

clear all
close all
clc
syms Ds  kp1 kp2 gamma x10 y10 x1t y1t Dt alpha Dinit DG1 DG2 kp DG 'real'
global Ds
assume(Dinit,'positive')                        % Initial distance between robots at t=0
assume(DG1,'positive')                          % Distance between robot 1 and its goal at t=0
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

u1cap    = -kp1*(p10-pd1);                      % Equation 20, prescribed nominal controller for robot 1
u2cap    = -kp2*(p20-pd2);                      % Equation 20, prescribed nominal controller for robot 2

a12      = -(p10-p20);
b12      = gamma*h(p10,p20)/4;
f        = (a12'*u1cap) - b12;                  % Equation 23, flag for constraint for robot 1
c        = subs(f,'Dinit',0);                   % computing c of the quadratic equation ax^2 + bx + c where x = Dinit
b        = subs(diff(f,Dinit),'Dinit',0);       % computing b of the quadratic equation ax^2 + bx + c where x = Dinit
a        = diff(diff(f,Dinit),Dinit)/2 ;        % computing a of the quadratic equation ax^2 + bx + c where x = Dinit    
p        = [a b c];                             % computing roots of the quadratic equation ax^2 + bx + c where x = Dinit
beta1    = simplify(roots(p));                  % Equation 24
beta1plus= beta1(1);                            % Equation 22

a21      = -(p20-p10);
b21      = gamma*h(p20,p10)/4;
f        = (a21'*u2cap) - b21;                  % flag for constraint for robot 2
c        = subs(f,'Dinit',0);                   % computing c of the quadratic equation ax^2 + bx + c where x = Dinit
b        = subs(diff(f,Dinit),'Dinit',0);       % computing b of the quadratic equation ax^2 + bx + c where x = Dinit
a        = diff(diff(f,Dinit),Dinit)/2 ;        % computing a of the quadratic equation ax^2 + bx + c where x = Dinit
p        = [a b c];                             % computing roots of the quadratic equation ax^2 + bx + c where x = Dinit
beta2    = simplify(roots(p));                  % Equation 24
beta2plus= beta2(1);                            % Equation 22


function hval = h(p1,p2)
global Ds
hval     = ((p1-p2)'*(p1-p2)) - Ds^2;
end