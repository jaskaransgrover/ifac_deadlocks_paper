% Analysis of Phase 0 of Three Robot Deadlock. This is the proof of Lemma 7
% in the paper.

clear all
close all
clc

syms Di Ds Dg kp gamma x y alpha 'real'


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


u1cap0 = -kp*(p10-pd1);                                     % Prescribed nominal control for robot 1
u2cap0 = -kp*(p20-pd2);                                     % Prescribed nominal control for robot 2
u3cap0 = -kp*(p30-pd3);                                     % Prescribed nominal control for robot 3


h_12   = norm(p10-p20)^2 - Ds^2;    
h_13   = norm(p10-p30)^2 - Ds^2;
h_23   = norm(p20-p30)^2 - Ds^2;


a12    = -(p10-p20);  b12 = (gamma*h_12/4);
a13    = -(p10-p30);  b13 = (gamma*h_13/4);
a23    = -(p20-p30);  b23 = (gamma*h_23/4);
a21    = -a12      ;  b21 = b12 ;
a31    = -a13      ;  b31 = b13 ; 
a32    = -a23      ;  b32 = b23 ; 


f12    = simplify((a12'*u1cap0) - b12);                     % Equation 58, f12 for robot 1, robot 1's constraint with robot 2
f13    = simplify((a13'*u1cap0) - b13);                     % Equation 58, f13 for robot 1, robot 1's constraint with robot 3

f21    = simplify((a21'*u2cap0) - b21);                     % f21 for robot 2, robot 2's constraint with robot 1
f23    = simplify((a23'*u2cap0) - b23);                     % f23 for robot 2, robot 2's constraint with robot 3

f31    = simplify((a31'*u3cap0) - b31);                     % f31 for robot 3, robot 3's constraint with robot 1           
f32    = simplify((a32'*u3cap0) - b32);                     % f32 for robot 3, robot 3's constraint with robot 2     


F      = [f12 f13 f21 f23 f31 f32];
expr   = F(1) ; 
c      = subs(expr,'Di',0);                                 % c in ax^2 + bx + c = 0, x = Dinit     
b      = subs(diff(expr,Di),'Di',0);                        % b in ax^2 + bx + c = 0, x = Dinit     
a      = diff(diff(expr,Di),Di)/2 ;                         % a in ax^2 + bx + c = 0, x = Dinit     
p      = [a b c];                                           % calculating roots of ax^2 + bx + c = 0 in Equation 58

beta  = roots(p);
beta1 = beta(1)                                             % Equation 57, Initial critical distance
