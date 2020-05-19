% Simulation to generate plots in Figure 2

clear all
close all
clc

global Dinit Ds DG1 DG2 kp1 kp2  gamma pd1 pd2 L12 L21
Dinit    = 7    ;                       % Initial separation between robots
Ds       = 1    ;                       % Allowed Safety Margin
kp1      = .5   ;                       % Proportional Gain of Robot 1
kp2      = .25  ;                       % Proportional Gain of Robot 2
gamma    = 5    ;                       % Barrier Function Parameter
DG1      = 10   ;                       % Initial distance of robot 1 from its goal
DG2      = 8    ;                       % Initial distance of robot 2 from its goal

alpha    = pi/4 ;                       % Line segment angle alpha, refer to Figure 1
e1       = [cos(alpha) ; sin(alpha)];   
eperp    = [-sin(alpha) ; cos(alpha)];

p10      = [1; 0];                      % initial position of robot 1
p20      = p10 + Dinit*e1;              % initial position of robot 2
pd1      = p10 + DG1*e1 ;               % goal position of robot 1
pd2      = p20 - DG2*e1 ;               % goal position of robot 2

beta1    = (2*DG1*kp1 + (4*DG1^2*kp1^2 + Ds^2*gamma^2)^(1/2))/gamma ; % Equation 22, initial critical distance robot 1
beta2    = (2*DG2*kp2 + (4*DG2^2*kp2^2 + Ds^2*gamma^2)^(1/2))/gamma ; % Equation 22, initial critical distance robot 2


disp('if Dinit is greater than beta1, then both multipliers are off because beta1 is greater than beta2')


tspan = 0:0.003:3;
beta1_tminus   = zeros(1,length(tspan));
beta1_tplus    = zeros(1,length(tspan));
beta2_tminus   = zeros(1,length(tspan));
beta2_tplus    = zeros(1,length(tspan));
Dt             = zeros(1,length(tspan));

% The loop below precomputes the critical distances beta1plus, beta2plus.
for i = 1 : length(tspan)
    t                 = tspan(i);
    eta1              = exp(-kp1*t);
    eta2              = exp(-kp2*t);
    beta1_tplus(i)    = (2*DG1*kp1*eta1 + (4*DG1^2*kp1^2*eta1^2 + Ds^2*gamma^2)^(1/2))/gamma;
    beta2_tplus(i)    = (2*DG2*kp2*eta2 + (4*DG2^2*kp2^2*eta2^2 + Ds^2*gamma^2)^(1/2))/gamma;
    Dt(i)             = -((DG1*(1-eta1)) + (DG2*(1-eta2)) -Dinit);
end


% lines 51-94 are for robot animation
close all
figure('units','normalized','outerposition',[0 0 1 1],'color','white')
plot(p10(1),p10(2),'*g'); hold on ; viscircles(p10',Ds/2); text(p10(1),p10(2),'$p_{1}$','interpreter','latex','fontsize',28)
plot(p20(1),p20(2),'*g'); viscircles(p20',Ds/2); text(p20(1),p20(2),'$p_{2}$','interpreter','latex','fontsize',28)
grid on
hold on
plot(pd1(1),pd1(2),'*r'); text(pd1(1),pd1(2),'$p_{d_1}$','interpreter','latex','fontsize',28)
plot(pd2(1),pd2(2),'*r'); text(pd2(1),pd2(2),'$p_{d_2}$','interpreter','latex','fontsize',28)
axis equal
xlim([-3 10])
ylim([-3 10])
zold  = [p10;p20];
hold off
dt       = tspan(2)-tspan(1);
distance = zeros(1,length(tspan));
Dtnew    = zeros(1,length(tspan));

for i = 1 : length(tspan)
    tspan(i)
    [zdot,lamb1,lamb2]       = cbf_control(zold);       % computes controls for robots 1 and 2 using control barrier functions
    L12(i)                   = lamb1;                   % Lagrange multiplier for robot 1
    L21(i)                   = lamb2;                   % Lagrange multiplier for robot 2
    znew                     = zold + (dt*zdot);        % Update robots' position by integrating velocity
    p1t                      = zold(1:2);               % Record old position of robot 1 for animation
    p2t                      = zold(3:4);               % Record old position of robot 2 for animation
    Dtnew(i)                 = (p2t-p1t)'*e1;           % Distance between robot 2 and robot 1
    zold                     = znew ;                   % Update
    
    plot(p1t(1),p1t(2),'*g'); hold on ; viscircles(p1t',Ds/2); text(p1t(1),p1t(2),'$p_{1}$','interpreter','latex','fontsize',28); hold on
    plot(p2t(1),p2t(2),'*g'); viscircles(p2t',Ds/2); text(p2t(1),p2t(2),'$p_{2}$','interpreter','latex','fontsize',28)
    plot(pd1(1),pd1(2),'*r'); text(pd1(1),pd1(2),'$p_{d_1}$','interpreter','latex','fontsize',28)
    plot(pd2(1),pd2(2),'*r'); text(pd2(1),pd2(2),'$p_{d_2}$','interpreter','latex','fontsize',28)
    line([pd1(1) pd2(1)],[pd1(2) pd2(2)])
    xlim([-3 10])
    ylim([-3 10])
    axis equal
    grid on
    set(gca,'fontsize',28)
    drawnow
    hold off
    clc
    
end

% lines 97 - 107 compute times t1, t2 according to equation 27 and 36 in
% the paper
tol     = 0.006; 
index1  = find(abs(beta1_tplus-Dtnew) < tol);
t_1     = tspan(index1(1));
yval1   = Dtnew(index1(1));
Y1      = -1:0.01:yval1;

index2  = find(abs(beta2_tplus-Dtnew) < tol);
t_2     = tspan(index2(1));
yval2   = Dtnew(index2(1));
Y2      = -1:0.01:yval2;

% Everything below this is graphics and plotting
close all
plot_graphs(t_1,t_2,Dtnew,tspan,beta1_tplus,beta2_tplus,Y1,Y2)





function [u,lamb1,lamb2] = cbf_control(z)

global kp1 kp2 pd1 pd2  gamma  D 

options                = optimoptions('quadprog','Display','Off');
p1                     = z(1:2);
p2                     = z(3:4);

u1cap                  = -kp1*(p1-pd1);
u2cap                  = -kp2*(p2-pd2);


A1                     = -(p1-p2)';
b1                     = 0.25*gamma*h_function(p1,p2);
[u1star,~,~,~,lambda1] = quadprog(2*eye(2,2),-2*u1cap,A1,b1,[],[],[],[],[],options); 
lamb1                  = lambda1.ineqlin(1);


A2                     = -(p2-p1)';
b2                     = 0.25*gamma*h_function(p2,p1);
[u2star,~,~,~,lambda2] = quadprog(2*eye(2,2),-2*u2cap,A2,b2,[],[],[],[],[],options); 
lamb2                  = lambda2.ineqlin(1);

u                      = [u1star;u2star];
D                      = [D,norm(p1-p2)];
end

function hval = h_function(pa,pb)

global Ds
hval = (norm(pa-pb))^2 - (Ds^2);

end

function plot_graphs(t_1,t_2,Dtnew,tspan,beta1_tplus,beta2_tplus,Y1,Y2)


global Dinit Ds DG1 DG2 kp1 kp2  gamma L12 L21

f4      = figure('units','normalized','outerposition',[0 0 1 1],'color','white');
K       = length(L12);

% Phase 1 patch
x       = [0 t_1 t_1 0];
y       = [0 0 7 7];
p       = patch(x,y,'r');
set(p,'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off'); hold on

% Phase 2 patch
x       = [t_1 t_2 t_2 t_1];
y       = [0 0 7 7];
p       = patch(x,y,'g');
set(p,'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');

% Phase 3 patch
x       = [t_2 3 3 t_2];
y       = [0 0 7 7];
p       = patch(x,y,'b');
set(p,'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');



plot(tspan(1:K),Ds*ones(1,K),'-g','linewidth',9);
plot(tspan(1:K),Dtnew(1:K),'color',[0,0.5,0],'linewidth',9)
hold on
plot(tspan(1:K),L12,'-m','linewidth',9);
plot(tspan(1:K),beta1_tplus(1:K),'-b','linewidth',9);
X1     = ones(1,length(Y1))*t_1;
plot(X1,Y1,'--k','linewidth',9,'HandleVisibility','off');
plot(X1(end),Y1(end),'*c','linewidth',14,'HandleVisibility','off');


plot(tspan(1:K),L21,'-r','linewidth',9);
plot(tspan(1:K),beta2_tplus(1:K),'color',[0.9290, 0.6940, 0.1250],'linewidth',9);
X2     = ones(1,length(Y2))*t_2;
plot(X2,Y2,'--k','linewidth',9,'HandleVisibility','off');
plot(X2(end),Y2(end),'*c','linewidth',14,'HandleVisibility','off');


h = legend('$D_s$','$D(t)$','$\mu_{12}(t)$','$\beta^1_{+}(t)$','$\mu_{21}(t)$','$\beta^2_{+}(t)$');
xlabel('Time [s]','interpreter','latex');
grid on
xlim([0 3])
dim = [0.8    0.57    0.0836    0.3328];
set(h,'interpreter','latex','fontsize',48,'Position',dim)

text(t_1-0.08,0.3,'$\boldmath{t_{1}}$','interpreter','latex','fontsize',48)
text(t_2-0.08,0.3,'$\boldmath{t_{2}}$','interpreter','latex','fontsize',48)
plot(t_1,0,'*c','linewidth',14,'HandleVisibility','off')
plot(t_2,0,'*c','linewidth',14,'HandleVisibility','off')
ylim([-0.01 7])



text(0.1,6.5,'Phase 1','interpreter','latex','fontsize',40)
text(t_1+0.00,6.5,'Phase 2','interpreter','latex','fontsize',40)
text(t_2+.5,6.5,'Phase 3','interpreter','latex','fontsize',40)
xx = [0:0.1:t_1];
yy = zeros(1,length(xx));
plot(xx,yy,'--m','linewidth',9,'HandleVisibility','off');
grid  on
set(gca,'fontsize',48)
ax = gca;
ax.GridAlpha = .2;



dim = [.79 .45 .105 .05];
stra = '$D_{G_1}=';
strb = num2str(DG1);
strc = '$';
str1 = strcat(stra,strb,strc);
annotation('textbox',dim,'String',str1,'FitBoxToText','off','interpreter','latex','fontsize',40,'EdgeColor','none','BackgroundColor','white');


dim = [.79 .40 .105 .05];
stra = '$D_{G_2}=';
strb = num2str(DG2);
strc = '$';
str1 = strcat(stra,strb,strc);
annotation('textbox',dim,'String',str1,'FitBoxToText','off','interpreter','latex','fontsize',40,'EdgeColor','none','BackgroundColor','white');


dim = [.79 .35 .105 .05];
stra = '$D_{init.}\hspace{-0.1cm}=';
strb = num2str(Dinit);
strc = '$';
str1 = strcat(stra,strb,strc);
annotation('textbox',dim,'String',str1,'FitBoxToText','off','interpreter','latex','fontsize',40,'EdgeColor','none','BackgroundColor','white');


dim = [.79 .30 .105 .05];
stra = '$k_{p_1}\hspace{0.18cm}=';
strb = num2str(kp1);
strc = '$';
str1 = strcat(stra,strb,strc);
annotation('textbox',dim,'String',str1,'FitBoxToText','off','interpreter','latex','fontsize',40,'EdgeColor','none','BackgroundColor','white');


dim = [.79 .25 .105 .05];
stra = '$k_{p_2}\hspace{0.18cm}=';
strb = num2str(kp2);
strc = '$';
str1 = strcat(stra,strb,strc);
annotation('textbox',dim,'String',str1,'FitBoxToText','off','interpreter','latex','fontsize',40,'EdgeColor','none','BackgroundColor','white');


dim = [.79 .20 .105 .05];
stra = '$\hspace{0.4cm}\gamma=';
strb = num2str(gamma);
strc = '$';
str1 = strcat(stra,strb,strc);
annotation('textbox',dim,'String',str1,'FitBoxToText','off','interpreter','latex','fontsize',40,'EdgeColor','none','BackgroundColor','white');

end


