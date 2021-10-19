% this script computes steady-state solution of linear oscillator subject
% to random excitation

function linear_oscillator_random_noise
clear;
close all;
clc;
global A_g tspan_g zeta_g wn_g omega_g phi_g distlag

tspan_g=linspace(0,100,1000);
wn_g=1;
zeta_g=0.05;
distlag=0;
omega_g=3.5;
phi_g=0;

A_g=5*rand(size(tspan_g));

phi=atan(2*zeta_g*omega_g*wn_g/(wn_g^2-omega_g^2));
x0=A_g(1)*cos(omega_g*tspan_g(1)-phi);
xd0=-A_g(1)*omega_g*sin(omega_g*tspan_g(1)-phi);
xv0=[x0,xd0];
[t,x]=ode45(@randdx,tspan_g,xv0);

figure(1)
plot(t,x(:,1))

end
% 
% function [t,x,A] = randamp(tspan,wn,zeta,distlag,omega)
%     % distlag = 0 for uniform and 1 for Gaussian distribution
%     global A_g tspan_g zeta_g wn_g omega_g phi_g
%     % these global parameters are needed by ODE definition
%     if nargin<1 tspan=linspace(0,30,1000); end % from 0-30 with 1000 steps
%     if nargin<2 wn=1; end
%     if nargin<3 zeta=0.05; end
%     if nargin<4 distlag=0; end
%     if nargin<5 omega=3.5; end
%     
%     % assign global values
%     tspan_g=tspan; 
%     zeta_g=zeta;
%     wn_g=wn; 
%     omega_g=omega; 
%     phi_g=0; % zero phase angle
% 
% %     m=max(size(tspan));
% %     if m<1000
% %         tspan=linspace(t(1),t(m),1000);
% %     end
% 
%     switch distlag
%         case 0
%             A=5*rand(size(tspan)); % A now has a range from 0 to 5.
%         case 1
%             A=2.5+randn(size(tspan)); % One standard deviation will span [1.5,3.5].
%             A=A.*sign(A); % Multiplying A by the signum function of itself.
%         otherwise
%             error('Distribution lag must be zero or one.');
%     end
%     A_g=A;
%     
%     % next: define initial condition to remove transients
%     phi=atan(2*zeta*omega*wn/(wn^2-omega^2));
%     x0=A(1)*cos(omega*tspan(1)-phi);
%     xd0=-A(1)*omega*sin(omega*tspan(1)-phi);
%     xv0=[x0,xd0];
% 
%     
% 
%     [t,x]=ode45('randdx',tspan,xv0);
% 
% 
% end
% 
% 

function dx = randdx(t,x)
    global A_g tspan_g zeta_g wn_g omega_g phi_g
    At=interp1(tspan_g,A_g,t);
    dx = [x(2); At*cos(omega_g*t-phi_g)-(2*zeta_g*wn_g*x(2)+wn_g^2*x(1))];
%     dx(1,1)=x(2,1);
%     dx(2,1)=At*cos(omega*t-phi)-(2*zeta*wn*x(2,1)+wn^2*x(1,1));
end


