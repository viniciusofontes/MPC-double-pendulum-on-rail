clc, clear, close all
set(0,'DefaultFigureWindowStyle','docked')
%% Mechanical Nonlinear Model of a double inverted pendulum on rails
% Professor:  Dr. Helon Vicente Hultmann Ayala (helon@puc-rio.br)

% Source: K. Furuta, T. Okutani, H. Sone, Computer control of a
% double inverted pendulum, Computers & Electrical Engineering,
% Volume 5, Issue 1, 1978, Pages 67-84, ISSN 0045-7906,
% https://doi.org/10.1016/0045-7906(78)90018-6.

% Conventions:
% x (scalar): horizontal position of the cart
% x (array): state vector
% x = [qp1 qp2 qp3 q1 q2 q3]' = [dxdt dtheta1dt dtheta2dt x theta1 theta2]'
% y = [qp1 qp2 qp3 qpp1 qpp2 qpp3] = ...

%% Parameters (in SI)
% l1 = 0.273;
% l2 = 0.245;
% L1 = 0.490;
% L2 = 0.490;
% M1 = 0.123;
% M2 = 0.092;
% J1 = 0.0362 - M1*l1^2- M2*l2^2;
% J2 = 0.00873-M2*l2^2;
% M0 = 0.818-M1-M2;
% F0 = 2.21;
% F1 = 0.00164;
% F2 = 0.00050;
% g = 9.8;

% PROP = [l1 l2 L1 L2 M1 M2 J1 J2 M0 F0 F1 F2 g];
%% Initial conditions

% Horizontal position of the cart (m)
x = 0;
% Angle between the lower pendulum and vertical axis (rad)
theta1 = deg2rad(0);
% Angle between the upper pendulum and vertical axis (rad)
theta2 = deg2rad(45);
% Horizontal velocity of the cart (m)
dxdt = 0;
% Angular velocity of the lower pendulum (rad/s)
dtheta1dt = 0;
% Angular velocity of the upper pendulum (rad/s)
dtheta2dt = 0;
% Initial state vector
x0 = [x theta1 theta2 dxdt dtheta1dt dtheta2dt]';
% Horizontal force applied on the cart (N)
F = 0;
% Initial and final time of the simulation (s)
tspan = [0 5];


%% Processing

[t,x] = ode45(@acc,tspan,x0);
L1 = 0.490;
L2 = 0.490;
% Plor cart
set(gcf,'color','white')

grid on
for i = 1:length(t)
  % Update states
  x_cart = x(i,1);
  x1 = x_cart + L1*sin(x(i,2));
  x2 = x1 + L2*sin(x(i,3));
  y_cart = 0;
  y1 = y_cart + L1*cos(x(i,2));
  y2 = y1 + L2*cos(x(i,3));
  clf
  hold on
  grid on
  axis((L1+L2)*[-1 +1 -1 +1])
  % Plot the cart
  plot(x_cart,y_cart,'ks','markersize',5)
  % Plot the lower pendulum
  plot([x_cart x1],[y_cart y1],'-g','linewidth',2)
  % Plot the upper pendulum
  plot([x1 x2],[y1 y2],'-r','linewidth',2)
  hold off
  legend({'Cart';'Lower pendulum';'Upper pendulum'},...
    'interpreter','latex','location','nw')
  pause(.001)
end