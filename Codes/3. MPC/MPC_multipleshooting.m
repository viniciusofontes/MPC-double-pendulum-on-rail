% first casadi test for mpc fpr mobile robots
clear all
close all
clc

% CasADi v3.4.5
% addpath('C:\Users\mehre\OneDrive\Desktop\CasADi\casadi-windows-matlabR2016a-v3.4.5')
% CasADi v3.5.5
addpath('D:\Casadi\casadi-windows-matlabR2016a-v3.5.5')

import casadi.*

%% Parâmetros do problema
l1 = 0.273;
l2 = 0.245;
L1 = 0.490;
L2 = 0.490;
M1 = 0.123;
M2 = 0.092;
J1 = 0.0362 - M1*l1^2- M2*l2^2;
J2 = 0.00873-M2*l2^2;
M0 = 0.818-M1-M2;
F0 = 2.21;
F1 = 0.00164;
F2 = 0.00050;
g = 9.8;

T = 0.02; %[s]
N = 19; % prediction horizon
% rob_diam = 0.3;

% Weight matrix for the associated with the state vector
Q = eye(6);
% Position and speed of the cart, respectively
Q(1,1) = .005; Q(4,4) = 1e-4;
% Angular positions of the lower and upper pendulums, respectively
Q(2,2) = 2; Q(3,3) = 1;
% Angular velocities of the lower and upper pendulums, respectively
Q(5,5) = 0*1e-4; Q(6,6) = 0*1e-4;
% Weight associated with the control vector
R = 1*1e-9; % Only 1 component since the only control variable is the force F
S = 1e-8;

%% Valores máximo e mínimo de F
F_max = 1e3; F_min = -1e3;

% v_max = 0.6; v_min = -v_max;
% omega_max = pi/4; omega_min = -omega_max;

%% Variáveis de estado
x = SX.sym('x'); dx = SX.sym('dx');
theta1 = SX.sym('theta1'); dtheta1 = SX.sym('dtheta1');
theta2 = SX.sym('theta2'); dtheta2 = SX.sym('dtheta2');
states = [x;theta1;theta2;dx;dtheta1;dtheta2]; n_states = length(states);

% x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
% states = [x;y;theta]; n_states = length(states);

%% Ações de controle (Apenas força sobre o carrinho)
F = SX.sym('F');
controls = [F]; n_controls = length(controls);
rhs = [
  dx;
  dtheta1;
  dtheta2;
  (2*F*J1*J2 + 2*F*J2*M1*l1^2 + 2*F*J1*M2*l2^2 - 2*F0*J1*J2*dx + F*L1^2*M2^2*l2^2 + 2*F*J2*L1^2*M2 - F*L1^2*M2^2*l2^2*cos(2*theta1 - 2*theta2) - F0*L1^2*M2^2*l2^2*dx - 2*F0*J2*L1^2*M2*dx - 2*F0*J2*M1*l1^2*dx - 2*F0*J1*M2*l2^2*dx - J2*L1^2*M2^2*g*sin(2*theta1) + 2*J2*L1^3*M2^2*dtheta1^2*sin(theta1) - J2*M1^2*g*l1^2*sin(2*theta1) - J1*M2^2*g*l2^2*sin(2*theta2) + 2*F*M1*M2*l1^2*l2^2 + 2*J2*M1^2*l1^3*dtheta1^2*sin(theta1) + 2*J1*M2^2*l2^3*dtheta2^2*sin(theta2) + J1*L1*M2^2*l2^2*dtheta1^2*sin(theta1) + J2*L1^2*M2^2*l2*dtheta2^2*sin(theta2) - 2*F0*M1*M2*l1^2*l2^2*dx + 2*J1*J2*L1*M2*dtheta1^2*sin(theta1) - F2*L1^2*M2^2*l2*cos(2*theta1 - theta2)*dtheta2 - J1*L1*M2^2*l2^2*sin(theta1 - 2*theta2)*dtheta1^2 + 2*J1*J2*M1*l1*dtheta1^2*sin(theta1) + 2*J1*J2*M2*l2*dtheta2^2*sin(theta2) - M1^2*M2*g*l1^2*l2^2*sin(2*theta1) - M1*M2^2*g*l1^2*l2^2*sin(2*theta2) + F0*L1^2*M2^2*l2^2*cos(2*theta1 - 2*theta2)*dx + 2*M1^2*M2*l1^3*l2^2*dtheta1^2*sin(theta1) + 2*M1*M2^2*l1^2*l2^3*dtheta2^2*sin(theta2) + F1*L1*M2^2*l2^2*dtheta1*cos(theta1) + F2*L1^2*M2^2*l2*dtheta2*cos(theta2) + J2*L1^2*M2^2*l2*dtheta2^2*sin(2*theta1 - theta2) + 2*F1*J2*L1*M2*dtheta1*cos(theta1) - F1*L1*M2^2*l2^2*cos(theta1 - 2*theta2)*dtheta1 + 2*F1*J2*M1*l1*dtheta1*cos(theta1) + 2*F2*J1*M2*l2*dtheta2*cos(theta2) + 2*J2*L1*M1*M2*l1^2*dtheta1^2*sin(theta1) + 2*J2*L1^2*M1*M2*l1*dtheta1^2*sin(theta1) + L1*M1*M2^2*l1*l2^3*dtheta2^2*sin(2*theta1 - theta2) - L1*M1*M2^2*l1^2*l2^2*sin(theta1 - 2*theta2)*dtheta1^2 + L1^2*M1*M2^2*l1*l2^2*sin(theta1 - 2*theta2)*dtheta1^2 + 2*J1*M1*M2*l1*l2^2*dtheta1^2*sin(theta1) + 2*J2*M1*M2*l1^2*l2*dtheta2^2*sin(theta2) - L1*M1*M2^2*g*l1*l2^2*sin(2*theta1) + L1*M1*M2^2*g*l1*l2^2*sin(2*theta2) - L1*M1*M2^2*l1*l2^3*dtheta2^2*sin(theta2) - 2*J2*L1*M1*M2*g*l1*sin(2*theta1) + 2*F1*M1*M2*l1*l2^2*dtheta1*cos(theta1) + 2*F2*M1*M2*l1^2*l2*dtheta2*cos(theta2) + L1*M1*M2^2*l1^2*l2^2*dtheta1^2*sin(theta1) + L1^2*M1*M2^2*l1*l2^2*dtheta1^2*sin(theta1) - J2*L1*M1*M2*l1*l2*dtheta2^2*sin(theta2) - F2*L1*M1*M2*l1*l2*cos(2*theta1 - theta2)*dtheta2 - F2*L1*M1*M2*l1*l2*dtheta2*cos(theta2) + J2*L1*M1*M2*l1*l2*dtheta2^2*sin(2*theta1 - theta2))/(J2*L1^2*M2^2 + J2*M1^2*l1^2 + J1*M2^2*l2^2 + 2*J1*J2*M0 + 2*J1*J2*M1 + 2*J1*J2*M2 + 2*J2*M0*M1*l1^2 + 2*J1*M0*M2*l2^2 + 2*J1*M1*M2*l2^2 + 2*J2*M1*M2*l1^2 - J2*L1^2*M2^2*cos(2*theta1) - J2*M1^2*l1^2*cos(2*theta1) - J1*M2^2*l2^2*cos(2*theta2) + L1^2*M0*M2^2*l2^2 + L1^2*M1*M2^2*l2^2 + 2*J2*L1^2*M0*M2 + 2*J2*L1^2*M1*M2 + M1*M2^2*l1^2*l2^2 + M1^2*M2*l1^2*l2^2 - M1^2*M2*l1^2*l2^2*cos(2*theta1) - M1*M2^2*l1^2*l2^2*cos(2*theta2) - L1^2*M0*M2^2*l2^2*cos(2*theta1 - 2*theta2) - L1^2*M1*M2^2*l2^2*cos(2*theta1 - 2*theta2) - L1*M1*M2^2*l1*l2^2 + 2*M0*M1*M2*l1^2*l2^2 - 2*J2*L1*M1*M2*l1 - L1*M1*M2^2*l1*l2^2*cos(2*theta1) + L1*M1*M2^2*l1*l2^2*cos(2*theta2) - 2*J2*L1*M1*M2*l1*cos(2*theta1) + L1*M1*M2^2*l1*l2^2*cos(2*theta1 - 2*theta2));
  -(F1*M2^2*l2^2*dtheta1 + 2*F1*J2*M0*dtheta1 + 2*F1*J2*M1*dtheta1 + 2*F1*J2*M2*dtheta1 - 2*J2*M1^2*g*l1*sin(theta1) + J2*L1^2*M2^2*sin(2*theta1)*dtheta1^2 + J2*M1^2*l1^2*sin(2*theta1)*dtheta1^2 + F*L1*M2^2*l2^2*cos(theta1) + 2*F*J2*L1*M2*cos(theta1) - F*L1*M2^2*l2^2*cos(theta1 - 2*theta2) + 2*F1*M0*M2*l2^2*dtheta1 + 2*F1*M1*M2*l2^2*dtheta1 + 2*F*J2*M1*l1*cos(theta1) - F1*M2^2*l2^2*cos(2*theta2)*dtheta1 - 2*J2*L1*M2^2*g*sin(theta1) - 2*J2*M0*M1*g*l1*sin(theta1) - 2*J2*M1*M2*g*l1*sin(theta1) + J2*L1*M2^2*l2*sin(theta1 - theta2)*dtheta2^2 + M1^2*M2*l1^2*l2^2*sin(2*theta1)*dtheta1^2 + M1*M2^2*g*l1*l2^2*sin(theta1 - 2*theta2) + L1^2*M0*M2^2*l2^2*dtheta1^2*sin(2*theta1 - 2*theta2) + L1^2*M1*M2^2*l2^2*dtheta1^2*sin(2*theta1 - 2*theta2) + M1*M2^2*l1*l2^3*dtheta2^2*sin(theta1 + theta2) + F2*L1*M2^2*l2*dtheta2*cos(theta1 + theta2) + 2*F*M1*M2*l1*l2^2*cos(theta1) + 2*L1*M0*M2^2*l2^3*sin(theta1 - theta2)*dtheta2^2 + 2*L1*M1*M2^2*l2^3*sin(theta1 - theta2)*dtheta2^2 - M1*M2^2*l1*l2^3*sin(theta1 - theta2)*dtheta2^2 - F2*L1*M2^2*l2*cos(theta1 - theta2)*dtheta2 + J2*L1*M2^2*l2*dtheta2^2*sin(theta1 + theta2) - F0*L1*M2^2*l2^2*dx*cos(theta1) - L1*M0*M2^2*g*l2^2*sin(theta1) - L1*M1*M2^2*g*l2^2*sin(theta1) - 2*F0*J2*L1*M2*dx*cos(theta1) - 2*J2*L1*M0*M2*g*sin(theta1) - 2*J2*L1*M1*M2*g*sin(theta1) - M1*M2^2*g*l1*l2^2*sin(theta1) - 2*M1^2*M2*g*l1*l2^2*sin(theta1) + F0*L1*M2^2*l2^2*cos(theta1 - 2*theta2)*dx - L1*M0*M2^2*g*l2^2*sin(theta1 - 2*theta2) - L1*M1*M2^2*g*l2^2*sin(theta1 - 2*theta2) - 2*F0*J2*M1*l1*dx*cos(theta1) - J2*M1*M2*l1*l2*sin(theta1 - theta2)*dtheta2^2 + 2*J2*L1*M1*M2*l1*sin(2*theta1)*dtheta1^2 - L1*M1*M2^2*l1*l2^2*dtheta1^2*sin(2*theta1 - 2*theta2) + F2*M1*M2*l1*l2*dtheta2*cos(theta1 + theta2) - 2*F2*L1*M0*M2*l2*cos(theta1 - theta2)*dtheta2 - 2*F2*L1*M1*M2*l2*cos(theta1 - theta2)*dtheta2 + F2*M1*M2*l1*l2*cos(theta1 - theta2)*dtheta2 + J2*M1*M2*l1*l2*dtheta2^2*sin(theta1 + theta2) - 2*F0*M1*M2*l1*l2^2*dx*cos(theta1) - 2*M0*M1*M2*g*l1*l2^2*sin(theta1) + 2*J2*L1*M0*M2*l2*sin(theta1 - theta2)*dtheta2^2 + 2*J2*L1*M1*M2*l2*sin(theta1 - theta2)*dtheta2^2 + L1*M1*M2^2*l1*l2^2*sin(2*theta1)*dtheta1^2)/(J2*L1^2*M2^2 + J2*M1^2*l1^2 + J1*M2^2*l2^2 + 2*J1*J2*M0 + 2*J1*J2*M1 + 2*J1*J2*M2 + 2*J2*M0*M1*l1^2 + 2*J1*M0*M2*l2^2 + 2*J1*M1*M2*l2^2 + 2*J2*M1*M2*l1^2 - J2*L1^2*M2^2*cos(2*theta1) - J2*M1^2*l1^2*cos(2*theta1) - J1*M2^2*l2^2*cos(2*theta2) + L1^2*M0*M2^2*l2^2 + L1^2*M1*M2^2*l2^2 + 2*J2*L1^2*M0*M2 + 2*J2*L1^2*M1*M2 + M1*M2^2*l1^2*l2^2 + M1^2*M2*l1^2*l2^2 - M1^2*M2*l1^2*l2^2*cos(2*theta1) - M1*M2^2*l1^2*l2^2*cos(2*theta2) - L1^2*M0*M2^2*l2^2*cos(2*theta1 - 2*theta2) - L1^2*M1*M2^2*l2^2*cos(2*theta1 - 2*theta2) - L1*M1*M2^2*l1*l2^2 + 2*M0*M1*M2*l1^2*l2^2 - 2*J2*L1*M1*M2*l1 - L1*M1*M2^2*l1*l2^2*cos(2*theta1) + L1*M1*M2^2*l1*l2^2*cos(2*theta2) - 2*J2*L1*M1*M2*l1*cos(2*theta1) + L1*M1*M2^2*l1*l2^2*cos(2*theta1 - 2*theta2));
  (2*J1*M2^2*g*l2*sin(theta2) - F2*M1^2*l1^2*dtheta2 - 2*F2*J1*M0*dtheta2 - 2*F2*J1*M1*dtheta2 - 2*F2*J1*M2*dtheta2 - F2*L1^2*M2^2*dtheta2 - J1*M2^2*l2^2*sin(2*theta2)*dtheta2^2 - F*L1^2*M2^2*l2*cos(theta2) - 2*F2*L1^2*M0*M2*dtheta2 - 2*F2*L1^2*M1*M2*dtheta2 - 2*F2*M0*M1*l1^2*dtheta2 - 2*F2*M1*M2*l1^2*dtheta2 - 2*F*J1*M2*l2*cos(theta2) + F2*L1^2*M2^2*cos(2*theta1)*dtheta2 + F2*M1^2*l1^2*cos(2*theta1)*dtheta2 + F*L1^2*M2^2*l2*cos(2*theta1 - theta2) + 2*J1*M0*M2*g*l2*sin(theta2) + 2*J1*M1*M2*g*l2*sin(theta2) + J1*L1*M2^2*l2*sin(theta1 - theta2)*dtheta1^2 - M1*M2^2*l1^2*l2^2*sin(2*theta2)*dtheta2^2 + L1^2*M0*M2^2*l2^2*dtheta2^2*sin(2*theta1 - 2*theta2) + L1^2*M1*M2^2*l2^2*dtheta2^2*sin(2*theta1 - 2*theta2) - M1^2*M2*l1^3*l2*dtheta1^2*sin(theta1 + theta2) - F1*L1*M2^2*l2*dtheta1*cos(theta1 + theta2) - 2*F*M1*M2*l1^2*l2*cos(theta2) - F0*L1^2*M2^2*l2*cos(2*theta1 - theta2)*dx - L1^2*M0*M2^2*g*l2*sin(2*theta1 - theta2) - L1^2*M1*M2^2*g*l2*sin(2*theta1 - theta2) + 2*L1^3*M0*M2^2*l2*sin(theta1 - theta2)*dtheta1^2 + 2*L1^3*M1*M2^2*l2*sin(theta1 - theta2)*dtheta1^2 + 2*F2*L1*M1*M2*l1*dtheta2 + M1^2*M2*g*l1^2*l2*sin(2*theta1 - theta2) - M1^2*M2*l1^3*l2*sin(theta1 - theta2)*dtheta1^2 + F1*L1*M2^2*l2*cos(theta1 - theta2)*dtheta1 - J1*L1*M2^2*l2*dtheta1^2*sin(theta1 + theta2) + F0*L1^2*M2^2*l2*dx*cos(theta2) + L1^2*M0*M2^2*g*l2*sin(theta2) + L1^2*M1*M2^2*g*l2*sin(theta2) + 2*M1*M2^2*g*l1^2*l2*sin(theta2) + M1^2*M2*g*l1^2*l2*sin(theta2) + 2*F0*J1*M2*l2*dx*cos(theta2) - J1*M1*M2*l1*l2*sin(theta1 - theta2)*dtheta1^2 - L1*M1*M2^2*l1*l2^2*dtheta2^2*sin(2*theta1 - 2*theta2) + F*L1*M1*M2*l1*l2*cos(theta2) - F1*M1*M2*l1*l2*dtheta1*cos(theta1 + theta2) + L1*M1*M2^2*g*l1*l2*sin(2*theta1 - theta2) - L1*M1^2*M2*g*l1*l2*sin(2*theta1 - theta2) - L1*M1*M2^2*l1^2*l2*dtheta1^2*sin(theta1 + theta2) + L1*M1^2*M2*l1^2*l2*dtheta1^2*sin(theta1 + theta2) + L1^2*M1*M2^2*l1*l2*dtheta1^2*sin(theta1 + theta2) + 2*F1*L1*M0*M2*l2*cos(theta1 - theta2)*dtheta1 + 2*F1*L1*M1*M2*l2*cos(theta1 - theta2)*dtheta1 - F1*M1*M2*l1*l2*cos(theta1 - theta2)*dtheta1 + 2*F2*L1*M1*M2*l1*cos(2*theta1)*dtheta2 - J1*M1*M2*l1*l2*dtheta1^2*sin(theta1 + theta2) + 2*F0*M1*M2*l1^2*l2*dx*cos(theta2) - 3*L1*M1*M2^2*g*l1*l2*sin(theta2) - L1*M1^2*M2*g*l1*l2*sin(theta2) + 2*M0*M1*M2*g*l1^2*l2*sin(theta2) + F*L1*M1*M2*l1*l2*cos(2*theta1 - theta2) + L1*M1*M2^2*l1^2*l2*sin(theta1 - theta2)*dtheta1^2 + L1*M1^2*M2*l1^2*l2*sin(theta1 - theta2)*dtheta1^2 - 3*L1^2*M1*M2^2*l1*l2*sin(theta1 - theta2)*dtheta1^2 + 2*J1*L1*M0*M2*l2*sin(theta1 - theta2)*dtheta1^2 + 2*J1*L1*M1*M2*l2*sin(theta1 - theta2)*dtheta1^2 + L1*M1*M2^2*l1*l2^2*sin(2*theta2)*dtheta2^2 - F0*L1*M1*M2*l1*l2*cos(2*theta1 - theta2)*dx - L1*M0*M1*M2*g*l1*l2*sin(2*theta1 - theta2) - F0*L1*M1*M2*l1*l2*dx*cos(theta2) - L1*M0*M1*M2*g*l1*l2*sin(theta2) + 2*L1*M0*M1*M2*l1^2*l2*sin(theta1 - theta2)*dtheta1^2)/(J2*L1^2*M2^2 + J2*M1^2*l1^2 + J1*M2^2*l2^2 + 2*J1*J2*M0 + 2*J1*J2*M1 + 2*J1*J2*M2 + 2*J2*M0*M1*l1^2 + 2*J1*M0*M2*l2^2 + 2*J1*M1*M2*l2^2 + 2*J2*M1*M2*l1^2 - J2*L1^2*M2^2*cos(2*theta1) - J2*M1^2*l1^2*cos(2*theta1) - J1*M2^2*l2^2*cos(2*theta2) + L1^2*M0*M2^2*l2^2 + L1^2*M1*M2^2*l2^2 + 2*J2*L1^2*M0*M2 + 2*J2*L1^2*M1*M2 + M1*M2^2*l1^2*l2^2 + M1^2*M2*l1^2*l2^2 - M1^2*M2*l1^2*l2^2*cos(2*theta1) - M1*M2^2*l1^2*l2^2*cos(2*theta2) - L1^2*M0*M2^2*l2^2*cos(2*theta1 - 2*theta2) - L1^2*M1*M2^2*l2^2*cos(2*theta1 - 2*theta2) - L1*M1*M2^2*l1*l2^2 + 2*M0*M1*M2*l1^2*l2^2 - 2*J2*L1*M1*M2*l1 - L1*M1*M2^2*l1*l2^2*cos(2*theta1) + L1*M1*M2^2*l1*l2^2*cos(2*theta2) - 2*J2*L1*M1*M2*l1*cos(2*theta1) + L1*M1*M2^2*l1*l2^2*cos(2*theta1 - 2*theta2))];

% v = SX.sym('v'); omega = SX.sym('omega');
% controls = [v;omega]; n_controls = length(controls);
% rhs = [v*cos(theta);v*sin(theta);omega]; % system r.h.s

%% Funções necessárias praotimização

% Cálcula derivadas do estado em i+1, dado estado atual e F em i (EDO)
f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
% Variáveis de PROJETO: somente as forças (Estados são adicionados dps)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
% Parâmetros: 6 Estados iniciais  (p/ CIs) e 6 estados finais (para calculo
% da func obj)
P = SX.sym('P',n_states + n_states);

% f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
% U = SX.sym('U',n_controls,N); % Decision variables (controls)
% P = SX.sym('P',n_states + n_states);
% parameters (which include the initial state and the reference state)

%% Matriz com variáveis de estados (linhas) em cada evento (colunas)
X = SX.sym('X',n_states,(N+1)); % prevejo N eventos (+1 é o inicial)
% A vector that represents the states over the optimization problem.

% AQUI NÃO SE FAZ PREVISÃO DO X USANDO EULER, POIS O ESTADO X NO HORIZONTE
% DE EVENTOS É CALCULADO COMO RESTRIÇÃO (E NÃO DE ANTEMÃO COMO NO SINGLESHOOTING)


%% Quantidades de otimização

% Função objetivo
obj = 0; % Objective function
% Vetor de restrições
g = [];  % constraints vector

% Matriz Q: penaliza estados X; R: penaliza o controle
% Q = eye(6);
% Q(1,1) = 0; Q(4,4) = 0;  % Posição e vel. do carrinho nem tanto
% Q(2,2) = 1; Q(3,3) = 1; % ângulos theta1 e theta2 são mais importantes
% Q(5,5) = 0; Q(6,6) = 0; % Velocidade dos pêndulos meio importantes
% R = 0; % Penalizar força aplicada



% Q = zeros(3,3); Q(1,1) = 1;Q(2,2) = 5;Q(3,3) = 0.1; % weighing matrices (states)
% R = zeros(2,2); R(1,1) = 0.5; R(2,2) = 0.05; % weighing matrices (controls)

st  = X(:,1); % initial state
g = [g;st-P(1:6)]; % initial condition constraints
h=T;
for k = 1:N
%   st = X(:,k);  con = U(:,k);
%   obj = obj+(st-P(7:12))'*Q*(st-P(7:12)) + con'*R*con;
%   %     obj = obj+(st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con; % calculate obj
%   st_next = X(:,k+1);
%   f_value = f(st,con);
%   st_next_euler = st+ (T*f_value);
%   g = [g;st_next-st_next_euler]; % compute constraints

    st = X(:,k);  con = U(:,k);
    obj = obj+(st-P(7:12))'*Q*(st-P(7:12)) + con'*R*con; % calculate obj
    st_next = X(:,k+1);
    k1 = f(st, con);   % new 
    k2 = f(st + h/2*k1, con); % new
    k3 = f(st + h/2*k2, con); % new
    k4 = f(st + h*k3, con); % new
    st_next_RK4=st +h/6*(k1 +2*k2 +2*k3 +k4); % new    
    % f_value = f(st,con);
    % st_next_euler = st+ (h*f_value);
    % g = [g;st_next-st_next_euler]; % compute constraints
    g = [g;st_next-st_next_RK4]; % compute constraints % new
end

% Penalize change of control variables
for k=2:N
  % Get the state X and control actions U for the each event in the horizon
  con0 = U(:,k-1); con = U(:,k);
  % Compute the objective function using the weigth matrices defined before
  obj = obj+ (con-con0)'*S*(con-con0); % calculate obj
end

% Restrição do carrinho
% for k = 1:N+1   % nos N instates (+1 é o inicial)
%   g = [g ; X(1,k)]; % Posição hor. do carrinho (X(1,~))
% %   g = [g ; X(2,k)];   %state theta1
% %   g = [g ; X(3,k)];   %state theta2
% end

%% Montar objeto de otimização
% make the decision variable one column  vector
% VARIAVEIS DE PROEJETO: ESTADOS (MULTISHOOTING PRECISA) E FORÇAS
OPT_variables = [reshape(X,n_states*(N+1),1);reshape(U,n_controls*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;

% As EDOS x - T*dxdt=0 estão nas restrições g, que devem valer zero (ou
% próximo de zero) numa solução viável. Por isso as restrições devem ter
% valor máximo=mínimo=0.
args.lbg(1:6*(N+1)) = 0;  % -1e-20  % Equality constraints
args.ubg(1:6*(N+1)) = 0;  % 1e-20   % Equality constraints
% posição do carrinho (DEIXA COMENTADO)
% args.lbg(6*(N+1)+1:6*(N+1)+N+1) = 0;  % Lim. esq. do carrinho
% args.ubg(6*(N+1)+1:6*(N+1)+N+1) = 2;  % Lim. dir. do carrinho

%% Limites do estado e controle

% Posição mínima (esq) e máxima (dir) do carrinho
args.lbx(1:6:6*(N+1),1) = -1; %state x lower bound
args.ubx(1:6:6*(N+1),1) = +1; %state x upper bound
% Ângulo 1 mínimo e máximo (dir)
args.lbx(2:6:6*(N+1),1) = -2*pi; %state y lower bound
args.ubx(2:6:6*(N+1),1) = +2*pi; %state y upper bound
% Ângulo 2 mínimo e máximo (dir)
args.lbx(3:6:6*(N+1),1) = -inf; %state y lower bound
args.ubx(3:6:6*(N+1),1) = +inf; %state y upper bound
% Velocidade mínima e máxima do carrinho
args.lbx(4:6:6*(N+1),1) = -inf; %state x lower bound
args.ubx(4:6:6*(N+1),1) = +inf; %state x upper bound
% Velocidade mínima e máxima do theta1
args.lbx(5:6:6*(N+1),1) = -inf; %state x lower bound
args.ubx(5:6:6*(N+1),1) = +inf; %state x upper bound
% Velocidade mínima e máxima do theta2
args.lbx(6:6:6*(N+1),1) = -inf; %state x lower bound
args.ubx(6:6:6*(N+1),1) = +inf; %state x upper bound

% Força mínima
args.lbx(6*(N+1)+1:1:6*(N+1)+1*N,1) = F_min; %v lower bound
args.ubx(6*(N+1)+1:1:6*(N+1)+1*N,1) = F_max; %v lower bound

% args.lbx(3*(N+1)+1:2:3*(N+1)+2*N,1) = v_min; %v lower bound
% args.ubx(3*(N+1)+1:2:3*(N+1)+2*N,1) = v_max; %v upper bound
% args.lbx(3*(N+1)+2:2:3*(N+1)+2*N,1) = omega_min; %omega lower bound
% args.ubx(3*(N+1)+2:2:3*(N+1)+2*N,1) = omega_max; %omega upper bound
%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = [0 ; 0 ; deg2rad(1); 0; 0; 0];    % initial condition.
xs = [0 ; 0 ; 0; 0; 0; 0]; % Reference posture.

% Todos os estados pro gráfico no final
xx(:,1) = x0; % xx contains the history of states
% Todos os instantes pro gráfico
t(1) = t0;

% Controles pra cada iteração
u0 = zeros(N,1);        % two control inputs for each robot
% u0 = zeros(N,2);        % two control inputs for each robot

% Estados de cada iteração (começa com CI em x0) PREVISAO INICIAL DE X
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

% Tempo de simulação (s_
sim_tim = 1; % Maximum simulation time

% Start MPC
mpciter = 0; % Número de iterações de MPC
xx1 = [];
u_cl=[];

% the main simulaton loop... it works as long as the error is greater
% than 10^-6 and the number of mpc steps is less than its maximum
% value.
tic
while(norm((x0-xs),2) > 1e-2 && mpciter < sim_tim / T)
  % args.p tem x0 (condições iniciais) e xs (condições de referencia des.)
  args.p   = [x0;xs]; % set the values of the parameters vector
  % Variável de projeto: vetor coluna com 6 estados (vezes N+1 eventos) e 1
  % controle (vezes N+1 eventos)
  args.x0  = [reshape(X0',6*(N+1),1);reshape(u0',1*N,1)];
  sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
    'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
  % u: N Ações de controle desta iteração do MPC
  u = reshape(full(sol.x(6*(N+1)+1:end))',1,N)'; % get controls only from the solution
  % xx1: eventos (linhas), estado (colunas) e iteração MPC (paginas)
  xx1(:,1:6,mpciter+1)= reshape(full(sol.x(1:6*(N+1)))',6,N+1)'; % get solution TRAJECTORY
  % u_c1: iteração MPC (linha) e ação de controle (coluna_
  u_cl= [u_cl ; u(1,:)];
  % Tempo (duh)
  t(mpciter+1) = t0;
  % Apply the control and shift the solution
  [t0, x0, u0] = shift(T, t0, x0, u,f);
  xx(:,mpciter+2) = x0;
  X0 = reshape(full(sol.x(1:6*(N+1)))',6,N+1)'; % get solution TRAJECTORY
  % Shift trajectory to initialize the next step
  X0 = [X0(2:end,:);X0(end,:)];
  mpciter
  mpciter = mpciter + 1;
end
toc

%% Animação
set(0,'DefaultFigureWindowStyle','docked')

% L1 = 0.490;
% L2 = 0.490;
% Plor cart
set(gcf,'color','white')

for i = 1:length(t)
  set(gcf,'Visible','on')
  set(gcf,'color','white')
  % Update states
  x_cart = xx(1,i);
  t1 = xx(2,i); % theta1
  t2 = xx(3,i); % theta2
  x1 = x_cart + L1*sin(t1);
  x2 = x1 + L2*sin(t2);
  y_cart = 0;
  y1 = y_cart + L1*cos(t1);
  y2 = y1 + L2*cos(t2);
  clf
  subplot(1,2,1)
  hold on
  grid on
  %   axis((L1+L2)*[-1 +1 -1 +1])
  ylim((L1+L2)*[-1 +1])
  axis equal
  % Plot the cart
  plot(x_cart,y_cart,'ks','markersize',5)
  % Plot the lower pendulum
  plot([x_cart x1],[y_cart y1],'-g','linewidth',2)
  % Plot the upper pendulum
  plot([x1 x2],[y1 y2],'-r','linewidth',2)
  hold off
  %   legend({'Cart';'Lower pendulum';'Upper pendulum'},...
  %     'interpreter','latex','location','nw')
  grid on
  titulo = sprintf('t = %g s',t(i));
  title(titulo,'interpreter','latex')
  xlabel('x $(m)$','Interpreter',"latex")
  ylabel('y $(m)$','Interpreter',"latex")
  xlim([-1,+1])


  % Plot control action
  subplot(1,2,2)
  hold on
  plot(t(1:i),u_cl(1:i))
  plot([0 t(end)],[F_min F_min],'--r')
  plot([0 t(end)],[F_max F_max],'--r')
  hold off
  %   ylim(1e3*[-6 4])
  axis tight
  xlim([0 t(end)])
  grid on
  xlabel('t $(s)$','Interpreter',"latex")
  ylabel('F $(N)$','Interpreter',"latex")

  pause(.001)
end
