% point stabilization + Single shooting
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

T = 0.01; % sampling time [s]
N = 20; % prediction horizon
% rob_diam = 0.3;

%% Valores máximo e mínimo de F
F_max = 1000; F_min = -1000;

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

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
P = SX.sym('P',n_states + n_states);
% parameters (which include the initial and the reference state of the robot)

%% Matriz com variáveis de estados (linhas) em cada evento (colunas
X = SX.sym('X',n_states,(N+1));
% A Matrix that represents the states over the optimization problem.

% compute solution symbolically
X(:,1) = P(1:6); % initial state
for k = 1:N
  st = X(:,k);  con = U(:,k);
  f_value  = f(st,con);
  st_next  = st+ (T*f_value);
  X(:,k+1) = st_next;
end

%% Função que retorna a matriz X para valores de U (controle) e P
% this function to get the optimal trajectory knowing the optimal solution
ff=Function('ff',{U,P},{X});

%% Quantidades de otimização

% Função objetivo
obj = 0; % Objective function
% Vetor de restrições
g = [];  % constraints vector

% Matriz Q: penaliza estados X; R: penaliza o controle
Q = eye(6);
Q(1,1) = 0.01; Q(4,4) = 1; 
Q(2,2) = 10; Q(3,3) = 10;
Q(5,5) = 5; Q(6,6) = 5;
R = 5;

% Q = zeros(3,3); Q(1,1) = 1;Q(2,2) = 5;Q(3,3) = 0.1; % weighing matrices (states)
% R = zeros(2,2); R(1,1) = 0.5; R(2,2) = 0.05; % weighing matrices (controls)

% Computar função objetivo e restrições simbólicamente
% compute objective
for k=1:N
  st = X(:,k);  con = U(:,k);
  obj = obj+(st-P(7:12))'*Q*(st-P(7:12)) + con'*R*con; % calculate obj
  
  %     obj = obj+(st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con; % calculate obj
end

% compute constraints
for k = 1:N+1   % box constraints due to the map margins
  g = [g ; X(2,k)];   %state theta1
  g = [g ; X(3,k)];   %state theta2
end

% make the decision variables one column vector

%% Montar objeto de otimização
OPT_variables = reshape(U,1*N,1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

%% struct de restrições
args = struct;

% Restrições dos estados
args.lbg = deg2rad(-1);  % Valor mínimo de theta1 e theta2
args.ubg = deg2rad(+1);   % Valor mínimo de theta1 e theta2

% inequality constraints (state constraints)
% args.lbg = -2;  % lower bound of the states x and y
% args.ubg = 2;   % upper bound of the states x and y

% Restrições da ação de controle (F: força no carrinho)
args.lbx(1:N,1) = F_min;
args.ubx(1:N,1) = F_max;

% % input constraints
% args.lbx(1:2:2*N-1,1) = v_min; args.lbx(2:2:2*N,1)   = omega_min;
% args.ubx(1:2:2*N-1,1) = v_max; args.ubx(2:2:2*N,1)   = omega_max;

%%
%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SETTING UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
%% Condições iniciais
t0 = 0;
x0 = [0 ; 0 ; deg2rad(90); 0; 0; 0];    % initial condition.
xs = [0 ; 0 ; 0; 0 ; 0 ; 0]; % Valor desejado.

% t0 = 0;
% x0 = [0 ; 0 ; 0];    % initial condition.
% xs = [1.5 ; 1.5 ; 0]; % Reference posture.

%% Histórico de eventos
xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

% Previsão de controle em cada iteração de MPC
u0 = zeros(N,1);

% u0 = zeros(N,2);  % two control inputs

% Tempo de simulação
sim_tim = 2 ; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];

%% Simulação
% the main simulaton loop... it works as long as the error is greater
% than 10^-2 and the number of mpc steps is less than its maximum
% value.
main_loop = tic;
while(norm((x0-xs),2) > 1e-2 && mpciter < sim_tim / T)
  % args.p: [estado atual, estado desejado no final]
  args.p   = [x0;xs]; % set the values of the parameters vector
  % ars.x0: AÇÃO DE CONTROLE atual (zero no começo, mas nas iterações
  % seguintes usa-se a solução anterior como chute inicial!)
  args.x0 = reshape(u0',1*N,1); % initial value of the optimization variables
  %tic
  sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
    'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
  %toc
  
  % A ação de controle nos N eventos é guardada em u
  u = reshape(full(sol.x)',1,N)'; % 1 ação de controle só
  %     u = reshape(full(sol.x)',2,N)';
  
  % A estimativa dos estados a partir de u e args.p é computada:
  ff_value = ff(u',args.p); % compute OPTIMAL solution TRAJECTORY
  % Mas apenas o próximo evento (de N) é adotado para simulação
  xx1(:,1:6,mpciter+1)= full(ff_value)';
  
  %     xx1(:,1:3,mpciter+1)= full(ff_value)';
  
  % Guardar apenas a primeira ação de controle no histórico, pois as
  % próximas "são jogadas fora" (na vdd são usadas como chute inicial da
  % próxima iteração do MPC)
  u_cl= [u_cl ; u(1,:)];
  
  % Guardar valor do instante de tempo
  t(mpciter+1) = t0;
  
  % Dar 1 passo pra frente T segundo no tempo, estado e controle!
  [t0, x0, u0] = shift(T, t0, x0, u,f); % get the initialization of the next optimization step
  
  xx(:,mpciter+2) = x0;
  mpciter
  mpciter = mpciter + 1;
end
main_loop_time = toc(main_loop)

ss_error = norm((x0-xs),2)

%% Animação
set(0,'DefaultFigureWindowStyle','docked')

% L1 = 0.490;
% L2 = 0.490;
% Plor cart
set(gcf,'color','white')

grid on
for i = 1:length(t)
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
  titulo = sprintf('t = %g s',t(i));
  title(titulo,'interpreter','latex')
  pause(.001)
end


