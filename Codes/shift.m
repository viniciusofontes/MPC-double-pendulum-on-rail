function [t0, x0, u0] = shift(T, t0, x0, u,f)
%% Esta função dá um passo pra frente no tempo, estado e controle

% O estado e controle atual são guardados
st = x0; % essa variável não é do tipo double, mas casadi.DM
con = u(1,:)';

% Computar as derivadas do estado (dxdt: velocidades e acelerações)
f_value = f(st,con);
% Estimar o estado seguinte com Euler explícito (ugh)
st = st+ (T*f_value);
% Transformar estado st em um vetor coluna do tipo double
x0 = full(st);

% Atualizar o instante de tempo
t0 = t0 + T;

% Atualizar
u0 = [u(2:size(u,1),:);u(size(u,1),:)];

% u0 = [u(2:size(u,1),:);u(size(u,1),:)];
end