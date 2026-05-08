clear;clc;

x = [1;0]; % [v;s]
u = 0.1; % F

dot_x = model(x,u);

fprintf('dot_x = \n');
disp(dot_x)

function dot_x = model(x, u)

v = x(1); % m/s
s = x(2); % m

F = u(1); % N

mass = 1; % kg

dot_v = F/mass;
dot_s = v;

dot_x = [dot_v;dot_s];

end