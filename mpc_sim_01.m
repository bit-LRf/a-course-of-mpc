clear;clc;close all

x_0 = [1;0]; % [v;s]

N = 10; % 仿真步数
h = 0.1; % 仿真步长 （s）
T = h:h:h*N; % 未来N步时间

% 初始化N个x
X = zeros(2,N); 

% 假设控制量始终是0.1，初始化N个u
U = 0.1*ones(1,N); 

% 设置x的初始值
x_k = x_0;
for k = 1:N
    % 提取当前控制量
    u_k = U(k);

    % 计算当前导数
    dot_x = model(x_k,u_k);

    % 前向欧拉
    x_new = x_k + dot_x*h;
    X(:,k) = x_new;
    x_k = x_new;
end

figure()
hold on
stairs(T,X(1,:))
stairs(T,X(2,:))
legend('v','s');

function dot_x = model(x, u)

v = x(1); % m/s
s = x(2); % m

F = u(1); % N

mass = 1; % kg

dot_v = F/mass;
dot_s = v;

dot_x = [dot_v;dot_s];

end