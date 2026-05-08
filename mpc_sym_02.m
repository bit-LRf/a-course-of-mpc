clear;clc;close all

% 初始值
x_0 = casadi.SX.sym('x_init',2,1); % [v;s]

N = 10; % 仿真步数
h = 0.1; % 仿真步长 （s）
T = h:h:h*N; % 未来N步时间

% 初始化N个x
X = casadi.SX.zeros(2,N); 

% 初始化N个u
U = casadi.SX.sym('u',N,1); 

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

% 声明函数
fun_X = casadi.Function('predict_X',{x_0,U},{X});

% 进行仿真（预测）
X_predict = full(fun_X([1;0],0.1*ones(1,N)));

figure()
hold on
stairs(T,X_predict(1,:))
stairs(T,X_predict(2,:))
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