clear;
clc;
close all

% 初始值
x_0 = casadi.SX.sym('x_init',2,1); % [v;s]

N = 10; % 仿真步数
h = 0.1; % 仿真步长 （s）
T = h:h:h*N; % 未来N步时间

% 初始化N个x
X = casadi.SX.sym('x',2,N); 

% 初始化N个u
U = casadi.SX.sym('u',1,N); 

% 初始化g
g = {};

v_ref = casadi.SX.sym('v_ref',1,N); % 期望速度
J = 0; % 代价
x_k = x_0; % 设置x的初始值
for k = 1:N
    % 提取当前控制量
    u_k = U(k);

    % 计算当前导数
    dot_x = model(x_k,u_k);

    % 前向欧拉
    x_new = x_k + dot_x*h;
    g = [g;{X(:,k) - x_new}];

    % 更新状态
    x_k = X(:,k);

    % 计算代价
    J = J + (x_k(1) - v_ref(k))^2 + 0.01*(u_k)^2;
end

% 优化变量记为r
r = [U(:);X(:)];

% 参数
p = [x_0;v_ref(:)];

% 约束
g = vertcat(g{:});

% 上下界
lbx = [0;-inf];
lbu = -inf;

ubx = [inf;inf];
ubu = inf;

lbr = [];
ubr = [];
for i = 1:N
    lbr = [lbr;lbu];
    ubr = [ubr;ubu];
end
for i = 1:N
    lbr = [lbr;lbx];
    ubr = [ubr;ubx];
end

nlp_prob = struct('f',J,'x',r,'g',g,'p',p);
solver = casadi.nlpsol('solver','ipopt',nlp_prob);

sol = solver('x0',zeros(size(r)),'lbx',lbr,'ubx',ubr,...
    'lbg',zeros(size(g)),'ubg',zeros(size(g)),'p',[[0;0];1.2*ones(N,1)]);

opt_r = full(sol.x);

opt_U = opt_r(1:N);
opt_X = opt_r(N + 1:end);
opt_X = reshape(opt_X,2,N);

figure()
hold on
stairs(T,opt_X(1,:))
stairs(T,opt_X(2,:))
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