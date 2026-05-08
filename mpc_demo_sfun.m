function [sys,x0,str,ts] = mpc_demo_sfun(t,x,u,flag)
switch flag
    case 0
        [sys,x0,str,ts] = mdlInitializeSizes;
        clear mdlOutputs
    case 2
        sys = mdlUpdates(t,x,u);
    case 3
        sys = mdlOutputs(t,x,u);
    case {1,4,9}
        sys = [];
    otherwise
        error(['unhandled flag = ' , num2str(flag)]);
end
end

function [sys,x0,str,ts] = mdlInitializeSizes
% 仿真参数设置
sizes = simsizes;
sizes.NumContStates = 0;
sizes.NumDiscStates = 0;
sizes.NumOutputs = 1;
sizes.NumInputs = 2;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
str = [];
ts = [0.05,0];
x0 = [];
end

function sys = mdlUpdates(t,x,u)
sys = x;

end

function sys = mdlOutputs(t,x,x_init)
N = 10; % 仿真步数
h = 0.1; % 仿真步长 （s）
T = t + h:h:t + h*N; % 未来N步时间

persistent prob
if isempty(prob)
    % 初始值
    x_0 = casadi.SX.sym('x_init',2,1); % [v;s]
    
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
    
    H_sym = hessian(J,r);
    c_sym = jacobian(J,r);
    
    fun_H = casadi.Function('hessian',{r,p},{H_sym});
    fun_c = casadi.Function('jacobian',{r,p},{c_sym});
    
    Ae_sym = jacobian(g,r);
    fun_Ae = casadi.Function('A_e',{r,p},{Ae_sym});
    fun_be = casadi.Function('b_e',{r,p},{-g});
    
    prob.fun_H = fun_H;
    prob.fun_c = fun_c;
    prob.fun_Ae = fun_Ae;
    prob.fun_be = fun_be;
    prob.n_r = size(r);
    prob.lbr = lbr;
    prob.ubr = ubr;
end

r0 = zeros(prob.n_r);
p = [x_init;1.2*ones(N,1)];
H = full(prob.fun_H(r0,p));
c = full(prob.fun_c(r0,p));
Ae = full(prob.fun_Ae(r0,p));
be = full(prob.fun_be(r0,p));

sol_qp = quadprog(H,c,[],[],Ae,be,prob.lbr,prob.ubr,r0);

opt_U = sol_qp(1:N);

sys = opt_U(1);

end

function dot_x = model(x, u)

v = x(1); % m/s
s = x(2); % m

F = u(1); % N

mass = 1; % kg

dot_v = F/mass;
dot_s = v;

dot_x = [dot_v;dot_s];

end