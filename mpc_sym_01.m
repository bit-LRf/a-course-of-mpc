clear;clc;

x = casadi.SX.sym('x',2,1); % 声明一个2×1大小的符号变量
disp(x)

y = x(1) + x(2)^2; % 建立符号表达式
disp(y)

fun_y = casadi.Function('fun_y',{x},{y}); % 声明一个由x到y的映射函数
disp(fun_y)

y_11 = fun_y([1;1]); % 假设x=[1;1]，计算y的值
disp(y_11)