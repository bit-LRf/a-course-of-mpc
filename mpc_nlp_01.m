x = casadi.SX.sym('x',2,1);

f = x(1) + x(2)^2;
g = x(1) + x(2) - 1;
lbx = [-1;-1];
ubx = [1;1];

nlp_prob = struct('f',f,'x',x,'g',g);
solver = casadi.nlpsol('solver','ipopt',nlp_prob);

sol = solver('x0',0,'lbx',lbx,'ubx',ubx,...
    'lbg',0,'ubg',0,'p',[]);