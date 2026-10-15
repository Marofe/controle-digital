%% Sintonia do compensador C(z) = K (z-a)/(z-b) por PSO (Aula 6 - SEL0359, 2026)
clc; clear all; close all;
rng(1);
G  = tf(1,[1 1 0]);  Ts = 0.1;  Gd = c2d(G, Ts);  z = tf('z', Ts);
a  = exp(-Ts);                                    % zero cancela o polo estável
lambda = 1e-3;                                    % peso do esforço de controle
J  = @(p) custo(feedback(p(1)*(z-a)/(z-p(2))*Gd, 1), ...
                feedback(p(1)*(z-a)/(z-p(2)), Gd), lambda);

nP = 20; nIt = 40;  lb = [0 -0.9];  ub = [40 0.95];
w = 0.7;  c1 = 1.5;  c2 = 1.5;
x  = lb + rand(nP,2).*(ub-lb);  v = zeros(nP,2);
pb = x;  Jpb = arrayfun(@(i) J(x(i,:)), (1:nP)');
[Jg, ig] = min(Jpb);  g = pb(ig,:);
for it = 1:nIt
    r1 = rand(nP,2);  r2 = rand(nP,2);
    v  = w*v + c1*r1.*(pb - x) + c2*r2.*(g - x);   % velocidades
    x  = min(max(x + v, lb), ub);                  % posições (com limites)
    for i = 1:nP
        Ji = J(x(i,:));
        if Ji < Jpb(i), Jpb(i) = Ji;  pb(i,:) = x(i,:); end
    end
    [Jm, im] = min(Jpb);
    if Jm < Jg, Jg = Jm;  g = pb(im,:); end
end
K = g(1),  b = g(2)
S = stepinfo(feedback(K*(z-a)/(z-b)*Gd, 1))

function Jv = custo(T, U, lambda)
    if any(abs(pole(T)) >= 1), Jv = 1e6; return; end     % instável
    [y, t] = step(T, 10);  u = step(U, 10);  Ts = t(2) - t(1);
    Jv = sum(t.*abs(1 - y))*Ts + lambda*sum(u.^2)*Ts ... % ITAE + esforço
         + 100*max(0, max(y) - 1.05);                    % penaliza Mp > 5 %
end
