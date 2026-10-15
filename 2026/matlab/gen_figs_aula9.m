%% Figuras da Aula 9 (SEL0359 - Controle Digital, 2026): Controle LQR
% Exemplo: A = [3 0; 5 4], B = [1; 2], Q = diag(2,3), R = 1
% Gera PNGs em ../images/   (Prof. Marcos R. Fernandes)

clc; clear all; close all;
out = fullfile(fileparts(mfilename('fullpath')), '..', 'images');
set(groot, 'defaultAxesFontSize', 13, 'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
cB = [0.165 0.494 0.878]; cR = [0.878 0.141 0.141]; cG = [0.122 0.620 0.122]; cO = [0.93 0.55 0.10];

A = [3 0; 5 4];  B = [1; 2];  Q = diag([2 3]);  R = 1;
x0 = [1; 1];

%% ---------- Horizonte finito: N = 20, P_N = 0.1 I ----------
N = 20;  P = 0.1*eye(2);  K = zeros(N, 2);
for k = N-1:-1:0                                   % recursão para trás
    Kk = (R + B'*P*B) \ (B'*P*A);
    P  = Q + A'*P*A - A'*P*B*Kk;
    K(k+1, :) = Kk;
end
P0 = P;
[Kinf, Pinf] = dlqr(A, B, Q, R);                   % horizonte infinito
x = zeros(2, N+1); u = zeros(1, N); x(:,1) = x0; J = 0;
for k = 1:N
    u(k) = -K(k,:)*x(:,k);
    J = J + x(:,k)'*Q*x(:,k) + u(k)'*R*u(k);
    x(:,k+1) = A*x(:,k) + B*u(k);
end
J = J + x(:,N+1)'*(0.1*eye(2))*x(:,N+1);
fprintf('Horizonte finito: K_19 = [%.4f %.4f], K_0 = [%.4f %.4f]\n', K(N,:), K(1,:));
fprintf('J_0->N = %.4f  =  x0''P0 x0 = %.4f\n', J, x0'*P0*x0);

f = figure('Color','w','Position',[100 100 1150 450]);
tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
kk = 0:N-1;
stairs(ax, kk, K(:,1), '-o', 'Color', cB, 'LineWidth', 2, 'MarkerFaceColor', cB, 'MarkerSize', 5);
stairs(ax, kk, K(:,2), '-s', 'Color', cR, 'LineWidth', 2, 'MarkerFaceColor', cR, 'MarkerSize', 5);
yline(ax, Kinf(1), '--', 'Color', cB, 'LineWidth', 1.2);
yline(ax, Kinf(2), '--', 'Color', cR, 'LineWidth', 1.2);
xlabel(ax, '$k$'); ylabel(ax, 'Ganho $K_k$');
legend(ax, {'$K_k(1)$', '$K_k(2)$', '$K_\infty(1)$', '$K_\infty(2)$'}, 'Location', 'southwest');
title(ax, 'Ganho variante no tempo (recurs\~ao de Riccati de $k=N-1$ at\''e $0$)', 'FontSize', 13);
xlim(ax, [0 N-1]); ylim(ax, [0 2.8]);
ax2 = nexttile; hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
stairs(ax2, 0:N, x(1,:), 'Color', cB, 'LineWidth', 2);
stairs(ax2, 0:N, x(2,:), 'Color', cR, 'LineWidth', 2);
stairs(ax2, 0:N-1, u, 'Color', cG, 'LineWidth', 2);
xlabel(ax2, '$k$'); legend(ax2, {'$x_{1,k}$', '$x_{2,k}$', '$u_k$'}, 'Location', 'southeast');
title(ax2, sprintf('Malha fechada, $x_0=[1\\;1]^{\\top}$: $J_{0\\to N}^*=x_0^{\\top} P_0x_0=%.2f$', J), 'FontSize', 13);
xlim(ax2, [0 10]);
exportgraphics(f, fullfile(out, 'aula9_lqr_horizonte_finito.png'), 'Resolution', 200);
close(f);

%% ---------- Horizonte infinito ----------
[P2, K2, polos] = idare(A, B, Q, R);              % mesma solução que dlqr
fprintf('P_inf = [%.4f %.4f; %.4f %.4f]\n', P2');
fprintf('K_inf = [%.4f %.4f], polos MF = %.4f, %.4f\n', K2, polos);
M = 15; x = zeros(2, M+1); u = zeros(1, M); x(:,1) = x0; Jc = zeros(1, M);
for k = 1:M
    u(k) = -Kinf*x(:,k);
    Jc(k) = x(:,k)'*Q*x(:,k) + u(k)^2*R;
    x(:,k+1) = A*x(:,k) + B*u(k);
end
Jinf = x0'*Pinf*x0;
fprintf('J_inf = x0''P x0 = %.4f, soma simulada = %.4f\n', Jinf, sum(Jc));

f = figure('Color','w','Position',[100 100 1150 450]);
tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
tc = linspace(0, 2*pi, 400); plot(ax, cos(tc), sin(tc), 'k-', 'LineWidth', 1);
pa = eig(A); pf = eig(A - B*Kinf);
plot(ax, real(pa), imag(pa), 'x', 'Color', cR, 'MarkerSize', 14, 'LineWidth', 2.5);
plot(ax, real(pf), imag(pf), 'o', 'Color', cB, 'MarkerSize', 11, 'LineWidth', 2, 'MarkerFaceColor', cB);
for i = 1:2
    text(ax, real(pa(i)), 0.35, sprintf('$%.0f$', pa(i)), 'HorizontalAlignment', 'center', 'Color', cR, 'FontSize', 13);
end
text(ax, 0.26, -0.35, sprintf('$%.3f,\\;%.3f$', sort(pf)), 'HorizontalAlignment', 'center', 'Color', cB, 'FontSize', 13);
axis(ax, 'equal'); xlim(ax, [-1.3 4.5]); ylim(ax, [-1.5 1.5]);
xlabel(ax, 'Re$\{z\}$'); ylabel(ax, 'Im$\{z\}$');
legend(ax, {'$|z|=1$', 'malha aberta: eig$(A)$', 'LQR: eig$(A-BK_\infty)$'}, 'Location', 'northeast');
title(ax, 'Polos: o LQR estabiliza o sistema inst\''avel', 'FontSize', 13);
ax2 = nexttile; hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
yyaxis(ax2, 'left');
stairs(ax2, 0:M, x(1,:), '-', 'Color', cB, 'LineWidth', 2);
stairs(ax2, 0:M, x(2,:), '-', 'Color', cR, 'LineWidth', 2);
stairs(ax2, 0:M-1, u, '-', 'Color', cG, 'LineWidth', 2);
ylabel(ax2, 'estados e entrada');
yyaxis(ax2, 'right');
plot(ax2, 0:M-1, cumsum(Jc), '-^', 'Color', cO, 'LineWidth', 1.6, 'MarkerFaceColor', cO, 'MarkerSize', 5);
yline(ax2, Jinf, ':', 'Color', cO, 'LineWidth', 1.5);
ylabel(ax2, 'custo acumulado'); ylim(ax2, [0 1.15*Jinf]);
ax2.YAxis(1).Color = [0 0 0]; ax2.YAxis(2).Color = cO;
xlabel(ax2, '$k$'); xlim(ax2, [0 M-1]);
legend(ax2, {'$x_{1,k}$', '$x_{2,k}$', '$u_k$', '$\sum_{j\le k}(x_j^{\top} Qx_j+u_j^{\top} Ru_j)$', '$x_0^{\top} P_\infty x_0$'}, ...
    'Location', 'east', 'FontSize', 11);
title(ax2, sprintf('$u_k=-K_\\infty x_k$, $x_0=[1\\;1]^{\\top}$: $J_\\infty^*=%.2f$', Jinf), 'FontSize', 13);
exportgraphics(f, fullfile(out, 'aula9_lqr_horizonte_infinito.png'), 'Resolution', 200);
close(f);
