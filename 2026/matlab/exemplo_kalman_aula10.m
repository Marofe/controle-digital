%% EXEMPLO_KALMAN_AULA10  Exemplo do filtro de Kalman (Aula 10, slide "Exemplo")
% Aula 10 (SEL0359 - Controle Digital, 2026) - Prof. Marcos R. Fernandes
%
%   x_{k+1} = A x_k + B u_k + w_k,   w_k ~ N(0, Q)
%   y_k     = H x_k + eps_k,         eps_k ~ N(0, R)
% com u_k = 1, Q = 0.1 I, R = 10. Só x1 é medido; x2 é estimado pelo filtro.
% Compara o filtro recursivo (K_k variável) com o de horizonte infinito
% (ganho constante obtido da equação de Riccati via idare).
% Gera: ../images/aula10_kalman_exemplo.png

clc; clear; close all;
rng(2);
A = [0.9677 -0.3226; 0 0.6452];  B = [1; 2];  H = [1 0];
Q = 0.1*eye(2);  R = 10;  u = 1;  N = 150;

%% Filtro de Kalman recursivo
x  = [0; 0];  xh = [0; 0];  P = 10*eye(2);                % x_0 e x_0|0 ~ N(xh, P)
X = zeros(2, N);  XH = X;  Y = zeros(1, N);  KK = X;  PP = zeros(2, N);
for k = 1:N
    % planta (simulação)
    x = A*x + B*u + sqrtm(Q)*randn(2,1);
    y = H*x + sqrt(R)*randn;
    % predição
    xh = A*xh + B*u;
    P  = A*P*A' + Q;
    % correção
    K  = P*H'/(R + H*P*H');
    xh = xh + K*(y - H*xh);
    P  = (eye(2) - K*H)*P;
    X(:,k) = x;  XH(:,k) = xh;  Y(k) = y;  KK(:,k) = K;  PP(:,k) = sqrt(diag(P));
end

%% Horizonte infinito (P_inf da equação de Riccati; dualidade com o LQR)
Pinf = idare(A', H', Q, R);                                % P_{k+1|k} em regime
Kinf = Pinf*H'/(R + H*Pinf*H');                            % K_f em regime
L    = A*Kinf;                                             % ganho do observador
fprintf('K_f(k=N) = [%.4f; %.4f]   K_inf = [%.4f; %.4f]   L = [%.4f; %.4f]\n', KK(:,end), Kinf, L);
fprintf('RMS erro: x1 medido = %.2f | x1 estimado = %.2f | x2 estimado = %.2f\n', ...
        rms(Y - X(1,:)), rms(XH(1,:) - X(1,:)), rms(XH(2,:) - X(2,:)));

%% Figura
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');
co = get(groot, 'defaultAxesColorOrder');
k = 1:N;
fig = figure('Color', 'w', 'Position', [50 50 760 860]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on'); ax.FontSize = 13;
plot(ax, k, Y, '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 9);
plot(ax, k, X(1,:), '-', 'Color', co(1,:), 'LineWidth', 2);
plot(ax, k, XH(1,:), '-', 'Color', co(2,:), 'LineWidth', 2);
ylabel(ax, '$x_1$');
legend(ax, {'medi\c{c}\~ao $y_k$', 'real', 'estimado $\hat{x}_{k|k}$'}, 'Location', 'northeast', 'FontSize', 12);
title(ax, '$x_1$ (medido, $R=10$)', 'FontSize', 15);

ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on'); ax.FontSize = 13;
fill(ax, [k fliplr(k)], [XH(2,:) + 3*PP(2,:) fliplr(XH(2,:) - 3*PP(2,:))], co(2,:), ...
     'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(ax, k, X(2,:), '-', 'Color', co(1,:), 'LineWidth', 2);
plot(ax, k, XH(2,:), '-', 'Color', co(2,:), 'LineWidth', 2);
ylabel(ax, '$x_2$');
legend(ax, {'$\pm 3\sigma$ (de $P_{k|k}$)', 'real', 'estimado $\hat{x}_{k|k}$'}, 'Location', 'southeast', 'FontSize', 12);
title(ax, '$x_2$ (n\~ao medido: estimado pelo filtro)', 'FontSize', 15);

ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on'); ax.FontSize = 13;
plot(ax, k, KK(1,:), '-', 'Color', co(4,:), 'LineWidth', 2);
plot(ax, k, KK(2,:), '-', 'Color', co(5,:), 'LineWidth', 2);
yline(ax, Kinf(1), '--', 'Color', co(4,:), 'LineWidth', 1.5);
yline(ax, Kinf(2), '--', 'Color', co(5,:), 'LineWidth', 1.5);
xlabel(ax, '$k$');  ylabel(ax, '$K_f$');
legend(ax, {'$K_f(1)$', '$K_f(2)$', 'horizonte infinito'}, 'Location', 'northeast', 'FontSize', 12);
title(ax, 'Ganho de Kalman: converge para o de horizonte infinito', 'FontSize', 15);

exportgraphics(fig, fullfile(fileparts(mfilename('fullpath')), '..', 'images', 'aula10_kalman_exemplo.png'), ...
               'Resolution', 150);
