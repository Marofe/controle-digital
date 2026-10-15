%% ANIM_LQR_PID_PENDULO  Comparativo LQR x PID no pêndulo invertido sobre carro (GIF)
% Aula 9 (SEL0359 - Controle Digital, 2026) - Prof. Marcos R. Fernandes
%
% Planta não linear (carro de massa M, pêndulo de massa m na ponta de uma haste l,
% força F no carro, atrito viscoso c; theta medido a partir da vertical para cima):
%   (M + m - m cos^2 th) xdd = F - c xd + m l thd^2 sin th - m g sin th cos th
%   thdd = (g sin th - xdd cos th)/l
% Estado: s = [x; xd; th; thd].  Controle digital com Ts = 20 ms (ZOH).
%   LQR: modelo linearizado discretizado (c2d) + dlqr,   F_k = -K s_k
%   PID: PID no ângulo + PD na posição do carro (sintonia propositalmente
%        pouco amortecida):  F_k = Kp e + Ki Ts sum(e) + Kd (e_k-e_{k-1})/Ts + kx x + kv xd
% Gera: ../images/pendulo_carro_lqr_vs_pid.gif

clc; clear all; close all;
M = 1.0;  m = 0.3;  l = 0.5;  c = 0.1;  g = 9.81;  Ts = 0.02;  T = 8;
s0 = [0; 0; 0.25; 0];                            % pêndulo inclinado 0.25 rad (~14 graus)

%% Projeto
A = [0 1 0 0; 0 -c/M -m*g/M 0; 0 0 0 1; 0 c/(M*l) (g/l)*(1 + m/M) 0];
B = [0; 1/M; 0; -1/(M*l)];
sysd = c2d(ss(A, B, eye(4), 0), Ts, 'zoh');
Q = diag([10 1 50 1]);  R = 0.1;
K = dlqr(sysd.A, sysd.B, Q, R);
Kp = 40;  Ki = 10;  Kd = 3;  kx = 2;  kv = 3;    % PID + PD (pouco amortecido)
fprintf('K_LQR = [%.2f %.2f %.2f %.2f]\n', K);

%% Simulação não linear (força constante entre amostras)
f  = @(s, F) cartpole(s, F, M, m, l, c, g);
op = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
N  = round(T/Ts);  tk = (0:N)*Ts;
SL = zeros(4, N+1);  SP = SL;  UL = zeros(1, N+1);  UP = UL;
SL(:,1) = s0;  SP(:,1) = s0;  I = 0;  ep = s0(3);
for k = 1:N
    UL(k) = -K*SL(:,k);
    e = SP(3,k);  I = I + e*Ts;
    UP(k) = Kp*e + Ki*I + Kd*(e - ep)/Ts + kx*SP(1,k) + kv*SP(2,k);  ep = e;
    [~, z] = ode45(@(t,s) f(s, UL(k)), [0 Ts], SL(:,k), op);  SL(:,k+1) = z(end,:)';
    [~, z] = ode45(@(t,s) f(s, UP(k)), [0 Ts], SP(:,k), op);  SP(:,k+1) = z(end,:)';
end
UL(end) = UL(end-1);  UP(end) = UP(end-1);
tsf = @(th) tk(find(abs(th) > 0.02, 1, 'last') + 1);
fprintf('LQR: ts(theta) = %.2f s, min theta = %.3f rad, max|F| = %.1f N\n', tsf(SL(3,:)), min(SL(3,:)), max(abs(UL)));
fprintf('PID: ts(theta) = %.2f s, min theta = %.3f rad, max|F| = %.1f N\n', tsf(SP(3,:)), min(SP(3,:)), max(abs(UP)));

%% Figura
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');
co = get(groot, 'defaultAxesColorOrder');        % cores padrão do MATLAB
cL = co(2,:);  cP = co(4,:);  cG = [0.6 0.6 0.6];  % laranja (LQR), roxo (PID)
fig = figure('Color', 'w', 'Position', [30 40 1400 800]);
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

wC = 0.28;  hC = 0.14;  rW = 0.035;               % dimensões do carro (desenho)
nm = {'LQR: $F_k=-Ks_k$', 'PID: $\theta$ (PID) + posi\c{c}\~ao do carro (PD)'};
cc = {cL, cP};  S = {SL, SP};
hCart = gobjects(1,2);  hW1 = hCart;  hW2 = hCart;  hRod = hCart;  hBob = hCart;  hTr = hCart;  hTxt = hCart;
for i = 1:2
    ax = nexttile(tl, i); hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off');
    xlim(ax, [-0.75 1.05]); ylim(ax, [-0.12 0.82]);
    plot(ax, [-0.75 1.05], [0 0], 'k-', 'LineWidth', 2);                       % trilho
    for xm = -0.5:0.25:1.0, plot(ax, [xm xm], [-0.03 0], 'k-', 'LineWidth', 1); end
    plot(ax, [0 0], [-0.07 -0.03], 'k-', 'LineWidth', 2.5);                    % referência x = 0
    text(ax, 0, -0.1, '$x=0$', 'FontSize', 13, 'HorizontalAlignment', 'center');
    hTr(i)   = plot(ax, nan, nan, '-', 'Color', [cc{i} 0.35], 'LineWidth', 2);
    hCart(i) = patch(ax, nan, nan, cc{i}, 'FaceAlpha', 0.85, 'EdgeColor', 'k', 'LineWidth', 1.5);
    hW1(i)   = plot(ax, nan, nan, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'k');
    hW2(i)   = plot(ax, nan, nan, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'k');
    hRod(i)  = plot(ax, nan, nan, 'k-', 'LineWidth', 4);
    hBob(i)  = plot(ax, nan, nan, 'o', 'MarkerSize', 24, 'MarkerFaceColor', cc{i}, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    hTxt(i)  = text(ax, 0.55, 0.74, '', 'FontSize', 16);
    title(ax, nm{i}, 'FontSize', 18, 'Color', cc{i});
end

ax3 = nexttile(tl, 3); hold(ax3, 'on'); box(ax3, 'on'); grid(ax3, 'on'); ax3.FontSize = 14;
yline(ax3, 0, 'k:');
h3L = plot(ax3, nan, nan, '-', 'Color', cL, 'LineWidth', 2.5);
h3P = plot(ax3, nan, nan, '-', 'Color', cP, 'LineWidth', 2.5);
xlim(ax3, [0 T]); ylim(ax3, [-0.3 0.3]);
xlabel(ax3, '$t$ [s]'); ylabel(ax3, '$\theta$ [rad]');
legend(ax3, [h3L h3P], {'LQR', 'PID'}, 'Location', 'northeast', 'FontSize', 15);
title(ax3, '\^Angulo do p\^endulo', 'FontSize', 17);

ax4 = nexttile(tl, 4); hold(ax4, 'on'); box(ax4, 'on'); grid(ax4, 'on'); ax4.FontSize = 14;
yline(ax4, 0, 'k:');
h4L = plot(ax4, nan, nan, '-', 'Color', cL, 'LineWidth', 2.5);
h4P = plot(ax4, nan, nan, '-', 'Color', cP, 'LineWidth', 2.5);
xlim(ax4, [0 T]); ylim(ax4, [-0.2 0.45]);
xlabel(ax4, '$t$ [s]'); ylabel(ax4, '$x$ [m]');
legend(ax4, [h4L h4P], {'LQR', 'PID'}, 'Location', 'northeast', 'FontSize', 15);
title(ax4, 'Posi\c{c}\~ao do carro', 'FontSize', 17);

%% Animação -> GIF (câmera lenta 2x)
outFile = fullfile(fileparts(mfilename('fullpath')), '..', 'images', 'pendulo_carro_lqr_vs_pid.gif');
dpi = 110;  fig.InvertHardcopy = 'off';
kk = 1:2:N+1;  szFr = [];  nTr = 30;
for n = 1:numel(kk)
    k = kk(n);
    for i = 1:2
        x = S{i}(1,k);  th = S{i}(3,k);
        yc = rW*2 + hC/2;                                  % centro do carro
        set(hCart(i), 'XData', x + wC/2*[-1 1 1 -1], 'YData', rW*2 + [0 0 hC hC]);
        set(hW1(i), 'XData', x - wC/3, 'YData', rW);  set(hW2(i), 'XData', x + wC/3, 'YData', rW);
        px = x + l*sin(th);  py = yc + l*cos(th);
        set(hRod(i), 'XData', [x px], 'YData', [yc py]);
        set(hBob(i), 'XData', px, 'YData', py);
        j0 = max(1, k - nTr);
        set(hTr(i), 'XData', S{i}(1,j0:k) + l*sin(S{i}(3,j0:k)), 'YData', yc + l*cos(S{i}(3,j0:k)));
        set(hTxt(i), 'String', sprintf('$\\theta=%+5.1f^\\circ$', th*180/pi));
    end
    set(h3L, 'XData', tk(1:k), 'YData', SL(3,1:k));  set(h3P, 'XData', tk(1:k), 'YData', SP(3,1:k));
    set(h4L, 'XData', tk(1:k), 'YData', SL(1,1:k));  set(h4P, 'XData', tk(1:k), 'YData', SP(1,1:k));
    drawnow;

    img = print(fig, '-RGBImage', sprintf('-r%d', dpi), '-opengl');
    if isempty(szFr), szFr = size(img, [1 2]); end
    tmp = 255*ones([szFr 3], 'uint8');
    h = min(szFr(1), size(img,1)); wd = min(szFr(2), size(img,2));
    tmp(1:h, 1:wd, :) = img(1:h, 1:wd, :);
    [Ai, map] = rgb2ind(tmp, 128, 'nodither');
    dly = 0.08 + 1.0*(n==1) + 3.0*(n==numel(kk));
    if n == 1
        imwrite(Ai, map, outFile, 'gif', 'LoopCount', Inf, 'DelayTime', dly);
    else
        imwrite(Ai, map, outFile, 'gif', 'WriteMode', 'append', 'DelayTime', dly);
    end
end
close(fig);
fprintf('%d quadros. GIF salvo em: %s\n', numel(kk), outFile);

function ds = cartpole(s, F, M, m, l, c, g)
    v = s(2);  th = s(3);  w = s(4);  S = sin(th);  C = cos(th);
    a  = (F - c*v + m*l*w^2*S - m*g*S*C) / (M + m - m*C^2);
    ds = [v; a; w; (g*S - a*C)/l];
end
