%% ANIM_CUSTO_QUADRATICO  GIF girando 360 graus da função custo quadrática J(u0,u1)
% Aula 9 (SEL0359 - Controle Digital, 2026) - Prof. Marcos R. Fernandes
%
% Sistema escalar  x_{k+1} = a x_k + b u_k,  horizonte N = 2:
%   J = x_2' P_N x_2 + sum_{k=0}^{1} ( x_k' Q x_k + u_k' R u_k )
% Substituindo x_1 = a x0 + b u0 e x_2 = a x_1 + b u1, J é uma função
% QUADRÁTICA de (u0,u1): um paraboloide convexo com mínimo único,
% que coincide com a solução LQR obtida pela recursão de Riccati.
% Gera: ../images/custo_quadratico_3d.gif

clc; clear all; close all;
a = 1.2;  b = 1;  x0 = 1;          % sistema (instável em malha aberta)
Q = 1;    R = 0.5;  PN = 1;  N = 2;

Jfun = @(u0,u1) Q*x0.^2 + R*u0.^2 + Q*(a*x0 + b*u0).^2 + R*u1.^2 ...
              + PN*(a*(a*x0 + b*u0) + b*u1).^2;

%% Solução LQR (Riccati para trás) -> mínimo de J
P = PN;  K = zeros(1,N);
for k = N-1:-1:0
    K(k+1) = (R + b*P*b) \ (b*P*a);
    P = Q + a*P*a - a*P*b*K(k+1);
end
u0s = -K(1)*x0;  x1 = a*x0 + b*u0s;  u1s = -K(2)*x1;
Js  = Jfun(u0s, u1s);
fprintf('u0* = %.4f, u1* = %.4f, J* = %.4f  (x0''P0x0 = %.4f)\n', u0s, u1s, Js, P*x0^2);

%% Superfície
[U0, U1] = meshgrid(linspace(u0s-2, u0s+2, 81), linspace(u1s-2, u1s+2, 81));
JJ = Jfun(U0, U1);
zf = -10;                                     % plano das curvas de nível (abaixo da superfície)

set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex');
fig = figure('Color', 'w', 'Position', [50 50 1000 720]);
tl  = tiledlayout(fig, 1, 1, 'Padding', 'compact');
ax  = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'off');
ax.FontSize = 15;
surf(ax, U0, U1, JJ, 'EdgeColor', 'none', 'FaceAlpha', 0.92, 'FaceLighting', 'gouraud');
mesh(ax, U0(1:5:end,1:5:end), U1(1:5:end,1:5:end), JJ(1:5:end,1:5:end), ...
     'EdgeColor', [0.25 0.25 0.25], 'EdgeAlpha', 0.35, 'FaceColor', 'none');
[~, hc] = contour(ax, U0, U1, JJ, 18, 'LineWidth', 1.1);
hc.ContourZLevel = zf;
colormap(ax, turbo);  clim(ax, [Js max(JJ(:))*0.8]);
% mínimo e projeção
plot3(ax, u0s, u1s, Js, 'p', 'MarkerSize', 22, 'MarkerFaceColor', [0.85 0.1 0.1], 'MarkerEdgeColor', 'k');
plot3(ax, [u0s u0s], [u1s u1s], [zf Js], 'k--', 'LineWidth', 1.5);
plot3(ax, u0s, u1s, zf, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.85 0.1 0.1], 'MarkerEdgeColor', 'k');
text(ax, u0s + 0.25, u1s, Js + 0.15*max(JJ(:)), sprintf('$J^*=%.3f$', Js), 'FontSize', 18, ...
     'HorizontalAlignment', 'left', 'Color', [0.7 0 0]);
light(ax, 'Position', [0 0 1], 'Style', 'infinite');  lighting(ax, 'gouraud');  material(ax, 'dull');
xlabel(ax, '$u_0$', 'FontSize', 18); ylabel(ax, '$u_1$', 'FontSize', 18);
zlabel(ax, '$J(u_0,u_1)$', 'FontSize', 18);
xlim(ax, [u0s-2 u0s+2]); ylim(ax, [u1s-2 u1s+2]); zlim(ax, [zf max(JJ(:))]);
% (sem título: o slide já traz o título)
zticks(ax, 0:5:30);
ax.Projection = 'orthographic';
pbaspect(ax, [1 1 0.8]);
view(ax, -37.5, 28);  axis(ax, 'vis3d');         % tamanho fixo durante a rotação
camva(ax, 1.03*camva(ax));                         % enquadramento (sem cortar durante a rotação)

%% Rotação 360 graus -> GIF
outFile = fullfile(fileparts(mfilename('fullpath')), '..', 'images', 'custo_quadratico_3d.gif');
dpi = 120;  fig.InvertHardcopy = 'off';
el  = 28;  azs = -37.5 + (0:2:358);         % passo de 2 graus
szFr = [];
for k = 1:numel(azs)
    view(ax, azs(k), el);
    drawnow;
    img = print(fig, '-RGBImage', sprintf('-r%d', dpi), '-opengl');
    if isempty(szFr), szFr = size(img, [1 2]); end
    tmp = 255*ones([szFr 3], 'uint8');
    h = min(szFr(1), size(img,1)); wd = min(szFr(2), size(img,2));
    tmp(1:h, 1:wd, :) = img(1:h, 1:wd, :);
    [A, map] = rgb2ind(tmp, 200, 'nodither');
    if k == 1
        imwrite(A, map, outFile, 'gif', 'LoopCount', Inf, 'DelayTime', 0.10);
    else
        imwrite(A, map, outFile, 'gif', 'WriteMode', 'append', 'DelayTime', 0.10);
    end
end
close(fig);
fprintf('%d quadros. GIF salvo em: %s\n', numel(azs), outFile);
