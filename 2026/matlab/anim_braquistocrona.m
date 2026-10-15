%% ANIM_BRAQUISTOCRONA  Animação (GIF) do problema da braquistócrona
% Aula 9 (SEL0359 - Controle Digital, 2026) - Prof. Marcos R. Fernandes
%
% Uma bolinha parte do repouso em A e desliza sem atrito até B sob a ação
% da gravidade. Compara trajetórias subótimas (reta, parábola, curva
% profunda) com a solução ótima: a cicloide (Johann Bernoulli, 1696).
%
% Tempo de percurso ao longo de uma curva:  T = int ds / v,  v = sqrt(2 g |y|)
% Gera: ../images/braquistocrona.gif

clc; clear all; close all;
g  = 9.81;                    % [m/s^2]
xB = 4;  yB = -1.5;           % A = (0,0), B = (4,-1.5) [m]

%% Cicloide por A e B:  x = r(th - sin th),  y = -r(1 - cos th)
thB = fzero(@(t) (t - sin(t))./(1 - cos(t)) - xB/(-yB), [0.5 2*pi-0.01]);
r   = -yB/(1 - cos(thB));
Tcic = sqrt(r/g)*thB;         % tempo exato na cicloide

%% Curvas candidatas (todas de A a B)
s  = linspace(0, 1, 6000).^2;                 % malha refinada perto de A (v = 0)
C(1).nome = 'Reta';            C(1).x = xB*s;  C(1).y = yB*s;
C(2).nome = 'Par\''abola';     C(2).x = xB*s;  C(2).y = yB*(2*s - s.^2);
C(3).nome = 'Curva profunda';  C(3).x = xB*s;  C(3).y = yB*s - 2.5*sin(pi*s);
th = thB*s;
C(4).nome = 'Cicloide (\''otima)'; C(4).x = r*(th - sin(th)); C(4).y = -r*(1 - cos(th));
cor = [0.50 0.50 0.50; 0.165 0.494 0.878; 0.93 0.55 0.10; 0.122 0.620 0.122];
lw  = [2 2 2 3.5];

for i = 1:numel(C)
    v  = sqrt(2*g*max(-C(i).y, 0));
    ds = hypot(diff(C(i).x), diff(C(i).y));
    vm = (v(1:end-1) + v(2:end))/2;
    dt = zeros(size(ds)); ok = vm > 0; dt(ok) = ds(ok)./vm(ok);
    C(i).t = [0 cumsum(dt)];                   % tempo acumulado ao longo da curva
    C(i).T = C(i).t(end);
    fprintf('%-22s T = %.3f s\n', C(i).nome, C(i).T);
end
fprintf('Cicloide (exato): T = sqrt(r/g)*thB = %.3f s  (r = %.4f m, thB = %.4f rad)\n', Tcic, r, thB);

%% Figura
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');
fig = figure('Color', 'w', 'Position', [50 50 1200 560]);
tl  = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
ax  = nexttile(tl, [1 3]); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
ax.FontSize = 14;
hC = gobjects(1, numel(C)); hB = hC;
for i = 1:numel(C)
    hC(i) = plot(ax, C(i).x, C(i).y, '-', 'Color', cor(i,:), 'LineWidth', lw(i));
end
plot(ax, [0 xB], [0 yB], 'ks', 'MarkerSize', 11, 'MarkerFaceColor', 'k');
text(ax, -0.08, 0.18, '$A$', 'FontSize', 20, 'HorizontalAlignment', 'right');
text(ax, xB + 0.08, yB + 0.2, '$B$', 'FontSize', 20);
for i = 1:numel(C)
    hB(i) = plot(ax, 0, 0, 'o', 'MarkerSize', 15, 'MarkerFaceColor', cor(i,:), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
end
quiver(ax, 4.45, -0.4, 0, -0.7, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.8);
text(ax, 4.55, -0.75, '$\vec g$', 'FontSize', 18);
axis(ax, 'equal'); xlim(ax, [-0.4 4.9]); ylim(ax, [-3.5 0.5]);
xlabel(ax, '$x$ [m]'); ylabel(ax, '$y$ [m]');
legend(ax, hC, {C.nome}, 'Location', 'southwest', 'FontSize', 14);
title(ax, 'Braquist\''ocrona: qual curva leva a bolinha de $A$ a $B$ no menor tempo?', 'FontSize', 16);
hT = text(ax, 2.0, 0.3, '', 'FontSize', 16, 'HorizontalAlignment', 'center');

% corrida: barras de tempo
ax2 = nexttile(tl); hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on'); ax2.FontSize = 13;
hBar = gobjects(1, numel(C)); hLab = hBar;
for i = 1:numel(C)
    yy = numel(C) + 1 - i;
    hBar(i) = patch(ax2, [0 0 0 0], yy + [-0.32 -0.32 0.32 0.32], cor(i,:), 'EdgeColor', 'none');
    hLab(i) = text(ax2, 0, yy, '', 'FontSize', 13, 'VerticalAlignment', 'middle');
end
xline(ax2, Tcic, '--', 'Color', cor(4,:), 'LineWidth', 1.5);
xlim(ax2, [0 2.4]); ylim(ax2, [0.4 numel(C) + 0.6]);
yticks(ax2, 1:numel(C)); yticklabels(ax2, fliplr({'Reta', 'Par\''abola', 'Profunda', 'Cicloide'}));
xlabel(ax2, 'tempo [s]'); title(ax2, 'Tempo de percurso', 'FontSize', 15);

%% Animação -> GIF
outFile = fullfile(fileparts(mfilename('fullpath')), '..', 'images', 'braquistocrona.gif');
dpi = 144;  fig.InvertHardcopy = 'off';
dtSim = 0.01;  dtGif = 0.04;                  % câmera lenta (4x)
Tmax  = max([C.T]);
tk    = [0:dtSim:Tmax, Tmax];
szFr  = [];
for k = 1:numel(tk)
    t = tk(k);
    for i = 1:numel(C)
        tc = min(t, C(i).T);
        [tu, iu] = unique(C(i).t);
        set(hB(i), 'XData', interp1(tu, C(i).x(iu), tc), 'YData', interp1(tu, C(i).y(iu), tc));
        set(hBar(i), 'XData', [0 tc tc 0]);
        if t >= C(i).T
            set(hLab(i), 'Position', [tc + 0.03, numel(C) + 1 - i], 'String', sprintf('$%.3f$ s', C(i).T));
        end
    end
    set(hT, 'String', sprintf('$t = %.3f$ s \\quad (c\\^amera lenta $4\\times$)', t));
    drawnow;

    img = print(fig, '-RGBImage', sprintf('-r%d', dpi), '-opengl');
    if isempty(szFr), szFr = size(img, [1 2]); end
    tmp = 255*ones([szFr 3], 'uint8');
    h = min(szFr(1), size(img,1)); wd = min(szFr(2), size(img,2));
    tmp(1:h, 1:wd, :) = img(1:h, 1:wd, :);
    [A, map] = rgb2ind(tmp, 128, 'nodither');
    dly = dtGif + 1.0*(k==1) + 3.0*(k==numel(tk));   % pausa no início e no fim
    if k == 1
        imwrite(A, map, outFile, 'gif', 'LoopCount', Inf, 'DelayTime', dly);
    else
        imwrite(A, map, outFile, 'gif', 'WriteMode', 'append', 'DelayTime', dly);
    end
end
close(fig);
fprintf('%d quadros. GIF salvo em: %s\n', numel(tk), outFile);
