function anim_pendulo(caso)
%ANIM_PENDULO Animação (GIF) do pêndulo simples (estável) e do pêndulo invertido
% em malha aberta (instável), com espaço de fase e função de energia.
% SEL0359 - Controle Digital - Aula 7 (Prof. Marcos R. Fernandes)
%
%   anim_pendulo('simples')    -> pêndulo simples com atrito: energia diminui (estável)
%   anim_pendulo('invertido')  -> pêndulo invertido em malha aberta: afasta-se do
%                                 equilíbrio e V cresce (instável)
%   anim_pendulo               -> gera os dois GIFs
%
% Modelo (theta medido a partir do equilíbrio analisado, x = [theta; omega]):
%   simples:    theta'' = -(g/l) sin(theta) - b theta'
%   invertido:  theta'' = +(g/l) sin(theta) - b theta'
% Simulação não linear (ode45), amostrada com período Ts: x_k = x(k Ts).
%
% Gera: ../images/pendulo_simples_estavel.gif e ../images/pendulo_invertido_instavel.gif

if nargin < 1
    anim_pendulo('simples');
    anim_pendulo('invertido');
    return
end

g = 9.81; l = 1;
switch lower(caso)
    case 'simples'
        s   = -1;  b = 0.4;               % sinal da gravidade no modelo, atrito
        x0  = [2.6; 0];                   % grande amplitude inicial
        Ts  = 0.05; T = 14;
        E   = @(th, w) 0.5*l^2*w.^2 + g*l*(1 - cos(th));   % energia mecânica (por massa)
        Vlab = '$V(x)=\frac{1}{2}\ell^2\omega^2+g\ell(1-\cos\theta)$';
        ttlV = '$V(x_{k+1})-V(x_k)\le 0$ \ (est\''avel)';
        thLim = [-pi pi]; wLim = [-7 7];
        arq  = 'pendulo_simples_estavel.gif';
        dtGif = Ts;                       % tempo real
    case 'invertido'
        s   = +1;  b = 0.05;
        x0  = [0.05; 0];                  % pequena perturbação em torno do topo
        Ts  = 0.02; T = 2.2;
        E   = @(th, w) 0.5*w.^2 + 0.5*(g/l)*th.^2;         % V = x'Px, P = diag(g/l, 1)/2
        Vlab = '$V(x)=\frac{1}{2}\omega^2+\frac{1}{2}\frac{g}{\ell}\theta^2$';
        ttlV = '$V(x_{k+1})-V(x_k)>0$ \ (inst\''avel)';
        thLim = [-0.5 3.4]; wLim = [-2 8];
        arq  = 'pendulo_invertido_instavel.gif';
        dtGif = 0.06;                     % câmera lenta (~3x)
    otherwise
        error('caso deve ser ''simples'' ou ''invertido''.');
end

f = @(t, x) [x(2); s*(g/l)*sin(x(1)) - b*x(2)];
tk = (0:Ts:T)';
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);
[tk, X] = ode45(f, tk, x0, opts);
if strcmp(caso, 'invertido')     % para quando o pêndulo chega embaixo (theta = pi)
    kEnd = find(X(:,1) >= pi, 1);
    if ~isempty(kEnd), tk = tk(1:kEnd); X = X(1:kEnd,:); end
end
th = X(:,1); w = X(:,2); Vk = E(th, w); N = numel(tk);

%% Cores
cMain = [0.878 0.141 0.141];     % vermelho
cBlue = [0.165 0.494 0.878];     % azul
cGray = [0.6 0.6 0.6];

%% Figura
fig = figure('Color', 'w', 'Position', [30 60 1500 560]);
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (1) Pêndulo
ax1 = nexttile(tl, 1); hold(ax1, 'on'); axis(ax1, 'equal'); axis(ax1, 'off');
xlim(ax1, [-1.35 1.35]*l); ylim(ax1, [-1.35 1.35]*l);
if s < 0
    bob = @(a) [l*sin(a), -l*cos(a)];
    plot(ax1, [0 0], [0 -1.25*l], ':', 'Color', cGray, 'LineWidth', 1.2);
    fill(ax1, [-0.3 0.3 0.3 -0.3], [0 0 0.08 0.08], [0.75 0.75 0.75], 'EdgeColor', 'k');   % teto
    title(ax1, 'P\^endulo simples', 'Interpreter', 'latex', 'FontSize', 18);
else
    bob = @(a) [l*sin(a), l*cos(a)];
    plot(ax1, [0 0], [0 1.25*l], ':', 'Color', cGray, 'LineWidth', 1.2);
    fill(ax1, [-0.3 0.3 0.3 -0.3], [-0.08 -0.08 0 0], [0.75 0.75 0.75], 'EdgeColor', 'k'); % base
    title(ax1, 'P\^endulo invertido (malha aberta)', 'Interpreter', 'latex', 'FontSize', 18);
end
tc = linspace(0, 2*pi, 200);
plot(ax1, l*cos(tc), l*sin(tc), '-', 'Color', [0.9 0.9 0.9]);
hTrail = plot(ax1, nan, nan, '-', 'Color', [cBlue 0.5], 'LineWidth', 1.5);
hRod   = plot(ax1, [0 nan], [0 nan], 'k-', 'LineWidth', 4);
hBob   = plot(ax1, nan, nan, 'o', 'MarkerSize', 26, 'MarkerFaceColor', cMain, 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
plot(ax1, 0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
hT = text(ax1, -1.3*l, -1.25*l, '', 'Interpreter', 'latex', 'FontSize', 16);

% (2) Espaço de fase
ax2 = nexttile(tl, 2); hold(ax2, 'on'); box(ax2, 'on');
[T1, W1] = meshgrid(linspace(thLim(1), thLim(2), 300), linspace(wLim(1), wLim(2), 300));
Eg = E(T1, W1);
if s < 0
    contour(ax2, T1, W1, Eg, linspace(0, max(Vk), 12), 'Color', cGray, 'LineWidth', 0.7);
else  % curvas de energia mecânica real (separatriz em torno do topo)
    Em = 0.5*W1.^2 + (g/l)*(cos(T1) - 1);
    contour(ax2, T1, W1, Em, linspace(-2*g/l, 2*g/l, 11), 'Color', [0.72 0.72 0.72], 'LineWidth', 0.7);
    contour(ax2, T1, W1, Em, [0 0], 'Color', [0.122 0.620 0.122], 'LineWidth', 1.6, 'LineStyle', '--');
    text(ax2, 1.9, 7.0, 'separatriz', 'Interpreter', 'latex', 'FontSize', 14, 'Color', [0.122 0.620 0.122]);
end
plot(ax2, 0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 1.8);
xlim(ax2, thLim); ylim(ax2, wLim);
xlabel(ax2, '$\theta$ [rad]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel(ax2, '$\omega=\dot\theta$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 15);
title(ax2, 'Espa\c{c}o de fase: $x_k=[\theta_k\;\;\omega_k]^\top$', 'Interpreter', 'latex', 'FontSize', 16);
hTr2 = plot(ax2, nan, nan, '.-', 'Color', cMain, 'LineWidth', 2, 'MarkerSize', 9);
hPt2 = plot(ax2, nan, nan, 'o', 'MarkerFaceColor', cMain, 'MarkerEdgeColor', 'k', 'MarkerSize', 9);

% (3) V(x_k)
ax3 = nexttile(tl, 3); hold(ax3, 'on'); grid(ax3, 'on'); box(ax3, 'on');
xlim(ax3, [0 tk(end)]); ylim(ax3, [0 1.1*max(Vk)]);
xlabel(ax3, '$t_k=kT_s$ [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel(ax3, Vlab, 'Interpreter', 'latex', 'FontSize', 15);
title(ax3, ttlV, 'Interpreter', 'latex', 'FontSize', 16);
hV = plot(ax3, nan, nan, '.-', 'Color', cMain, 'LineWidth', 1.5, 'MarkerSize', 9);

%% GIF (alta resolução, cores sem dithering)
outFile = fullfile(fileparts(mfilename('fullpath')), '..', 'images', arq);
dpi = 144;                                  % 1.5x a resolução da tela
fig.InvertHardcopy = 'off';
szFr = [];
nTrail = round(1.0/Ts);                     % rastro de ~1 s
for k = 1:N
    p = bob(th(k));
    set(hRod, 'XData', [0 p(1)], 'YData', [0 p(2)]);
    set(hBob, 'XData', p(1), 'YData', p(2));
    i0 = max(1, k - nTrail); pt = bob(th(i0:k));
    set(hTrail, 'XData', pt(:,1), 'YData', pt(:,2));
    set(hT, 'String', sprintf('$t=%.2f$ s, $\\theta=%.2f$ rad', tk(k), th(k)));
    set(hTr2, 'XData', th(1:k), 'YData', w(1:k));
    set(hPt2, 'XData', th(k), 'YData', w(k));
    set(hV, 'XData', tk(1:k), 'YData', Vk(1:k));
    drawnow;

    img = print(fig, '-RGBImage', sprintf('-r%d', dpi), '-opengl');
    if isempty(szFr), szFr = size(img, [1 2]); end
    tmp = 255*ones([szFr 3], 'uint8');
    h = min(szFr(1), size(img,1)); wd = min(szFr(2), size(img,2));
    tmp(1:h, 1:wd, :) = img(1:h, 1:wd, :);
    [A, map] = rgb2ind(tmp, 128, 'nodither');

    dly = dtGif + 1.0*(k==1) + 2.0*(k==N);  % pausa no início e no fim
    if k == 1
        imwrite(A, map, outFile, 'gif', 'LoopCount', Inf, 'DelayTime', dly);
    else
        imwrite(A, map, outFile, 'gif', 'WriteMode', 'append', 'DelayTime', dly);
    end
end
close(fig);
fprintf('[%s] %d quadros, V: %.2f -> %.2f. GIF salvo em: %s\n', caso, N, Vk(1), Vk(end), outFile);
end
