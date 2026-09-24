function anim_lyapunov_energia(caso)
%ANIM_LYAPUNOV_ENERGIA Animação da função de energia (Lyapunov) e do espaço de fase
% SEL0359 - Controle Digital - Aula 7
%
%   anim_lyapunov_energia('estavel')      -> V(x_{k+1}) - V(x_k) = 0   (estável, conservativo)
%   anim_lyapunov_energia('dissipativo')  -> V(x_{k+1}) - V(x_k) <= 0  (estável, V -> V_inf > 0)
%   anim_lyapunov_energia('assintotico')  -> V(x_{k+1}) - V(x_k) < 0   (assint. estável)
%   anim_lyapunov_energia('instavel')     -> V(x_{k+1}) - V(x_k) > 0   (instável, ganha energia)
%   anim_lyapunov_energia                 -> gera os quatro vídeos
%
% Sistemas (Ts = 0.1 s, x_{k+1} = Ad x_k, Ad = expm(A Ts)):
%   estável:       A = [0 1; -4 0]   (oscilador sem amortecimento, |lambda(Ad)| = 1)
%   dissipativo:   A = [0 1; 0 -1]   (massa com atrito: posição x1, velocidade x2;
%                                     lambda(Ad) = {1, e^-0.1} -> para em x1 ~= 0)
%   assintótico:   A = [-1 1; 0 -2]  (exemplo da aula, |lambda(Ad)| < 1)
% Função de energia: V(x) = x' P x
%
% Gera: ../images/lyapunov_{estavel,dissipativo,assintotico}.mp4

if nargin < 1
    anim_lyapunov_energia('estavel');
    anim_lyapunov_energia('dissipativo');
    anim_lyapunov_energia('assintotico');
    anim_lyapunov_energia('instavel');
    return
end

Ts = 0.1;
eqLine = false;                                   % reta de equilíbrios (caso dissipativo)
switch lower(caso)
    case 'estavel'
        A  = [0 1; -4 0];
        Ad = expm(A*Ts);
        P  = diag([4 1]);                        % Ad'*P*Ad - P = 0
        X0 = [3 4; 1.5 -3; -1 1.5]';             % condições iniciais
        L  = 8;                                   % limites dos eixos
        N  = 64;                                  % ~2 voltas
        ttlV = '$V(x_{k+1})-V(x_k)=0$ \ (conservativo)';
        arq  = 'lyapunov_estavel.mp4';
    case 'dissipativo'
        A  = [0 1; 0 -1];
        Ad = expm(A*Ts);
        P  = [1 1; 1 2];                          % V = (x1+x2)^2 + x2^2
        % Ad'*P*Ad - P = -(1 - e^{-2Ts}) e2 e2'  <= 0  (semidefinida)
        X0 = [-3 5; 2 3; 0 -4; 4 -2]';
        N  = 50;
        L  = 12;
        eqLine = true;
        ttlV = '$V(x_{k+1})-V(x_k)\le 0$ \ ($V\to V_\infty>0$)';
        arq  = 'lyapunov_dissipativo.mp4';
    case 'assintotico'
        A  = [-1 1; 0 -2];
        Ad = expm(A*Ts);
        Q  = eye(2);                              % P - Ad'*P*Ad = Q
        P  = reshape((eye(4) - kron(Ad', Ad')) \ Q(:), 2, 2);
        th = linspace(0, 2*pi, 9); th(end) = [];
        X0 = 10*[cos(th + pi/8); sin(th + pi/8)];
        X0(:,1) = [10; 10];                       % exemplo da aula
        N  = 50;
        L  = 13;
        ttlV = '$V(x_{k+1})-V(x_k)<0$ \ (assint. est\''avel)';
        arq  = 'lyapunov_assintotico.mp4';
    case 'instavel'
        A  = [0.25 2; -2 0.25];                   % espiral divergente
        Ad = expm(A*Ts);                          % |lambda(Ad)| = e^{0.025} > 1
        P  = eye(2);                              % V = x'x
        % Ad'*P*Ad - P = (e^{0.05} - 1) I > 0  -> energia cresce a cada passo
        X0 = [2 0; 0 -1.5; -1 1]';
        N  = 60;                                  % ~2 voltas
        L  = 11;
        ttlV = '$V(x_{k+1})-V(x_k)>0$ \ (inst\''avel)';
        arq  = 'lyapunov_instavel.mp4';
    otherwise
        error('caso deve ser ''estavel'', ''dissipativo'', ''assintotico'' ou ''instavel''.');
end
P = (P + P')/2;
fprintf('[%s] |eig(Ad)| = %s,  eig(P) = %s,  max|Ad''PAd-P| = %.2e\n', caso, ...
    mat2str(abs(eig(Ad)).', 4), mat2str(eig(P).', 4), max(max(abs(Ad'*P*Ad - P))));

V = @(x) sum(x .* (P*x), 1);                     % V(x) = x'Px (colunas)

%% Trajetórias
nIC = size(X0, 2);
X = zeros(2, N+1, nIC);
for i = 1:nIC
    X(:,1,i) = X0(:,i);
    for k = 1:N
        X(:,k+1,i) = Ad*X(:,k,i);
    end
end
Vk = reshape(V(reshape(X, 2, [])), N+1, nIC);

%% Superfície V(x)
[g1, g2] = meshgrid(linspace(-L, L, 201));
Vg   = reshape(V([g1(:) g2(:)]'), size(g1));
Vmax = 1.15*max(Vk(:));                          % (máx. da trajetória: cobre o caso instável)
% "tigela" em grade polar/elíptica: x = sqrt(v) P^{-1/2} [cos t; sin t]
Ph = sqrtm(inv(P));
[rr, tt] = meshgrid(linspace(0, sqrt(Vmax), 60), linspace(0, 2*pi, 121));
e1 = Ph(1,1)*rr.*cos(tt) + Ph(1,2)*rr.*sin(tt);
e2 = Ph(2,1)*rr.*cos(tt) + Ph(2,2)*rr.*sin(tt);
Vs = rr.^2;
lev  = sort(max(Vk(:)) * (0.7.^(0:14)));

cMain = [0.878 0.141 0.141];     % vermelho (#E02424)
cOth  = [0.165 0.494 0.878];     % azul     (#2A7EE0)
cLev  = [0.6 0.6 0.6];

%% Figura
fig = figure('Color', 'w', 'Position', [30 60 1500 560]);
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (1) Superfície de energia 3D
ax1 = nexttile(tl, 1); hold(ax1, 'on');
surf(ax1, e1, e2, Vs, 'EdgeColor', 'none', 'FaceAlpha', 0.55);
tc = linspace(0, 2*pi, 200);
for lv = lev(lev <= Vmax)
    ec = sqrt(lv)*Ph*[cos(tc); sin(tc)];
    plot3(ax1, ec(1,:), ec(2,:), lv*ones(size(tc)), 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5);
end
colormap(ax1, parula); grid(ax1, 'on');
xlabel(ax1, '$x_1$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(ax1, '$x_2$', 'Interpreter', 'latex', 'FontSize', 16);
zlabel(ax1, '$V(x)=x^\top P x$', 'Interpreter', 'latex', 'FontSize', 16);
title(ax1, 'Fun\c{c}\~ao de energia $V(x)$', 'Interpreter', 'latex', 'FontSize', 17);
xlim(ax1, [-L L]); ylim(ax1, [-L L]); zlim(ax1, [0 Vmax]);

% (2) Espaço de fase
ax2 = nexttile(tl, 2); hold(ax2, 'on'); axis(ax2, 'equal'); box(ax2, 'on');
contour(ax2, g1, g2, Vg, lev, 'Color', cLev, 'LineWidth', 0.6);
[q1, q2] = meshgrid(linspace(-L, L, 15));
dq = Ad*[q1(:) q2(:)]' - [q1(:) q2(:)]';
quiver(ax2, q1(:), q2(:), dq(1,:)', dq(2,:)', 1.2, 'Color', [0.78 0.78 0.78]);
plot(ax2, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
if eqLine   % conjunto de equilíbrios: x2 = 0 (em 2D e sobre a superfície)
    plot(ax2, [-L L], [0 0], '--', 'Color', [0.122 0.620 0.122], 'LineWidth', 1.8);
    text(ax2, L-0.5, 0.9, 'equil\''ibrios ($x_2=0$)', 'Interpreter', 'latex', ...
         'FontSize', 13, 'Color', [0.122 0.620 0.122], 'HorizontalAlignment', 'right');
    xe = linspace(-sqrt(Vmax/P(1,1)), sqrt(Vmax/P(1,1)), 100);
    plot3(ax1, xe, 0*xe, P(1,1)*xe.^2 + 0.01*Vmax, '--', 'Color', [0.122 0.620 0.122], 'LineWidth', 1.8);
end
xlabel(ax2, '$x_1$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel(ax2, '$x_2$', 'Interpreter', 'latex', 'FontSize', 15);
title(ax2, 'Espa\c{c}o de fase: $x_{k+1}=A_d x_k$', 'Interpreter', 'latex', 'FontSize', 16);
xlim(ax2, [-L L]); ylim(ax2, [-L L]);

% (3) V(x_k) ao longo de k
ax3 = nexttile(tl, 3); hold(ax3, 'on'); grid(ax3, 'on'); box(ax3, 'on');
xlabel(ax3, '$k$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel(ax3, '$V(x_k)$', 'Interpreter', 'latex', 'FontSize', 15);
title(ax3, ttlV, 'Interpreter', 'latex', 'FontSize', 16);
xlim(ax3, [0 N]); ylim(ax3, [0 1.1*max(Vk(:))]);
if eqLine   % nível de energia final da trajetória destacada
    plot(ax3, [0 N], Vk(end,1)*[1 1], '--', 'Color', [0.122 0.620 0.122], 'LineWidth', 1.6);
    text(ax3, N, Vk(end,1), sprintf('$V_\\infty=%.1f$ ', Vk(end,1)), 'Interpreter', 'latex', ...
         'FontSize', 14, 'Color', [0.122 0.620 0.122], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end

% Objetos animados
[hTr3, hPt3, hTr2, hPt2, hV] = deal(gobjects(nIC,1));
for i = 1:nIC
    if i == 1, c = cMain; lw = 2.6; ms = 9; else, c = cOth; lw = 1.4; ms = 6; end
    hTr3(i) = plot3(ax1, nan, nan, nan, '-', 'Color', c, 'LineWidth', lw);
    hPt3(i) = plot3(ax1, nan, nan, nan, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'MarkerSize', ms);
    hTr2(i) = plot(ax2, nan, nan, '.-', 'Color', c, 'LineWidth', lw, 'MarkerSize', 8);
    hPt2(i) = plot(ax2, nan, nan, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'MarkerSize', ms);
    if i == 1
        hV(i) = stem(ax3, nan, nan, 'filled', 'Color', c, 'LineWidth', 1.2, 'MarkerSize', 4);
    else
        hV(i) = plot(ax3, nan, nan, '-', 'Color', [cOth 0.6], 'LineWidth', 1.2);
    end
end
hK = text(ax2, -L+0.8, -L+1.6, '', 'FontSize', 14, 'Interpreter', 'latex', ...
          'BackgroundColor', 'w', 'Margin', 2);

%% Vídeo
outFile = fullfile(fileparts(mfilename('fullpath')), '..', 'images', arq);
vw = VideoWriter(outFile, 'MPEG-4');
vw.FrameRate = 8; vw.Quality = 100;
open(vw);
% Alta resolução: renderiza cada quadro com print (DPI definido) em vez de getframe
dpi   = 192;                                      % 192 dpi -> ~2x a resolução da tela
fig.InvertHardcopy = 'off';                       % mantém fundo branco da figura
szFr  = [];                                       % tamanho fixo (par) dos quadros
dz = 0.01*Vmax;                                   % leve offset p/ ficar sobre a superfície
for k = 0:N
    idx = 1:k+1;
    for i = 1:nIC
        x1 = X(1,idx,i); x2 = X(2,idx,i); vv = Vk(idx,i)';
        set(hTr3(i), 'XData', x1, 'YData', x2, 'ZData', vv + dz);
        set(hPt3(i), 'XData', x1(end), 'YData', x2(end), 'ZData', vv(end) + dz);
        set(hTr2(i), 'XData', x1, 'YData', x2);
        set(hPt2(i), 'XData', x1(end), 'YData', x2(end));
        set(hV(i),   'XData', idx-1, 'YData', vv);
    end
    set(hK, 'String', sprintf('$k=%d,\\;V(x_k)=%.1f$', k, Vk(k+1,1)));
    view(ax1, -35 + 0.4*k, 30);                   % leve rotação da superfície
    drawnow;
    img = print(fig, '-RGBImage', sprintf('-r%d', dpi), '-opengl');
    if isempty(szFr), szFr = 2*floor(size(img, [1 2])/2); end   % H.264 exige dimensões pares
    tmp = 255*ones([szFr 3], 'uint8');            % garante tamanho idêntico em todos os quadros
    h = min(szFr(1), size(img,1)); w = min(szFr(2), size(img,2));
    tmp(1:h, 1:w, :) = img(1:h, 1:w, :);  img = tmp;
    fr  = im2frame(img);
    nRep = 1 + 12*(k==0) + 20*(k==N);             % pausa no início e no fim
    for r = 1:nRep, writeVideo(vw, fr); end
end
close(vw); close(fig);
fprintf('Vídeo salvo em: %s\n', outFile);
end
