%% ANIM_KALMAN_APOLLO  Filtro de Kalman na navegação Terra -> Lua (estilo Apollo) - GIF
% Aula 10 (SEL0359 - Controle Digital, 2026) - Prof. Marcos R. Fernandes
%
% Ilustração DIDÁTICA (fora de escala) da navegação de uma nave Apollo:
%   1) órbita de estacionamento na Terra;
%   2) TLI (Trans-Lunar Injection): queima que envia a nave para a Lua;
%   3) viagem translunar (coast);
%   4) MCC (Midcourse Correction): correção de meio curso;
%   5) LOI (Lunar Orbit Insertion): queima atrás da Lua para entrar em órbita;
%   6) órbita lunar.
%
% Modelo (plano, cônicas "remendadas"/patched conics, Terra e Lua fixas):
%   estado     x = [r; v] em R^4  (posição e velocidade no plano)
%   dinâmica   x_{k+1} = f(x_k) + B u_k + w_k ,  w_k ~ N(0, Q)
%              f: gravidade da Terra OU da Lua (esfera de influência da Lua)
%              u_k: queimas dos motores (TLI, MCC e LOI), conhecidas pelo computador
%   medição    y_k = H x_k + eps_k ,  H = [I 0],  eps_k ~ N(0, R_k)
%              R_k cresce com a distância à Terra (rastreio por radar);
%              sem medição atrás da Lua (perda de sinal) e no "descanso".
% Estimadores comparados:
%   - Filtro de Kalman (estendido: f não linear, linearizada a cada passo);
%   - Só predição (dead reckoning): integra o modelo sem usar medições.
% As queimas TLI/MCC/LOI são calculadas pelo computador de bordo a partir do
% estado ESTIMADO pelo filtro de Kalman.
%
% Imagens reais: Terra - foto "Blue Marble" (Apollo 17, NASA, AS17-148-22727);
%                Lua - mapa de cores LRO/LROC (NASA/GSFC/Arizona State Univ.).
% Gera: ../images/kalman_apollo_terra_lua.gif

clc; clear; close all;
rng(3);
imgDir  = fullfile(fileparts(mfilename('fullpath')), '..', 'images');
makeGif = true;

%% Parâmetros (unidades didáticas: raio da Terra = 1)
muE = 1.0;  muM = 0.08;  D = 10;  Mpos = [D; 0];  Rsoi = 2.5;    % Lua fixa em (D,0)
RE  = 1.0;  RMd = 0.45;                                          % raios desenhados
r0  = 1.6;  tTLI = 10.0;  dvTLI = 0.235;  rLLO = 0.87;           % órbitas e TLI
dt  = 0.05;  T = 75;  N = round(T/dt);  t = (0:N)*dt;
tSleep = [18 24];  tMCC = 28;                                                % tripulação dormindo: sem rastreio

f   = @(x) dyn(x, muE, muM, Mpos, Rsoi);
Bd  = [dt^2/2*eye(2); dt*eye(2)];                                % u_k = aceleração do motor
Sw  = 2e-6;                                                      % PSD do ruído de aceleração
Q   = Sw*[dt^3/3*eye(2) dt^2/2*eye(2); dt^2/2*eye(2) dt*eye(2)];
H   = [eye(2) zeros(2)];
sigR = @(r) 0.03 + 0.012*norm(r);                                % desvio-padrão da medição
dvErr = 0.02;                                                    % erro de execução das queimas (2%)

%% Condição inicial e trajetória nominal (planejada em solo, sem ruídos)
th0 = -pi/2;
x0  = [r0*cos(th0); r0*sin(th0); -sin(th0)*sqrt(muE/r0); cos(th0)*sqrt(muE/r0)];
xn = x0;  Xn = zeros(4, N+1);  Xn(:,1) = x0;
for k = 1:N
    if abs(t(k) - tTLI) < dt/2, xn(3:4) = xn(3:4)*(1 + dvTLI/norm(xn(3:4))); end
    xn = rk4(f, xn, dt);  Xn(:,k+1) = xn;
end
dM = vecnorm(Xn(1:2,:) - Mpos);  dM(t < tTLI + 5) = inf;
[~, kPer] = min(dM);  rTgt = Xn(1:2,kPer);              % alvo: periápside atrás da Lua
P0  = diag([0.04 0.04 0.01 0.01].^2);

X  = zeros(4, N+1);  XH = X;  XD = X;  PP = zeros(4, 4, N+1);
Y  = nan(2, N+1);  Kg = nan(1, N+1);  burn = zeros(1, N+1);  hasY = false(1, N+1);  uDir = zeros(2, N+1);
X(:,1)  = x0;
XH(:,1) = x0 + chol(P0, 'lower')*randn(4,1);
XD(:,1) = XH(:,1);  P = P0;  PP(:,:,1) = P;  loiDone = false;  tLOI = nan;

%% Simulação
for k = 1:N
    xh = XH(:,k);
    % ---- guiagem: queimas calculadas com o estado ESTIMADO
    u = [0; 0];  Qb = zeros(4);
    if abs(t(k) - tTLI) < dt/2                                   % TLI
        u = dvTLI*xh(3:4)/norm(xh(3:4))/dt;  burn(k) = 1;
    end
    if abs(t(k) - tMCC) < dt/2                                   % MCC (mira no periápside planejado)
        u = mccSolve(f, xh, kPer - k, rTgt, dt)/dt;  burn(k) = 3;
    end
    d = xh(1:2) - Mpos;  dn = norm(d);  vr = d'*xh(3:4)/dn;
    if ~loiDone && t(k) > tTLI + 5 && dn < Rsoi && vr > 0       % LOI (após o periápside)
        vt = xh(3:4) - vr*d/dn;
        u  = (vt/norm(vt)*sqrt(muM/dn) - xh(3:4))/dt;  burn(k) = 2;
        loiDone = true;  tLOI = t(k);
    end
    if any(u), Qb = Bd*(dvErr*norm(u))^2*eye(2)*Bd';  uDir(:,k) = u/norm(u); end   % incerteza da queima
    % ---- planta real (queima com erro de execução + perturbações)
    uR = u + dvErr*norm(u)*randn(2,1);
    X(:,k+1) = rk4(f, X(:,k), dt) + Bd*uR + chol(Q + 1e-14*eye(4), 'lower')*randn(4,1);
    % ---- só predição (dead reckoning)
    XD(:,k+1) = rk4(f, XD(:,k), dt) + Bd*u;
    % ---- filtro de Kalman: predição
    xp = rk4(f, xh, dt) + Bd*u;
    A  = jacF(xh, muE, muM, Mpos, Rsoi, dt);
    P  = A*P*A' + Q + Qb;
    % ---- filtro de Kalman: correção (rastreio da Terra a cada 0.5 u.t.)
    rt = X(1:2,k+1);
    vis = ~behindMoon(rt, Mpos, RMd) && ~(t(k+1) > tSleep(1) && t(k+1) < tSleep(2));
    if mod(k, 10) == 0 && vis
        sR = sigR(rt);  Rk = sR^2*eye(2);
        y  = rt + sR*randn(2,1);
        K  = P*H'/(H*P*H' + Rk);
        xp = xp + K*(y - H*xp);
        P  = (eye(4) - K*H)*P*(eye(4) - K*H)' + K*Rk*K';
        Y(:,k+1) = y;  Kg(k+1) = trace(K(1:2,1:2))/2;  hasY(k+1) = true;
    end
    XH(:,k+1) = xp;  PP(:,:,k+1) = P;
end

eKF = vecnorm(XH(1:2,:) - X(1:2,:));
eDR = vecnorm(XD(1:2,:) - X(1:2,:));
s3  = 3*sqrt(squeeze(PP(1,1,:) + PP(2,2,:)))';
fprintf('TLI em t = %.2f,  LOI em t = %.2f,  %d medições\n', tTLI, tLOI, nnz(hasY));
fprintf('Erro final:  KF = %.3f   |  só predição = %.3f\n', eKF(end), eDR(end));
fprintf('Raio da órbita lunar real (min/max): %.3f / %.3f\n', ...
        min(vecnorm(X(1:2,t > tLOI + 1) - Mpos)), max(vecnorm(X(1:2,t > tLOI + 1) - Mpos)));
fprintf('Fração do tempo com erro KF < 3 sigma: %.1f %%\n', 100*mean(eKF <= s3));

%% Figura
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');
bg = [0.02 0.02 0.06];  fg = [0.93 0.93 0.93];
cT = [0.85 0.85 0.85];  cK = [1.00 0.62 0.10];  cD = [0.95 0.25 0.25];  cY = [0.30 0.95 0.45];
fig = figure('Color', bg, 'Position', [30 40 1400 800], 'InvertHardcopy', 'off');

% ---- cena
ax1 = axes(fig, 'Position', [0.01 0.30 0.98 0.62]);  hold(ax1, 'on');
axis(ax1, 'equal');  axis(ax1, 'off');  set(ax1, 'Color', bg);
xlim(ax1, [-2.4 12.0]);  ylim(ax1, [-3.6 2.3]);
ns = 350;  plot(ax1, -2.4 + 14.4*rand(1,ns), -3.6 + 5.9*rand(1,ns), '.', 'Color', [0.75 0.75 0.85], 'MarkerSize', 3);
[Ie, Ae] = earthDisk(fullfile(imgDir, 'terra_apollo17_blue_marble.jpg'));
image(ax1, 'XData', [-RE RE], 'YData', [RE -RE], 'CData', Ie, 'AlphaData', Ae);
[Im, Am] = moonDisk(fullfile(imgDir, 'moon_lroc_color_2k.jpg'), 500);
image(ax1, 'XData', D + [-RMd RMd], 'YData', [RMd -RMd], 'CData', Im, 'AlphaData', Am);
set(ax1, 'YDir', 'normal');
text(ax1, 0, -1.25, 'Terra', 'Color', fg, 'FontSize', 15, 'HorizontalAlignment', 'center');
text(ax1, D, -1.2, 'Lua', 'Color', fg, 'FontSize', 15, 'HorizontalAlignment', 'center');
text(ax1, 11.95, 2.2, '(fora de escala)', 'Color', [0.6 0.6 0.6], 'FontSize', 12, 'HorizontalAlignment', 'right');
% sombra de rádio atrás da Lua (região sem sinal da Terra)
patch(ax1, [D D 12 12], [RMd -RMd -RMd*12/D RMd*12/D], [1 1 1], 'FaceAlpha', 0.06, 'EdgeColor', 'none');

plot(ax1, X(1,:), X(2,:), '-', 'Color', [cT 0.18], 'LineWidth', 1);
hTrT = plot(ax1, nan, nan, '-',  'Color', cT, 'LineWidth', 1.5);
hTrK = plot(ax1, nan, nan, '-',  'Color', cK, 'LineWidth', 2);
hYm  = plot(ax1, nan, nan, 'x',  'Color', cY, 'MarkerSize', 7, 'LineWidth', 1.5);
hYn  = plot(ax1, nan, nan, 'o',  'Color', cY, 'MarkerSize', 12, 'LineWidth', 2.5);
hEl  = patch(ax1, nan, nan, cK, 'FaceAlpha', 0.20, 'EdgeColor', cK, 'LineWidth', 1.5);
hKF  = plot(ax1, nan, nan, '+',  'Color', cK, 'MarkerSize', 12, 'LineWidth', 2.5);
hFl  = patch(ax1, nan, nan, [1 0.55 0.1], 'EdgeColor', [1 0.9 0.3], 'LineWidth', 1);
[csmP, csmC] = csmShape();  sc = 0.04;  hCSM = gobjects(1, numel(csmP));
for i = 1:numel(csmP), hCSM(i) = patch(ax1, nan, nan, csmC{i}, 'EdgeColor', [0.2 0.2 0.2]); end
legend(ax1, [hTrT hTrK hEl hYm], {'trajet\''oria real', 'estimativa do filtro de Kalman', ...
       'incerteza $3\sigma$ (de $P_{k|k}$)', 'medi\c{c}\~oes $y_k$'}, ...
       'FontSize', 12, 'TextColor', fg, 'Color', bg, 'EdgeColor', [0.4 0.4 0.4], 'Position', [0.795 0.35 0.19 0.125]);
title(ax1, 'Apollo: navega\c{c}\~ao Terra $\rightarrow$ Lua com filtro de Kalman', 'Color', fg, 'FontSize', 20);
hPh = text(ax1, -2.3, 2.15, '', 'Color', fg, 'FontSize', 16, 'VerticalAlignment', 'top');
hSt = text(ax1, -2.3, 1.80, '', 'Color', cY, 'FontSize', 14, 'VerticalAlignment', 'top');

% ---- erro de posição
ax2 = axes(fig, 'Position', [0.06 0.07 0.42 0.20]);  hold(ax2, 'on');  styleAx(ax2, bg, fg);
set(ax2, 'YScale', 'log');
shadeGaps(ax2, t, hasY, tSleep, [1e-3 1e2]);
hS3 = plot(ax2, nan, nan, ':',  'Color', cK, 'LineWidth', 1.8);
hEK = plot(ax2, nan, nan, '-',  'Color', cK, 'LineWidth', 2.5);
hED = plot(ax2, nan, nan, '--', 'Color', cD, 'LineWidth', 2.5);
xlim(ax2, [0 T]);  ylim(ax2, [5e-3 20]);
xlabel(ax2, '$t$');  ylabel(ax2, '$\|r_k-\hat r_k\|$');
title(ax2, 'Erro de posi\c{c}\~ao', 'Color', fg, 'FontSize', 15);
legend(ax2, [hEK hS3 hED], {'Kalman', 'Kalman: $3\sigma$', 's\''o predi\c{c}\~ao'}, 'Location', 'northwest', ...
       'FontSize', 11, 'TextColor', fg, 'Color', bg, 'EdgeColor', [0.4 0.4 0.4], 'NumColumns', 3);

% ---- ganho de Kalman
ax3 = axes(fig, 'Position', [0.56 0.07 0.42 0.20]);  hold(ax3, 'on');  styleAx(ax3, bg, fg);
shadeGaps(ax3, t, hasY, tSleep, [0 1]);
hKg = plot(ax3, nan, nan, 'o-', 'Color', cY, 'MarkerSize', 3.5, 'MarkerFaceColor', cY, 'LineWidth', 1);
xlim(ax3, [0 T]);  ylim(ax3, [0 1]);
xlabel(ax3, '$t$');  ylabel(ax3, '$K_f$ (posi\c{c}\~ao)');
title(ax3, 'Ganho de Kalman: $\uparrow$ ap\''os ficar sem medi\c{c}\~oes ($P$ grande)', 'Color', fg, 'FontSize', 15);

%% Animação -> GIF
outFile = fullfile(imgDir, 'kalman_apollo_terra_lua.gif');
dpi = 75;  kk = 1:6:N+1;  if kk(end) ~= N+1, kk = [kk N+1]; end
szFr = [];  burnK = find(burn);
for n = 1:numel(kk)
    k = kk(n);  r = X(1:2,k);  v = X(3:4,k);
    % fase
    if t(k) < tTLI,                   ph = '1) \''Orbita de estacionamento na Terra';
    elseif t(k) < tTLI + 1.2,         ph = '2) TLI: queima que envia a nave \`a Lua ($u_k \neq 0$)';
    elseif t(k) > tMCC - 0.3 && t(k) < tMCC + 1.2, ph = '4) MCC: corre\c{c}\~ao de meio curso ($u_k$ calculado com $\hat x_k$)';
    elseif isnan(tLOI) || t(k) < tLOI, ph = '3) Viagem translunar';
    elseif t(k) < tLOI + 1.2,         ph = '5) LOI: queima atr\''as da Lua ($u_k$ calculado com $\hat x_k$)';
    else,                             ph = '6) \''Orbita lunar';
    end
    set(hPh, 'String', sprintf('$t=%5.1f$ \\quad %s', t(k), ph));
    jl = find(hasY(1:k), 1, 'last');
    if behindMoon(r, Mpos, RMd)
        st = 'Atr\''as da Lua: sem sinal da Terra $\Rightarrow$ s\''o predi\c{c}\~ao, $P$ cresce';  sc0 = [1 0.6 0.5];
    elseif t(k) > tSleep(1) && t(k) < tSleep(2)
        st = 'Tripula\c{c}\~ao dormindo: sem medi\c{c}\~oes $\Rightarrow$ s\''o predi\c{c}\~ao, $P$ cresce';  sc0 = [1 0.6 0.5];
    elseif ~isempty(jl) && k - jl < 5
        st = 'Medi\c{c}\~ao recebida $\Rightarrow$ corre\c{c}\~ao: $\hat x_{k|k}=\hat x_{k|k-1}+K_f(y_k-H\hat x_{k|k-1})$';  sc0 = cY;
    else
        st = 'Predi\c{c}\~ao: $\hat x_{k|k-1}=f(\hat x_{k-1|k-1})+Bu_{k-1}$';  sc0 = fg;
    end
    set(hSt, 'String', st, 'Color', sc0);
    % trilhas
    set(hTrT, 'XData', X(1,1:k),  'YData', X(2,1:k));
    set(hTrK, 'XData', XH(1,1:k), 'YData', XH(2,1:k));
    set(hKF,  'XData', XH(1,k), 'YData', XH(2,k));
    iy = find(hasY(max(1,k-150):k)) + max(1,k-150) - 1;
    set(hYm, 'XData', Y(1,iy), 'YData', Y(2,iy));
    if ~isempty(jl) && k - jl < 5, set(hYn, 'XData', Y(1,jl), 'YData', Y(2,jl));
    else, set(hYn, 'XData', nan, 'YData', nan); end
    % elipse 3 sigma
    [Ve, De] = eig(PP(1:2,1:2,k));  th = linspace(0, 2*pi, 60);
    E = XH(1:2,k) + 3*Ve*sqrt(max(De,0))*[cos(th); sin(th)];
    set(hEl, 'XData', E(1,:), 'YData', E(2,:));
    % nave (orientada pela velocidade; vira na LOI para frear)
    e1 = v/norm(v);  kb = burnK(k - burnK >= 0 & k - burnK < 20);  firing = ~isempty(kb);
    if firing, e1 = uDir(:,kb(end)); end                         % aponta o motor na direção da queima
    Rb = [e1 [-e1(2); e1(1)]];
    for i = 1:numel(csmP)
        W = r + Rb*(sc*csmP{i});  set(hCSM(i), 'XData', W(1,:), 'YData', W(2,:));
    end
    if firing
        fl = r + Rb*(sc*[-10.3 -15-3*rand -10.3; 1.1 0 -1.1]);  set(hFl, 'XData', fl(1,:), 'YData', fl(2,:));
    else
        set(hFl, 'XData', nan, 'YData', nan);
    end
    % gráficos
    set(hEK, 'XData', t(1:k), 'YData', max(eKF(1:k), 1e-3));
    set(hED, 'XData', t(1:k), 'YData', max(eDR(1:k), 1e-3));
    set(hS3, 'XData', t(1:k), 'YData', s3(1:k));
    ik = find(hasY(1:k));  set(hKg, 'XData', t(ik), 'YData', Kg(ik));
    drawnow;

    if makeGif
        img = print(fig, '-RGBImage', sprintf('-r%d', dpi), '-opengl');
        if isempty(szFr), szFr = size(img, [1 2]); end
        tmp = zeros([szFr 3], 'uint8');
        h = min(szFr(1), size(img,1));  wd = min(szFr(2), size(img,2));
        tmp(1:h, 1:wd, :) = img(1:h, 1:wd, :);
        [Ai, map] = rgb2ind(tmp, 160, 'nodither');
        dly = 0.06 + 1.0*(n==1) + 4.0*(n==numel(kk));
        if n == 1
            imwrite(Ai, map, outFile, 'gif', 'LoopCount', Inf, 'DelayTime', dly);
        else
            imwrite(Ai, map, outFile, 'gif', 'WriteMode', 'append', 'DelayTime', dly);
        end
    end
end
if makeGif, fprintf('%d quadros. GIF salvo em: %s\n', numel(kk), outFile); end

%% Funções auxiliares
function dx = dyn(x, muE, muM, Mpos, Rsoi)          % cônicas remendadas
    r = x(1:2);  d = r - Mpos;
    if norm(d) < Rsoi, a = -muM*d/norm(d)^3; else, a = -muE*r/norm(r)^3; end
    dx = [x(3:4); a];
end

function x = rk4(f, x, dt)
    k1 = f(x);  k2 = f(x + dt/2*k1);  k3 = f(x + dt/2*k2);  k4 = f(x + dt*k3);
    x = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

function A = jacF(x, muE, muM, Mpos, Rsoi, dt)      % A_k = df/dx (linearização)
    r = x(1:2);  d = r - Mpos;
    if norm(d) < Rsoi, mu = muM; s = d; else, mu = muE; s = r; end
    ns = norm(s);  G = mu*(3*(s*s')/ns^5 - eye(2)/ns^3);
    F = [zeros(2) eye(2); G zeros(2)];
    A = eye(4) + F*dt + (F*dt)^2/2;
end

function dv = mccSolve(f, xh, nSteps, rTgt, dt)   % tiro (shooting): Newton em dv
    dv = [0; 0];  g = @(dv) propR(f, xh + [0; 0; dv], nSteps, dt) - rTgt;
    for it = 1:8
        g0 = g(dv);  J = zeros(2);  e = 1e-6;
        for i = 1:2, di = zeros(2,1); di(i) = e; J(:,i) = (g(dv + di) - g0)/e; end
        dv = dv - J\g0;
    end
end

function r = propR(f, x, n, dt)
    for i = 1:n, x = rk4(f, x, dt); end
    r = x(1:2);
end

function b = behindMoon(r, Mpos, Rm)                % Lua bloqueia a linha Terra -> nave?
    if r(1) <= Mpos(1), b = false; return; end
    s = Mpos'*r/(r'*r);  b = norm(s*r - Mpos) < Rm;
end

function shadeGaps(ax, t, hasY, tSleep, yl)
    ty = t(hasY);  gaps = find(diff(ty) > 1.0);
    for g = gaps
        patch(ax, [ty(g) ty(g+1) ty(g+1) ty(g)], [yl(1) yl(1) yl(2) yl(2)], [1 1 1], ...
              'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    if yl(1) > 0, yt = 0.1; else, yt = mean(yl); end
    text(ax, mean(tSleep), yt, 'sem medi\c{c}\~oes', 'Color', [0.7 0.7 0.7], 'FontSize', 10, 'HorizontalAlignment', 'center', 'Rotation', 90);
end

function styleAx(ax, bg, fg)
    set(ax, 'Color', bg*1.8, 'XColor', fg, 'YColor', fg, 'GridColor', [0.7 0.7 0.7], ...
        'GridAlpha', 0.25, 'FontSize', 12, 'Box', 'on');
    grid(ax, 'on');
end

function [I, A] = earthDisk(file)                    % recorta o disco da foto da Terra
    I0 = imread(file);  g = mean(double(I0), 3);  m = g > 30;
    [rr, cc] = find(m);  c = [mean(cc) mean(rr)];  R = sqrt(nnz(m)/pi);
    L = round(R);  [Xg, Yg] = meshgrid(-L:L);
    xi = min(max(round(c(1) + Xg), 1), size(I0,2));  yi = min(max(round(c(2) + Yg), 1), size(I0,1));
    I = zeros(2*L+1, 2*L+1, 3, 'uint8');
    for ch = 1:3, Ic = I0(:,:,ch); I(:,:,ch) = Ic(sub2ind(size(Ic), yi, xi)); end
    A = double(hypot(Xg, Yg) <= L);
end

function [I, A] = moonDisk(file, n)                  % projeção ortográfica do mapa LROC (lado visível)
    T = imread(file);  [h, w, ~] = size(T);
    [u, v] = meshgrid(linspace(-1, 1, n));  v = -v;  m = u.^2 + v.^2 <= 1;
    z = sqrt(max(0, 1 - u.^2 - v.^2));  lat = asin(v);  lon = atan2(u, z);
    col = min(max(round((lon + pi)/(2*pi)*w), 1), w);  row = min(max(round((pi/2 - lat)/pi*h), 1), h);
    shade = 0.35 + 0.65*z.^0.5;
    I = zeros(n, n, 3, 'uint8');
    for ch = 1:3, Tc = double(T(:,:,ch)); I(:,:,ch) = uint8(Tc(sub2ind([h w], row, col)).*shade); end
    A = double(m);
end

function [P, C] = csmShape()                         % silhueta do CSM (metros, nariz em +x)
    P = { [-7.4 0 0 -7.4; -1.96 -1.96 1.96 1.96], ...             % módulo de serviço
          [0 3.4 3.4 0; -1.96 -0.35 0.35 1.96], ...                % módulo de comando
          [-10.3 -7.4 -7.4 -10.3; -1.25 -0.45 0.45 1.25], ...      % tubeira do SPS
          [-2.6 -1.8 -1.8 -2.6; -1.96 -1.96 1.96 1.96] };          % faixa de radiadores
    C = {[0.72 0.74 0.78], [0.86 0.87 0.90], [0.38 0.38 0.42], [0.97 0.97 0.97]};
end
