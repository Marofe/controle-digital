%% ANIM_KALMAN_GPS  Filtro de Kalman aplicado ao GPS - antena em São Carlos (GIF)
% Aula 10 (SEL0359 - Controle Digital, 2026) - Prof. Marcos R. Fernandes
%
% Constelação GPS simplificada: 24 satélites em 6 planos orbitais (inclinação 55 graus,
% raio 26 560 km, período ~12 h). Uma antena fixa em São Carlos-SP mede a
% pseudodistância a cada satélite visível (elevação > 10 graus):
%       rho_i = || r_sat,i - r || + b + m_i + eps_i ,   eps_i ~ N(0, sigma^2), sigma = 5 m
% m_i: erro lento (atmosfera/multicaminho, Gauss-Markov 1.5 m, tau = 10 min) não modelado no filtro
% em que b = c*dt é o erro do relógio do receptor (em metros).
% Estado do EKF:  x = [r (ECEF, 3); b; bdot]   (antena parada, relógio com deriva)
%       x_{k+1} = A x_k + w_k ,   A = blkdiag(I_3, [1 dt; 0 1])
%       y_k     = h(x_k) + eps_k  (não linear -> H_k = dh/dx: vetores unitários de visada)
% Comparação: solução instantânea por mínimos quadrados (uma época por vez) x EKF.
% Imagem da Terra: NASA Blue Marble (Visible Earth, NASA/GSFC).
% Gera: ../images/kalman_gps_sao_carlos.gif

clc; clear; close all;
rng(5);
imgDir  = fullfile(fileparts(mfilename('fullpath')), '..', 'images');
makeGif = true;

%% Constantes (m, s)
RE = 6371e3;  aG = 26560e3;  mu = 3.986004418e14;  wE = 7.2921159e-5;
nG = sqrt(mu/aG^3);  inc = 55*pi/180;
dt = 10;  T = 6000;  t = 0:dt:T;  Nt = numel(t);
sig = 5;  elMin = 10;

% constelação: 6 planos x 4 satélites
[pl, sl] = meshgrid(0:5, 0:3);  pl = pl(:)';  sl = sl(:)';  Ns = numel(pl);
RAAN = pl*60*pi/180;  u0 = (sl*90 + pl*15)*pi/180;
satECI = @(tt) aG*[cos(RAAN).*cos(u0 + nG*tt) - sin(RAAN).*sin(u0 + nG*tt)*cos(inc);
                   sin(RAAN).*cos(u0 + nG*tt) + cos(RAAN).*sin(u0 + nG*tt)*cos(inc);
                   sin(u0 + nG*tt)*sin(inc)];
Rz = @(a) [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
satECEF = @(tt) Rz(-wE*tt)*satECI(tt);

% antena em São Carlos (Terra esférica para simplificar)
lat = -22.0087;  lon = -47.8909;  alt = 850;
up = [cosd(lat)*cosd(lon); cosd(lat)*sind(lon); sind(lat)];
ea = [-sind(lon); cosd(lon); 0];  no = cross(up, ea);  Renu = [ea no up]';
rTrue = (RE + alt)*up;

%% Simulação: medições, mínimos quadrados e EKF
bT = 120;  dT = 0.3;                                       % relógio real: 120 m, deriva 0.3 m/s
Sf = 0.05;  Sg = 0.005;                                    % PSD do relógio
Qc = [Sf*dt + Sg*dt^3/3, Sg*dt^2/2; Sg*dt^2/2, Sg*dt];
Qd = blkdiag((1e-3)^2*eye(3), Qc);  A = blkdiag(eye(3), [1 dt; 0 1]);
vis = false(Ns, Nt);  rho = nan(Ns, Nt);
sM = 1.5;  tauM = 600;  aM = exp(-dt/tauM);  mE = sM*randn(Ns, 1);   % erro lento por satélite
XL = nan(4, Nt);  XK = nan(5, Nt);  PK = nan(5, 5, Nt);
xl = [0; 0; 0; 0];  xk = [];  P = [];
for k = 1:Nt
    S = satECEF(t(k));
    L = S - rTrue;  Ln = vecnorm(L);  el = asind((up'*L)./Ln);
    v = el > elMin;  vis(:,k) = v;
    rho(v,k) = Ln(v)' + bT + mE(v) + sig*randn(nnz(v), 1);
    mE = aM*mE + sM*sqrt(1 - aM^2)*randn(Ns, 1);
    Sv = S(:,v);  yv = rho(v,k);
    % ---- mínimos quadrados (Gauss-Newton), só esta época
    for it = 1:8
        d = Sv - xl(1:3);  dn = vecnorm(d)';
        Hl = [-(d./dn')' ones(numel(dn),1)];
        xl = xl + Hl\(yv - dn - xl(4));
    end
    XL(:,k) = xl;
    % ---- EKF
    if k == 1
        xk = [xl; 0];  P = diag([30 30 30 60 5].^2);
    else
        xk = A*xk;  P = A*P*A' + Qd;                       % predição
        d  = Sv - xk(1:3);  dn = vecnorm(d)';
        H  = [-(d./dn')' ones(numel(dn),1) zeros(numel(dn),1)];
        Rk = sig^2*eye(numel(dn));
        K  = P*H'/(H*P*H' + Rk);                           % correção
        xk = xk + K*(yv - dn - xk(4));
        P  = (eye(5) - K*H)*P*(eye(5) - K*H)' + K*Rk*K';
    end
    XK(:,k) = xk;  PK(:,:,k) = P;
    bT = bT + dT*dt + sqrt(Sf*dt)*randn;  dT = dT + sqrt(Sg*dt)*randn;
end
eL = Renu*(XL(1:3,:) - rTrue);  eK = Renu*(XK(1:3,:) - rTrue);
hL = vecnorm(eL(1:2,:));  hK = vecnorm(eK(1:2,:));
s3 = zeros(1, Nt);
for k = 1:Nt, Pe = Renu*PK(1:3,1:3,k)*Renu'; s3(k) = 3*sqrt(Pe(1,1) + Pe(2,2)); end
fprintf('Satélites visíveis: %d a %d\n', min(sum(vis)), max(sum(vis)));
fprintf('Erro horizontal RMS:  mínimos quadrados = %.2f m  |  EKF (após 10 min) = %.2f m\n', ...
        rms(hL), rms(hK(t > 600)));

%% Figura
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');
bg = [0.01 0.01 0.04];  fg = [0.93 0.93 0.93];
cK = [1.00 0.62 0.10];  cL = [0.65 0.65 0.70];  cV = [0.35 1.00 0.50];
fig = figure('Color', bg, 'Position', [30 40 1400 800], 'InvertHardcopy', 'off');

% ---- Terra 3D (NASA Blue Marble, nuvens parcialmente transparentes)
ax1 = axes(fig, 'Position', [0.0 0.0 0.60 0.93]);  hold(ax1, 'on');
I1 = double(imread(fullfile(imgDir, 'terra_blue_marble_nuvens_2048.jpg')));
I0 = double(imread(fullfile(imgDir, 'terra_blue_marble_2048.jpg')));
Ie = uint8(0.65*I1 + 0.35*I0);
[SX, SY, SZ] = sphere(200);
surf(ax1, RE*SX, RE*SY, RE*SZ, 'CData', flipud(Ie), 'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
     'AmbientStrength', 0.35, 'DiffuseStrength', 0.9, 'SpecularStrength', 0.15, 'SpecularExponent', 8);
surf(ax1, 1.03*RE*SX, 1.03*RE*SY, 1.03*RE*SZ, 'FaceColor', [0.45 0.70 1.0], 'FaceAlpha', 0.12, ...
     'EdgeColor', 'none', 'FaceLighting', 'none');                                        % atmosfera
sunDir = Rz(35*pi/180)*up + [0; 0; 0.25];
light(ax1, 'Position', sunDir', 'Style', 'infinite');  lighting(ax1, 'gouraud');
axis(ax1, 'equal', 'off', 'vis3d');  set(ax1, 'Color', bg, 'Clipping', 'off');
st = randn(3, 500);  st = 1e9*st./vecnorm(st);
plot3(ax1, st(1,:), st(2,:), st(3,:), '.', 'Color', [0.8 0.8 0.9], 'MarkerSize', 3);

% antena em São Carlos
pA = 1.012*rTrue;
plot3(ax1, pA(1), pA(2), pA(3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [1 0.2 0.2], 'MarkerEdgeColor', 'w');
pT = 1.10*rTrue - 0.10*RE*ea;
text(ax1, pT(1), pT(2), pT(3), 'S\~ao Carlos', 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'right');

% órbitas, satélites e linhas de visada
th = linspace(0, 2*pi, 120);
hOrb = gobjects(1, 6);
for p = 1:6, hOrb(p) = plot3(ax1, nan, nan, nan, '-', 'Color', [0.55 0.65 0.85 0.35], 'LineWidth', 1); end
[bV, bF] = boxMesh([1 1 1]);  [pV, pF] = boxMesh([0.15 7 1.6]);       % corpo e painéis (unid. de escala)
scS = 0.55e6;
hBody = gobjects(1, Ns);  hPan = hBody;  hLos = hBody;  hTxt = hBody;
for i = 1:Ns
    hLos(i)  = plot3(ax1, nan, nan, nan, '-', 'Color', [cV 0.75], 'LineWidth', 1.4);
    hBody(i) = patch(ax1, 'Vertices', bV, 'Faces', bF, 'FaceColor', [0.85 0.70 0.30], 'EdgeColor', 'none');
    hPan(i)  = patch(ax1, 'Vertices', pV, 'Faces', pF, 'FaceColor', [0.15 0.25 0.60], 'EdgeColor', [0.6 0.6 0.7]);
    hTxt(i)  = text(ax1, nan, nan, nan, '', 'Color', cV, 'FontSize', 10);
end
camproj(ax1, 'perspective');  camva(ax1, 26);
camDir = Rz(-25*pi/180)*up + [0; 0; 0.15];  camDir = camDir/norm(camDir);
campos(ax1, 1.1e8*camDir');  camtarget(ax1, (0.25*RE*up)');  camup(ax1, [0 0 1]);
title(ax1, 'GPS: pseudodist\^ancias at\''e uma antena em S\~ao Carlos', 'Color', fg, 'FontSize', 19);
hInfo = annotation(fig, 'textbox', [0.01 0.01 0.58 0.08], 'String', '', 'Color', fg, 'EdgeColor', 'none', ...
                   'FontSize', 15, 'Interpreter', 'latex');

% ---- posição estimada (plano Leste x Norte)
ax2 = axes(fig, 'Position', [0.64 0.50 0.33 0.40]);  hold(ax2, 'on');  styleAx(ax2, bg, fg);
hLS = plot(ax2, nan, nan, '.', 'Color', cL, 'MarkerSize', 8);
hEl = patch(ax2, nan, nan, cK, 'FaceAlpha', 0.15, 'EdgeColor', cK, 'LineWidth', 1.2);
hKt = plot(ax2, nan, nan, '-', 'Color', cK, 'LineWidth', 2);
hKp = plot(ax2, nan, nan, 'o', 'MarkerSize', 9, 'MarkerFaceColor', cK, 'MarkerEdgeColor', 'k');
plot(ax2, 0, 0, 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
axis(ax2, 'equal');  xlim(ax2, [-16 16]);  ylim(ax2, [-16 16]);
xlabel(ax2, 'Leste [m]');  ylabel(ax2, 'Norte [m]');
title(ax2, 'Posi\c{c}\~ao da antena (erro)', 'Color', fg, 'FontSize', 16);
legend(ax2, [hLS hKt hEl], {'m\''inimos quadrados (cada \''epoca)', 'filtro de Kalman', '$3\sigma$ de $P_{k|k}$'}, ...
       'Location', 'southoutside', 'NumColumns', 3, 'FontSize', 11, 'TextColor', fg, 'Color', bg, 'EdgeColor', [0.4 0.4 0.4]);

% ---- erro horizontal no tempo
ax3 = axes(fig, 'Position', [0.64 0.08 0.33 0.26]);  hold(ax3, 'on');  styleAx(ax3, bg, fg);
hHL = plot(ax3, nan, nan, '-', 'Color', cL, 'LineWidth', 1);
hS3 = plot(ax3, nan, nan, ':', 'Color', cK, 'LineWidth', 1.8);
hHK = plot(ax3, nan, nan, '-', 'Color', cK, 'LineWidth', 2.5);
xlim(ax3, [0 T/60]);  ylim(ax3, [0 15]);
xlabel(ax3, '$t$ [min]');  ylabel(ax3, 'erro horizontal [m]');
title(ax3, '$\rho_i=\|r_{\mathrm{sat},i}-r\|+c\,\delta t+\epsilon_i$', 'Color', fg, 'FontSize', 15);

%% Animação -> GIF
outFile = fullfile(imgDir, 'kalman_gps_sao_carlos.gif');
dpi = 80;  kk = 1:3:Nt;  if kk(end) ~= Nt, kk = [kk Nt]; end
szFr = [];  tm = t/60;
for n = 1:numel(kk)
    k = kk(n);  S = satECEF(t(k));  Ri = Rz(-wE*t(k));
    for p = 1:6                                             % órbitas (anel inercial -> ECEF)
        Om = (p-1)*60*pi/180;
        ring = aG*[cos(Om)*cos(th) - sin(Om)*sin(th)*cos(inc); sin(Om)*cos(th) + cos(Om)*sin(th)*cos(inc); sin(th)*sin(inc)];
        ring = Ri*ring;  set(hOrb(p), 'XData', ring(1,:), 'YData', ring(2,:), 'ZData', ring(3,:));
    end
    nv = 0;
    for i = 1:Ns
        s = S(:,i);  ez = -s/norm(s);                        % corpo aponta para a Terra
        hN = Ri*[sin(RAAN(i))*sin(inc); -cos(RAAN(i))*sin(inc); cos(inc)];   % normal da órbita
        ey = hN - (hN'*ez)*ez;  ey = ey/norm(ey);  ex = cross(ey, ez);
        Rb = [ex ey ez];
        set(hBody(i), 'Vertices', (s + Rb*(scS*bV'))');
        set(hPan(i),  'Vertices', (s + Rb*(scS*pV'))');
        if vis(i,k)
            nv = nv + 1;
            set(hLos(i), 'XData', [s(1) pA(1)], 'YData', [s(2) pA(2)], 'ZData', [s(3) pA(3)]);
            set(hTxt(i), 'Position', (s + 0.9e6*ex)', 'String', sprintf('%.0f km', norm(s - rTrue)/1e3));
            set(hBody(i), 'FaceColor', [0.95 0.80 0.35]);
        else
            set(hLos(i), 'XData', nan, 'YData', nan, 'ZData', nan);  set(hTxt(i), 'String', '');
            set(hBody(i), 'FaceColor', [0.55 0.50 0.40]);
        end
    end
    set(hInfo, 'String', sprintf(['$t = %4.1f$ min \\quad sat\\''elites vis\\''iveis: %d \\quad ' ...
        '(linhas verdes: dist\\^ancias medidas, $\\sigma = %d$ m)'], tm(k), nv, sig));
    % painéis
    set(hLS, 'XData', eL(1,1:k), 'YData', eL(2,1:k));
    set(hKt, 'XData', eK(1,1:k), 'YData', eK(2,1:k));  set(hKp, 'XData', eK(1,k), 'YData', eK(2,k));
    Pe = Renu*PK(1:3,1:3,k)*Renu';  [Ve, De] = eig(Pe(1:2,1:2));  a2 = linspace(0, 2*pi, 60);
    E = eK(1:2,k) + 3*Ve*sqrt(max(De, 0))*[cos(a2); sin(a2)];
    set(hEl, 'XData', E(1,:), 'YData', E(2,:));
    set(hHL, 'XData', tm(1:k), 'YData', hL(1:k));
    set(hHK, 'XData', tm(1:k), 'YData', hK(1:k));
    set(hS3, 'XData', tm(1:k), 'YData', s3(1:k));
    drawnow;

    if makeGif
        img = print(fig, '-RGBImage', sprintf('-r%d', dpi), '-opengl');
        if isempty(szFr), szFr = size(img, [1 2]); end
        tmp = zeros([szFr 3], 'uint8');
        h = min(szFr(1), size(img,1));  wd = min(szFr(2), size(img,2));
        tmp(1:h, 1:wd, :) = img(1:h, 1:wd, :);
        [Ai, map] = rgb2ind(tmp, 255, 'nodither');
        dly = 0.07 + 1.0*(n==1) + 4.0*(n==numel(kk));
        if n == 1
            imwrite(Ai, map, outFile, 'gif', 'LoopCount', Inf, 'DelayTime', dly);
        else
            imwrite(Ai, map, outFile, 'gif', 'WriteMode', 'append', 'DelayTime', dly);
        end
    end
end
if makeGif, fprintf('%d quadros. GIF salvo em: %s\n', numel(kk), outFile); end

%% Funções auxiliares
function [V, F] = boxMesh(sz)                       % paralelepípedo centrado na origem
    [X, Y, Z] = ndgrid([-1 1]/2, [-1 1]/2, [-1 1]/2);
    V = [X(:)*sz(1) Y(:)*sz(2) Z(:)*sz(3)];
    F = [1 3 7 5; 2 4 8 6; 1 2 6 5; 3 4 8 7; 1 2 4 3; 5 6 8 7];
end

function styleAx(ax, bg, fg)
    set(ax, 'Color', bg*2.5, 'XColor', fg, 'YColor', fg, 'GridColor', [0.7 0.7 0.7], ...
        'GridAlpha', 0.25, 'FontSize', 12, 'Box', 'on');
    grid(ax, 'on');
end
