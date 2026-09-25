%% Quadricóptero rastreando uma trajetória em "8" com realimentação de estados + integrador
% Aula 8 - SEL0359 Controle Digital (2026) - Prof. Marcos R. Fernandes
%
% Modelo linearizado em torno do voo pairado (hover), eixos desacoplados:
%   eixo x: [p; v; theta; omega],  p'' = g*theta + w,   theta'' = tau_theta
%   eixo y: [p; v; phi;   omega],  p'' = g*phi   + w,   phi''   = tau_phi
%   eixo z: [p; v],                p'' = u_z     + w
% Discretização ZOH (Ts = 0.02 s) e controle
%   u_k = -K (x_k - x_ref,k) + u_ff,k + K_i m_k,   m_{k+1} = m_k + (r_k - y_k)
% com ganhos obtidos por dlqr do sistema aumentado [x; m].
% Comparação: mesma realimentação SEM integrador (u = -K0 x + N r).
% Uma rajada de vento constante atua a partir de t = 8 s. Ambos usam feedforward
% da referência de estados; só o integrador elimina o erro devido ao vento.
%
% Gera: ../images/aula8_quadricoptero_8.gif

clear; close all;
out = fullfile(fileparts(mfilename('fullpath')), '..', 'images');
set(groot, 'defaultTextInterpreter','latex', 'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex', 'defaultAxesFontSize', 12);
cB = [0.165 0.494 0.878]; cR = [0.878 0.141 0.141]; cK = [0.15 0.15 0.15]; cGr = [0.6 0.6 0.6];

g = 9.81; Ts = 0.02; T = 40; N = round(T/Ts);

%% Modelos discretos
Ah = [0 1 0 0; 0 0 g 0; 0 0 0 1; 0 0 0 0];  Bh = [0; 0; 0; 1];  Hh = [1 0 0 0];
Az = [0 1; 0 0];                             Bz = [0; 1];        Hz = [1 0];
[Ahd, Bhd] = c2d(Ah, Bh, Ts);  [Azd, Bzd] = c2d(Az, Bz, Ts);
Ewh = c2d(ss(Ah, [0;1;0;0], eye(4), 0), Ts); Ewh = Ewh.B;     % entrada de perturbação

%% Projeto: sistema aumentado + LQR discreto
R = 0.1;
[Kh, Kih] = projeta(Ahd, Bhd, Hh, diag([20 2 1 0.1 5]), R);
[Kz, Kiz] = projeta(Azd, Bzd, Hz, diag([20 2 2]), R);
% sem integrador (mesmas ponderações nos estados)
K0h = dlqr(Ahd, Bhd, diag([20 2 1 0.1]), R);
K0z = dlqr(Azd, Bzd, diag([20 2]), R);

%% Referência: lemniscata (figura em "8") + decolagem suave
% Para seguir uma trajetória variante no tempo usa-se também a referência de
% estados (feedforward): x_ref = [p; dp; d2p/g; d3p/g] e u_ff = d4p/g.
% O integrador elimina o erro causado pela perturbação (vento) desconhecida.
t  = (0:N)*Ts;  a = 2;  wr = 2*pi/20;  b = a/2;  w2 = 2*wr;
X = [a*sin(wr*t); a*wr*cos(wr*t); -a*wr^2*sin(wr*t); -a*wr^3*cos(wr*t); a*wr^4*sin(wr*t)];
Y = [b*sin(w2*t); b*w2*cos(w2*t); -b*w2^2*sin(w2*t); -b*w2^3*cos(w2*t); b*w2^4*sin(w2*t)];
tau = 1.2;  E = exp(-t/tau);
Z = [1.5*(1-E); 1.5/tau*E; -1.5/tau^2*E];
xref = {[X(1:2,:); X(3:4,:)/g], [Y(1:2,:); Y(3:4,:)/g]};
uff  = {X(5,:)/g, Y(5,:)/g};
rx = X(1,:); ry = Y(1,:); rz = Z(1,:);
w  = zeros(2, N+1);  w(1, t >= 8) = 1.5;  w(2, t >= 8) = -1.05;   % vento [m/s^2]

%% Simulação (com e sem integrador)
Xi = zeros(4,N+1,2); Zi = zeros(2,N+1); mi = zeros(3,1);   % com integrador (eixos x,y)
Xo = zeros(4,N+1,2); Zo = zeros(2,N+1);                    % sem integrador
for k = 1:N
    for ax = 1:2
        xr = xref{ax}(:,k);
        % com integrador:  u = -K (x - x_ref) + u_ff + Ki m
        u = -Kh*(Xi(:,k,ax) - xr) + uff{ax}(k) + Kih*mi(ax);
        mi(ax) = mi(ax) + (xr(1) - Hh*Xi(:,k,ax));
        Xi(:,k+1,ax) = Ahd*Xi(:,k,ax) + Bhd*u + Ewh*w(ax,k);
        % sem integrador:  u = -K0 (x - x_ref) + u_ff
        u0 = -K0h*(Xo(:,k,ax) - xr) + uff{ax}(k);
        Xo(:,k+1,ax) = Ahd*Xo(:,k,ax) + Bhd*u0 + Ewh*w(ax,k);
    end
    uz = -Kz*(Zi(:,k) - Z(1:2,k)) + Z(3,k) + Kiz*mi(3);  mi(3) = mi(3) + (rz(k) - Hz*Zi(:,k));
    Zi(:,k+1) = Azd*Zi(:,k) + Bzd*uz;
    uz0 = -K0z*(Zo(:,k) - Z(1:2,k)) + Z(3,k);
    Zo(:,k+1) = Azd*Zo(:,k) + Bzd*uz0;
end
Pi = [squeeze(Xi(1,:,1)); squeeze(Xi(1,:,2)); Zi(1,:)];
Po = [squeeze(Xo(1,:,1)); squeeze(Xo(1,:,2)); Zo(1,:)];
Ri = [rx; ry; rz];
ei = vecnorm(Pi - Ri);  eo = vecnorm(Po - Ri);
fprintf('Erro RMS (t > 15 s): com integrador %.3f m, sem integrador %.3f m\n', ...
    sqrt(mean(ei(t>15).^2)), sqrt(mean(eo(t>15).^2)));

%% Animação
f = figure('Color','w','Position',[40 40 1300 620]);
tl = tiledlayout(f, 2, 3, 'TileSpacing','compact', 'Padding','compact');
ax3 = nexttile(tl, 1, [2 2]); axT = nexttile(tl, 3); axE = nexttile(tl, 6);
gif = fullfile(out, 'aula8_quadricoptero_8.gif');
passo = 10; first = true; dpi = 110;
idxs = [ones(1,6) 1:passo:N+1 (N+1)*ones(1,15)];
for q = 1:numel(idxs)
    k = idxs(q);
    % ---------- 3D ----------
    cla(ax3); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
    plot3(ax3, rx, ry, rz(end)*ones(size(rx)), '--', 'Color', [0 0 0 0.35], 'LineWidth', 1.2);
    plot3(ax3, Po(1,1:k), Po(2,1:k), Po(3,1:k), '-', 'Color', [cR 0.55], 'LineWidth', 1.5);
    plot3(ax3, Pi(1,1:k), Pi(2,1:k), Pi(3,1:k), '-', 'Color', cB, 'LineWidth', 2.2);
    plot3(ax3, Ri(1,k), Ri(2,k), Ri(3,k), 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    quad(ax3, Po(:,k), Xo(3,k,1), Xo(3,k,2), cR, 0.5);
    quad(ax3, Pi(:,k), Xi(3,k,1), Xi(3,k,2), cB, 1);
    % sombras no chão
    plot3(ax3, Pi(1,1:k), Pi(2,1:k), 0*Pi(3,1:k), '-', 'Color', [cB 0.15], 'LineWidth', 1);
    if t(k) >= 8
        quiver3(ax3, -3.2, 2.2, 2.6, 0.9, -0.6, 0, 0, 'Color', [0.3 0.3 0.3], 'LineWidth', 2.5, 'MaxHeadSize', 0.8);
        text(ax3, -3.2, 2.2, 2.95, 'vento', 'FontSize', 13);
    end
    xlim(ax3, [-3.5 3.5]); ylim(ax3, [-2.5 2.5]); zlim(ax3, [0 3]);
    daspect(ax3, [1 1 1]); view(ax3, -35 + 0.08*k*Ts*10, 28);
    xlabel(ax3, '$x$ [m]'); ylabel(ax3, '$y$ [m]'); zlabel(ax3, '$z$ [m]');
    title(ax3, sprintf('Trajet\\''oria em ``8'''': $t=%.1f$ s', t(k)), 'FontSize', 15);
    % ---------- vista superior ----------
    cla(axT); hold(axT,'on'); grid(axT,'on'); box(axT,'on'); axis(axT,'equal');
    plot(axT, rx, ry, '--', 'Color', [0 0 0 0.4], 'LineWidth', 1.2);
    plot(axT, Po(1,1:k), Po(2,1:k), '-', 'Color', cR, 'LineWidth', 1.4);
    plot(axT, Pi(1,1:k), Pi(2,1:k), '-', 'Color', cB, 'LineWidth', 2);
    plot(axT, Pi(1,k), Pi(2,k), 'o', 'MarkerFaceColor', cB, 'MarkerEdgeColor', 'k');
    plot(axT, Po(1,k), Po(2,k), 'o', 'MarkerFaceColor', cR, 'MarkerEdgeColor', 'k');
    xlim(axT, [-2.8 2.8]); ylim(axT, [-1.6 1.6]);
    title(axT, 'Vista superior ($xy$)', 'FontSize', 13);
    % ---------- erro ----------
    cla(axE); hold(axE,'on'); grid(axE,'on'); box(axE,'on');
    plot(axE, t(1:k), eo(1:k), '-', 'Color', cR, 'LineWidth', 1.6);
    plot(axE, t(1:k), ei(1:k), '-', 'Color', cB, 'LineWidth', 2);
    xline(axE, 8, ':', 'Color', cK);
    xlim(axE, [0 T]); ylim(axE, [0 0.8]);
    xlabel(axE, '$t$ [s]'); ylabel(axE, '$\|r_k-y_k\|$ [m]');
    legend(axE, {'sem integrador','com integrador'}, 'Location','northwest', 'FontSize', 11);
    title(axE, 'Erro de rastreamento', 'FontSize', 13);
    drawnow;
    % ---------- grava ----------
    img = print(f, '-RGBImage', sprintf('-r%d', dpi), '-opengl');
    if first, sz = size(img, [1 2]); end
    tmp = 255*ones([sz 3], 'uint8'); h = min(sz(1), size(img,1)); wd = min(sz(2), size(img,2));
    tmp(1:h,1:wd,:) = img(1:h,1:wd,:);
    [I, map] = rgb2ind(tmp, 128, 'nodither');
    dly = 0.06 + 1.5*(q == numel(idxs));
    if first
        imwrite(I, map, gif, 'gif', 'LoopCount', Inf, 'DelayTime', dly); first = false;
    else
        imwrite(I, map, gif, 'gif', 'WriteMode', 'append', 'DelayTime', dly);
    end
end
close(f);
fprintf('K_x = [%s], Ki_x = %.4f\n', num2str(Kh, '%.3f '), Kih);

%% ---------------- funções ----------------
function [K, Ki] = projeta(A, B, H, Q, R)
    n = size(A,1);
    Aa = [A zeros(n,1); -H 1];  Ba = [B; 0];
    Ka = dlqr(Aa, Ba, Q, R);            % Ka = [K  -Ki]
    K = Ka(1:n);  Ki = -Ka(n+1);
end

function quad(ax, p, th, ph, col, alfa)
    % desenha o quadricóptero na posição p com arfagem th (eixo y) e rolagem ph (eixo x)
    L = 0.35; rr = 0.13;
    Ry = [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];
    Rx = [1 0 0; 0 cos(-ph) -sin(-ph); 0 sin(-ph) cos(-ph)];
    R = Ry*Rx;
    arms = L*[1 -1 0 0; 0 0 1 -1; 0 0 0 0];
    Wd = R*arms + p;
    plot3(ax, Wd(1,1:2), Wd(2,1:2), Wd(3,1:2), '-', 'Color', [0.1 0.1 0.1 alfa], 'LineWidth', 3);
    plot3(ax, Wd(1,3:4), Wd(2,3:4), Wd(3,3:4), '-', 'Color', [0.1 0.1 0.1 alfa], 'LineWidth', 3);
    tc = linspace(0, 2*pi, 30);
    for i = 1:4
        c = R*(arms(:,i) + rr*[cos(tc); sin(tc); 0*tc] + [0;0;0.03]) + p;
        fill3(ax, c(1,:), c(2,:), c(3,:), col, 'FaceAlpha', 0.45*alfa, 'EdgeColor', col, 'LineWidth', 1.2);
    end
    plot3(ax, p(1), p(2), p(3), 's', 'MarkerSize', 7, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k');
end
