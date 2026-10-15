%% Figuras da apostila da Aula 6 (SEL0359 - Controle Digital, 2026)
% Projeto direto em z pelo lugar geométrico das raízes (LGR)
% Gera PNGs em ../images/   (Prof. Marcos R. Fernandes)

clc; clear all; close all;
out = fullfile(fileparts(mfilename('fullpath')), '..', 'images');
set(groot, 'defaultAxesFontSize', 13, 'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
cB = [0.165 0.494 0.878]; cR = [0.878 0.141 0.141]; cG = [0.122 0.620 0.122];

casos = {struct('G', tf(1,[1 1 0]), 'Ts', 0.1, 'xi', 0.7, 'wn', 5, 'nome', 'aula6_exemplo'), ...
         struct('G', tf(1,[1 5 0]), 'Ts', 0.2, 'xi', -log(0.1)/sqrt(pi^2+log(0.1)^2), 'wn', 5, 'nome', 'aula6_exercicio')};

for c = 1:2
    P = casos{c}; Ts = P.Ts;
    Gd = c2d(P.G, Ts, 'zoh');                       % 1o passo
    s1 = -P.xi*P.wn + 1j*P.wn*sqrt(1-P.xi^2);       % 2o passo
    z1 = exp(s1*Ts);                                % 3o passo
    [zG, pG, kG] = zpkdata(Gd, 'v');
    a  = pG(abs(pG) < 0.999);                       % polo estável a cancelar
    th = angle(z1 - zG);  ph = angle(z1 - 1);
    al = th - ph - pi;                              % condição de fase
    b  = real(z1) - imag(z1)/tan(al);
    K  = abs((z1-b)*(z1-1)) / abs(kG*(z1-zG));      % condição de módulo
    C  = zpk(a, b, K, Ts);
    T  = feedback(C*Gd, 1);   U = feedback(C, Gd);
    S  = stepinfo(T, 'SettlingTimeThreshold', 0.02);
    fprintf('[%s] a=%.4f b=%.4f K=%.4f  Mp=%.2f%%  ts=%.2f s\n', P.nome, a, b, K, S.Overshoot, S.SettlingTime);

    % ---------- LGR ----------
    f = figure('Color','w','Position',[100 100 1100 470]);
    tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    % LGR calculado diretamente: raízes de D(z) + k N(z), k >= 0
    [Nz, Dz] = tfdata(minreal(C/K*Gd), 'v');
    Nz = Nz(find(abs(Nz) > 1e-12, 1):end);
    kk = [0 logspace(-3, 3, 4000)];
    R = zeros(numel(Dz)-1, numel(kk));
    for i = 1:numel(kk), R(:,i) = roots(Dz + kk(i)*[zeros(1,numel(Dz)-numel(Nz)) Nz]); end
    plot(ax, real(R(:)), imag(R(:)), '.', 'Color', cB, 'MarkerSize', 3);
    tc = linspace(0, 2*pi, 400); plot(ax, cos(tc), sin(tc), 'k-', 'LineWidth', 1);
    pG2 = roots(Dz); zG2 = roots(Nz);
    plot(ax, real(pG2), imag(pG2), 'kx', 'MarkerSize', 11, 'LineWidth', 2);
    plot(ax, real(zG2), imag(zG2), 'ko', 'MarkerSize', 9, 'LineWidth', 1.5);
    plot(ax, real(z1), [imag(z1) -imag(z1)], 'p', 'MarkerSize', 16, 'MarkerFaceColor', cG, 'MarkerEdgeColor', 'k');
    pcl = pole(T);
    plot(ax, real(pcl), imag(pcl), 's', 'MarkerSize', 11, 'MarkerEdgeColor', cR, 'LineWidth', 1.8);
    axis(ax, 'equal'); xlim(ax, [-2.8 1.4]); ylim(ax, [-1.7 1.7]);
    xlabel(ax, 'Re$\{z\}$'); ylabel(ax, 'Im$\{z\}$');
    title(ax, sprintf('LGR de $C(z)G_d(z)$ com $C(z)=K\\frac{z-%.4f}{z-%.3f}$', a, b), 'FontSize', 14);
    ax2 = nexttile;
    [y, t] = step(T, 3); [u, tu] = step(U, 3);
    yyaxis(ax2, 'left');  stairs(ax2, t, y, 'Color', cB, 'LineWidth', 2); ylabel(ax2, '$y_k$');
    yline(ax2, 1, ':k'); yline(ax2, 1 + S.Overshoot/100, '--', 'Color', [0.5 0.5 0.5]);
    yyaxis(ax2, 'right'); stairs(ax2, tu, u, 'Color', cR, 'LineWidth', 1.5); ylabel(ax2, '$u_k$');
    ax2.YAxis(1).Color = cB; ax2.YAxis(2).Color = cR; grid(ax2, 'on');
    xlabel(ax2, '$t=kT_s$ [s]');
    title(ax2, sprintf('Degrau em malha fechada: $M_p=%.1f\\%%$, $t_s=%.2f$ s', S.Overshoot, S.SettlingTime), 'FontSize', 14);
    exportgraphics(f, fullfile(out, [P.nome '_lgr_resposta.png']), 'Resolution', 200);
    close(f);
end
