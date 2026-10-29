%% ANIM_VEROSSIMILHANCA  Verossimilhança p(y_k|x_k) no filtro de Kalman (GIF)
% Aula 10 (SEL0359 - Controle Digital, 2026) - Prof. Marcos R. Fernandes
%
% Exemplo: sensor de temperatura com relação linear 1:2 (H = 2) que leu y_k = 80.
%   x_k: temperatura real [C],  y_k: leitura do sensor (y = 2x)
% Modelo de medição:  y_k = H x_k + eps_k,  H = 2,  eps_k ~ N(0, R),  sqrt(R) = 2
%   => p(y_k | x_k) = N(H x_k, R)   (função de DUAS variáveis)
% Etapas da animação:
%   1) x_k fixo  -> p(y_k|x_k) é a distribuição da medição (modelo do sensor);
%   2) y_k medido -> l(x_k) = p(y_k|x_k) vista como função de x_k: a
%      verossimilhança diz quais estados "explicam" bem a medição recebida;
%   3) efeito de R: sensor preciso (R pequeno) -> verossimilhança estreita.
% Gera: ../images/kalman_verossimilhanca.gif

clc; clear; close all;
h = 2;  sR0 = 2;  yObs = 80;                              % H, desvio do sensor, leitura
sMin = 0.9;  pk = 1/(sqrt(2*pi)*sMin);                   % maior pico mostrado (escala dos eixos)
N = @(z, m, v) exp(-(z - m).^2/(2*v))/sqrt(2*pi*v);
x = linspace(34, 46, 301);  y = linspace(68, 92, 301);
[Xg, Yg] = meshgrid(x, y);

%% Figura
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');
co = get(groot, 'defaultAxesColorOrder');
cS = co(1,:);  cY = co(2,:);  cL = co(4,:);
fig = figure('Color', 'w', 'Position', [30 40 1400 800]);
annotation(fig, 'textbox', [0.02 0.87 0.96 0.12], 'String', ...
    ['Sensor de temperatura: \ $y_k = Hx_k + \epsilon_k$, \ $H=2$, \ $\epsilon_k\sim\mathcal{N}(0,R)$' ...
     ' \ $\Rightarrow$ \ $p(y_k|x_k)=\mathcal{N}(Hx_k,\,R)$'], ...
    'Interpreter', 'latex', 'FontSize', 21, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

% ---- painel principal: p(y|x) no plano (x, y)
axM = axes(fig, 'Position', [0.12 0.30 0.46 0.52]);  hold(axM, 'on');
lv   = @(z) round(z*7)/7;
cBl  = reshape([0.10 0.45 0.80], 1, 1, 3);
band = @(sR) 1 - (1 - cBl).*(0.85*lv(N(Yg, h*Xg, sR^2)/N(0, 0, sR^2)));
hImg = image(axM, x, y, band(sR0));
set(axM, 'YDir', 'normal', 'FontSize', 13, 'Box', 'on', 'Layer', 'top', 'XTickLabel', []);
hLn = plot(axM, x, h*x, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.5);
angS = atand(h*(range(x)/range(y))*(0.52*800)/(0.46*1400));
text(axM, 34.5, 2*34.5 + 1.6, '$y_k = 2 x_k$', 'FontSize', 15, 'Rotation', angS, 'Color', [0.2 0.2 0.2]);
hV  = xline(axM, nan, '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 2.5);
hH  = yline(axM, nan, '-', 'Color', cY, 'LineWidth', 2.5);
hPt = plot(axM, nan, nan, 'o', 'MarkerSize', 11, 'MarkerFaceColor', cL, 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
xlim(axM, x([1 end]));  ylim(axM, y([1 end]));
ylabel(axM, '$y_k$: leitura do sensor', 'FontSize', 16);
hTit = title(axM, '', 'FontSize', 17);

% ---- painel direito: p(y|x) para x fixo (eixo vertical = y)
axR = axes(fig, 'Position', [0.60 0.30 0.24 0.52]);  hold(axR, 'on');
hRF = fill(axR, nan, nan, cS, 'FaceAlpha', 0.20, 'EdgeColor', 'none');
hR  = plot(axR, nan, nan, '-', 'Color', cS, 'LineWidth', 2.5);
hRy = yline(axR, nan, '-', 'Color', cY, 'LineWidth', 2.5);
hRp = plot(axR, nan, nan, 'o', 'MarkerSize', 11, 'MarkerFaceColor', cL, 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
hRt = text(axR, 0, 0, '', 'FontSize', 14, 'Color', cY);
xlim(axR, [0 1.15*pk]);  ylim(axR, y([1 end]));
set(axR, 'FontSize', 13, 'Box', 'on', 'XTick', [], 'YAxisLocation', 'right');
hRtit = title(axR, '', 'FontSize', 16, 'Color', cS);

% ---- painel inferior: verossimilhança l(x) = p(y_obs|x)
axB = axes(fig, 'Position', [0.12 0.08 0.46 0.20]);  hold(axB, 'on');
hBF = fill(axB, nan, nan, cL, 'FaceAlpha', 0.20, 'EdgeColor', 'none');
hB  = plot(axB, nan, nan, '-', 'Color', cL, 'LineWidth', 2.5);
hBs = stem(axB, nan, nan, 'filled', 'Color', cL, 'LineWidth', 2, 'MarkerSize', 7, 'MarkerEdgeColor', 'k');
hBt = text(axB, 34.2, 0, '', 'FontSize', 14, 'Color', cL, 'VerticalAlignment', 'top');
hBx = xline(axB, nan, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
xlim(axB, x([1 end]));  ylim(axB, [0 1.3*pk]);
set(axB, 'FontSize', 13, 'Box', 'on', 'YTick', []);
xlabel(axB, '$x_k$: temperatura real [$^{\circ}$C]', 'FontSize', 16);

%% Roteiro de quadros: [etapa, parâmetro]
xsA = linspace(35, 45, 50);                            % etapa 1: varre x
xsB = linspace(x(1), x(end), 70);                         % etapa 2: varre x com y fixo
sRC = sR0 + [linspace(0, -1.1, 18) linspace(-1.1, 2, 28) linspace(2, 0, 18)];   % etapa 3
plan = [ones(1, numel(xsA)) 2*ones(1, numel(xsB)) 3*ones(1, numel(sRC))];
par  = [xsA xsB sRC];
xML  = yObs/h;

outFile = fullfile(fileparts(mfilename('fullpath')), '..', 'images', 'kalman_verossimilhanca.gif');
dpi = 90;  szFr = [];
for n = 1:numel(plan)
    st = plan(n);  p = par(n);
    if st == 1                                            % ---- modelo do sensor
        xs = p;  py = N(y, h*xs, sR0^2);
        set(hTit, 'String', 'Para cada temperatura real $x_k$: leitura $y_k \sim \mathcal{N}(2x_k, R)$');
        set(hV, 'Value', xs);
        set(hR, 'XData', py, 'YData', y);  set(hRF, 'XData', [py 0*py], 'YData', [y fliplr(y)]);
        set(hRtit, 'String', sprintf('$p(y_k|x_k=%.0f\\,^{\\circ}\\mathrm{C})$', xs));
    elseif st == 2                                        % ---- verossimilhança
        xs = p;  py = N(y, h*xs, sR0^2);  L = N(yObs, h*x, sR0^2);  m = x <= xs;
        set(hTit, 'String', 'O sensor leu $y_k=80$: $\ell(x_k)=p(y_k|x_k)$ como fun\c{c}\~ao de $x_k$');
        set(hV, 'Value', xs);  set(hH, 'Value', yObs);  set(hRy, 'Value', yObs);
        set(hPt, 'XData', xs, 'YData', yObs);
        set(hR, 'XData', py, 'YData', y);  set(hRF, 'XData', [py 0*py], 'YData', [y fliplr(y)]);
        set(hRp, 'XData', N(yObs, h*xs, sR0^2), 'YData', yObs);
        set(hRt, 'Position', [0.45*pk yObs - 1.2 0], 'String', sprintf('$y_k = %.0f$ (leitura)', yObs));
        set(hRtit, 'String', sprintf('$p(y_k|x_k=%.0f\\,^{\\circ}\\mathrm{C})$', xs));
        set(hB, 'XData', x(m), 'YData', L(m));  set(hBF, 'XData', [x(m) fliplr(x(m))], 'YData', [L(m) 0*L(m)]);
        set(hBs, 'XData', xs, 'YData', N(yObs, h*xs, sR0^2));
        set(hBt, 'Position', [34.2 1.25*pk 0], 'String', 'Quais temperaturas reais $x_k$ explicam bem a leitura $y_k=80$?');
        if n == find(plan == 2, 1, 'last')
            set(hBx, 'Value', xML);  set([hBs hPt hRp], 'XData', nan, 'YData', nan);  set(hV, 'Value', xML);
            set(hBt, 'String', {'$\ell(x_k)$: m\''axima em $x_k = y_k/H = 80/2 = 40\,^{\circ}\mathrm{C}$', 'temperaturas pr\''oximas de $40\,^{\circ}\mathrm{C}$ s\~ao as mais plaus\''iveis'});
            set(hR, 'XData', N(y, yObs, sR0^2));  set(hRF, 'XData', [N(y, yObs, sR0^2) 0*y]);
            set(hRtit, 'String', sprintf('$p(y_k|x_k=%.0f\\,^{\\circ}\\mathrm{C})$', xML));
        end
    else                                                  % ---- efeito de R
        sR = p;  L = N(yObs, h*x, sR^2);  py = N(y, yObs, sR^2);
        set(hImg, 'CData', band(sR));
        set(hTit, 'String', sprintf('Efeito de $R$: \\ $\\sqrt{R} = %.1f$ \\ (%s)', sR, ...
            ternary(sR < sR0 - 0.05, 'sensor preciso $\Rightarrow$ $\ell$ estreita', ...
            ternary(sR > sR0 + 0.05, 'sensor ruidoso $\Rightarrow$ $\ell$ larga', 'valor nominal'))));
        set(hR, 'XData', py, 'YData', y);  set(hRF, 'XData', [py 0*py], 'YData', [y fliplr(y)]);
        set(hB, 'XData', x, 'YData', L);  set(hBF, 'XData', [x fliplr(x)], 'YData', [L 0*L]);
        set(hBt, 'String', '$\ell(x_k)=p(y_k|x_k)$');
    end
    drawnow;

    img = print(fig, '-RGBImage', sprintf('-r%d', dpi), '-opengl');
    if isempty(szFr), szFr = size(img, [1 2]); end
    tmp = 255*ones([szFr 3], 'uint8');
    hh = min(szFr(1), size(img,1));  wd = min(szFr(2), size(img,2));
    tmp(1:hh, 1:wd, :) = img(1:hh, 1:wd, :);
    [Ai, map] = rgb2ind(tmp, 255, 'nodither');
    last = n == numel(plan) || plan(min(n+1, end)) ~= st;   % pausa no fim de cada etapa
    dly = 0.07 + 0.8*(n == 1) + 3.0*last + 1.5*(n == numel(plan));
    if n == 1
        imwrite(Ai, map, outFile, 'gif', 'LoopCount', Inf, 'DelayTime', dly);
    else
        imwrite(Ai, map, outFile, 'gif', 'WriteMode', 'append', 'DelayTime', dly);
    end
end
close(fig);
fprintf('Máxima verossimilhança em x = y/H = %.1f C\n', xML);
fprintf('%d quadros. GIF salvo em: %s\n', numel(plan), outFile);

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
