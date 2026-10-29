%% ANIM_CHAPMAN_KOLMOGOROV  Passo de predição do filtro de Kalman (GIF)
% Aula 10 (SEL0359 - Controle Digital, 2026) - Prof. Marcos R. Fernandes
%
% Equação de Chapman-Kolmogorov (caso escalar):
%   p(x_k | y_{1:k-1}) = INT p(x_k | x_{k-1}) p(x_{k-1} | y_{1:k-1}) dx_{k-1}
% com
%   p(x_{k-1} | y_{1:k-1}) = N(xh, P)                     (estimativa anterior)
%   p(x_k | x_{k-1})       = N(a x_{k-1} + b u, q)         (modelo de transição)
% A integral é visualizada como uma "soma" de Gaussianas de transição, cada uma
% centrada em a x_{k-1} + b u e ponderada por p(x_{k-1}|y_{1:k-1}). O resultado
% é a Gaussiana N(a xh + b u, a^2 P + q): a média é propagada pelo modelo e a
% variância aumenta (a^2 P + q).
% Gera: ../images/kalman_chapman_kolmogorov.gif

clc; clear; close all;
a = 1.1;  bu = 1.0;  xh = 0;  P = 0.5^2;  q = 0.45^2;     % parâmetros
mP = a*xh + bu;  PP = a^2*P + q;                           % predição analítica
N  = @(x, m, v) exp(-(x - m).^2/(2*v))/sqrt(2*pi*v);

xs = linspace(-2.0, 2.4, 301);                             % x_{k-1}
xk = linspace(-1.6, 3.6, 301);                             % x_k
ds = xs(2) - xs(1);
pPrev = N(xs, xh, P);
[S, Xk] = meshgrid(xs, xk);
Tr = N(Xk, a*S + bu, q);                                   % p(x_k|x_{k-1})
J  = Tr.*N(S, xh, P);                                      % p(x_k, x_{k-1}|y)
pPred = N(xk, mP, PP);

%% Figura
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');
co = get(groot, 'defaultAxesColorOrder');
cA = co(1,:);  cT = co(2,:);  cP = co(5,:);  cG = [0.55 0.55 0.55];
fig = figure('Color', 'w', 'Position', [30 40 1400 800]);

% ---- painel principal: transição e conjunta
axM = axes(fig, 'Position', [0.12 0.30 0.46 0.52]);  hold(axM, 'on');
% imagem RGB composta: transição (cinza) + conjunta revelada (laranja)
lv   = @(z) round(z*7)/7;                                % níveis (GIF com poucas cores)
Tg   = repmat(1 - 0.40*lv(Tr/max(Tr(:))), 1, 1, 3);
cOr  = reshape([0.85 0.30 0.08], 1, 1, 3);
comp = @(Aj) Tg.*(1 - Aj) + cOr.*Aj;
hJ   = image(axM, xs, xk, comp(zeros(size(J))));
set(axM, 'YDir', 'normal', 'FontSize', 13, 'Box', 'on', 'Layer', 'top');
plot(axM, xs, a*xs + bu, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
angS = atand(a*(range(xs)/range(xk))*(0.52*800)/(0.46*1400));               % inclinação na tela
text(axM, -1.85, a*(-1.85) + bu + 0.35, '$x_k = a\,x_{k-1} + b\,u_{k-1}$', 'FontSize', 15, ...
     'Rotation', angS, 'Color', [0.25 0.25 0.25]);
hV  = xline(axM, nan, '-', 'Color', cA, 'LineWidth', 2.5);
hC  = plot(axM, nan, nan, '-', 'Color', cT, 'LineWidth', 2.5);
xlim(axM, xs([1 end]));  ylim(axM, xk([1 end]));
set(axM, 'XTickLabel', [], 'YTickLabel', []);
title(axM, '$p(x_k|x_{k-1})$ \ (modelo de transi\c{c}\~ao)', 'FontSize', 17);
hJt = text(axM, 2.3, -1.45, '', 'FontSize', 14, 'HorizontalAlignment', 'right', 'Color', [0.6 0.15 0.05]);

% ---- painel inferior: p(x_{k-1}|y_{1:k-1})
axB = axes(fig, 'Position', [0.12 0.08 0.46 0.20]);  hold(axB, 'on');
fill(axB, [xs fliplr(xs)], [pPrev zeros(size(xs))], cA, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(axB, xs, pPrev, '-', 'Color', cA, 'LineWidth', 2.5);
hSt = stem(axB, nan, nan, 'filled', 'Color', cA, 'LineWidth', 2.5, 'MarkerSize', 7);
xlim(axB, xs([1 end]));  ylim(axB, [0 1.35*max(pPrev)]);  set(axB, 'FontSize', 13, 'Box', 'on', 'YTick', []);
xlabel(axB, '$x_{k-1}$', 'FontSize', 17);
text(axB, -1.95, 1.2*max(pPrev), '$p(x_{k-1}|y_{1:k-1})=\mathcal{N}(\hat x_{k-1|k-1},P_{k-1|k-1})$', ...
     'FontSize', 13, 'Color', cA);

% ---- painel direito: p(x_k|y_{1:k-1}) (eixo vertical = x_k)
axR = axes(fig, 'Position', [0.60 0.30 0.26 0.52]);  hold(axR, 'on');
hAccF = fill(axR, nan, nan, cP, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
hAcc  = plot(axR, nan, nan, '-', 'Color', cP, 'LineWidth', 3);
hW    = plot(axR, nan, nan, '-', 'Color', cT, 'LineWidth', 2);
hAn   = plot(axR, nan, nan, 'k--', 'LineWidth', 1.8);
xmaxR = 1.25*max(max(pPred), max(max(Tr.*N(S, xh, P))));
xlim(axR, [0 xmaxR]);  ylim(axR, xk([1 end]));
set(axR, 'FontSize', 13, 'Box', 'on', 'XTick', [], 'YAxisLocation', 'right');
ylabel(axR, '$x_k$', 'FontSize', 17);
title(axR, '$p(x_k|y_{1:k-1})$', 'FontSize', 17, 'Color', cP);
legend(axR, [hW hAcc], {'$p(x_k|x_{k-1})\,p(x_{k-1}|y_{1:k-1})$', '$\sum (\cdot)\,\Delta x_{k-1}$ acumulada'}, ...
       'Location', 'southoutside', 'FontSize', 12);

% ---- textos
annotation(fig, 'textbox', [0.02 0.86 0.96 0.13], 'String', ...
    ['Chapman--Kolmogorov: \ $p(x_k|y_{1:k-1})=\displaystyle\int ' ...
     '\underbrace{p(x_k|x_{k-1})}_{\mbox{transi\c{c}\~ao}}\,' ...
     '\underbrace{p(x_{k-1}|y_{1:k-1})}_{\mbox{estimativa anterior}}\,dx_{k-1}$'], ...
    'Interpreter', 'latex', 'FontSize', 20, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Animação
outFile = fullfile(fileparts(mfilename('fullpath')), '..', 'images', 'kalman_chapman_kolmogorov.gif');
dpi = 90;  szFr = [];  frame = 0;
iSweep = round(linspace(1, numel(xs), 110));
nIntro = 3;  nTot = nIntro + numel(iSweep) + 1;
acc = zeros(size(xk));  iPrev = 0;
for n = 1:nTot
    if n <= nIntro                                         % introdução
    elseif n <= nIntro + numel(iSweep)                     % varredura (integral como soma)
        i = iSweep(n - nIntro);
        acc = acc + sum(Tr(:, iPrev+1:i).*pPrev(iPrev+1:i), 2)'*ds;  iPrev = i;
        s = xs(i);  w = pPrev(i);  cond = N(xk, a*s + bu, q);
        Aj = 0.9*lv(J/max(J(:)));  Aj(:, i+1:end) = 0;
        set(hJ, 'CData', comp(Aj));
        set(hV, 'Value', s);
        set(hC, 'XData', s + 0.0*xk, 'YData', xk);
        set(hSt, 'XData', s, 'YData', w);
        set(hW, 'XData', w*cond, 'YData', xk);
        set(hAcc, 'XData', acc, 'YData', xk);
        set(hAccF, 'XData', [acc zeros(size(acc))], 'YData', [xk fliplr(xk)]);
        set(hJt, 'String', 'conjunta $p(x_k,x_{k-1}|y_{1:k-1})$');
    else                                                   % resultado final
        set([hW hC], 'XData', nan, 'YData', nan);  set(hSt, 'XData', nan, 'YData', nan);
        set(hV, 'Value', nan);
        set(hAn, 'XData', pPred, 'YData', xk);
        legend(axR, [hAcc hAn], {'integral (soma completa)', '$\mathcal{N}(a\hat x+bu,\ a^2P+Q)$'}, ...
               'Location', 'southoutside', 'FontSize', 12);
    end
    drawnow;

    img = print(fig, '-RGBImage', sprintf('-r%d', dpi), '-opengl');
    if isempty(szFr), szFr = size(img, [1 2]); end
    tmp = 255*ones([szFr 3], 'uint8');
    h = min(szFr(1), size(img,1));  wd = min(szFr(2), size(img,2));
    tmp(1:h, 1:wd, :) = img(1:h, 1:wd, :);
    [Ai, map] = rgb2ind(tmp, 255, 'nodither');
    dly = 0.07 + 0.35*(n <= nIntro) + 6.0*(n == nTot);
    frame = frame + 1;
    if frame == 1
        imwrite(Ai, map, outFile, 'gif', 'LoopCount', Inf, 'DelayTime', dly);
    else
        imwrite(Ai, map, outFile, 'gif', 'WriteMode', 'append', 'DelayTime', dly);
    end
end
close(fig);
fprintf('Predição: média %.3f -> %.3f,  variância %.3f -> %.3f\n', xh, mP, P, PP);
fprintf('Erro da integral numérica: %.2e\n', max(abs(acc - pPred)));
fprintf('%d quadros. GIF salvo em: %s\n', frame, outFile);

