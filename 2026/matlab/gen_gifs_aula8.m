%% GIFs ilustrativos: controlabilidade e observabilidade (Aula 8 - SEL0359, 2026)
% Gera ../images/aula8_controlabilidade.gif e ../images/aula8_observabilidade.gif
% Prof. Marcos R. Fernandes

clear; close all;
out = fullfile(fileparts(mfilename('fullpath')), '..', 'images');
set(groot, 'defaultTextInterpreter','latex', 'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex', 'defaultAxesFontSize', 13);
cor.B = [0.165 0.494 0.878]; cor.R = [0.878 0.141 0.141]; cor.G = [0.122 0.620 0.122];
cor.Gr = [0.55 0.55 0.55];   cor.K = [0.15 0.15 0.15];
dpi = 130;

%% =================== 1) CONTROLABILIDADE ===================
x0 = [-2.5; -1];  xf = [1; 2.5];
cs(1).A = 0.9*[0 -1; 1 0];  cs(1).B = [1; 0];   cs(1).nome = 'Control\''avel';
cs(2).A = 0.8*eye(2);       cs(2).B = [1; 1];   cs(2).nome = 'N\~ao control\''avel';
for c = 1:2
    A = cs(c).A; B = cs(c).B;
    C = [B A*B];
    v = pinv(C)*(xf - A^2*x0);          % v = [u1; u0] (mínima norma / mínimos quadrados)
    u = flipud(v);
    X = zeros(2,3); X(:,1) = x0;
    for k = 1:2, X(:,k+1) = A*X(:,k) + B*u(k); end
    cs(c).u = u; cs(c).X = X; cs(c).r = rank(C);
end

f = figure('Color','w','Position',[50 50 1200 560]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax = [nexttile(tl,1) nexttile(tl,2)];
gif = fullfile(out,'aula8_controlabilidade.gif'); first = true;
nF = 12;                                       % quadros por fase
seq = [];                                      % [passo k, fase, t]
for k = 1:2
    for ph = 1:3
        for t = linspace(0,1,nF), seq(end+1,:) = [k ph t]; end %#ok<AGROW>
    end
end
seq(end+1:end+25,:) = repmat([3 0 1],25,1);    % pausa final com conclusão
seq = [repmat([1 0 0],10,1); seq];              % pausa inicial

for q = 1:size(seq,1)
    k = seq(q,1); ph = seq(q,2); t = seq(q,3);
    for c = 1:2
        desenha_ctrb(ax(c), cs(c), x0, xf, k, ph, t, cor);
    end
    drawnow;
    first = grava(f, gif, first, dpi, 0.06 + 1.2*(q==size(seq,1)));
end
close(f);

%% =================== 2) OBSERVABILIDADE ===================
x0 = [1.5; -1];  H = [1 0];
os(1).A = [2 1; 5 -3];  os(1).nome = 'Observ\''avel';
os(2).A = 0.8*eye(2);   os(2).nome = 'N\~ao observ\''avel';
for c = 1:2
    os(c).O = [H; H*os(c).A]; os(c).y = os(c).O*x0; os(c).r = rank(os(c).O);
end
f = figure('Color','w','Position',[50 50 1200 560]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax = [nexttile(tl,1) nexttile(tl,2)];
gif = fullfile(out,'aula8_observabilidade.gif'); first = true;
seq = [repmat([0 0],8,1)];
for ph = 1:3, for t = linspace(0,1,14), seq(end+1,:) = [ph t]; end, end   %#ok<AGROW>
for t = linspace(0,1,40), seq(end+1,:) = [4 t]; end                     %#ok<AGROW>
for t = linspace(1,0,40), seq(end+1,:) = [4 t]; end                     %#ok<AGROW>
seq(end+1:end+20,:) = repmat([5 1],20,1);
for q = 1:size(seq,1)
    for c = 1:2
        desenha_obsv(ax(c), os(c), x0, H, seq(q,1), seq(q,2), cor);
    end
    drawnow;
    first = grava(f, gif, first, dpi, 0.07 + 1.2*(q==size(seq,1)));
end
close(f);
disp('GIFs gerados.');

%% =================== funções ===================
function first = grava(f, arq, first, dpi, dly)
    img = print(f, '-RGBImage', sprintf('-r%d', dpi), '-opengl');
    persistent sz
    if first || isempty(sz), sz = size(img,[1 2]); end
    tmp = 255*ones([sz 3],'uint8'); h = min(sz(1),size(img,1)); w = min(sz(2),size(img,2));
    tmp(1:h,1:w,:) = img(1:h,1:w,:);
    [I,map] = rgb2ind(tmp, 64, 'nodither');
    if first
        imwrite(I,map,arq,'gif','LoopCount',Inf,'DelayTime',dly); first = false;
    else
        imwrite(I,map,arq,'gif','WriteMode','append','DelayTime',dly);
    end
end

function seta(ax, p, d, col, lw, estilo)
    if norm(d) < 1e-6, return; end
    quiver(ax, p(1), p(2), d(1), d(2), 0, 'Color', col, 'LineWidth', lw, ...
        'MaxHeadSize', min(0.6, 0.35/norm(d)+0.08), 'LineStyle', estilo);
end

function desenha_ctrb(ax, s, x0, xf, k, ph, t, cor)
    cla(ax); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
    xlim(ax,[-4 4]); ylim(ax,[-4 4]);
    A = s.A; B = s.B; X = s.X; u = s.u;
    % direções de controle disponíveis: span{B, AB}
    if s.r == 1
        L = A^2*x0 + B*linspace(-6,6,2)/norm(B);
        plot(ax, L(1,:), L(2,:), '--', 'Color', [cor.G 0.8], 'LineWidth', 2);
        text(ax, -3.8, -3.3, '$x_2\in A^2x_0+\mathrm{span}\{B\}$', 'FontSize', 14, 'Color', cor.G);
    end
    plot(ax, xf(1), xf(2), 'p', 'MarkerSize', 22, 'MarkerFaceColor', cor.G, 'MarkerEdgeColor', 'k');
    text(ax, xf(1)+0.3, xf(2)+0.35, '$x_f$', 'FontSize', 17, 'Color', cor.G);
    % trajetória já percorrida
    kk = min(k,3);
    plot(ax, X(1,1:kk), X(2,1:kk), '-o', 'Color', cor.B, 'LineWidth', 2.2, 'MarkerFaceColor', cor.B, 'MarkerSize', 8);
    text(ax, x0(1)-0.2, x0(2)-0.45, '$x_0$', 'FontSize', 16);
    if k <= 2 && ph >= 1
        xk = X(:,k); Ax = A*xk; Bu = B*u(k);
        tt1 = (ph==1)*t + (ph>1);
        seta(ax, xk, tt1*(Ax - xk), cor.Gr, 2, '--');
        if ph >= 2
            tt2 = (ph==2)*t + (ph>2);
            seta(ax, Ax, tt2*Bu, cor.R, 2.8, '-');
            % direção de B (reta tracejada fina)
            dB = B/norm(B);
            plot(ax, Ax(1)+[-1 1]*5*dB(1), Ax(2)+[-1 1]*5*dB(2), ':', 'Color', [cor.R 0.5], 'LineWidth', 1.2);
        end
        if ph == 3
            xn = Ax + Bu;
            plot(ax, xn(1), xn(2), 'o', 'MarkerSize', 8 + 4*t, 'MarkerFaceColor', cor.B, 'MarkerEdgeColor', 'k');
        end
        texto = sprintf('passo %d: $x_{%d}=\\underbrace{Ax_{%d}}_{\\mathrm{deriva}}+\\underbrace{Bu_{%d}}_{\\mathrm{controle}}$, $u_{%d}=%.2f$', k, k, k-1, k-1, k-1, u(k));
    else
        texto = '';
    end
    if k > 2
        err = norm(X(:,3) - xf);
        if err < 1e-6
            texto = '$x_2=x_f$: $B$ e $AB$ geram o plano $\Rightarrow$ qualquer $x_f$';
        else
            plot(ax, [X(1,3) xf(1)], [X(2,3) xf(2)], ':', 'Color', 'k', 'LineWidth', 1.5);
            texto = sprintf('$x_2\\neq x_f$ (erro $%.2f$): $B$ e $AB$ paralelos', err);
        end
    end
    xlabel(ax,'$x_1$'); ylabel(ax,'$x_2$');
    title(ax, {sprintf('%s: $\\mathrm{rank}\\,\\mathcal{C}=%d$', s.nome, s.r), texto}, 'FontSize', 14);
    % legenda manual
    text(ax, 1.0, -3.0, '-- deriva $Ax_k$', 'FontSize', 13, 'Color', cor.Gr);
    text(ax, 1.0, -3.6, '$\rightarrow$ controle $Bu_k$', 'FontSize', 13, 'Color', cor.R);
end

function desenha_obsv(ax, s, x0, H, ph, t, cor)
    cla(ax); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
    xlim(ax,[-4 4]); ylim(ax,[-4 4]);
    A = s.A; y = s.y;
    r1 = H;  r2 = H*A;                       % restrições: r_i * x0 = y_i
    lin = @(r, yy, tt) linha(r, yy, tt);
    if ph >= 1
        tt = (ph==1)*t + (ph>1);
        P = lin(r1, y(1), tt); plot(ax, P(1,:), P(2,:), '-', 'Color', cor.B, 'LineWidth', 3);
        text(ax, y(1)+0.15, 3.4, '$Hx_0=y_0$', 'FontSize', 15, 'Color', cor.B);
    end
    if ph >= 2
        tt = (ph==2)*t + (ph>2);
        P = lin(r2, y(2), tt);
        if s.r == 2
            plot(ax, P(1,:), P(2,:), '-', 'Color', cor.R, 'LineWidth', 3);
            text(ax, -3.7, -1.3, '$HAx_0=y_1$', 'FontSize', 15, 'Color', cor.R);
        else
            plot(ax, P(1,:), P(2,:), '--', 'Color', cor.R, 'LineWidth', 3);
            text(ax, y(1)+0.15, 2.8, '$HAx_0=y_1$ (mesma reta!)', 'FontSize', 15, 'Color', cor.R);
        end
    end
    if ph >= 3
        if s.r == 2
            x0h = s.O \ y;
            plot(ax, x0h(1), x0h(2), 'p', 'MarkerSize', 12 + 12*min(1,(ph==3)*t + (ph>3)), 'MarkerFaceColor', cor.G, 'MarkerEdgeColor', 'k');
            text(ax, x0h(1)+0.35, x0h(2)-0.45, '$x_0=\mathcal{O}^{-1}\Delta y$', 'FontSize', 15, 'Color', cor.G);
        else
            plot(ax, [y(1) y(1)], [-4 4], '-', 'Color', [cor.G 0.25], 'LineWidth', 14*min(1,(ph==3)*t + (ph>3)) + 0.1);
        end
    end
    txt = '';
    if ph == 4 || ph == 5
        % candidato a condição inicial deslizando sobre a reta y0
        cand = [y(1); -3.5 + 7*t];
        yh1 = H*A*cand;
        ok = abs(yh1 - y(2)) < 0.08;
        col = cor.R; if ok, col = cor.G; end
        plot(ax, cand(1), cand(2), 'o', 'MarkerSize', 13, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
        if ok
            txt = sprintf('candidato $\\hat x_0$: $\\hat y_1=%.2f=y_1$ (ok)', yh1);
        else
            txt = sprintf('candidato $\\hat x_0$: $\\hat y_1=%.2f\\neq y_1=%.2f$', yh1, y(2));
        end
    end
    if ph == 5
        if s.r == 2
            txt = 'Um \''unico $x_0$ explica as medidas $y_0,y_1$';
        else
            txt = 'Infinitos $x_0$ explicam as mesmas medidas';
        end
    end
    xlabel(ax,'$x_{0,1}$'); ylabel(ax,'$x_{0,2}$');
    title(ax, {sprintf('%s: $\\mathrm{rank}\\,\\mathcal{O}=%d$', s.nome, s.r), txt}, 'FontSize', 14);
end

function P = linha(r, yy, tt)
    % segmento da reta r*x = yy dentro de [-4,4]^2, desenhado até a fração tt
    if abs(r(2)) < 1e-9
        P = [yy/r(1) yy/r(1); -4 -4+8*tt];
    else
        xs = linspace(-4, 4, 200);
        ys = (yy - r(1)*xs)/r(2);
        m = ys >= -4 & ys <= 4; xs = xs(m); ys = ys(m);
        n = max(2, round(tt*numel(xs)));
        P = [xs(1:n); ys(1:n)];
    end
end
