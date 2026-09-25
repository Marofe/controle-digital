%% Figuras ilustrativas da Aula 8 (SEL0359 - Controle Digital, 2026)
% Controlabilidade, observabilidade, alocação de polos, observador e rastreamento.
% Gera PNGs em ../images/  (Prof. Marcos R. Fernandes)

clear; close all;
out = fullfile(fileparts(mfilename('fullpath')), '..', 'images');
cB = [0.165 0.494 0.878];  cR = [0.878 0.141 0.141];  cG = [0.122 0.620 0.122];
cK = [0.2 0.2 0.2];        cGr = [0.6 0.6 0.6];
fs = 13;
salva = @(f, nome) exportgraphics(f, fullfile(out, nome), 'Resolution', 200);
set(groot, 'defaultAxesFontSize', fs, 'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');

%% 1) Controlabilidade: levar x0 a xf em n = 2 passos
A = [2 1; 5 -3];  B = [0; 1];
C = [B A*B];                         % C = [B AB]
x0 = [1; -1];  xf = [2; 1];
dx = xf - A^2*x0;                    % Delta x = x_2 - A^2 x_0
v  = C \ dx;                         % v = [u1; u0]  (ordem invertida!)
u  = flipud(v);                      % u = [u0; u1]
x1 = A*x0 + B*u(1);  x2 = A*x1 + B*u(2);
xl = A*x0;  xl2 = A*xl;              % evolução livre (u = 0)

f = figure('Color','w','Position',[100 100 1100 440]);
tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
plot(ax,[x0(1) xl(1) xl2(1)],[x0(2) xl(2) xl2(2)],'--o','Color',cGr,'LineWidth',1.5,'MarkerFaceColor',cGr);
plot(ax,[x0(1) x1(1) x2(1)],[x0(2) x1(2) x2(2)],'-o','Color',cB,'LineWidth',2.5,'MarkerFaceColor',cB,'MarkerSize',8);
plot(ax,xf(1),xf(2),'p','MarkerSize',20,'MarkerFaceColor',cG,'MarkerEdgeColor','k');
text(ax,x0(1)+0.2,x0(2)-0.35,'$x_0$','FontSize',16);
text(ax,x1(1)-0.55,x1(2)+0.1,'$x_1$','FontSize',16,'Color',cB);
text(ax,xf(1)+0.3,xf(2)+0.4,'$x_2=x_f$','FontSize',16,'Color',cG);
text(ax,1.25,3.3,'livre ($u=0$): diverge','FontSize',14,'Color',cGr);
xlim(ax,[-2 4]); ylim(ax,[-2.5 4]);
xlabel(ax,'$x_1$'); ylabel(ax,'$x_2$');
title(ax,'Espa\c{c}o de estados: $x_0\to x_f$ em $n=2$ passos','FontSize',15);
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
stem(ax2,[0 1],u,'filled','Color',cR,'LineWidth',2.5,'MarkerSize',8);
xlim(ax2,[-0.5 1.5]); xticks(ax2,[0 1]); xticklabels(ax2,{'$u_0$','$u_1$'});
for i = 1:2
    text(ax2,i-1+0.08,u(i),sprintf('$%.2f$',u(i)),'FontSize',15,'Color',cR);
end
yline(ax2,0,'k');
title(ax2,'$\left[\begin{array}{c}u_1\\u_0\end{array}\right]=\mathcal{C}^{-1}(x_f-A^2x_0)$','FontSize',15);
salva(f,'aula8_controlabilidade_passos.png'); close(f);
fprintf('u0 = %.4f, u1 = %.4f, x2 = [%g %g]\n', u(1), u(2), x2);

%% 2) Controlável x não controlável: estados alcançáveis a partir da origem
rng(3); N = 2; nTraj = 150;
casos = {struct('A',[0.6 0.8;-0.8 0.6],'B',[1;0],'tit','Control\''avel: $\mathrm{rank}\,\mathcal{C}=2$'), ...
         struct('A',0.8*eye(2),     'B',[1;1],'tit','N\~ao control\''avel: $\mathrm{rank}\,\mathcal{C}=1$')};
f = figure('Color','w','Position',[100 100 1100 470]);
tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
for c = 1:2
    Ac = casos{c}.A; Bc = casos{c}.B;
    ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
    for t = 1:nTraj
        x = zeros(2,N+1); uu = 2*(2*rand(1,N)-1);
        for k = 1:N, x(:,k+1) = Ac*x(:,k) + Bc*uu(k); end
        plot(ax,x(1,:),x(2,:),'-','Color',[cB 0.25],'LineWidth',1);
        plot(ax,x(1,end),x(2,end),'o','MarkerFaceColor',cB,'MarkerEdgeColor','none','MarkerSize',5);
    end
    if c == 2
        plot(ax,[-4 4],[-4 4],'--','Color',cG,'LineWidth',2);
        text(ax,1.2,2.9,'$\mathcal{R}\{\mathcal{C}\}=\mathrm{span}\{B\}$','FontSize',15,'Color',cG);
    end
    plot(ax,0,0,'k+','MarkerSize',14,'LineWidth',2);
    xlim(ax,[-4 4]); ylim(ax,[-4 4]);
    xlabel(ax,'$x_1$'); ylabel(ax,'$x_2$');
    title(ax,casos{c}.tit,'FontSize',15);
end
salva(f,'aula8_controlavel_vs_nao.png'); close(f);

%% 3) Observabilidade: cada medida y_k define uma reta de restrição para x0
A = [2 1; 5 -3]; H = [1 0]; x0 = [1.5; -1];
O = [H; H*A]; y = O*x0;                  % y0, y1 (u = 0)
x0h = O \ y;
f = figure('Color','w','Position',[100 100 1100 470]);
tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
s = linspace(-4,4,2);
xline(ax,y(1),'-','Color',cB,'LineWidth',2.5);                         % x1 = y0
plot(ax,s,(y(2)-2*s)/1,'-','Color',cR,'LineWidth',2.5);                 % 2x1 + x2 = y1
plot(ax,x0h(1),x0h(2),'p','MarkerSize',20,'MarkerFaceColor',cG,'MarkerEdgeColor','k');
text(ax,y(1)+0.15,3.3,'$Hx_0=y_0$','FontSize',15,'Color',cB);
text(ax,-3.6,-1.2,'$HAx_0=y_1$','FontSize',15,'Color',cR);
text(ax,x0h(1)+0.3,x0h(2)-0.5,'$x_0=\mathcal{O}^{-1}\Delta y$','FontSize',15,'Color',cG);
xlim(ax,[-4 4]); ylim(ax,[-4 4]); xlabel(ax,'$x_{0,1}$'); ylabel(ax,'$x_{0,2}$');
title(ax,'Observ\''avel: retas se cruzam (solu\c{c}\~ao \''unica)','FontSize',15);
% não observável: A = 0.8 I, H = [1 0]
A2 = 0.8*eye(2); O2 = [H; H*A2]; y2 = O2*x0;
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
xline(ax,y2(1),'-','Color',cB,'LineWidth',5);
xline(ax,y2(2)/0.8,'--','Color',cR,'LineWidth',2.5);
plot(ax,[x0(1) x0(1)],[-4 4],'-','Color',[cG 0.25],'LineWidth',12);
text(ax,x0(1)+0.25,3.3,'$y_0$ e $y_1$: mesma reta','FontSize',15);
text(ax,x0(1)+0.25,-3.3,'infinitas solu\c{c}\~oes para $x_{0,2}$','FontSize',14,'Color',cG);
xlim(ax,[-4 4]); ylim(ax,[-4 4]); xlabel(ax,'$x_{0,1}$'); ylabel(ax,'$x_{0,2}$');
title(ax,'N\~ao observ\''avel: $\mathrm{rank}\,\mathcal{O}=1$','FontSize',15);
salva(f,'aula8_observabilidade_retas.png'); close(f);

%% 4) Alocação de polos: exemplo A = [0 1; -0.96 -2], polos desejados 0.8 e 0.9
A = [0 1; -0.96 -2]; B = [0; 1];
K = acker(A,B,[0.8 0.9]);
pol = eig(A); pcl = eig(A - B*K);
N = 40; x0 = [1; 0];
xa = zeros(2,N+1); xf_ = xa; xa(:,1) = x0; xf_(:,1) = x0;
for k = 1:N
    xa(:,k+1)  = A*xa(:,k);
    xf_(:,k+1) = (A - B*K)*xf_(:,k);
end
f = figure('Color','w','Position',[100 100 1150 460]);
tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
th = linspace(0,2*pi,300);
fill(ax,cos(th),sin(th),cB,'FaceAlpha',0.07,'EdgeColor',cK,'LineWidth',1.2);
xline(ax,0,'Color',cGr); yline(ax,0,'Color',cGr);
plot(ax,real(pol),imag(pol),'x','Color',cR,'MarkerSize',16,'LineWidth',3);
plot(ax,real(pcl),imag(pcl),'x','Color',cB,'MarkerSize',16,'LineWidth',3);
text(ax,-1.35,0.25,'malha aberta','FontSize',14,'Color',cR);
text(ax,0.55,-0.3,'malha fechada','FontSize',14,'Color',cB);
xlim(ax,[-1.5 1.3]); ylim(ax,[-1.2 1.2]); xlabel(ax,'Re$\{z\}$'); ylabel(ax,'Im$\{z\}$');
title(ax,'$\lambda\{A\}=\{-0{,}8;\,-1{,}2\}\ \to\ \lambda\{A-BK\}=\{0{,}8;\,0{,}9\}$','FontSize',14);
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax,0:8,xa(1,1:9),'-o','Color',cR,'LineWidth',1.8,'MarkerFaceColor',cR,'MarkerSize',5);
text(ax,4.5,2.6,'malha aberta diverge ($|\lambda|=1{,}2>1$)','FontSize',13,'Color',cR);
stairs(ax,0:N,xf_(1,:),'Color',cB,'LineWidth',2.5);
ylim(ax,[-3 3]); xlim(ax,[0 N]);
legend(ax,{'malha aberta ($u_k=0$)','malha fechada ($u_k=-Kx_k$)'},'Location','southeast','FontSize',13);
xlabel(ax,'$k$'); ylabel(ax,'$x_{1,k}$');
title(ax,sprintf('$K=[%.2f\\;\\;%.2f]$, \\ $x_0=[1\\;\\;0]^\\top$',K),'FontSize',14);
salva(f,'aula8_alocacao_polos.png'); close(f);

%% 5) Ackermann: exemplo A = [-1.2 1; 0 -0.8], B = [0.5; 1]
A = [-1.2 1; 0 -0.8]; B = [0.5; 1];
K = acker(A,B,[0.8 0.9]);
N = 50; x = zeros(2,N+1); x(:,1) = [1; -1]; uu = zeros(1,N);
for k = 1:N, uu(k) = -K*x(:,k); x(:,k+1) = A*x(:,k) + B*uu(k); end
f = figure('Color','w','Position',[100 100 1000 430]);
tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
stairs(ax,0:N,x(1,:),'Color',cB,'LineWidth',2.2); stairs(ax,0:N,x(2,:),'Color',cR,'LineWidth',2.2);
legend(ax,{'$x_{1,k}$','$x_{2,k}$'},'FontSize',14); xlabel(ax,'$k$'); title(ax,'Estados em malha fechada','FontSize',15);
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
stairs(ax,0:N-1,uu,'Color',cG,'LineWidth',2.2); xlabel(ax,'$k$'); ylabel(ax,'$u_k$');
title(ax,sprintf('$u_k=-Kx_k,\\ K=[%.3f\\;\\;%.3f]$',K),'FontSize',15);
salva(f,'aula8_ackermann_exemplo.png'); close(f);

%% 6) Observador de estados (Luenberger) + realimentação do estado estimado
A = [0 1; -0.96 -2]; B = [0; 1]; H = [1 0];
K = acker(A,B,[0.8 0.9]);
L = acker(A',H',[0.4 0.5])';          % observador ~2x mais rápido (dualidade)
N = 40; x = zeros(2,N+1); xh = x; x(:,1) = [1; -0.5]; xh(:,1) = [0; 0];
for k = 1:N
    y = H*x(:,k); yh = H*xh(:,k);
    u = -K*xh(:,k);
    x(:,k+1)  = A*x(:,k) + B*u;
    xh(:,k+1) = A*xh(:,k) + B*u + L*(y - yh);
end
e = x - xh;
f = figure('Color','w','Position',[100 100 1150 430]);
tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
for i = 1:2
    ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    stairs(ax,0:N,x(i,:),'Color',cB,'LineWidth',2.4);
    stairs(ax,0:N,xh(i,:),'--','Color',cR,'LineWidth',2.2);
    legend(ax,{sprintf('$x_{%d,k}$',i),sprintf('$\\hat x_{%d,k}$',i)},'FontSize',14);
    xlabel(ax,'$k$'); title(ax,sprintf('Estado $x_%d$: real $\\times$ estimado',i),'FontSize',14);
end
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
stairs(ax,0:N,vecnorm(e),'Color',cG,'LineWidth',2.4);
xlabel(ax,'$k$'); ylabel(ax,'$\|\varepsilon_k\|$');
title(ax,'$\varepsilon_{k+1}=(A-LH)\varepsilon_k$','FontSize',14);
salva(f,'aula8_observador_simulacao.png'); close(f);
fprintf('L = [%.4f %.4f]''\n', L);

%% 7) Rastreamento com ação integral x realimentação pura
A = [0 1; -0.96 -2]; B = [0; 1]; H = [1 0];
Aa = [A zeros(2,1); -H 1]; Ba = [B; 0];
Ka = acker(Aa,Ba,[0.5 0.55 0.6]); K = Ka(1:2); Ki = -Ka(3);
K0 = acker(A,B,[0.5 0.6]);
N = 60; r = ones(1,N+1); w = zeros(1,N+1); w(31:end) = 0.05;   % perturbação na entrada em k = 60
x1 = zeros(2,N+1); x2 = x1; m = 0; y1 = zeros(1,N+1); y2 = y1;
for k = 1:N+1
    y1(k) = H*x1(:,k); y2(k) = H*x2(:,k);
    if k > N, break; end
    u1 = -K0*x1(:,k);                        % só regulação
    u2 = -K*x2(:,k) + Ki*m;                  % com integrador
    x1(:,k+1) = A*x1(:,k) + B*(u1 + w(k));
    x2(:,k+1) = A*x2(:,k) + B*(u2 + w(k));
    m = m + (r(k) - y2(k));
end
f = figure('Color','w','Position',[100 100 1000 430]);
ax = axes(f); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
stairs(ax,0:N,r,'k--','LineWidth',1.5);
stairs(ax,0:N,y1,'Color',cR,'LineWidth',2.2);
stairs(ax,0:N,y2,'Color',cB,'LineWidth',2.6);
xline(ax,30,':','perturba\c{c}\~ao $w=0{,}05$','Interpreter','latex','FontSize',13,'LabelVerticalAlignment','middle','LabelOrientation','horizontal');
legend(ax,{'refer\^encia $r_k$','$u_k=-Kx_k$ (regula\c{c}\~ao)','$u_k=-Kx_k+K_im_k$ (integrador)'},'Location','north','FontSize',13);
xlabel(ax,'$k$'); ylabel(ax,'$y_k$'); ylim(ax,[-0.2 1.5]);
title(ax,sprintf('Rastreamento de degrau: $K=[%.3f\\;\\;%.3f]$, $K_i=%.4f$',K,Ki),'FontSize',15);
salva(f,'aula8_rastreamento_integrador.png'); close(f);
fprintf('K = [%.4f %.4f], Ki = %.4f\n', K, Ki);

function annotation_arrow(ax, xa, xb, i)
    % seta curva entre polo de malha aberta e polo de malha fechada
    t = linspace(0,1,50);
    xm = (xa(1)+xb(1))/2; h = 0.35*(-1)^i;
    xx = (1-t).^2*xa(1) + 2*(1-t).*t*xm + t.^2*xb(1);
    yy = 2*(1-t).*t*h;
    plot(ax,xx,yy,'-','Color',[0.4 0.4 0.4],'LineWidth',1.2);
    plot(ax,xx(end-1:end),yy(end-1:end),'-','Color',[0.4 0.4 0.4]);
    quiver(ax,xx(end-3),yy(end-3),xx(end)-xx(end-3),yy(end)-yy(end-3),0,'Color',[0.4 0.4 0.4],'MaxHeadSize',3,'LineWidth',1.2);
end
