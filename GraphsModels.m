% =========================================================================
% NAME:  GraphsModels
% DESCRIPTION:  Computational framework for analyzing unperturbed tumor
%               growth through mKJMA, MBBG, and Richards models, focusing
%               on fractal dimension, heterogeneous nucleation, dynamic
%               transformations, cooperativity, and phase transitions
%
% See Article:  Close relationship and role of fractal dimension,
%               heterogeneous nucleation, dynamic transformations
%               and cooperativity mechanisms in unperturbed tumor growth
%
% Article authors: Luis Enrique Bergues Cabrales, J. Roberto Romero-Arias,
%                  Guillermo Ramrez-Santiago, Jose Luis Gracia and 
%                  Juan Ignacio Montijano Torcal.
%
% AUTHOR:   J. Roberto Romero-Arias
% ADDRESS:  IIMAS, UNAM. Mexico.
% DATE:     July 2026
% =========================================================================

clear all
tic
aVec = 0.05:0.005:1.2;  %alpha variations
bdfVec = 0.05:0.005:3.0;%beta variations

% Sarcoma
au2   = 0.576;
b    = 0.252;
Vo    = 0.5;
Vobs  = 0.148;
cc    = 1.0;
Df    = 1.256;

dataalphasar = modelTau(Vo,Vobs,b,cc,aVec);
datafractalsar = modelTauFractal(au2,Vo,Vobs,Df,cc,bdfVec);

% Ehrlich
au2   = 0.687;
b    = 0.261;
Vo    = 0.5;
Vobs  = 0.190;
cc    = 1.0;
Df    = 1.023;
dataalphaEhr = modelTau(Vo,Vobs,b,cc,aVec);
datafractalEhr = modelTauFractal(au2,Vo,Vobs,Df,cc,bdfVec);

%% PLOTS
So.alpha = dataalphasar(:,1);   S1.df = datafractalsar(:,1);
So.l1 = dataalphasar(:,2);      S1.l1 = datafractalsar(:,2);
So.l2 = dataalphasar(:,3);      S1.l2 = datafractalsar(:,3);
So.n1 = dataalphasar(:,4);      S1.n1 = datafractalsar(:,4);
So.n2 = dataalphasar(:,5);      S1.n2 = datafractalsar(:,5);
So.K1 = dataalphasar(:,6);      S1.K1 = datafractalsar(:,6);
So.K2 = dataalphasar(:,7);      S1.K2 = datafractalsar(:,7);
So.KR = dataalphasar(:,8);      S1.KR = datafractalsar(:,8);
So.m = dataalphasar(:,9);       S1.m = datafractalsar(:,9);
So.r = dataalphasar(:,10);      S1.r = datafractalsar(:,10);
S1.alpha = datafractalsar(:,11);
S1.beta = datafractalsar(:,12);

E1.alpha = dataalphaEhr(:,1);   E3.df = datafractalEhr(:,1);
E1.l1 = dataalphaEhr(:,2);      E3.l1 = datafractalEhr(:,2);
E1.l2 = dataalphaEhr(:,3);      E3.l2 = datafractalEhr(:,3);
E1.n1 = dataalphaEhr(:,4);      E3.n1 = datafractalEhr(:,4);
E1.n2 = dataalphaEhr(:,5);      E3.n2 = datafractalEhr(:,5);
E1.K1 = dataalphaEhr(:,6);      E3.K1 = datafractalEhr(:,6);
E1.K2 = dataalphaEhr(:,7);      E3.K2 = datafractalEhr(:,7);
E1.KR = dataalphaEhr(:,8);      E3.KR = datafractalEhr(:,8);
E1.m = dataalphaEhr(:,9);       E3.m = datafractalEhr(:,9);
E1.r = dataalphaEhr(:,10);      E3.r = datafractalEhr(:,10);
E3.alpha = datafractalEhr(:,11);
E3.beta = datafractalEhr(:,12);

%% Figure 1a
ppin = 4;
pp = 4;
ppend = length(E1.alpha)-80;

figure;
clf
hold on
plot(E1.alpha(ppin:pp:ppend), E1.l2(ppin:pp:ppend),'ko','LineWidth',1,'MarkerSize',6)
plot(E1.alpha(ppin:pp:ppend), real(E1.n2(ppin:pp:ppend)),'b^','LineWidth',1,'MarkerSize',6)
plot(E1.alpha(ppin:pp:ppend), (1)*real(E1.K2(ppin:pp:ppend)),'ks','LineWidth',1,'MarkerSize',6)
plot(E1.alpha(ppin:pp:ppend), (10)*E1.m(ppin:pp:ppend),'kp','LineWidth',1,'MarkerSize',6)
plot(E1.alpha(ppin:pp:ppend), (0.1)*E1.r(ppin:pp:ppend),'.','Color',[0.6,0,0],'LineWidth',2,'MarkerSize',16)
plot(E1.alpha(ppin:pp:ppend), ((1)*E1.KR(ppin:pp:ppend)),'kd','LineWidth',1,'MarkerSize',6)

axis([0.3 1.01 -1 5.1])
box on
legend('\lambda','n','10K','10m','0.1r','K_{R}','Box','off','FontSize',18,'FontWeight','bold','Location','Northwest','NumColumns',2)
xlabel('\alpha (days^{-1})','FontWeight','bold','FontSize',18)
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',22,'FontWeight','bold')
xtickformat('%.1f');

%%  Figure 1b
ppin = 4;
pp = 4;
ppend = length(So.alpha)-94;

figure;
clf
hold on
plot(So.alpha(ppin:pp:ppend), So.l2(ppin:pp:ppend),'ko','LineWidth',1,'MarkerSize',6)
plot(So.alpha(ppin:pp:ppend), So.n2(ppin:pp:ppend),'b^','LineWidth',1,'MarkerSize',6)
plot(So.alpha(ppin:pp:ppend), So.K2(ppin:pp:ppend),'ks','LineWidth',1,'MarkerSize',6)
plot(So.alpha(ppin:pp:ppend), (10)*So.m(ppin:pp:ppend),'kp','LineWidth',1,'MarkerSize',6)
plot(So.alpha(ppin:pp:ppend), (0.1)*So.r(ppin:pp:ppend),'.','Color',[0.6,0,0],'LineWidth',1,'MarkerSize',16)
plot(So.alpha(ppin:pp:ppend), (1)*So.KR(ppin:pp:ppend),'kd','LineWidth',1,'MarkerSize',6)

axis([0.3 1.01 -1 5.1])
box on
legend('\lambda','n','10K','10m','0.1r','K_{R}','Box','off','FontSize',18,'FontWeight','bold','Location','Northwest','NumColumns',2)
xlabel('\alpha (days^{-1})','FontWeight','bold','FontSize',18)
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',22,'FontWeight','bold')
xtickformat('%.1f');

%% Figure Erhlich fractal
ppin = 6;
pp = 6;
ppend = length(E3.df);

figure; %  Volume V_0 = 0.5 cm^3 Ehrlich Fractal
clf
hold on
plot(E3.df(ppin:pp:ppend), E3.l2(ppin:pp:ppend),'ko','LineWidth',1,'MarkerSize',6)
plot(E3.df(ppin:pp:ppend), real(E3.n2(ppin:pp:ppend)),'b^','LineWidth',1,'MarkerSize',10)
plot(E3.df(ppin:pp:ppend), real(0.1*E3.K2(ppin:pp:ppend)),'ks','LineWidth',1,'MarkerSize',10)
plot(E3.df(ppin:pp:ppend), E3.m(ppin:pp:ppend),'k*','LineWidth',1,'MarkerSize',8)
plot(E3.df(ppin:pp:ppend), E3.r(ppin:pp:ppend),'.','Color',[0.6,0,0],'LineWidth',1,'MarkerSize',30)
plot(E3.df(ppin:pp:ppend), E3.KR(ppin:pp:ppend),'kd','LineWidth',1,'MarkerSize',8)
plot(E3.df(ppin:pp:ppend), E3.alpha(ppin:pp:ppend),'rp','LineWidth',1,'MarkerSize',10)
plot(E3.df(ppin:pp:ppend), E3.beta(ppin:pp:ppend),'m+','LineWidth',1,'MarkerSize',8)

axis([0.0 3.01 -4.22 10.22])
box on
legend('\lambda','n','K','m','r','K_{R}','\alpha','\beta','Box','off','FontSize',18,'FontWeight','bold','Location','Northeast','NumColumns',2)
xlabel('d_f','FontWeight','bold','FontSize',18)
txt = sprintf('Ehrlich');
text(0.1, 9.7, txt,'FontName','Helvetica','FontSize',18,'FontWeight','bold');
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',18,'FontWeight','bold')

%% Figure Fibrosarcoma fractal
ppin = 6;
pp = 6;
ppend = length(S1.df);

figure; %  Volume V_0 = 0.5 cm^3 Sarcoma Fractal
clf
hold on
plot(S1.df(ppin:pp:ppend), S1.l2(ppin:pp:ppend),'ko','LineWidth',1,'MarkerSize',6)
plot(S1.df(ppin:pp:ppend), real(S1.n2(ppin:pp:ppend)),'b^','LineWidth',1,'MarkerSize',10)
plot(S1.df(ppin:pp:ppend), real(0.1*S1.K2(ppin:pp:ppend)),'ks','LineWidth',1,'MarkerSize',10)
plot(S1.df(ppin:pp:ppend), S1.m(ppin:pp:ppend),'k*','LineWidth',1,'MarkerSize',8)
plot(S1.df(ppin:pp:ppend), S1.r(ppin:pp:ppend),'.','Color',[0.6,0,0],'LineWidth',1,'MarkerSize',30)
plot(S1.df(ppin:pp:ppend), S1.KR(ppin:pp:ppend),'kd','LineWidth',1,'MarkerSize',8)
plot(S1.df(ppin:pp:ppend), S1.alpha(ppin:pp:ppend),'rp','LineWidth',1,'MarkerSize',10)
plot(S1.df(ppin:pp:ppend), S1.beta(ppin:pp:ppend),'m+','LineWidth',1,'MarkerSize',6)

axis([0.0 3.01 -4.18 10.22])
box on
legend('\lambda','n','K','m','r','K_{R}','\alpha','\beta','Box','off','FontSize',18,'FontWeight','bold','Location','Northeast','NumColumns',2)
xlabel('d_f','FontWeight','bold','FontSize',18)
txt = sprintf('Fibrosarcoma Sa-37');
text(0.7, 9.7, txt,'FontName','Helvetica','FontSize',18,'FontWeight','bold');
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',18,'FontWeight','bold')

%% Figure 1c
ppin = 6;
ppin1 = 9;
pp = 12;
pp1 = 16;
ppend = length(S1.df);

figure; %  Volume V_0 = 0.5 cm^3 Sarcoma Fractal
clf
hold on

plot(E3.df(ppin:pp:ppend), real(E3.n2(ppin:pp:ppend)),'k^','LineWidth',1,'MarkerSize',6)
plot(E3.df(ppin:pp:ppend), real(1*E3.K2(ppin:pp:ppend)),'ks','LineWidth',1,'MarkerSize',6)
plot(E3.df(ppin:pp:ppend), 10*E3.m(ppin:pp:ppend),'kp','LineWidth',1,'MarkerSize',6)
plot(E3.df(ppin:pp:ppend), 0.1*E3.r(ppin:pp:ppend),'ko','LineWidth',1,'MarkerSize',6)

plot(S1.df(ppin:pp:ppend), real(S1.n2(ppin:pp:ppend)),'b^','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','b')
plot(S1.df(ppin:pp:ppend), real(1*S1.K2(ppin:pp:ppend)),'bs','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','b')
plot(S1.df(ppin:pp:ppend), 10*S1.m(ppin:pp:ppend),'bp','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','b')
plot(S1.df(ppin:pp:ppend), 0.1*S1.r(ppin:pp:ppend),'bo','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','b')

axis([0.0 3.11 -0.18 5.22])
box on
legend('n_E','10K_E','10m_E','0.1r_E','n_S','10K_S','10m_S','0.1r_S','Box','off','FontSize',14,'FontWeight','bold','Location','Northeast','NumColumns',2)
xlabel('d_f','FontWeight','bold','FontSize',22)
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',22,'FontWeight','bold')
%ytickformat('%.2f');
xtickformat('%.1f');

%% Figure 1d
ppin = 74;
pp = 3;
ppend = length(S1.df)-420;

figure; %
clf
hold on

plot(E3.df(ppin:pp:ppend), real(E3.n2(ppin:pp:ppend)),'k^','LineWidth',1,'MarkerSize',6)
plot(E3.df(ppin:pp:ppend), real(1*E3.K2(ppin:pp:ppend)),'ks','LineWidth',1,'MarkerSize',6)
plot(E3.df(ppin:pp:ppend), 10*E3.m(ppin:pp:ppend),'kp','LineWidth',1,'MarkerSize',6)
plot(E3.df(ppin:pp:ppend), 0.1*E3.r(ppin:pp:ppend),'ko','LineWidth',1,'MarkerSize',6)

plot(S1.df(ppin:pp:ppend), real(S1.n2(ppin:pp:ppend)),'b^','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','b')
plot(S1.df(ppin:pp:ppend), real(1*S1.K2(ppin:pp:ppend)),'bs','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','b')
plot(S1.df(ppin:pp:ppend), 10*S1.m(ppin:pp:ppend),'bp','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','b')
plot(S1.df(ppin:pp:ppend), 0.1*S1.r(ppin:pp:ppend),'bo','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','b')

xline(0.611,'--k','LineWidth',1)
xline(0.688,'--b','LineWidth',1)

axis([0.4 0.91 -0.58 5.22])
box on
legend('n_E','10K_E','10m_E','0.1r_E','n_S','10K_S','10m_S','0.1r_S','Box','off','FontSize',14,'FontWeight','bold','Location','Northeast','NumColumns',2)

xlabel('d_f','FontWeight','bold','FontSize',22)
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',22,'FontWeight','bold')
%ytickformat('%.2f');
xtickformat('%.1f');

%%
% Maps
dfVec = 0.4:0.02:0.8;  % df variations
DfVec = 1.0:0.02:1.45; % Df variations

% Generic tumors
au2   = 0.620;
b    = 0.252;
Vo    = 0.5;
Vobs  = 0.15;
cc    = 1.0;
Df    = 1.1;

alphasurf = NaN(length(dfVec),length(DfVec));
betasurf = alphasurf;
nsurf = alphasurf;
lsurf = alphasurf;
Ksurf = alphasurf;
msurf = alphasurf;
KRsurf = alphasurf;
easurf = alphasurf;

for j = 1:length(DfVec)

    Df = DfVec(j);
    data = modelTauFractal(au2,Vo,Vobs,Df,cc,bdfVec);
    df_model = S1.df;
    alpha_model = (data(:,11));
    beta_model = (data(:,12));
    n_model = (data(:,5));
    l_model = (data(:,3));
    K_model = (data(:,7));
    m_model = (data(:,9));
    KR_model =(data(:,8));

    for i=1:length(df_model)
        idx = find(abs(dfVec-df_model(i))<1e-4);
        if ~isempty(idx)
            alphasurf(idx,j) = alpha_model(i);
            betasurf(idx,j) = beta_model(i);
            nsurf(idx,j) = n_model(i);
            lsurf(idx,j) = l_model(i);
            Ksurf(idx,j) = K_model(i);
            msurf(idx,j) = m_model(i);
            KRsurf(idx,j) = KR_model(i);
            easurf(idx,j) = -2568.85*log(0.1*Ksurf(idx,j)/0.109);
        end
    end

end

[Dfgrid,dfgrid] = meshgrid(DfVec,dfVec);

%% Figure 2a
Kcrit = 1;
phase = zeros(size(alphasurf));
phase(alphasurf > Kcrit) = 1;

figure;
clf
hold on
contourf(Dfgrid,dfgrid,alphasurf,20)
contour(Dfgrid,dfgrid,alphasurf,[Kcrit Kcrit],'k','LineWidth',4)
plot(1.256, 0.688,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(1.023, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')

ccp = colorbar;
ccp.YTick = [0.72 0.8 0.9 1.0 1.1];
ccp.YTickLabel = [0.7 0.8 0.9 1.0 1.1];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([1.0 1.4 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%% Figure 2b

figure;
clf
hold on
contourf(Dfgrid,dfgrid,betasurf,20)
plot(1.256, 0.688,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(1.023, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')

ccp = colorbar;
ccp.YTick = [0.13 0.22 0.32 0.41];
ccp.YTickLabel = [0.12 0.22 0.32 0.42];
box on
axis([1.0 1.4 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%% Figure 2c
Kcrit = 1;
phase = zeros(size(nsurf));
phase(nsurf > Kcrit) = 1;

figure;
clf
hold on
contourf(Dfgrid,dfgrid,nsurf,20)
contour(Dfgrid,dfgrid,nsurf,[Kcrit Kcrit],'k','LineWidth',4)
plot(1.256, 0.688,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(1.023, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')

ccp = colorbar;
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([1.0 1.4 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%% Figure 2d
figure;
clf
hold on
contourf(Dfgrid,dfgrid,0.1*real(Ksurf),20)
plot(1.256, 0.688,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(1.023, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')

ccp = colorbar;
ccp.YTick = 0.1*[0.58 0.95 1.35 1.75 2.15 2.54];
ccp.YTickLabel = 0.1*[0.6 1.0 1.4 1.8 2.2 2.6];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([1.0 1.4 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%% Figure 2e
figure;
clf
hold on
contourf(Dfgrid,dfgrid,msurf,20)
plot(1.256, 0.688,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(1.023, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')

ccp = colorbar
ccp.YTick = [0.002 0.01 0.02 0.03 0.04 0.05 0.06];
ccp.YTickLabel = [0 0.01 0.02 0.03 0.04 0.05 0.06];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([1.0 1.4 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%% Figure 2f
KRsurf1 = KRsurf;
KRsurf1(find(KRsurf>10.3))=10;

Kcrit = 9;
phase = zeros(size(KRsurf1));
phase(KRsurf1 > Kcrit) = 9;

figure;
clf
hold on
contourf(Dfgrid,dfgrid,KRsurf1,20)
contour(Dfgrid,dfgrid,KRsurf,[Kcrit Kcrit],'k','LineWidth',3)
plot(1.256, 0.688,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(1.023, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')

ccp = colorbar;
ccp.YTick = [1.33 3 4.8 6.5 8.2 9.6];
ccp.YTickLabel = [1 3 5 7 8 40]; %K_R
box on
axis([1.0 1.4 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%%
% Maps for Vobs
clear KRsurf
dfVec = 0.4:0.02:0.8;
DVobs = 0.1:0.01:0.3;

% Generic tumor
au2   = 0.620;
b    = 0.252;
Vo    = 0.5;
Vobs  = 0.15;
cc    = 1.0;
Df    = 1.1;

alphasurf = NaN(length(dfVec),length(DVobs));
betasurf = alphasurf;
Ksurf = alphasurf;
msurf = alphasurf;

for j = 1:length(DVobs)

    Vobs = DVobs(j);
    data = modelTauFractal(au2,Vo,Vobs,Df,cc,bdfVec);
    df_model = S1.df;
    alpha_model = real(data(:,11));
    beta_model = real(data(:,12));
    K_model = real(data(:,7));
    m_model = real(data(:,9));

    for i=1:length(df_model)
        idx = find(abs(dfVec-df_model(i))<1e-4);
        if ~isempty(idx)
            alphasurf(idx,j) = alpha_model(i);
            betasurf(idx,j) = beta_model(i);
            msurf(idx,j) = m_model(i);
            Ksurf(idx,j) = K_model(i);
        end
    end

end

[Dfgrid,dfgrid] = meshgrid(DVobs,dfVec);

%% Figure 3a
figure;
clf
hold on
contourf(Dfgrid,dfgrid,alphasurf,20)

ccp = colorbar;
ccp.YTick = [0.6 0.7 0.8 0.9 1.0];
ccp.YTickLabel = [0.6 0.7 0.8 0.9 1.0];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0.1 0.3 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.2f');

%% Figure 3b
figure;
clf
hold on
contourf(Dfgrid,dfgrid,betasurf,20)

ccp = colorbar;
ccp.YTick = [0.17 0.21 0.25 0.29 0.33 0.37];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0.1 0.3 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.2f');

%% Figure 3c
figure;
clf
hold on
contourf(Dfgrid,dfgrid,msurf,20)

ccp = colorbar;
ccp.YTick = [-0.1 -0.05 0.00 0.05 0.10];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0.1 0.3 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.2f');

%% Figure 3d
figure;
clf
hold on
contourf(Dfgrid,dfgrid,0.1*Ksurf,20)

ccp = colorbar;
ccp.YTick = 0.1*[0.8 1.1 1.4 1.7 2.0 2.3];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0.1 0.3 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.2f');

%%
clear RR AV BG

dfVec = E3.df(70:150);
Time = 0:0.1:15.0;
[dfgrid,Timegrid] = meshgrid(dfVec',Time);
Vobs = 0.15;

for  t= 1:length(Time)

    for j = 1:length(dfVec)

        RR(t,j) = (Richards(E3.KR(j),E3.m(j),E3.r(j),Vobs,Time(t)) - Vobs)./(E3.KR(j) - Vobs);
        AV(t,j) = real(Avrami(0.1*E3.K2(j),E3.n2(j),(E3.l2(j)),Time(t)));
        BG(t,j) = (Gompertz(E3.alpha(j),E3.beta(j),Vobs,Vo,Time(t)-1.8)...
            - Vobs )./(Vobs*exp(E3.alpha(j)./E3.beta(j))-Vobs);
    end
end

AV = smoothdata2(AV);

%% Figure 4a
figure;
clf
surf(dfgrid,Timegrid,(RR),'FaceColor','blue','FaceAlpha',.5)
hold on
surf(dfgrid,Timegrid,(BG),'FaceColor','red','FaceAlpha',.7)
hold off
view(-128,14)
grid off

xlabel('d_f')
ylabel('Days','Rotation',20)
zlabel('Normalized models')
axis([0.4 0.8 0 15 0 1])
box on
legend('Richards','MBBG','Box','off','FontSize',18,'FontWeight','bold','Location','best','NumColumns',1)
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',22,'FontWeight','bold')
ztickformat('%.1f');
ytickformat('%.0f');
xtickformat('%.1f');


%% Figure 4b
figure;
clf
hold on
surf(dfgrid,Timegrid,BG,'FaceColor','red','FaceAlpha',.7)
mesh(dfgrid,Timegrid,AV,'FaceColor','cyan','FaceAlpha',.1)
hold off
view(174,14)
grid off

xlabel('d_f')
ylabel('Days','Rotation',-60)
zlabel('Normalized models')
axis([0.4 0.8 0 15 0 1])
box on
legend('MBBG','mKJMA','Box','off','FontSize',18,'FontWeight','bold','Location','best','NumColumns',1)
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',22,'FontWeight','bold')
ztickformat('%.1f');
ytickformat('%.0f');
xtickformat('%.1f');

%%
% Maps for Vobs
clear KRsurf
dfVec = 0.4:0.02:0.8;
DVobs = 0.1:0.01:0.3;

% Generic tumor
au2   = 0.620;
b    = 0.252;
Vo    = 0.5;
Vobs  = 0.15;
cc    = 1.0;
Df    = 1.1; % change for complete figue 4f

Gomsurf = NaN(length(dfVec),length(DVobs));
Richsurf = Gomsurf;
AVsurf = Gomsurf;

t = 10;  %Time in days

for j = 1:length(DVobs)

    Vobs = DVobs(j);
    data = modelTauFractal(au2,Vo,Vobs,Df,cc,bdfVec);
    df_model = S1.df;
    Soo.alpha = real(data(:,11));
    Soo.beta = real(data(:,12));
    Soo.K2 = real(data(:,7));
    Soo.n2 = real(data(:,5));
    Soo.l2 = real(data(:,3));
    Soo.KR = real(data(:,8));
    Soo.m = real(data(:,9));
    Soo.r = real(data(:,10));

    for i=1:length(df_model)
        idx = find(abs(dfVec-df_model(i))<1e-4);
        if ~isempty(idx)
            Gomsurf(idx,j) = (Gompertz(Soo.alpha(i),Soo.beta(i),Vobs,Vo,t-1.7)...
                - Vobs )./(Vobs*exp(Soo.alpha(i)./Soo.beta(i))-Vobs);
            AVsurf(idx,j) = Avrami(0.1*Soo.K2(i),Soo.n2(i),(Soo.l2(i)),t);
            Richsurf(idx,j) = (Richards(Soo.KR(i),Soo.m(i),Soo.r(i),Vobs,t) - Vobs)./(Soo.KR(i) - Vobs);
        end
    end

end

[Dfgrid,dfgrid] = meshgrid(DVobs,dfVec);

%% Figure 4c
figure;
clf
hold on
contourf(Dfgrid,dfgrid,real(Gomsurf),20)
plot(0.15, 0.688,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(0.19, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')

ccp = colorbar;
clim([0, 1])
ccp.YTick = [0. 0.2 0.4 0.6 0.8 1.0];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0.1 0.3 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.2f');

%% Figure 4d
figure;
clf
hold on
contourf(Dfgrid,dfgrid,real(Richsurf),20)
plot(0.15, 0.688,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(0.19, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')

ccp = colorbar
clim([0, 1])
ccp.YTick = [0. 0.2 0.4 0.6 0.8 1.0];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0.1 0.3 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.2f');

%% Figure 4e
figure;
clf
hold on
contourf(Dfgrid,dfgrid,real(AVsurf),20)
plot(0.15, 0.688,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(0.19, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')

colormap(parula);
ccp = colorbar;
clim([0, 1])
ccp.YTick = [0. 0.2 0.4 0.6 0.8 1.0];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0.1 0.3 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.2f');

%% Figure 4f
FAVsurf =real(AVsurf);
FAVsurf(find(FAVsurf>=0.76))=1;
FAVsurf(find(FAVsurf<0.76))=0;

figure;
clf
hold on
contourf(Dfgrid,dfgrid,FAVsurf,30)

colormap(gray);
box on
axis([0.1 0.3 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.2f');

%%
% Maps for Vobs
%clear nrlocal, lrlocal, Krlocal
deltat = 0.1;
dfVec = 0.4:0.01:0.8;
TIMET = deltat:deltat:21;
% Chose the correct tumor
%
% Generic tumor
au2   = 0.620;
b    = 0.252;
Vo    = 0.5;
Vobs  = 0.15;
cc    = 1.0;
Df    = 1.1;

% % F3II
% au2   = 0.252;
% b    = 0.094;
% Vo    = 0.5;
% Vobs  = 0.16;
% cc    = 1.0;
% Df    = 1.1;

% % Sarcoma
% au2   = 0.576;
% b    = 0.252;
% Vo    = 0.5;
% Vobs  = 0.15;
% cc    = 1.0;
% Df    = 1.256;

% % Ehrlich
% au2   = 0.687;
% b    = 0.261;
% Vo    = 0.5;
% Vobs  = 0.190;
% cc    = 1.0;
% Df    = 1.023;

nrlocal = NaN(length(dfVec),length(TIMET));
lrlocal = NaN(length(dfVec),length(TIMET));
Krlocal = NaN(length(dfVec),length(TIMET));
AVlocal = NaN(length(dfVec),length(TIMET));
AVTlocal = NaN(length(dfVec),length(TIMET));
AVMlocal = NaN(length(dfVec),length(TIMET));
AVNlocal = NaN(length(dfVec),length(TIMET));
AVNNlocal = NaN(length(dfVec),length(TIMET));
Gomlocal = NaN(length(dfVec),length(TIMET));
Richlocal = NaN(length(dfVec),length(TIMET));
nnrlocal = NaN(length(dfVec),length(TIMET));
nKrlocal = NaN(length(dfVec),length(TIMET));
earlocal = NaN(length(dfVec),length(TIMET));
AVTTlocal = NaN(length(dfVec),length(TIMET));

data = modelTauFractal(au2,Vo,Vobs,Df,cc,bdfVec);

df_model = S1.df;
Soo.alpha = real(data(:,11));
Soo.beta = real(data(:,12));
Soo.K2 = real(data(:,7));
Soo.n2 = real(data(:,5));
Soo.l2 = real(data(:,3));
Soo.KR = real(data(:,8));
Soo.m = real(data(:,9));
Soo.r = real(data(:,10));

% mean values
SSoo.n2 = Soo.n2;
SSoo.n2 = SSoo.n2(~isnan(SSoo.n2));
SSoo.n2(find(SSoo.n2 ==0))=[];
Soon2 = mean(SSoo.n2);

SSoo.l2 = Soo.l2;
SSoo.l2 = SSoo.l2(~isnan(SSoo.l2));
SSoo.l2(find(SSoo.l2 ==0))=[];
Sool2 = mean(SSoo.l2);

SSoo.K2 = Soo.K2;
SSoo.K2 = SSoo.K2(~isnan(SSoo.K2));
SSoo.K2(find(SSoo.K2 == 0))=[];
SSoo.K2(find(SSoo.K2 < 0))=[];
SooK2 = mean(SSoo.K2);


for j = 1:length(TIMET)

    t = TIMET(j);

    for i=1:length(df_model)
        idx = find(abs(dfVec-df_model(i))<1e-4);
        if ~isempty(idx)

            nrlocal(idx,j) = (nlocal(Soo.K2(i),Soo.n2(i),Soo.l2(i),t+deltat)- nlocal(Soo.K2(i),Soo.n2(i),Soo.l2(i),t-deltat))/(log(t+deltat)-log(t-deltat));

            lrlocal(idx,j) =  (1-qTau(Soo.alpha(i),Soo.beta(i),Vo,Vobs)).*(nrlocal(idx,j) -1)./ (Tau(Soo.alpha(i),Soo.beta(i),Vo,Vobs)*qpTau(Soo.alpha(i),Soo.beta(i),Vo,Vobs));

            Krlocal(idx,j) = (( (nrlocal(idx,j) -1)./(nrlocal(idx,j) + lrlocal(idx,j) -1) ).^(1./nrlocal(idx,j)) )./Tau(Soo.alpha(i),Soo.beta(i),Vo,Vobs);

            AVlocal(idx,j) =  Avrami(Krlocal(idx,j),nrlocal(idx,j),lrlocal(idx,j),t);

            Gomlocal(idx,j) = (Gompertz(Soo.alpha(i),Soo.beta(i),Vobs,Vo,t)...
                - Vo )./(Vobs*exp(Soo.alpha(i)./Soo.beta(i))-Vo);

            Richlocal(idx,j) = (Richards(Soo.KR(i),Soo.m(i),Soo.r(i),Vo,t) - Vo)./(Soo.KR(i) - Vo);

            AVTlocal(idx,j) = Avrami(0.1*Soo.K2(i),Soo.n2(i),Soo.l2(i)+0.000001,t);

            AVMlocal(idx,j) = Avrami(0.1*SooK2,Soon2,Sool2,t); %Avrami average

            auxn2 = 1./( 1 + log( 1 - qTau(Soo.alpha(i),Soo.beta(i),Vo,Vobs) )  );

            auxK2 = (1./Tau(Soo.alpha(i),Soo.beta(i),Vo,Vobs))*((auxn2-1)/auxn2 ).^((auxn2-1)/auxn2);

            nnrlocal(idx,j) = (nnlocal(auxK2,auxn2,t+deltat)- nnlocal(auxK2,auxn2,t-deltat))/(log(t+deltat)-log(t-deltat));

            nKrlocal(idx,j) = (1./Tau(Soo.alpha(i),Soo.beta(i),Vo,Vobs))*((nnrlocal(idx,j)-1)/nnrlocal(idx,j) ).^((nnrlocal(idx,j)-1)/nnrlocal(idx,j));

            AVNlocal(idx,j) = 1 - exp(-(nKrlocal(idx,j)*t).^nnrlocal(idx,j));

            AVNNlocal(idx,j) = Avrami(0.1*Soo.K2(i),Soo.n2(i),1.01,t);

            ealocal(idx,j) = -2568.85*log(Krlocal(idx,j)/0.109);

            nealocal(idx,j) = -2568.85*log(Krlocal(idx,j)/0.109);

            Mealocal(idx,j) = -2568.85*log(Krlocal(idx,j)/0.13); 

            AVTTlocal(idx,j) = Soo.l2(idx);

        end
    end

end

[Dfgrid,dfgrid] = meshgrid(TIMET,dfVec);

%% Figure 6a
tstep = ones(length(TIMET),1);

nrlocal1 = nrlocal;
nrlocal1(find(nrlocal>3.3))=3;
nrlocal1(find(nrlocal<0))=0;

figure;
clf
hold on
contourf(log(Dfgrid),dfgrid,real(nrlocal1),13)

plot(log(TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.688 + 0.003), 'sr','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r') %sarcoma
plot(log(TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.611 + 0.005), 'or','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r')
xline( log((1/0.252)*log((0.833/0.252 - log(Vo/Vobs)))),'-r','LineWidth',3) %sarcoma
xline( log((1/0.261)*log((0.719/0.261 - log(Vo/Vobs)))),':r','LineWidth',3)
plot(log(TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.655 + 0.003), 'r^','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r') % table 2

ccp =colorbar ;
ccp.YTick = [0. 0.5 1.0 1.5 2.0 2.5 3.0];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0 3. 0.40 0.8])
%axis([2 2.3 0.6 0.7])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%% Figure 6b
tstep = ones(length(TIMET),1);

nrlocal1 = nrlocal;
nrlocal1(find(nrlocal>3.3))=3;
nrlocal1(find(nrlocal<0))=0;

figure;
clf
hold on
contourf(log(Dfgrid),dfgrid,real(nrlocal1),13)

plot(log(TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.688 + 0.003), 'sr','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r') %sarcoma
plot(log(TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.611 + 0.005), 'or','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r')
xline( log((1/0.252)*log((0.833/0.252 - log(Vo/Vobs)))),'-r','LineWidth',3) %sarcoma
xline( log((1/0.261)*log((0.719/0.261 - log(Vo/Vobs)))),':r','LineWidth',3)
plot(log(TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.655 + 0.003), 'r^','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r') % table 2

ccp =colorbar ;
ccp.YTick = [0. 0.5 1.0 1.5 2.0 2.5 3.0];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([2 2.3 0.60 0.7])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%% Figure 6c

lrlocal1 = lrlocal;
lrlocal1(find(lrlocal > 6.5))=6;
lrlocal1(find(lrlocal < 0))=0;

figure;
clf
hold on
contourf(log(Dfgrid),dfgrid,real(lrlocal1),13)

plot(log(TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.688 + 0.003), 'sr','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r') %sarcoma
plot(log(TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.611 + 0.005), 'or','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r')
xline( log((1/0.252)*log((0.833/0.252 - log(Vo/Vobs)))),'-r','LineWidth',3) %sarcoma
xline( log((1/0.261)*log((0.719/0.261 - log(Vo/Vobs)))),':r','LineWidth',3)
plot(log(TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.655 + 0.003), 'r^','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r') % table 2

colorbar;
box on
axis([0 3 0.40 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');


%% Figure 6d

Krlocal1 = Krlocal;
Krlocal1(find(Krlocal>0.65))=0.65;
Krlocal1(find(Krlocal<0))=0;

figure;
clf
hold on
contourf(log(Dfgrid),dfgrid,real(Krlocal1),13)

plot(log(TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.688 + 0.003), 'sr','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r') %sarcoma
plot(log(TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.611 + 0.005), 'or','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r')
xline( log((1/0.252)*log((0.833/0.252 - log(Vo/Vobs)))),'-r','LineWidth',3) %sarcoma
xline( log((1/0.261)*log((0.719/0.261 - log(Vo/Vobs)))),':r','LineWidth',3)
plot(log(TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.655 + 0.003), 'r^','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r') % table 2

ccp =colorbar ;
box on
axis([0 3 0.40 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%% Figura 7a

AVlocal1 = AVTlocal;
AVlocal1(find(AVTlocal>1.05))=1.05;

figure;
clf
hold on
contour((Dfgrid),dfgrid,real(AVlocal1),300)

plot((TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.688 + 0.003), 'sr','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r') %sarcoma
plot((TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.611 + 0.005), 'or','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r')
xline( ((1/0.252)*log((0.833/0.252 - log(Vo/Vobs)))),'-r','LineWidth',2) %sarcoma
xline( ((1/0.261)*log((0.719/0.261 - log(Vo/Vobs)))),':r','LineWidth',2)
plot((TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.655 + 0.003), 'r^','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r') %sarcoma

ccp =colorbar ;
ccp.YTick = [0.01 0.5 1];
ccp.YTickLabel = [0.0 0.5 1];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0 20.3 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.0f');

%% Figure 7b

AVlocal1 = AVlocal;
AVlocal1(find(AVlocal>1.55))=1.5;
figure;
clf
hold on
contour((Dfgrid),dfgrid,real(AVlocal1),200)

plot((TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.688 + 0.003), 'sr','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r') %sarcoma
plot((TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.611 + 0.005), 'or','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r')
xline( ((1/0.252)*log((0.833/0.252 - log(Vo/Vobs)))),'-r','LineWidth',2) %sarcoma
xline( ((1/0.261)*log((0.719/0.261 - log(Vo/Vobs)))),':r','LineWidth',2)
plot((TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.655 + 0.003), 'r^','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r') %sarcoma

ccp =colorbar ;
ccp.YTick = [-0.5 0.0 0.5 1 1.5];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0 20.3 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.0f');

%% Figure 7c
%29, 22  Ehrlich
ppin = 1;
pp = 2;
ppend = length(TIMET)-14;
ppinn = 90;

figure;
clf
hold on
plot(TIMET(ppin:pp:ppend), Gomlocal(22, ppin:pp:ppend), '-k','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r')
plot(TIMET(ppin:pp:ppend), Richlocal(22, ppin:pp:ppend), ':k','LineWidth',2,'MarkerSize',9,'MarkerFaceColor','r')
plot(TIMET(ppin:pp:ppend), AVTlocal(22, ppin:pp:ppend), '-or','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','r')
plot(TIMET(ppin:pp:ppend-110), real(AVlocal(22,ppin:pp:ppend-110)), '-ob','LineWidth',1,'MarkerSize',8)
plot(TIMET(ppin:pp:ppend-10), AVNlocal(22, ppin:pp:ppend-10), '-.g','LineWidth',1,'MarkerSize',8)
plot(TIMET(ppinn:pp:ppend), AVNNlocal(22, ppinn:pp:ppend), '-ob','LineWidth',1,'MarkerSize',8)

box on
axis([0 20 -0.01 1.05])
%axis([3 5 0.35 0.6 ])
legend('MBBG','Richards','mKJMA','mKJMA_{loc}','KJMA','Box','off','FontSize',18,'FontWeight','bold','Location','best','NumColumns',1)
%txt = sprintf('Ehrlich');
%text(3, 0.97, txt,'FontName','Helvetica','FontSize',18,'FontWeight','bold');
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.1f');
%xtickformat('%.1f');

%% Figure 7d

%29, 22  Sarcoma
ppin = 1;
pp = 2;
ppend = length(TIMET)-14;
ppinn = 90;

figure;
clf
hold on
plot(TIMET(ppin:pp:ppend)+0.0, Gomlocal(29, ppin:pp:ppend), '-k','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r')
plot(TIMET(ppin:pp:ppend), Richlocal(29, ppin:pp:ppend), ':k','LineWidth',2,'MarkerSize',9,'MarkerFaceColor','r')
plot(TIMET(ppin:pp:ppend), AVTlocal(29, ppin:pp:ppend), '-sr','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','r') % sarcoma
plot(TIMET(ppin:pp+1:ppend-110), real(AVlocal(29, ppin:pp+1:ppend-110)), '-sb','LineWidth',1,'MarkerSize',8) % sarcoma
plot(TIMET(ppin:pp:ppend), AVNlocal(29, ppin:pp:ppend), '-.g','LineWidth',1,'MarkerSize',8)
plot(TIMET(ppinn:pp+3:ppend), AVNNlocal(29, ppinn:pp+3:ppend), '-sb','LineWidth',1,'MarkerSize',8)
plot( ((1/0.252)*log((0.833/0.252 - log(Vo/Vobs)))),0.32, 'ok','LineWidth',2,'MarkerSize',14,'MarkerFaceColor','k')

box on
axis([0 20 -0.01 1.05])
%axis([3 5 0.35 0.6 ])
legend('MBBG','Richards','mKJMA','mKJMA_{loc}','KJMA','Box','off','FontSize',16,'FontWeight','bold','Location','best','NumColumns',1)
%txt = sprintf('Fibrosarcoma Sa-37');
%text(2.53, 0.97, txt,'FontName','Helvetica','FontSize',18,'FontWeight','bold');
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');


%% Figure 7e
%29, 22  F3II
AVlocal1 = AVlocal;
AVlocal1(find(AVlocal>1.55))=1.5;

figure;
clf
hold on
contour((Dfgrid),dfgrid,real(AVlocal1),200)

xline( ((1/0.252)*log((0.833/0.252 - log(Vo/Vobs)))),'-r','LineWidth',2) %sarcoma
xline( ((1/0.261)*log((0.719/0.261 - log(Vo/Vobs)))),':r','LineWidth',2)
plot((TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.650 + 0.003), 'r<','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','r') % table 2
plot((TIMET(2:9:end-10)), tstep(2:9:end-10)*(0.67 + 0.003), 'rp','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','r') % table 2

ccp =colorbar ;
ccp.YTick = [-0.5 0.0 0.5 1 1.5];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0 20.3 0.4 0.8])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.0f');

%% Figure 7f
%29, 22  F3II
ppin = 1;
pp = 2;
ppend = length(TIMET)-14;
ppinn = 140;

figure;
clf
hold on

plot(TIMET(ppin:pp:ppend), Gomlocal(22, ppin:pp:ppend), '-k','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r')
plot(TIMET(ppin:pp:ppend), Richlocal(22, ppin:pp:ppend), ':k','LineWidth',2,'MarkerSize',9,'MarkerFaceColor','r')
plot(TIMET(ppin:pp:ppend), AVTlocal(22, ppin:pp:ppend), '-pr','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','r')
plot(TIMET(ppin:pp:ppend-60), real(AVlocal(22,ppin:pp:ppend-60)), '-pb','LineWidth',1,'MarkerSize',8)
plot(TIMET(ppin:pp:ppend-10), AVNlocal(22, ppin:pp:ppend-10), '-.g','LineWidth',1,'MarkerSize',8)
plot(TIMET(ppinn:pp:ppend), AVNNlocal(22, ppinn:pp:ppend), '-pb','LineWidth',1,'MarkerSize',8)

box on
axis([0 20 -0.01 1.05])
%axis([3 5 0.35 0.6 ])
legend('MBBG','Richards','mKJMA','mKJMA_{loc}','KJMA','Box','off','FontSize',18,'FontWeight','bold','Location','best','NumColumns',1)
%txt = sprintf('Ehrlich');
%text(3, 0.97, txt,'FontName','Helvetica','FontSize',18,'FontWeight','bold');
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.1f');
%xtickformat('%.1f');

%% Fugure 8a  %Sarcoma
texp = [0.3 4.5 6.9 7.9 9.0 10. 11. 12. 14.0 16.5 21.5];
Pexp = [0.01 0.01 0.1 0.2 0.30 0.37 0.48 0.58 0.69 0.82 0.960];
errexp = [0.01 0.01 0.02 0.03 0.035 0.06 0.05 0.05 0.055 0.07 0.02];

ppin = 1;
pp = 3;
ppend = length(TIMET)-110;
ppinn = 89;
ffa = -0.0005 ;
esc = 0.0;

figure;
clf
hold on
plot(TIMET(ppin:pp:ppend), Gomlocal(29, ppin:pp:ppend), '-k','LineWidth',2,'MarkerSize',5,'MarkerFaceColor','r')
plot(TIMET(ppin:pp:ppend-14), real(AVlocal(29, ppin:pp:ppend-14)), '-sb','LineWidth',1,'MarkerSize',8) % sarcoma
plot(TIMET(ppin:pp:ppend), ffa*( real(ealocal(29, ppin:pp:ppend)) ) + esc, '^k','LineWidth',1,'MarkerSize',9,'MarkerFaceColor','k') % sarcoma
plot(texp-6.6,Pexp-0.03,'om','MarkerSize',9,'LineWidth',1.5)
errorbar(texp-6.6,Pexp-0.03,errexp,'om','MarkerSize',9,'LineWidth',1.5)
plot(TIMET(ppinn:pp:ppend), ffa*(real(nealocal(29, ppinn:pp:ppend)) )+ esc  , '-^k','LineWidth',1,'MarkerSize',9,'MarkerFaceColor','k') % sarcoma
plot(TIMET(ppinn:pp:ppend), AVNNlocal(29, ppinn:pp:ppend) , '-sb','LineWidth',1,'MarkerSize',8)

box on
axis([0 10.1 -0.01 .85])
legend('MBBG','mKJMA_{loc}','^{*}Ea','Data','Box','off','FontSize',16,'FontWeight','bold','Location','best','NumColumns',1)
%txt = sprintf('Fibrosarcoma Sa-37');
%text(2.53, 0.97, txt,'FontName','Helvetica','FontSize',18,'FontWeight','bold');
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.1f');
%xtickformat('%.1f');

%%
%% Figure 8b  Ehrlich

texpl = [1 2 3 4 5 10 15];
% Volumen tumoral (mm^3)
Pextl = [700/3000 900/3000 1200/3000 1600/3000 2000/3000 2500/3000 3000/3000];
% Incertidumbre (mm^3)
errexpl = [80/3000 90/3000 120/3000 150/3000 150/3000 180/3000 200/3000];

ppin = 1;
pp = 2;
ppend = length(TIMET)-110;
ppinn = 89;
ffa = -0.0005 ;
esc =0;% 2.31;

figure;
clf
hold on
plot(TIMET(ppin:pp:ppend)+0.05, Gomlocal(22, ppin:pp:ppend)-0.0, '-k','LineWidth',1.6,'MarkerSize',5,'MarkerFaceColor','r')
plot(TIMET(ppin:pp:ppend-14), real(AVlocal(22, ppin:pp:ppend-14)), '-ob','LineWidth',1,'MarkerSize',8) % sarcoma
plot(TIMET(ppin:pp:ppend), ffa*( real(Mealocal(22, ppin:pp:ppend)) ) + esc, '^k','LineWidth',1,'MarkerSize',9,'MarkerFaceColor','k') % sarcoma
plot(texpl-0.9,Pextl-0.215,'om','MarkerSize',9,'LineWidth',1.5)
errorbar(texpl-0.9,Pextl-0.215,errexpl,'om','MarkerSize',9,'LineWidth',1.5)
plot(TIMET(ppinn:pp:ppend), AVNNlocal(22, ppinn:pp:ppend) , '-sb','LineWidth',1,'MarkerSize',8)

box on
axis([0 10.1 -0.01 .85])
legend('MBBG','mKJMA_{loc}','^{*}Ea','Data','Box','off','FontSize',16,'FontWeight','bold','Location','best','NumColumns',1)
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',24,'FontWeight','bold')
ytickformat('%.1f');
%xtickformat('%.1f');

%% Figure 9a time qx days
clear lambda_vals n_vals CHI X
% Parameters
K0 = 0.109;
Ea = 610;
R = 8.314;
T = 310;
t = linspace(0.01,6,30); % move the second parameter for qx time in days
% Sarcoma
au2   = 0.576;
b    = 0.252;
Vo    = 0.5;
Vobs  = 0.15;
cc    = 1.0;
Df    = 1.256;

data = modelTauFractal(au2,Vo,Vobs,Df,cc,bdfVec);

df_model = data(:,1);
Soo.alpha = real(data(:,11));
Soo.beta = real(data(:,12));
Soo.K2 = real(data(:,7));
Soo.n2 = real(data(:,5));
Soo.l2 = real(data(:,3));
Soo.KR = real(data(:,8));
Soo.m = real(data(:,9));
Soo.r = real(data(:,10));

n_vals = sort(Soo.n2(~isnan(Soo.n2)));
lambda_vals = sort(Soo.l2(~isnan(Soo.l2)));
CHI = zeros(length(n_vals),length(lambda_vals));

for i = 1:length(n_vals)
    for j = 1:length(lambda_vals)
        n = n_vals(i);
        lambda = lambda_vals(j);
        K = K0*exp(-Ea/(R*T));
        xi = (K*t).^n;
        X = 1 - (1 + (lambda-1).*xi).^(-1/(lambda-1));
        % efective variable
        Phi = log(xi + 1e-10);
        % susceptibility
        dXdPhi = gradient(X, Phi);
        CHI(i,j) = real(max(dXdPhi));
    end
end

figure;
clf
hold on
h = imagesc(lambda_vals,n_vals,CHI);

set(gca,'YDir','normal')
colorbar;
%xlabel('\lambda');
%ylabel('n');
clim([0, 18])
box on
plot([1 1],[min(n_vals) max(n_vals)],'w--','LineWidth',2)
plot([min(lambda_vals) max(lambda_vals)],[1 1],'w--','LineWidth',2)
axis([0 2.5 0.0 5.05])
txt = sprintf('I');
text(2.35, 1.5, txt,'FontName','Helvetica','FontSize',26,'FontWeight','bold','Color','White');
txt = sprintf('IV');
text(2.3, 0.7, txt,'FontName','Helvetica','FontSize',26,'FontWeight','bold','Color','White');
txt = sprintf('II');
text(0.1, 1.5, txt,'FontName','Helvetica','FontSize',26,'FontWeight','bold','Color','White');
txt = sprintf('III');
text(0.1, 0.7, txt,'FontName','Helvetica','FontSize',26,'FontWeight','bold','Color','White');

ytickformat('%.0f');
xtickformat('%.1f');
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',18,'FontWeight','bold')

%% Figure 9b time qx days

clear lambda_vals n_vals CHI X
% Parameters
K0 = 0.109;
Ea = 1400;
R = 8.314;
T = 310;
t = linspace(0.01,6,30);% move the second parameter for qx time in days
% Ehrlich
au2   = 0.687;
b    = 0.261;
Vo    = 0.5;
Vobs  = 0.190;
cc    = 1.0;
Df    = 1.023;

data = modelTauFractal(au2,Vo,Vobs,Df,cc,bdfVec);

df_model = S1.df;
Soo.alpha = real(data(:,11));
Soo.beta = real(data(:,12));
Soo.K2 = real(data(:,7));
Soo.n2 = real(data(:,5));
Soo.l2 = real(data(:,3));
Soo.KR = real(data(:,8));
Soo.m = real(data(:,9));
Soo.r = real(data(:,10));

n_vals = sort(Soo.n2(~isnan(Soo.n2)));
lambda_vals = sort(Soo.l2(~isnan(Soo.l2)));

CHI = zeros(length(n_vals),length(lambda_vals));

for i = 1:length(n_vals)
    for j = 1:length(lambda_vals)
        n = n_vals(i);
        lambda = lambda_vals(j);
        K = K0*exp(-Ea/(R*T));
        xi = (K*t).^n;
        X = 1 - (1 + (lambda-1).*xi).^(-1/(lambda-1));
        Phi = log(xi + 1e-10);
        dXdPhi = gradient(X, Phi);
        CHI(i,j) = real(max(dXdPhi));
    end
end

figure;
clf
hold on
h = imagesc(lambda_vals,n_vals,CHI);
set(gca,'YDir','normal')

colorbar;
clim([0, 18])
box on
plot([1 1],[min(n_vals) max(n_vals)],'w--','LineWidth',2)
plot([min(lambda_vals) max(lambda_vals)],[1 1],'w--','LineWidth',2)
axis([0 2.5 0.0 5.05])
txt = sprintf('I');
text(2.35, 1.5, txt,'FontName','Helvetica','FontSize',26,'FontWeight','bold','Color','White');
txt = sprintf('IV');
text(2.3, 0.7, txt,'FontName','Helvetica','FontSize',26,'FontWeight','bold','Color','White');
txt = sprintf('II');
text(0.1, 1.5, txt,'FontName','Helvetica','FontSize',26,'FontWeight','bold','Color','White');
txt = sprintf('III');
text(0.1, 0.7, txt,'FontName','Helvetica','FontSize',26,'FontWeight','bold','Color','White');
ytickformat('%.0f');
xtickformat('%.1f');
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',18,'FontWeight','bold')

%%

% Maps
dfVec = 0.60:0.01:0.67;  %%
DfVec = 0.05:0.02:1.4; %%
% Generic
au2   = 0.12;
b    = 0.06;
Vo    = 0.5;
Vobs  = 0.16;
cc    = 1.0;
Df    = 1.1;

alphasurf = NaN(length(dfVec),length(DfVec));
betasurf = alphasurf;
nsurf = alphasurf;
lsurf = alphasurf;
Ksurf = alphasurf;
msurf = alphasurf;
KRsurf = alphasurf;
easurf = alphasurf;

for j = 1:length(DfVec)
    au2 = DfVec(j);
    data = modelTauFractal(au2,Vo,Vobs,Df,cc,bdfVec);
    df_model = S1.df;
    alpha_model = (data(:,11));
    beta_model = (data(:,12));
    n_model = (data(:,5));
    l_model = (data(:,3));
    K_model = (data(:,7));
    m_model = (data(:,9));
    KR_model =(data(:,8));

    for i=1:length(df_model)
        idx = find(abs(dfVec-df_model(i))<1e-4);
        if ~isempty(idx)
            alphasurf(idx,j) = alpha_model(i);
            betasurf(idx,j) = beta_model(i);
            nsurf(idx,j) = n_model(i);
            lsurf(idx,j) = l_model(i);
            Ksurf(idx,j) = K_model(i);
            msurf(idx,j) = m_model(i);
            KRsurf(idx,j) = KR_model(i);

        end
    end

end

[Dfgrid,dfgrid] = meshgrid(DfVec,dfVec);

%% Figure 5a
Kcrit = 1;
phase = zeros(size(alphasurf));
phase(alphasurf > Kcrit) = 1;

figure;
clf
hold on
contourf(Dfgrid,dfgrid,alphasurf,20)

plot(0.687, 0.669,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(0.576, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(0.252, 0.66,'rp','LineWidth',1,'MarkerSize',12,'MarkerFaceColor','r')

plot(0.76, 0.646,'k>','LineWidth',2,'MarkerSize',10)
text(0.79, 0.646, '1','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.576, 0.643,'k>','LineWidth',2,'MarkerSize',10)
text(0.606, 0.643, '8','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.614, 0.630,'k>','LineWidth',2,'MarkerSize',10)
text(0.644, 0.630, '9','FontSize',14,'FontName','Helvetica','FontWeight','bold')
%
plot(0.934, 0.653,'k^','LineWidth',2,'MarkerSize',10)
text(0.964, 0.653, '2','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.997, 0.658,'k^','LineWidth',2,'MarkerSize',10)
text(1.03, 0.658, '3','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(1.181, 0.645,'k^','LineWidth',2,'MarkerSize',10)
text(1.211, 0.645, '4','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(1.255, 0.655,'k^','LineWidth',2,'MarkerSize',10)
text(1.276, 0.655, '5','FontSize',14,'FontName','Helvetica','FontWeight','bold')

plot(0.139, 0.652,'k<','LineWidth',2,'MarkerSize',10)
text(0.159, 0.652, '6','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.065, 0.635,'k<','LineWidth',2,'MarkerSize',10)
text(0.085, 0.635, '7','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.047, 0.656,'k<','LineWidth',2,'MarkerSize',10)
text(0.067, 0.656, '10','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.278, 0.641,'k<','LineWidth',2,'MarkerSize',10)
text(0.298, 0.641, '11','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.044, 0.668,'k<','LineWidth',2,'MarkerSize',10)
text(0.064, 0.668, '12','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.164, 0.646,'k<','LineWidth',2,'MarkerSize',10)
text(0.184, 0.646, '13','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.083, 0.649,'k<','LineWidth',2,'MarkerSize',10)
text(0.103, 0.649, '14','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.143, 0.6435,'k<','LineWidth',2,'MarkerSize',10)
text(0.163, 0.6435, '15','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.078, 0.646,'k<','LineWidth',2,'MarkerSize',10)
text(0.092, 0.646, '16','FontSize',14,'FontName','Helvetica','FontWeight','bold')

ccp = colorbar;
ccp.YTick = [0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6];
ccp.YTickLabel = [0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0.0 1.3 0.6 0.67])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',22,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%% Figure 5b

figure;
clf
hold on
contourf(Dfgrid,dfgrid,betasurf,20)

plot(0.687, 0.669,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(0.576, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(0.252, 0.66,'rp','LineWidth',1,'MarkerSize',12,'MarkerFaceColor','r')

plot(0.76, 0.646,'k>','LineWidth',2,'MarkerSize',10)
text(0.79, 0.646, '1','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.576, 0.643,'k>','LineWidth',2,'MarkerSize',10)
text(0.606, 0.643, '8','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.614, 0.630,'k>','LineWidth',2,'MarkerSize',10)
text(0.644, 0.630, '9','FontSize',14,'FontName','Helvetica','FontWeight','bold')

plot(0.934, 0.653,'k^','LineWidth',2,'MarkerSize',10)
text(0.964, 0.653, '2','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.997, 0.658,'k^','LineWidth',2,'MarkerSize',10)
text(1.03, 0.658, '3','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(1.181, 0.645,'k^','LineWidth',2,'MarkerSize',10)
text(1.211, 0.645, '4','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(1.255, 0.655,'k^','LineWidth',2,'MarkerSize',10)
text(1.276, 0.655, '5','FontSize',14,'FontName','Helvetica','FontWeight','bold')

plot(0.139, 0.652,'k<','LineWidth',2,'MarkerSize',10)
text(0.159, 0.652, '6','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.065, 0.635,'k<','LineWidth',2,'MarkerSize',10)
text(0.085, 0.635, '7','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.047, 0.656,'k<','LineWidth',2,'MarkerSize',10)
text(0.067, 0.656, '10','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.278, 0.641,'k<','LineWidth',2,'MarkerSize',10)
text(0.298, 0.641, '11','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.044, 0.668,'k<','LineWidth',2,'MarkerSize',10)
text(0.064, 0.668, '12','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.164, 0.646,'k<','LineWidth',2,'MarkerSize',10)
text(0.184, 0.646, '13','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.083, 0.649,'k<','LineWidth',2,'MarkerSize',10)
text(0.103, 0.649, '14','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.143, 0.6435,'k<','LineWidth',2,'MarkerSize',10)
text(0.163, 0.6435, '15','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.078, 0.646,'k<','LineWidth',2,'MarkerSize',10)
text(0.092, 0.646, '16','FontSize',14,'FontName','Helvetica','FontWeight','bold')

ccp = colorbar;
ccp.YTick = [0.02 0.1 0.2 0.3 0.4 0.5 0.6];
ccp.YTickLabel = [.0 0.1 0.2 0.3 0.4 0.5 0.6];
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0.0 1.3 0.6 0.67])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',18,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%% Figure 5c
Kcrit = 1;
phase = zeros(size(nsurf));
phase(nsurf > Kcrit) = 1;

figure;
clf
hold on
contourf(Dfgrid,dfgrid,nsurf,20)
contour(Dfgrid,dfgrid,nsurf,[Kcrit Kcrit],'k','LineWidth',4)

plot(0.687, 0.669,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(0.576, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(0.252, 0.66,'rp','LineWidth',1,'MarkerSize',12,'MarkerFaceColor','r')

plot(0.76, 0.646,'k>','LineWidth',2,'MarkerSize',10)
text(0.79, 0.646, '1','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.576, 0.643,'k>','LineWidth',2,'MarkerSize',10)
text(0.606, 0.643, '8','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.614, 0.630,'k>','LineWidth',2,'MarkerSize',10)
text(0.644, 0.630, '9','FontSize',14,'FontName','Helvetica','FontWeight','bold')

plot(0.934, 0.653,'k^','LineWidth',2,'MarkerSize',10)
text(0.964, 0.653, '2','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.997, 0.658,'k^','LineWidth',2,'MarkerSize',10)
text(1.03, 0.658, '3','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(1.181, 0.645,'k^','LineWidth',2,'MarkerSize',10)
text(1.211, 0.645, '4','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(1.255, 0.655,'k^','LineWidth',2,'MarkerSize',10)
text(1.276, 0.655, '5','FontSize',14,'FontName','Helvetica','FontWeight','bold')

plot(0.139, 0.652,'k<','LineWidth',2,'MarkerSize',10)
text(0.159, 0.652, '6','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.065, 0.635,'k<','LineWidth',2,'MarkerSize',10)
text(0.085, 0.635, '7','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.047, 0.656,'k<','LineWidth',2,'MarkerSize',10)
text(0.067, 0.656, '10','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.278, 0.641,'k<','LineWidth',2,'MarkerSize',10)
text(0.298, 0.641, '11','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.044, 0.668,'k<','LineWidth',2,'MarkerSize',10)
text(0.064, 0.668, '12','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.164, 0.646,'k<','LineWidth',2,'MarkerSize',10)
text(0.184, 0.646, '13','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.083, 0.649,'k<','LineWidth',2,'MarkerSize',10)
text(0.103, 0.649, '14','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.143, 0.6435,'k<','LineWidth',2,'MarkerSize',10)
text(0.163, 0.6435, '15','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.078, 0.646,'k<','LineWidth',2,'MarkerSize',10)
text(0.092, 0.646, '16','FontSize',14,'FontName','Helvetica','FontWeight','bold')

ccp = colorbar;
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0.0 1.3 0.6 0.67])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',18,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%% Figure 5d
figure;
clf
hold on
contourf(Dfgrid,dfgrid,0.1*real(Ksurf),20)

plot(0.687, 0.669,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(0.576, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(0.252, 0.66,'rp','LineWidth',1,'MarkerSize',12,'MarkerFaceColor','r')

plot(0.76, 0.646,'k>','LineWidth',2,'MarkerSize',10)
text(0.79, 0.646, '1','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.576, 0.643,'k>','LineWidth',2,'MarkerSize',10)
text(0.606, 0.643, '8','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.614, 0.630,'k>','LineWidth',2,'MarkerSize',10)
text(0.644, 0.630, '9','FontSize',14,'FontName','Helvetica','FontWeight','bold')

plot(0.934, 0.653,'k^','LineWidth',2,'MarkerSize',10)
text(0.964, 0.653, '2','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.997, 0.658,'k^','LineWidth',2,'MarkerSize',10)
text(1.03, 0.658, '3','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(1.181, 0.645,'k^','LineWidth',2,'MarkerSize',10)
text(1.211, 0.645, '4','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(1.255, 0.655,'k^','LineWidth',2,'MarkerSize',10)
text(1.276, 0.655, '5','FontSize',14,'FontName','Helvetica','FontWeight','bold')

plot(0.139, 0.652,'k<','LineWidth',2,'MarkerSize',10)
text(0.159, 0.652, '6','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.065, 0.635,'k<','LineWidth',2,'MarkerSize',10)
text(0.085, 0.635, '7','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.047, 0.656,'k<','LineWidth',2,'MarkerSize',10)
text(0.067, 0.656, '10','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.278, 0.641,'k<','LineWidth',2,'MarkerSize',10)
text(0.298, 0.641, '11','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.044, 0.668,'k<','LineWidth',2,'MarkerSize',10)
text(0.064, 0.668, '12','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.164, 0.646,'k<','LineWidth',2,'MarkerSize',10)
text(0.184, 0.646, '13','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.083, 0.649,'k<','LineWidth',2,'MarkerSize',10)
text(0.103, 0.649, '14','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.143, 0.6435,'k<','LineWidth',2,'MarkerSize',10)
text(0.163, 0.6435, '15','FontSize',14,'FontName','Helvetica','FontWeight','bold')
plot(0.078, 0.646,'k<','LineWidth',2,'MarkerSize',10)
text(0.092, 0.646, '16','FontSize',14,'FontName','Helvetica','FontWeight','bold')

ccp = colorbar;
ticks = get(ccp, 'Ticks');
tickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'un', 0);
set(ccp, 'TickLabels', tickLabels);
box on
axis([0.0 1.3 0.6 0.67])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',18,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%%
%% Definition of functions
%%
function val = Richards(KR,m,r,Vobs,t)
val = KR./((1 + ((KR/Vobs).^(m) -1).*exp(-r.*m*t)).^(1/m));
end

function val = Avrami(K,n,lambda,t)
val = 1 - (1 + (lambda-1).*(K*t).^(n)).^(-1/(lambda-1));
end

function val = Gompertz(alpha,beta,Vobs,Vo,t)
val = (Vobs*(Vo/Vobs).^(exp(-beta*t))).*exp((alpha./beta).*(1-exp(-beta*t)));
end

function val = nlocal(K,n,lambda,t)
pt = Avrami(K,n,lambda,t);
val = log(  ((1-pt).^(1-lambda) -1)./(lambda-1) ) ;
end

function val = nnlocal(K,n,t)
pt = 1 - exp(-(K*t).^n);
val = log(- log( 1-pt ) ) ;
end

function val = qTau(a,b,Vo,Vobs)
val = (-Vobs + ...
    exp(a/b + a./(-a + b*log(Vo/Vobs))) .* ...
    Vobs .* ((Vo/Vobs).^(b./(a - b*log(Vo/Vobs))))) ...
    ./ (Vobs*exp(a/b) - Vobs);
end

function val = qpTau(a,b,Vo,Vobs)
val = (b* ...
    exp(a/b + a./(-a + b*log(Vo/Vobs))) .* ...
    Vobs .* ((Vo/Vobs).^(b./(a - b*log(Vo/Vobs))))) ...
    ./ (Vobs*exp(a/b) - Vobs);
end

function val = Tau(a,b,Vo,Vobs)
val = (1/b)*log((a/b - log(Vo/Vobs)));
end
%%
toc
