% =========================================================================
% NAME:  Hausdorff Distance
% DESCRIPTION:  Computational framework for analyzing unperturbed tumor
%               growth through mKJMA, MBBG, and Richards models, focusing
%               on Hausdorff distance
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
clear all;
tic
%% ===============================
% Domine
%% ===============================
t = linspace(0.01,30,200);
%% ===============================
% Models
%% ===============================
mbbg_model     = @(t,a,b,Vobs,Vo) (Vobs*((Vo/Vobs).^(exp(-b*t))).*exp((a./b).*(1 -exp(-b*t))) - Vo)./(Vobs*exp(a./b) - Vobs);
richards_model = @(t,KR,m,r,Vobs) (KR./((1 + ((KR/(Vobs)).^(m) -1).*exp(-r.*m*t)).^(1/m)) - Vobs)./(KR - Vobs);
function q = q_model_safe(t,lambda,K,n)

if lambda <= 0 || K <= 0 || n <= 0
    q = nan(size(t));
    return;
end

xi = (K*t).^n;

arg = 1 + (lambda-1)*xi;

if any(arg <= 0)
    q = nan(size(t));
    return;
end

%  lambda = 1 (limit)
if abs(lambda - 1) < 1e-3
    q = 1 - exp(-xi);
else
    q = 1 - arg.^(-1/(lambda-1));
end

q(~isfinite(q)) = nan;
end

%% ===============================
% Hausdorff
%% ===============================
function dH = hausdorff_distance(P,Q)

if isempty(P) || isempty(Q)
    dH = 1e6; return;
end

if any(~isfinite(P(:))) || any(~isfinite(Q(:)))
    dH = 1e6; return;
end

try
    D = pdist2(P,Q);
catch
    dH = 1e6; return;
end

if any(~isfinite(D(:)))
    dH = 1e6; return;
end

d1 = max(min(D,[],2));
d2 = max(min(D,[],1));
dH = max(d1,d2);
end
%% ===============================
% Susceptibility
%% ===============================
function chi = susceptibility(t,X,K,n)
dXdt = gradient(X,t);
Phi = log((K*t).^n + 1e-12);
dPhidt = gradient(Phi,t);
chi = dXdt ./ (dPhidt + 1e-12);
chi(~isfinite(chi)) = 0;
end

%% ===============================
% Objetive function
%% ===============================
gamma = 0; 

obj_fun = @(x) combined_distance(t,p,q_model_safe,x,gamma);

function D = combined_distance(t,p,q_model,x,gamma)

lambda = x(1);
K = x(2);
n = x(3);

if lambda <= 0 || K <= 0 || n <= 0
    D = 1e6; return;
end

q = q_model(t,lambda,K,n);

if any(~isfinite(q))
    D = 1e6; return;
end

% --- curves ---
P = [t(:), p(:)];
Q = [t(:), q(:)];

dH_curve = hausdorff_distance(P,Q);

Phi = log((K*t).^n + 1e-12);

chi_p = gradient(p, Phi);
chi_q = gradient(q, Phi);

if any(~isfinite(chi_p)) || any(~isfinite(chi_q))
    D = 1e6; return;
end

Pchi = [t(:), chi_p(:)];
Qchi = [t(:), chi_q(:)];

dH_chi = hausdorff_distance(Pchi,Qchi);

D = dH_curve + gamma*dH_chi;
end

%% ===============================
%  MBBG parameters
%% ===============================
alpha = 0.4;
beta  = 0.1;
Vobs  = 0.15;
Vo    = 0.5;
p = mbbg_model(t,alpha,beta,Vobs,Vo);

%% ===============================
%  Richards parameters
%% ===============================
KR = 9;
m  = 0.026;
r = 5;
Vobs  = 0.15;
p2 = richards_model(t,KR,m,r,Vobs);

%% ===============================
% Optimization
%% ===============================

options = optimset('Display','iter','TolX',1e-4);

%% ===============================
% Visualization
%% ===============================
%
% Maps
bdfVec = 0.05:0.005:3.0;
dfVec = 0.4:0.02:0.8;
DfVec = 1.0:0.02:1.45;

% Generic tumor
au2   = 0.620;
b    = 0.252;
Vo    = 0.5;
Vobs  = 0.15;
cc    = 1.0;
Df    = 1.1;

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

HaussurfG = NaN(length(dfVec),length(DfVec));
HaussurfR = NaN(length(dfVec),length(DfVec));

for j = 1:length(DfVec)

    Df = DfVec(j);
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

    for i=1:length(df_model)
        idx = find(abs(dfVec-df_model(i))<1e-4);
        if ~isempty(idx)

            [x_opt, Dmin] = fminsearch(@(x) ...
                combined_distance(t,mbbg_model(t,Soo.alpha(i),Soo.beta(i),Vobs,Vo),@q_model_safe,x,1), ...
                [Soo.l2(i),0.1*Soo.K2(i),Soo.n2(i)], options);

            HaussurfG(idx,j) = Dmin;

            [x_opt, Dmin] = fminsearch(@(x) ...
                combined_distance(t,richards_model(t,Soo.KR(i),Soo.m(i),Soo.r(i),Vobs),@q_model_safe,x,1), ...
                [Soo.l2(i),0.1*Soo.K2(i),Soo.n2(i)], options);

            HaussurfR(idx,j) = Dmin;

        end
    end

end

[Dfgrid,dfgrid] = meshgrid(DfVec,dfVec);

%% Figure 8c
HaussurfG1 = HaussurfG;
HaussurfG1(find(HaussurfG>1.3))=NaN;
HaussurfG1(find(HaussurfG<0.0228))=NaN;

figure;
clf
hold on
contourf(Dfgrid,dfgrid,HaussurfG1,20)
plot(1.256, 0.688,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(1.023, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')

ccp = colorbar;
ccp.Ruler.Exponent = -2;
box on
axis([1.0 1.4 0.4 0.75])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',22,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%% Figure 8d
HaussurfR1 = HaussurfR;
HaussurfR1(find(HaussurfR>0.039))=NaN;
HaussurfR1(find(HaussurfR<0.00918))=NaN;

figure;
clf
hold on
contourf(Dfgrid,dfgrid,HaussurfR1,20)
plot(1.256, 0.688,'rs','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')
plot(1.023, 0.611,'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r')

ccp = colorbar;
ccp.Ruler.Exponent = -3;
box on
axis([1.0 1.4 0.4 0.75])
set(gca,'linewidth',2,'FontName','Helvetica','FontSize',22,'FontWeight','bold')
ytickformat('%.2f');
xtickformat('%.1f');

%%
toc