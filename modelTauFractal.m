% =========================================================================
% NAME:  Model Tau fractal
% DESCRIPTION:  Computational framework for analyzing unperturbed tumor
%               growth through mKJMA, MBBG, and Richards models.
%               Relationship and role of fractal dimension.
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

function out = modelTauFractal(au2,Vo,Vobs,Df,cc,bdfVec)

data = zeros(length(bdfVec),12);

kvalid = 1;

for k = 1:length(bdfVec)
    
    df = bdfVec(k);

    a = AlphaF(au2,df,Df,Vobs);
    b = BetaF(au2,df,Df);

    if ~isreal(a) || ~isreal(b) || b==0
    continue
    end

    if (a/b - log(Vo/Vobs)) <= 0
    continue
    end

% ========= ROOT lambda ==========
    froot = @(x) ...
    (-1)/(1 - qTau(a,b,Vo,Vobs) + Tau(a,b,Vo,Vobs)*qpTau(a,b,Vo,Vobs)) ...
    + (-Tau(a,b,Vo,Vobs)*qpTau(a,b,Vo,Vobs).*x) ./ ...
      ((1-qTau(a,b,Vo,Vobs))*(1-qTau(a,b,Vo,Vobs) ...
      + Tau(a,b,Vo,Vobs)*qpTau(a,b,Vo,Vobs))) ...
    + (1/(1-qTau(a,b,Vo,Vobs))).^x;

    if ~isfinite(froot(0.5))
    continue
    end

    try
    r1 = fzero(froot,0.5);
    r2 = fzero(froot,2);
    catch
    continue
    end

    if abs(r1-1) > 1e-6
        lambda = r1;
    else
        lambda = r2;
    end

% ========= ROOT m ==========
    fm = @(m) ...
    (m+1)*Vobs*exp(a/b) ...
    - (Vobs*exp(a/b)-Vo)*b ./ ...
    (a - b*log(((qTau(a,b,Vo,Vobs))*(Vobs*exp(a/b)-Vo)+Vo)/Vobs));

    try
    mm = fsolve(fm,1,optimoptions('fsolve','Display','off'));
    catch
    continue
    end

% ========= DATA ==========
    data(kvalid,:) = [ ...
        df , ...
        cc , ...
        lambda , ...
        nTau(a,b,Vo,Vobs,cc) , ...
        nTau(a,b,Vo,Vobs,lambda) , ...
        10*KTau(a,b,Vo,Vobs,cc) , ...
        10*KTau(a,b,Vo,Vobs,lambda) , ...
        Vobs*exp(a/b) , ...
        mm , ...
        b/mm , ...
        a , ...
        b ];

        kvalid = kvalid + 1;
end

writematrix(data,'datosfractalsar.csv');

out = data;

end

%%% Functions

function val = BetaF(u2,df,Df)
val = u2*(1 - df/Df);
end

function val = AlphaF(u2,df,Df,Vobs)
val = u2*log((2/3*df - 1)/(df - 1)) ...
    - BetaF(u2,df,Df)*log(Vobs);
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

function val = nTau(a,b,Vo,Vobs,lambda)
val = 1 + lambda*(Tau(a,b,Vo,Vobs)*qpTau(a,b,Vo,Vobs)) ...
    /(1 - qTau(a,b,Vo,Vobs));
end

function val = KTau(a,b,Vo,Vobs,lambda)
nt = nTau(a,b,Vo,Vobs,lambda);
val = (1/Tau(a,b,Vo,Vobs))*((nt-1)/(nt+lambda-1)).^(1/nt);
end