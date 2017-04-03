function [] = summary0(ID) % ID = 1:36: 3 levels, 3 models, 4 variables
global X Y N p W

ch = str2double(num2str(ID));
nlevel = ceil(ch/12);
ch = ch - (nlevel-1)*12;
nmodel = ceil(ch/4);
ch = ch - (nmodel-1)*4;
nvar = ch; %ceil(ch/5);
fprintf('level=%d, model=%d, var=%d\n',[nlevel,nmodel,nvar])

% nlevel = 3; nmodel = 1; nvar = 1;

% Allmat = []; Alldic = zeros(3,12); Allr = [];

% for nlevel = 1:3
load(strcat('fulldat',num2str(nlevel),'.mat'))
Y0 = Y; X = zscore(X); X = [ones(size(X,1),1), X]; [N,p] = size(X);
% mats = [];
%     for nvar = 1:4
Y = Y0(:,nvar);
%         for nmodel = 1:3
% fprintf('level=%d, model=%d, var=%d\n',[nlevel,nmodel,nvar])
M = diag(sum(W,1));
if nmodel==3
    M = eye(N);
end
% [a,b] = regress(Y,X); disp([a,b])
chs = [1,2,3,4,5]; nch = length(chs);
niter = 6e3; burn = 5e3; thin = 1;
nsample = (niter-burn)/thin; tot = nch*nsample; npara = p+2;
MCSamples = nan(nch, npara, nsample);
mat = nan(npara,4);
for ch = 1:nch
    load(strcat('out',num2str(nlevel),num2str(nmodel),num2str(nvar),num2str(ch),'.mat'))
    MCSamples(ch,:,:) = xall((burn+1):end,:)';
end
R = psrf(MCSamples);
% if nmodel > 1
%    R = R(1:(length(R)-1));
% end
% display([R(1:2),mean(R(3:end)), std(R(3:end))])
% Allr = [Allr,R'];

for j = 1:npara
    vec = reshape(MCSamples(:,j,:),[tot,1]);
    % [lb, ub] = FindHPDset(vec, 0.95, []);
    qs = quantile(vec, [0.025, 0.975]);
    mat(j,:) = [mean(vec), std(vec), qs(1), qs(2)];
end
if nmodel > 1
    mat(end,:) = zeros(1,4);
end
disp(mat)

% compute DIC
betas = zeros(tot,p);
for j = 1:p
    betas(:,j) = reshape(MCSamples(:,j,:),[tot,1]);
end
tau2s = reshape(MCSamples(:,p+1,:),[tot,1]);
gammas = reshape(MCSamples(:,p+2,:),[tot,1]);
D1 = 0;
for i = 1:tot
    V = M-gammas(i)*W; L = chol(V,'lower');
    delta = Y - X*betas(i,:)'; delta = L'*delta;
    D1 = D1 - 0.5*N*log(2*pi) - 0.5*N*log(tau2s(i)) + sum(log(diag(L)))...
        - 0.5/tau2s(i).*sum(delta.^2);
end
D1 = -2*D1/tot;
tau2hat = mean(tau2s);
gammahat = mean(gammas);
betahat = mean(betas, 1);
V = M-gammahat*W; L = chol(V,'lower');
delta = Y - X*betahat'; delta = L'*delta; %residual
D2 = - 0.5*N*log(2*pi) - 0.5*N*log(tau2hat) + sum(log(diag(L)))...
    - 0.5/tau2hat.*sum(delta.^2);  %log-likelihood
D2 = -2*D2; % bug found 2013-12-02
DIC = 2*D1 - D2;
disp(DIC)

%             Alldic(nlevel, (nvar-1)*3+nmodel) = DIC;
%             mats = [mats, mat];
%         end
%     end
%     Allmat = [Allmat; mats];
% end
% save('outAll.mat','Allmat','Alldic','Allr')
save(strcat('summa',num2str(nlevel),num2str(nmodel),num2str(nvar),'.mat'),'mat','DIC','R')
end

function [] = collectAll()
Allmat = []; Alldic = zeros(3,12); Allr = [];
for nlevel = 1:3
    mats = [];
    for nvar = 1:4
        for nmodel = 1:3
            load(strcat('summa',num2str(nlevel),num2str(nmodel),num2str(nvar),'.mat'))
            Alldic(nlevel, (nvar-1)*3+nmodel) = DIC;
            mats = [mats, mat]; Allr = [Allr, R'];
        end
    end 
    Allmat = [Allmat; mats];
end
save('outAll.mat','Allmat','Alldic','Allr')
end

function [LBout,UBout] = FindHPDset(Samples,p,npoints)
%Function to find the 100p % HPD set based on Samples
if isempty(npoints)
    npoints=200;
end

[f,x] = ksdensity(Samples,'npoints',npoints); N = size(Samples,1); maxf = max(f);
step = maxf/npoints;
HPDdensity = (step:step:maxf);
NHPD = size(HPDdensity,2);
LB = cell(1,NHPD); UB = cell(1,NHPD); Prob = zeros(1,NHPD);
for i=1:NHPD
    indices0 = find(HPDdensity(NHPD-i+1) < f);
    if ~isempty(indices0)
        indices1 = find(diff(indices0)> 1);
        if isempty(indices1)
            LB{i} = x(indices0(1)); UB{i} = x(indices0(end));
        elseif (size(indices1,1)==1)
            LB{i} = [x(indices0(1)) x(indices0(indices1(1)+1))];
            UB{i} = [x(indices0(indices1(1))) x(indices0(end))];
        else
            LB{i} = x(indices0(1)); UB{i} = [];
            for j=1:(size(indices1,2)-1)
                LB{i} = [LB{i} x(indices0(indices1(j)+1))];
                UB{i} = [UB{i} x(indices0(indices1(j)))];
            end
            UB{i} =[UB{i} x(indices0(end))];
        end
    end
    Ns = size(LB{i},2);
    count = zeros(1,Ns);
    for j=1:Ns
        count(j) = sum((LB{i}(j) <= Samples).*(Samples <= UB{i}(j)));
    end
    Prob(i) = sum(count)/N;
end
[minval indexmin] = min(abs(Prob - p));
LBout = LB{indexmin};
UBout = UB{indexmin};
end

function [R,neff,V,W,B] = psrf(varargin)
%PSRF Potential Scale Reduction Factor
%
%   [R,NEFF,V,W,B] = PSRF(X) or
%   [R,NEFF,V,W,B] = PSRF(x1,x2,...,xs)
%   returns "Potential Scale Reduction Factor" (PSRF) for
%   collection of MCMC-simulations. X is a NxDxM matrix
%   which contains M MCMC simulations of length N, each with
%   dimension D. MCMC-simulations can be given as separate
%   arguments x1,x2,... which should have the same length.
%
%   Returns
%     R     PSRF in a row vector of length D
%     neff  estimated effective number of samples M*N*V/B
%     V     estimated mixture-of-sequences variances
%     W     estimated within sequence variances
%     B     estimated between sequence variances
%
%   The idea of the PSRF is that if R is not near 1 (below 1.1 for
%   example) one may conclude that the tested samples were not from
%   the same distribution (chain might not have been converged
%   yet).
%
%   If only one simulation is given, the factor is calculated
%   between first and last third of the chain. Note that use of
%   only one chain will produce over-optimistic result.
%
%   Method is from:
%      Brooks, S.P. and Gelman, A. (1998) General methods for
%      monitoring convergence of iterative simulations. Journal of
%      Computational and Graphical Statistics. 7, 434-455. Note that
%      this function returns square-root definiton of R (see Gelman
%      et al (2003), Bayesian Data Analsyis p. 297).
%
%   See also
%     CPSRF, MPSRF, IPSRF

% Copyright (C) 1999 Simo Särkk?
% Copyright (C) 2003 Aki Vehtari
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

% 2004-01-22 Aki.Vehtari@hut.fi Added neff, R^2->R, and cleaning

% In case of one argument split to two halves (first and last thirds)
onechain=0;
if nargin==1
    X = varargin{1};
    if size(X,3)==1
        n = floor(size(X,1)/3);
        x = zeros([n size(X,2) 2]);
        x(:,:,1) = X(1:n,:);
        x(:,:,2) = X((end-n+1):end,:);
        X = x;
        onechain=1;
    end
elseif nargin==0
    error('Cannot calculate PSRF of scalar');
else
    X = zeros([size(varargin{1}) nargin]);
    for i=1:nargin
        X(:,:,i) = varargin{i};
    end
end

[N,D,M]=size(X);

if N<1
    error('Too few samples');
end

% Calculate means W of the variances
W = zeros(1,D);
for n=1:M
    x = X(:,:,n) - repmat(mean(X(:,:,n)),N,1);
    W = W + sum(x.*x);
end
W = W / ((N-1) * M);

% Calculate variances B (in fact B/n) of the means.
Bpn = zeros(1,D);
m = mean(reshape(mean(X),D,M)');
for n=1:M
    x = mean(X(:,:,n)) - m;
    Bpn = Bpn + x.*x;
end
Bpn = Bpn / (M-1);

% Calculate reduction factors
S = (N-1)/N * W + Bpn;
R = (M+1)/M * S ./ W - (N-1)/M/N;
V = R .* W;
R = sqrt(R);
B = Bpn*N;
neff = min(M*N*V./B,M*N);
if onechain & (nargout>1)
    neff=neff*3/2;
end
end
