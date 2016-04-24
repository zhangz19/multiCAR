%% Multi-resolution CAR model for spatial areal data
% * The code is developed for modeling spatial areal data with
% multi-resolution nested spatial random effects. 
% For each resolution (group level), it has options: nonrandom, spatial (CAR), WLS (spatially weighted
% Least-Square) and OLS (Ordinary Least-Square with i.i.d. assumption). 
% It returns DIC4 for model comparisons with different model ID.
% * For detail:
% Daneshvar, F., Nejadhashemi, A. P., Zhang, Z., Herman, H. R., Shortridge, A.
% and Marquart-Pyatt, S. (2016). Evaluating Stream Health Based Environmental Justice Model
% Performance at Different Spatial Scales. _Journal of Hydrology_. To appear.
% * Please contact the authors if there are any questions or implementation issues:
%    Zhen Zhang, zhangz19@galton.uchicago.edu. 
%           Date coded: 2013-12-1

function [] = multiCAR(ID)
% main function
global X Z N p modind Ws Ms nlevel vec Ns Inds gapss gammass loglike0s npara ...
    verbose computeDIC dicLen %#ok<*NUSED>

% MCMC setup: and also see runModel() function below
nch = 3; % number of chains
nvar = 4; % number of response variables for separate analysis
niter = 1e3;  % number of total posterior samples
burn = 0;  % numrber of burn-in period of MCMC
thin = 1; % number of thin for posterior samples

verbose = 0; % show MCMC iteration id for tracking
computeDIC = 1; % compute DIC or not
dicLen = 50; %number of subsequent runs for calculating DIC4

id = str2double(num2str(ID)); % for one specific model id code

%translate code
if ~exist('modind.mat','file')
    getModind(); % pre-run, need input file grpind.mat. see getModind() function below
else
    load('modind.mat');
end
vec = modind(id,:); %vec = [county, cesus tract, block group]
% vec indicate model code for each group level: 
% 0 = non-random, 1 = spatial, 2 = weighted Least-square (LS), 3 = ordinary LS
disp(vec)
nlevel = find(vec>0,1,'last');
Ws = cell(1,nlevel); Ms = cell(1,nlevel); 
gapss = cell(1,nlevel); gammass = cell(1,nlevel); loglike0s = cell(1,nlevel);
for j = 1:nlevel
    load(strcat('fulldat',num2str(j),'.mat'))
    Ws{j} = W;
    Ms{j} = eye(size(X,1));
    if any([1,2] == vec(j)) % spatial or WLS
        if ~exist('EV.mat','file') % pre-calculation to greatly save computational efforts for CAR 
            M = diag(sum(W,1)); invM = inv(M);
            eigs = eig(sqrt(invM)*W*sqrt(invM));
            lgamma = max(1/min(eigs),-1); ugamma = 1/max(eigs);
            gaps = 1e-3; gammas = (lgamma+gaps):gaps:(ugamma-gaps); len = length(gammas);
            loglike0 = zeros(1,len);
            for i = 1:len
                gamma = gammas(i);
                B = (M-gamma*W);
                loglike0(i) = 0.5*sum(log(eig(B)));
            end
            save(strcat('EV',num2str(j),'.mat'),'M','gaps','gammas','loglike0')
        end
        load(strcat('EV',num2str(j),'.mat'))
        Ms{j} = M; gapss{j} = gaps; gammass{j} = gammas; loglike0s{j} = loglike0;
    end
end
X = zscore(X);
X = [ones(size(X,1),1), X];
[N,p] = size(X);

nsample = (niter-burn)/thin; tot = nch*nsample; npara = p+2*3;
% Rs = cell(1,nvar); % for MCMC convergence diagnostics
matparas = cell(1,nvar);

% set random seeds for reproducibility of results
rng('default'); rng(id*10); 

dics = nan(1,nvar); Yhat = nan(N,nvar);
for yind = 1:nvar
    fprintf('nvar = %d:\n', yind)
    MCSamples = nan(nch, npara, nsample);
    matparaAll = nan(tot, npara); EsAll = nan(tot, 2);
    for ch = 1:nch
        [xall, Es, Us] = runModel(Y, yind, ch);
        MCSamples(ch, :, :) = xall(burn+(1:nsample)*thin,:)';
        matparaAll((ch-1)*nsample+(1:nsample),:) = xall(burn+(1:nsample)*thin,:);
        EsAll((ch-1)*nsample+(1:nsample),:) = Es(burn+(1:nsample)*thin,:);
    end
    % % check convergence, need psrf.m file
    % R = psrf(MCSamples); disp(R); Rs{yind} = R; 
    matparas{yind} = matparaAll;
    meanPara = mean(matparaAll,1);
    Yhat(:,yind) = X*meanPara(1:p)' + sum(Us,2)/tot;
    meanEs = mean(EsAll,1);
    dics(yind) = -4*meanEs(1) + 2*meanEs(2);
end

save(strcat('paras',num2str(id),'.mat'),'matparas','meanEs','dics','Yhat') %'Rs',
end

%%
function [] = getModind()
% need input file: grpind.mat which contains:
% igrp2_1: integer vector of group level 1 for observations at group level 2
% igrp3_1: integer vector of group level 1 for observations at group level 3
% igrp3_2: integer vector of group level 2 for observations at group level 3
load('grpind.mat')
Ns = [max(igrp2_1), max(igrp3_2), length(igrp3_2)];
Inds = cell(2); Inds{1,1} = igrp2_1; Inds{2,1} = igrp3_1; Inds{2,2} = igrp3_2;

Z = cell(2);
for j = 1:2
    for k = 1:j
        Z{j,k} = zeros(Ns(j+1), Ns(k));
        for i = 1:Ns(k)
            Z{j,k}(Inds{j,k}==i,i) = 1;
        end
    end
end

% construct model index matrix
modind = zeros(6,3); modind(1:3,1) = 1:3; modind(4:6,2) = 1:3;
tmp = zeros(3*3,3); tmp(:,1:2) = [kron([1:3]',ones(3,1)), kron(ones(3,1),[1:3]')];
modind = [modind; tmp];
tmp = [kron([0:3]',ones(4,1)), kron(ones(4,1),[0:3]')];
tmp = [repmat(tmp,[3,1]), kron([1:3]', ones(16,1))];
modind = [modind; tmp];
disp(size(unique(modind, 'rows'),1))
% modind = [kron(modind, ones(nch,1)), kron(ones(size(modind,1),1), (1:nch)')];
save('modind.mat', 'modind','Ns','Inds','Z')
end

%%
function [xall, Es, Us] = runModel(Y0, yind, ch)
global X Z Ms Ws nlevel vec p Ns Inds gapss gammass loglike0s npara verbose computeDIC  a_tau2 invb_tau2

niter = 6e3; % number of chains 
burn = 5e3; % numrber of burn-in period of MCMC

Y = Y0(:,yind);
% initial values from pre-run
meantau2 = .01; vartau2 = 1e2;
a_tau2 = 2+meantau2^2/vartau2;
invb_tau2 = meantau2*(a_tau2-1);

% set initial values
x.beta = (X'*X)\(X'*Y).*(1 + (ch-3)/5*normrnd(0,1,[p,1]));
x.tau2 = zeros(1,length(vec)) + (vec~=0); x.gamma = zeros(1,length(vec));
x.tau2(nlevel) = sum((diag(Ms{nlevel}).*(Y-X*x.beta)).^2)/(Ns(nlevel)-p);
err = cell(1,nlevel); %random effect, last cell = error term
for j = 1:nlevel
    err{j} = zeros(Ns(j),1);
end

% run the chain, store the results
xall = nan(niter-burn, npara);
Es = nan(niter-burn, 2);
Us = zeros(Ns(1),1);
if nlevel > 1
    Us = zeros(Ns(nlevel), nlevel-1);
end

tic
for i = 1:niter
    
    [x, err, V] = updateBeta(x, err, Y);  Delta = Y - X*x.beta;
    
    % Update random effects, Tau2 and Gamma
    err{nlevel} = Delta;
    for j = 1:nlevel
        if vec(j)~=0
            V1 = V;
            if j < nlevel % update random effects
                Delta0 = Delta; ind0 = find((1:(nlevel-1))~=j);
                for k = ind0
                    Delta0 = Delta0 - err{k}(Inds{nlevel-1,k});
                end
                Sigma = Z{nlevel-1,j}'*V;  Mu = Sigma*Delta0;
                V1 = Ms{j} - x.gamma(j)*Ws{j};
                Sigma = Sigma*Z{nlevel-1,j} + V1*(x.tau2(nlevel)/x.tau2(j));
                Lo = chol(Sigma, 'lower'); Mu = Lo\Mu;
                Mu = Mu + normrnd(0, sqrt(x.tau2(nlevel)), [length(Mu),1]);
                err{j} = Lo'\Mu;
                err{nlevel} = err{nlevel} - err{j}(Inds{nlevel-1,j});
            end
            x.tau2(j) = 1/gamrnd(a_tau2+0.5*Ns(j), 1/(invb_tau2+0.5*err{j}'*V1*err{j}));
            if vec(j) == 1 % spatial model only
                x = updateGamma(x, err, j);
            end
        end
    end
    
    if i > burn
        xall(i-burn,:) = [x.beta', x.tau2, x.gamma];
        if computeDIC == 1
            Es(i-burn,:) = getDIC(x, err, Y);
        end
        if nlevel > 1
            for j = 1:(nlevel-1)
                Us(:,j) = Us(:,j) + err{j}(Inds{nlevel-1,j});
            end
        end
    end
    
    if verbose == 1
        fprintf('%6d', i)
        if(~mod(i,20))
            fprintf('\n')
        end
    end
end
runtime = toc/3600;
fprintf('\n%d iterations are done with elapsed time %.4f hours.\n', niter, runtime)
end

%%
function [x, err, V] = updateBeta(x, err, Y)
% sample fixed-effects
global nlevel Inds Ms Ws X
Delta = Y;
for j = 1:(nlevel-1)
    Delta = Delta - err{j}(Inds{nlevel-1,j});
end
V = Ms{nlevel} - x.gamma(nlevel)*Ws{nlevel};
Sigma = X'*V; Mu = Sigma*Delta; Sigma = Sigma*X;
Lo = chol(Sigma, 'lower'); Mu = Lo\Mu;
Mu = Mu + normrnd(0, sqrt(x.tau2(nlevel)), [length(Mu),1]);
x.beta = Lo'\Mu;
err{nlevel} = Delta - X*x.beta;
end

%%
function [x] = updateGamma(x, err, j)
% sample spatial dependence parameter for CAR model
global Ws loglike0s gammass gapss
tmp = err{j}'*Ws{j}*err{j}/(2*x.tau2(j));
loglike = loglike0s{j} + gammass{j}*tmp; MaxLogLike = max(loglike);
P = exp(loglike-MaxLogLike)/sum(exp(loglike-MaxLogLike));
u2 = rand(1);
cump = cumsum([0 P(1:(end-1))]); i0 = sum(u2 > cump);
x.gamma(j) = gammass{j}(1);
if(i0>1)
    x.gamma(j) = gammass{j}(i0-1) + gapss{j}/P(i0)*(u2-cump(i0));
end
end

%%
function [Es] = getDIC(x, err, Y)
% compute DIC4 for mixed-effects model
global X Ws Ms dicLen nlevel p vec a_tau2 invb_tau2 Ns Inds
log1 = 0;
for j = 1:nlevel
    if vec(j)~=0
        V = Ms{j}-x.gamma(j)*Ws{j}; L = chol(V,'lower'); eps = L'*err{j};
        log1 =log1 - 0.5/x.tau2(j)*sum(eps.^2) + sum(log(diag(L))) - 0.5*Ns(j)*log(2*pi*x.tau2(j));
    end
end
betabar = zeros(p,1); tau2bar = zeros(1,length(vec)); gammabar = zeros(1,length(vec));
for i = 1:dicLen
    [x, err, V] = updateBeta(x, err, Y); betabar = betabar + x.beta;
    for j = 1:nlevel
        if vec(j)~=0
            V1 = V;
            if j < nlevel
                V1 = Ms{j} - x.gamma(j)*Ws{j};
            end
            x.tau2(j) = 1/gamrnd(a_tau2+0.5*Ns(j), 1/(invb_tau2+0.5*err{j}'*V1*err{j}));
            if vec(j) == 1
                x = updateGamma(x, err, j);
            end
        end
    end
    tau2bar = tau2bar + x.tau2; gammabar = gammabar + x.gamma;
end
betabar = betabar/dicLen; tau2bar = tau2bar/dicLen; gammabar = gammabar/dicLen;
log2 = 0; err{nlevel} = Y - X*betabar;
for j = 1:nlevel
    if vec(j)~=0
        V = Ms{j}-gammabar(j)*Ws{j}; L = chol(V,'lower'); eps = L'*err{j};
        log2 =log2 - 0.5/tau2bar(j)*sum(eps.^2) + sum(log(diag(L))) - 0.5*Ns(j)*log(2*pi*tau2bar(j));
        if nlevel > 1 && j < nlevel
            err{nlevel} = err{nlevel} - err{j}(Inds{nlevel-1,j});
        end
    end
end
Es = [log1, log2];
end
