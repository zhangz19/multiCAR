function [] = mainPro0(ID) %ID = 1:180: 3 levels, 5 chains, 3 models, 4 variables
% spatial CAR model: revisited by Zhen, 2013-11-22
global W M X Y N p 

tot = 6e3;  simu = 0;
ch = str2double(num2str(ID));

nlevel = ceil(ch/60);
ch = ch - (nlevel-1)*60;
nmodel = ceil(ch/20);
ch = ch - (nmodel-1)*20;
nvar = ceil(ch/5);
ch = ch - (nvar-1)*5;
fprintf('level=%d, model=%d, var=%d, chain=%d\n',[nlevel,nmodel,nvar,ch])

switch nmodel
    case 1
        %modelnam = 'CAR';
        nonspatial = 0; simpleReg = 0;
    case 2
        %modelnam = 'weighReg';
        nonspatial = 1; simpleReg = 0;
    case 3
        %modelnam = 'simpleReg';
        nonspatial = 1; simpleReg = 1;
end

load(strcat('fulldat',num2str(nlevel),'.mat'))
Y = Y(:,nvar); X = zscore(X);
X = [ones(size(X,1),1), X];
[N,p] = size(X); M = diag(sum(W,1)); 

if ~nonspatial && ~simpleReg
    if ~exist(strcat('EV',num2str(nlevel),'.mat'),'file')
        invM = inv(M);
        eigs = eig(sqrt(invM)*W*sqrt(invM));
        lgamma = max(1/min(eigs),-1); ugamma = 1/max(eigs);
        gaps = 1e-3; gammas = (lgamma+gaps):gaps:(ugamma-gaps); len = length(gammas);
        loglike0 = zeros(1,len);
        for i = 1:len
            gamma = gammas(i);
            B = (M-gamma*W);
            loglike0(i) = 0.5*sum(log(eig(B)));
        end
        save(strcat('EV',num2str(nlevel),'.mat'),'M','gaps','gammas','loglike0')
    end
    load(strcat('EV',num2str(nlevel),'.mat'))
end

rng(10); %RandStream.setGlobalStream(RandStream('mt19937ar','seed',id0*200));

truePara = [];
if simu
    truePara.beta = (X'*X)\(X'*Y).*(1 + (6-1)/5*normrnd(0,1,[p,1]));
    truePara.tau2 = 1; 
    truePara.gamma = 0.0;
    u = mvnrnd(zeros(N,1), truePara.tau2.*(eye(N) - truePara.gamma.*invM*W)\invM);
    Y = X*truePara.beta + u';
end

meantau2 = .01; vartau2 = 1e2;
a_tau2 = 2+meantau2^2/vartau2;
invb_tau2 = meantau2*(a_tau2-1);

% set initial values
x.beta = (X'*X)\(X'*Y).*(1 + (ch-3)/5*normrnd(0,1,[p,1]));
x.tau2 = sum((diag(M).*(Y-X*x.beta)).^2)/(N-p);
x.gamma = 0 + (ch-3)*.45;
if simu == 1 
    x.beta = truePara.beta;
    x.tau2 = truePara.tau2;
    x.gamma = truePara.gamma;
end

if nonspatial
    x.gamma = 0;
end
if simpleReg
    x.gamma = 0; M = eye(N);
end

% run the chain, store the results 
xall = nan(tot, p+2);
rng(2^nlevel + 241*ch + 127*nvar + 8*nmodel);
u2 = rand(1,tot);

tic
for i = 1:tot
    
    % Update Beta
    V = M - x.gamma*W;
    Sigma = X'*V; Mu = Sigma*Y; Sigma = Sigma*X;
    Lo = chol(Sigma, 'lower'); Mu = Lo\Mu;
    Mu = Mu + normrnd(0, sqrt(x.tau2), [length(Mu),1]);
    x.beta = Lo'\Mu;
    err = Y - X*x.beta;
    
    % Update Tau2
    x.tau2 = 1/gamrnd(a_tau2+0.5*N, 1/(invb_tau2+0.5*err'*V*err));
    
    if ~nonspatial
        % update Gamma
        tmp = err'*W*err./(2*x.tau2);
        loglike = loglike0 + gammas*tmp; MaxLogLike = max(loglike);
        P = exp(loglike-MaxLogLike)/sum(exp(loglike-MaxLogLike));
        cump = cumsum([0 P(1:(end-1))]); i0 = sum(u2(i) > cump);
        gammanew = gammas(1);
        if(i0>1)
            gammanew = gammas(i0-1) + gaps/P(i0)*(u2(i)-cump(i0));
        end
        x.gamma = gammanew;
    end
    xall(i,:) = [x.beta', x.tau2, x.gamma];
    %     fprintf('%6d', i)
    %     if(~mod(i,20))
    %         fprintf('\n')
    %     end
end
toc
runtime = toc/3600;
% display(runtime) % elapsed time in hour
save(strcat('out',num2str(nlevel),num2str(nmodel),num2str(nvar),num2str(ch),'.mat'),'xall','runtime','truePara')
