function [] = sumPro()
load('modind.mat'); nm = size(modind,1); nvar = 4;
dicMat = zeros(nm, nvar);
for id = 1:nm
    load(strcat('paras',num2str(id),'.mat'))
    dicMat(id,:) = dics;
end
save('dicMat.mat', 'dicMat','modind')
% analyze dicMat
load('dicMat.mat');
dicMat1 = [dicMat(1:3,:), modind(1:3,:)];
dicMat2 = [dicMat(3+(1:12),:), modind(3+(1:12),:)];
dicMat3 = [dicMat(15+(1:48),:), modind(15+(1:48),:)];
optmat = zeros(3,4);
inds ={1:3, 3+(1:12), 15+(1:48)};
for nlevel = 1:3
    for yind = 1:4
        vec = dicMat(inds{nlevel}, yind);
        optmat(nlevel, yind) = find(vec == min(vec)) + inds{nlevel}(1) - 1;
    end
end
load('paras1.mat'); npara = size(matparas{1},2);
optpara = zeros(npara, 3*16); yhats = cell(1,3);
for nlevel = 1:3
    yhats{nlevel} = zeros(Ns(nlevel),4);
    for yind = 1:4
        load(strcat('paras',num2str(optmat(nlevel, yind)),'.mat'))
        optpara(:, (nlevel-1)*16+(yind-1)*4 + 1) = mean(matparas{yind},1)';
        optpara(:, (nlevel-1)*16+(yind-1)*4 + 2) = std(matparas{yind},[],1)';
        optpara(:, (nlevel-1)*16+(yind-1)*4 + (3:4)) = quantile(matparas{yind},[0.025,0.975],1)';
        yhats{nlevel}(:,yind) = Yhat(:,yind);
    end
end
save('final.mat','yhats','optpara','optmat')
end