function [Hstar,Sigma,WP,obj] = myLateFusionMVCminmaxCD(KH,numclass,option)

numker = size(KH,3);
Sigma = ones(numker,1)/numker;
%--------------------------------------------------------------------------------
% Options used in subroutines
%--------------------------------------------------------------------------------
if ~isfield(option,'goldensearch_deltmax')
    option.goldensearch_deltmax=5e-2;
end
if ~isfield(option,'goldensearchmax')
    optiongoldensearchmax=1e-8;
end
if ~isfield(option,'firstbasevariable')
    option.firstbasevariable='first';
end
%%-------------------------
%% Initializing Missing Elememts of KHs
%%-------------------------
[HP,WP] = myInitialiHp(KH,numclass);

flag =1;
iter = 1;
while flag
    %% MKKM with imputed kernel matrix
    [Hstar,Sigma,obj1] = minmaxLateFusionMVC(HP,WP,Sigma,option);
    obj(iter) = obj1(end);
    %% update kernel matrix with Hstar and Sigma
    for p = 1 : numker
        Wmatrix = HP(:,:,p)'*Hstar;
        [UH,~,VH] = svd(Wmatrix,'econ');
        WP(:,:,p) = UH*VH';
    end
    %     %% Calculate Objective
    %     Kmatrix = sumKbeta(KH,(Sigma.*Sigma));
    %     obj(iter) = trace(Hstar'*Kmatrix*Hstar);
    if iter >=2 &&  (abs(obj(iter)- obj(iter-1))/obj(iter)<1e-4)
        flag =0;
    else
        iter = iter + 1;
    end
end