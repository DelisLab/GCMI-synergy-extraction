function [Threshold,MI]=perm_mi(X,Y,Ym,percentile,iterations)
%%Input
%X=Source variable
%Y=Target variable (Range of values must be [0 Ym-1])
%Ym=Maximum value in Y
%Threshold=Empirical threshold
%iterations= Number of permutations to carry out

%%Functions used here include gcmi_mixture_cd, copnorm and mi_mixture_gd
%%See GCMI toolbox for details
MI=gcmi_mixture_cd(X,Y,Ym);

data_norm=copnorm(X);
MI_rand=[];
for iter=1:iterations
    trg_rand=Y(randperm(length(Y)));
    mi_rand=mi_mixture_gd(data_norm,trg_rand,Ym);
    MI_rand=[MI_rand;mi_rand];
end

MI_rand=reshape(MI_rand,[iterations,1]);
Threshold=prctile(MI_rand,percentile);
end