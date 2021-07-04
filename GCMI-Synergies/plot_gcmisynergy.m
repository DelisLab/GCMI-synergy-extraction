function [FIG1,FIG2]=plot_gcmisynergy(W,H,Type,Opt_rank,dim)

FIG1=figure;
for i=1:Opt_rank
   subplot(Opt_rank,1,i);
   component=W(:,i)';
   B = tril(ones(dim(2)));
   B = tril(ones(dim(2)-1,dim(2)));
   B(B==1) = component;
   zeros=repelem(0,dim(2));
   B=[zeros;B];
   B=tril(B)+tril(B,-1)';
   imagesc(B);
   colorbar
end

if Type=='ST'
    FIG2=figure;
    for i=1:Opt_rank
        subplot(Opt_rank,1,i);
        component=H(i,dim(1)+1:end);
        v=H(i,1:dim(1));
        B = tril(ones(dim(1)));
        B = tril(ones(dim(1)-1,dim(1)));
        B(B==1) = component;
        zeros=repelem(0,dim(1));
        B=[zeros;B];
        B=tril(B)+tril(B,-1)';
        B = B - diag(diag(B)) + diag(v);
        imagesc(B);
        colorbar
    end
else
    FIG2=figure;
    for i=1:Opt_rank
        subplot(Opt_rank,1,i);
        bar(H(i,:),'k');
    end
end

