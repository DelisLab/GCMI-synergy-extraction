function [FIG1,FIG2]=Plotter(W,H,Type,Opt_rank,dim)
         %A function to plot the extracted spatial, temporal or space time
         %synergies.
         
         %%Input
                %W: Synergy weights
                %H: For spatial/temporal synergies these are the activation
                    %coefficients, for space-time these are the temporal
                    %synergies.
                %Type: Spatial='S', Temporal='T', Space-Time='ST'.
                %Opt_rank: The rank of the extracted synergies
                %dim: For spatial synergies, the number of muscles. For
                    %temporal synergies, the number of timepoints per trial.,
                    %For space-time, a vector including both of the above in
                    %the following format: [No. of timepoints, No. of Muscles].
FIG1=figure;
if ~strcmp(Type,'T')
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
else
    for i=1:Opt_rank
        subplot(Opt_rank,1,i);
        component=W(:,i)';
        B = tril(ones(dim(1)));
        B = tril(ones(dim(1)-1,dim(1)));
        B(B==1) = component;
        zeros=repelem(0,dim(1));
        B=[zeros;B];
        B=tril(B)+tril(B,-1)';
        imagesc(B);
        colorbar
    end
end

if strcmp(Type,'ST')
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