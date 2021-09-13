function TMI=Temporal_GCMI(X,Trials)
%%Temporal synergy computation

%%Input
    %X= Input matrix of shape [Timepoints (rows) x Muscles (columns) x
        %Trials]
        %Note: All trials in the matrix must be of equal length
    %Trials= Number of trials in the matrix
%%Output
    %TMI= GCMI values [Rows: Trials x Columns: Unique pairings]



Mi_Temporal=[];
for i=1:Trials
    Mi_Temporal.(sprintf('Trial_%d',i))=[];
end

trial_len=length(X)/Trials;
x=mat2tiles(X,[trial_len,size(X,2)]);
combos=nchoosek(1:size(x{1},1),2);
x=permute(cat(3,x{:}),[2,1,3]);

for i=1:Trials
    trial=copnorm(x(:,:,i));
    Is=[];
    for s=1:length(combos)
        try
           Is=[Is;mi_gg(trial(:,combos(s,1)), trial(:,combos(s,2)))];
        catch message
              Is=[Is;ent_g(trial(:,combos(s,1)))];
        end
     end
     Mi_Temporal.(sprintf('Trial_%d',i))=[Mi_Temporal.(sprintf('Trial_%d',i));Is'];
end



TMI=cell2mat(struct2cell(Mi_Temporal));

