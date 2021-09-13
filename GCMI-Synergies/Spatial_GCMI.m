function SMI=Spatial_GCMI(X,Trials)
%%Spatial synergy computation

%%Input
    %X= Input matrix of shape [Timepoints (rows) x Muscles (columns) x
        %Trials]
        %Note: All trials in the matrix must be of equal length
    %Trials= Number of trials in the matrix
%%Output
    %SMI= GCMI values [Rows: Trials x Columns: Unique pairings]

for i=1:Trials
    Mi_Spatial.(sprintf('Trial_%d',i))=[];
end

trial_len=length(X)/Trials;
x=mat2tiles(X,[trial_len,size(X,2)]);
combos=nchoosek(1:size(x{1},2),2);
x=cat(3,x{:});

combos=nchoosek(1:size(X,2),2);
for i=1:Trials
    trial=reshape(x(:,:,i),[trial_len,size(X,2)]);
    for ii=1:length(combos)
        I=gcmi_cc(trial(:,combos(ii,1)), trial(:,combos(ii,2)));
        Mi_Spatial.(sprintf('Trial_%d',i))=[Mi_Spatial.(sprintf('Trial_%d',i));I];
    end
    Mi_Spatial.(sprintf('Trial_%d',i))=Mi_Spatial.(sprintf('Trial_%d',i))';
end

SMI=cell2mat(struct2cell(Mi_Spatial));