function [SP_COND]=SpaceTime_Cond(X,Y,Trials)
%%Space-Time computation conditioned on a discrete task variable
%%Input:
%X=Input matrix (Timepoints (Rows) x Muscles (Columns) x Trials)...
    %Each trial must have same number of timepoints
    
%Y= Task variable (ranging in values from [0 Ym-1]
%Ym= Maximum value in Y
%Trials= Number of trials

%%Output:
%SP_COND: a matrix of CMI values with the columns [1 length(X)] the CMI
%between muscles at the same timepoint across trials making up the diagonal of the adjacency matrix and the
%remaining columns the CMI between muscles at timepoint A and B


%%Note
%Output may have negative values due to bias correction. The output should
%be thresholded with these negative values included or if no thresholding
%required can be set to zero

%%For information on the gccmi_ccd function used here, please see:
%https://github.com/robince/gcmi

len=length(X)/Trials;
x=mat2tiles(X,[len,size(X,2)]);


combos_time=nchoosek(1:length(x{1}),2);
combos=nchoosek(1:size(x{1},2),2);
x=cat(3,x{:});

MIs=[];
for i=1:length(combos_time)
    for ii=1:length(combos)
        vars1=x(combos_time(i,1),combos(ii,1),:);
        vars2=x(combos_time(i,2),combos(ii,2),:);
        [CMI I]=gccmi_ccd(vars1(:),vars2(:),Y,max(Y)+1);
        MIs=[MIs;CMI];
    end
end
sp=reshape(MIs,[length(combos),length(combos_time)]);


MIs_linear=[];
for i=1:50
    for ii=1:length(combos)
        vars1=x(i,combos(ii,1),:);
        vars2=x(i,combos(ii,2),:);
        [CMI I]=gccmi_ccd(vars1(:),vars2(:),Y,max(Y)+1);
        MIs_linear=[MIs_linear;CMI];
    end
end
sp2=reshape(MIs_linear,[length(combos),length(x(:,:,1))]);


SP_COND=[sp2,sp];