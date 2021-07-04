function [SP_S,SP_R]=SpaceTime_SR(X,Y,Trials)
%%Space-Time computation of task-relevant synergies
%%Input:
%X=Input matrix (Timepoints (Rows) x Muscles (Columns) x Trials)...
    %Each trial must have same number of timepoints
        
%Y= Discrete task variable (ranging in values from [0 Ym-1]
%Ym= Maximum value in Y
%Trials= Number of trials


%%Output:
%SP_S: a matrix of Synergistic Interaction information values with the columns [1 length(X)] the MI
%between muscles at the same timepoint across trials making up the diagonal of the adjacency matrix and the
%remaining columns the MI between muscles at timepoint A and B

%SP_R: a matrix of Redundant Interaction information values (reverted to positive values)
%with the columns [1 length(X)] the MI between muscles at the same timepoint 
%across trials making up the diagonal of the adjacency matrix and the
%remaining columns the MI between muscles at timepoint A and B


%%For information on the gcmi_cc and gccmi_ccd function used here, please see:
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
        I=gcmi_cc(vars1(:),vars2(:));
        MIs=[MIs;I];
    end
end
sp=reshape(MIs,[length(combos),length(combos_time)]);


MIs_linear=[];
for i=1:50
    for ii=1:length(combos)
        vars1=x(i,combos(ii,1),:);
        vars2=x(i,combos(ii,2),:);
        I=gcmi_cc(vars1(:),vars2(:));
        MIs_linear=[MIs_linear;I];
    end
end
sp2=reshape(MIs_linear,[length(combos),length(x(:,:,1))]);
SP=[sp2,sp];
SP(SP<0)=0;

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
SP_COND(SP_COND<0)=0;

SP_I=SP-SP_COND;
SP_I=SP_I*-1;
SP_I(SP_I<0)=0;
SP_S=SP_I;

SP_I=SP-SP_COND;
SP_I=SP_I*-1;
SP_I(SP_I>0)=0;
SP_R=abs(SP_I);
