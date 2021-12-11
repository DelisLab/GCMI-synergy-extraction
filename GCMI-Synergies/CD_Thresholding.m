function [Mat_thresholded, Thresholds, Opt_rank, Q, S]=CD_Thresholding(X,Type,num)
    %%A function to compute the optimal model rank via a generalised
    %%Louvain algorithm and to threshold the input matrix using a modified
    %%percolation analysis
    
    %%Input
        %X: The 2D input matrix [rows: trials, columns: unique pairings]
            %Note: in the case of Space-Time matrix, rows=Spatial
            %components, and columns=Temporal components
        %Type: Spatial='S', Temporal='T', Space-Time='ST'
        %num: For Spatial or Temporal matrices, a scalar value describing
            %the number of Muscles or timepoints analyzed in a trial
            %respectively. For Space-Time, input the following= [No. of Timepoints, No. of Muscles].
    
    %%Output
        %Mat_thresholded: Thresholded version of the input matrix
        %Mat_unthresholded: Original Input matrix
        %Opt_rank: Optimal number of communities to extract using a
            %dimensionality reduction method
        %PNMF: Components extracted using projective-NMF.
        %Thresholds: A vector consisting of specific threshold values for
            %each layer of the multiplex network
        %Q: The normalised Q-statistic for multiplex modularity
        %S: A matrix of equal dimension to the input matrix specifying the
            %community affiliation of each datapoint
        
if contains('S',Type)
    NET={};
    for i=1:size(X,1)
        trial=X(i,:);
        B = tril(ones(num-1,num));
        B(B==1) = trial;
        i=repelem(0,num);
        B=[i;B];
        B=tril(B)+tril(B,-1)';
        NET=cat(2,[NET,B]);
    end
    
    NET_T={};
    Thresholds=[];
    for i=1:length(NET)
        A=NET{i};
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        A(A<0)=0;
        NET_T{i}=A;
        Thresholds=[Thresholds;threshold];
    end
    
    Mat_thresholded=[];
    for i=1:length(NET_T)
        A=NET_T{i};
        mask = tril(true(size(NET_T{1})),-1);
        Mat_thresholded = [Mat_thresholded,A(mask)];
    end
    
    
    N=length(NET_T{1});
    T=length(NET_T);
    ii=[]; jj=[]; vv=[];
    twomu=0;
    for s=1:T
        indx=[1:N]'+(s-1)*N;
        [i,j,v]=find(NET_T{s});
        ii=[ii;indx(i)]; jj=[jj;indx(j)]; vv=[vv;v];
        k=sum(NET_T{s});
        kv=zeros(N*T,1);
        twom=sum(k);
        twomu=twomu+twom;
        kv(indx)=k/twom;
        kcell{s}=kv;
    end
    gamma=1;
    omega=1;
    AA = sparse(ii,jj,vv,N*T,N*T);
    clear ii jj vv
    kvec = full(sum(AA));
    all2all = N*[(-T+1):-1,1:(T-1)];
    AA = AA + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
    twomu=twomu+T*omega*N*(T-1);
    B = @(i) AA(:,i) - gamma*kcell{ceil(i/(N+eps))}*kvec(i);
    [S,Q] = genlouvain(B);
    S = reshape(S,N,T);
    Q=Q/twomu;
    Opt_rank=max(max(S));
    
elseif contains('T',Type)
    NET={};
    for i=1:size(X,1)
        trial=X(i,:);
        B = tril(ones(num-1,num));
        B(B==1) = trial;
        i=repelem(0,num);
        B=[i;B];
        B=tril(B)+tril(B,-1)';
        NET=cat(2,[NET,B]);
    end
    
    NET_T={};
    Thresholds=[];
    for i=1:length(NET)
        A=NET{i};
        [threshold] = threshold_by_giant_component(A);
        A(A<threshold)=0;
        A(A<0)=0;
        NET_T{i}=A;
        Thresholds=[Thresholds;threshold];
    end
    
    Mat_thresholded=[];
    for i=1:length(NET_T)
        A=NET_T{i};
        mask = tril(true(size(NET_T{1})),-1);
        Mat_thresholded = [Mat_thresholded,A(mask)];
    end
    
    N=length(NET_T{1});
    T=length(NET_T);
    ii=[]; jj=[]; vv=[];
    twomu=0;
    for s=1:T
        indx=[1:N]'+(s-1)*N;
        [i,j,v]=find(NET_T{s});
        ii=[ii;indx(i)]; jj=[jj;indx(j)]; vv=[vv;v];
        k=sum(NET_T{s});
        kv=zeros(N*T,1);
        twom=sum(k);
        twomu=twomu+twom;
        kv(indx)=k/twom;
        kcell{s}=kv;
    end
    gamma=1;
    omega=1;
    AA = sparse(ii,jj,vv,N*T,N*T);
    clear ii jj vv
    kvec = full(sum(AA));
    all2all = N*[(-T+1):-1,1:(T-1)];
    AA = AA + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
    twomu=twomu+T*omega*N*(T-1);
    B = @(i) AA(:,i) - gamma*kcell{ceil(i/(N+eps))}*kvec(i);
    [S,Q] = genlouvain(B);
    S = reshape(S,N,T);
    Q=Q/twomu;
    Opt_rank=max(max(S));
    
elseif contains('ST',Type)
    NET={};
    for i=1:size(X,1)
        trial=X(i,num(1)+1:size(X,2));
        v=X(i,1:num(1));
        B = tril(ones(num(1)-1,num(1)));
        B(B==1) = trial;
        i=repelem(0,num(1));
        B=[i;B];
        B=tril(B)+tril(B,-1)';
        M = B - diag(diag(B)) + diag(v);
        NET=cat(2,[NET,M]);
    end
    
    NET_T={};
    Thresholds=[];
    for i=1:length(NET)
        A=NET{i};
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        A(A<0)=0;
        NET_T{i}=A;
        Thresholds=[Thresholds;threshold];
    end
    
    Mat_thresholded=[];
    vs=[];
    for i=1:length(NET_T)
        A=NET_T{i};
        v=diag(A);
        vs=[vs,v];
        mask = tril(true(size(NET_T{1})),-1);
        Mat_thresholded = [Mat_thresholded,A(mask)];
    end
    Mat_thresholded=[vs;Mat_thresholded]';
    
    N=length(NET_T{1});
    T=length(NET_T);
    ii=[]; jj=[]; vv=[];
    twomu=0;
    for s=1:T
        indx=[1:N]'+(s-1)*N;
        [i,j,v]=find(NET_T{s});
        ii=[ii;indx(i)]; jj=[jj;indx(j)]; vv=[vv;v];
        k=sum(NET_T{s});
        kv=zeros(N*T,1);
        twom=sum(k);
        twomu=twomu+twom;
        kv(indx)=k/twom;
        kcell{s}=kv;
    end
    gamma=1;
    omega=1;
    AA = sparse(ii,jj,vv,N*T,N*T);
    clear ii jj vv
    kvec = full(sum(AA));
    all2all = N*[(-T+1):-1,1:(T-1)];
    AA = AA + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
    twomu=twomu+T*omega*N*(T-1);
    B = @(i) AA(:,i) - gamma*kcell{ceil(i/(N+eps))}*kvec(i);
    [S,Q] = iterated_genlouvain(B);
    S = reshape(S,N,T);
    Q=Q/twomu;
    Opt_rank=max(max(S));
end
