function [W,H]=GNMFDMA(inputdata,params)

%  Wang M.N. et al. Combining Nonnegative Matrix Factorization with Graph Laplacian Regularization    
%   for Predicting Drug-MiRNA Associations based on Multi-layer Networks.
%  2021/12/20

K=params.K; 
a=params.a;
p=params.p;
Y_mat=inputdata.Y_mat=Y_mat; % drug-miRNA association adjacency matrix
SD_mat=inputdata.SD_mat; % drug similarity matrix
SM_mat=inputdata.SM_mat; % miRNA similarity matrix


Y = WKNN( Y_mat, SD_mat, SM_mat, K, r )
M=Y;
graph_D_mat = Graph( SD_mat , p );
graph_M_mat = Graph( SM_mat , p ); 
S_D = SD_mat.* graph_D_mat ;
S_M = SM_mat.* graph_M_mat ;  

clear K a p;
 
k=params.k;
lmada=params.lmada; 
beta=params.beta;
iterate=params.iterate; 
beta > 0  && lamda >0;
fprintf('k=%d  maxiter=%d  beta=%d  lamda=%d\n', k, iterate, beta, lamda ); 

[rows,cols] = size(M);
W=abs(rand(k,rows));        
H=abs(rand(k,cols));

D_d = diag(sum(S_D,2));
D_m = diag(sum(S_M,2));
L_d=D_d-S_D;
L_m=D_m-S_M;


fid = fopen( 'RunResult.txt','wt+');
for step=1:iterate
        U1=W.*((H*M'+lmada*W*S_D)./(H*H'*W+lmada*W*D_d+beta*W));
        V1=H.*((U1*M+lmada*H*S_M)./(U1*U1'*H+lmada*H*D_m+beta*H));
         
        ULU = sum(diag((U1*L_d)*U1'));
        VLV = sum(diag((V1*L_m)*V1'));
        obj = sum(sum((Y-U1'*V1).^2))+beta*(sum(sum(U1.^2)) )+beta*(sum(sum(V1.^2)))+lmada*ULU+lmada*VLV; 
        
        error=max([max(sum((U1-W).^2)),max(sum((V1-H).^2))]);      
        
        fprintf(fid,'%s\n',[sprintf('step = \t'),int2str(step),...
            sprintf('\t obj = \t'),num2str(obj),...
		    sprintf('\t error = \t'), num2str(error)]);
        fprintf('step=%d  obj=%d  error=%d\n',step, obj, error);   
        if error< 10^(-4)
            fprintf('step=%d\n',step);
            break;
        end
        
        W=U1; 
        H=V1;
        
end
fclose(fid);

end

 
