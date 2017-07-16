function [final_Rscore]=DMTHNDM(D_M,M_T,lamda)
% The best predictive results achieved by DMTHNDM when lamda=0.7
% Adjacency matrix D_M: known experimental verification of  diseases-miRNAs associations
% matrix D_M is 250*209:250 rows represent 250 diseases and 209 columns represent 209 miRNAs
% Adjacency matrix M_T: known experimental verification of  miRNAs-targets associations
% matrix M_L is 209*11001: 209 rows represent 209 miRNAs and 11001 columns represent 11001 targets

% nd:the number of diseases in D_M
% nm:the number of miRNAs in D_M
% nt:the number of targets in M_T
[nd,nm] = size(D_M);
[nm,nt] = size(M_T);

% Initialize matrix Rsocre1,Rsocre2 and weight matrix W
Rscore1=zeros(nd,nm);
Rsocre2=zeros(nm,nt);
W=[];

%calculate the corresponding weight matrix W_d in D_M
for i=1:nd
        q=bsxfun(@rdivide,repmat(D_M(i,:),nd,1).*D_M,sum(D_M));
        W_d(i,1:nd)=1./sum(D_M,2).*sum(q,2);
end
%calculate the corresponding weight matrix W_m in M_T
for i=1:nm
        q=bsxfun(@rdivide,repmat(M_T(i,:),nm,1).*M_T,sum(M_T));
        W_m(i,1:nm)=(1./sum(M_T,2)).*sum(q,2);
end
%calculate the level of consistency between the contribution of resource moved in both directions 
W=W_d';
W=W./(repmat(sum(W),nd,1));
W_d = (W_d+W);
%calculate the level of consistency between the contribution of resource moved in both directions 
W=W_m';
W=W./(repmat(sum(W),nm,1));
W_m = (W_m+W);

%obtain the first level of resource score about disease_miRNA_associations
Rscore1= W_d*D_M;
%obtain the second level of resource score about disease_miRNA_association by using miRNA-related targets as collaborative label 
Rsocre2=D_M*W_m;

%%%%%%%%%%%%%%%%%%construct Disease_MiRNA_Target Heterogeneous network%%%%%%%%%%%%%%%%%%
%calculate the final relevance resource score(Rscore) to infer potential disease-related miRNAs
final_Rscore=plus(lamda*(Rscore1),(1-lamda)*(Rsocre2));
%save final_Rscore;
pre_label_score = final_Rscore(:);

%load the matrix for the correspondence columns of disease names: 
disease_250=importdata('disease_250.txt','\n');
%load the matrix for the correspondence rows of miRNA names: 
miRNA_209=importdata('miRNA_209.txt');



%obtain corresponding ranks(descend) of predictions computed by DMTHNDM
[rank_score,previous_site]=sort(pre_label_score,'descend');


count=1;  
%match the corresponding site of prediction result in disease_250 and miRNA_209
for i=1:length(pre_label_score)
    [disease,miRNA]=ind2sub(size(final_Rscore),previous_site(i));
    rank_corr_position(i,:)=[disease,miRNA]; 

%remove the known interaction in prediction results, obtain the predicted potential disease-related miRNAs 
    if D_M(disease,miRNA) ~=1
         rank_ans_site(count,:)=[disease_250(disease),miRNA_209(miRNA)];
		 %the corresponding  relevance score of candidate miRNA
		 rank_candidate_score(count,3)=rank_score(i,:);
         count=count+1;
    end
end

%ranking predicted results of potential disease-related miRNAs in descending order
 predicted_results=rank_ans_site;
 %all potential candidate pairs are written in corresponding excel table(final prediction candidate pairs.xls)
 xlswrite('final prediction candidate pairs.xls',predicted_results(:,1),'A1:A48607');
 xlswrite('final prediction candidate pairs.xls',predicted_results(:,2),'B1:B48607');
 xlswrite('final prediction candidate pairs.xls',rank_candidate_score(:,3),'C1:C48607');
end


