function S = cov_gene_hereditarybc
%% 3226 genes, 22 samples 
%% size 22x3226 log2 scaled 

%% prepared by Yeung et al.(2005) 
% http://expression.washington.edu/publications/kayee/bma/
% Datasets/Hereditary breast cancer data/pre-processed data (text file)
% tab-delimited text file, after log2 transform. 

%% first analyzed by Hedenfalk et al. (2001) 
% New England Journal of Medicine, 244:539-548. (2001) 
% http://www.nejm.org/general/content/supplemental/hedenfalk/index.html
% http://research.nhgri.nih.gov/microarray/selected_publications.shtml

%% also analyzed by Storey and Tibshirani (2003)  
% Statistical significance for genomewide studies

%% loading and scaling
load hereditarybc_data.mat data textdata;   
% data size n by p, n: sample size; p: dimension
% textdata = textread('textdata.txt','%s');

fprintf('hereditary data 3226 genes 22 samples\n');

for i=1:(size(data,2))
    a = data(:,i);
    data(:,i) = (a-mean(a))/std(a);
end;


%% dimension reduction
 
 [S,textdata] = dim_red(data,textdata);


%% small pertubation to ensure eig(S)>=0
d = eig(S);
min_eig = min(d);
rank = sum(d>1e-8); 
% for small n large p sample, d usually contains many tiny negative eigenvalues 
if rank < length(d)
    fprintf('Sample covariance matrix is singular. ');
end
fprintf('min eig = %1.4e, rank = %d\n',min_eig,rank);
if min_eig<0
    S = S + max(0,-min_eig)*speye(length(d));
end

save hereditarybc S textdata;     