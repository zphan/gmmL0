function S = cov_gene_ER

%% 7027 genes 158 samples
%% size 158x7027 log-scaled and normalized

%% obtained from Dobra BMSS package
% http://www.stat.washington.edu/adobra/software/bmss/

%% first prepared by Pittman et al. (2005) 
% http://www.pnas.org/content/101/22/8431.full

%% loading and scaling
data = load('codex.158x7027.dat');  % size n by p, n: sample size; p: dimension                      
% log-scaled and normalized

textdata = textread('description.txt','%s'); 

fprintf('Estrogen receptor 1 data 7027 genes 158 samples\n');

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

save ER  S textdata;     