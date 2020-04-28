function S = cov_gene_Arabidopsis
%% 39+795 genes,  118 samples 
%% size 118x834, 

%% prepared by Wille et al. (2005) 
% http://genomebiology.com/2004/5/11/R92#IDAXXBEB
% Additional data files
% Additional data file 1. The gene expression values of the 40 isoprenoid genes
% Additional data file 2. The gene expression values of the 795 genes from other pathways

%% used by H. Li, and J. Gui (2006): 
% Gradient directed regularization for sparse Gaussian concentration graphs, 
% with applications to inference of genetic networks. Biostatistics, 7: 302-317. 

%% loading and scaling
load Arabidopsis_data.mat data795 textdata795; 

% data = data795; % n by p, sample size by dimension 
% textdata = textdata795;

load Arabidopsis_data.mat data39 textdata39;

data = [data39, data795];
textdata = [textdata39, textdata795];

data = log10(data); % all data are positive

fprintf('Arabidopsis data 834 genes 118 samples\n');

for i=1:size(data,2)
    a = data(:,i);
    data(:,i) = (a-mean(a))/std(a);
end;

%% dimension reduction
 
S = cov(data);
S = 0.5*(S+S');
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


save Arabidopsis S data textdata;
