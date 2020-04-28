function S = cov_gene_Lymph

%% 4514 genes 148 samples (100 low-risk, 48 high-risk)
% size 148x4515, last colum contains binary class info

%% obtained from Dobra BMSS package
% http://www.stat.washington.edu/adobra/software/bmss/

%%  analyzed by Hans et al.
% http://ftp.isds.duke.edu/WorkingPapers/05-10.html

%% first prepared by Pittman et al. (2005) 
% http://www.pnas.org/content/101/22/8431.full

%% loading and scaling
data = load('Lymph.dat');   % size n by p+1, n: sample size; p: dimension
                                              % last column is binary class info                             

fprintf('Lymph data 4514 genes 148 samples\n');
for i=1:(size(data,2)-1)
    a = data(:,i);
    data(:,i) = (a-mean(a))/std(a);
end;

%% dimension reduction FDR multiple testing

fprintf('Dimension reduction using FDR multiple testing\n');
mykendall = zeros(1,size(data,2)-1);
pval = zeros(1,size(data,2)-1);
t = 0;
for i=1:length(mykendall)    
    if (mod(t,100)==0); 
        fprintf('\n working [%d]',i); 
    else
        fprintf('.');
    end
    t = t+1;
  [mykendall(i),pval(i)] =...
corr(data(:,i),data(:,size(data,2)),'type','kendall');
end

[~, qvalues] = mafdr(pval);
%[pFDR, qvalues] = mafdr(pval, 'showplot', true);

idx = find(qvalues <= 0.05); % significance level

fprintf('sample size %d,\t dimension %d\n',size(data,1),length(idx));
new_data = data(:,idx); %size n by idx, n: sample size; idx: dimension
S = cov(new_data);
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


data = new_data;
save Lymph.mat S data;     