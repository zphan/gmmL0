function S = cov_gene_Leukemia
%% 3051 genes, 72 samples (38 training plus 34 testing) 
%% size 72x3052, log10 scaled, 2.0000 for log10(-ve)
%% the last column is binary class info 0: ALL; 1: AML

%% prepared by Yeung et al.(2005) 
% http://expression.washington.edu/publications/kayee/bma/
% Datasets/Leukemia data/pre-processed training set (text file)
% .../pre-processed test set (text file)
% .../class file for the training set (text file)
% .../class file for the test set (text file) 
% A tab-delimited text file with rows as samples and genes as
% columns, and genes sorted in descending order of the BSS/WSS ratio. 
% The data is already log10 scaled from the originial data of Golub
% using 2.0000 for log10(-ve)

%% also available at the homepage of Dobra A.
% http://www.stat.washington.edu/adobra/software/bmss/
% BMSS/Archive/Examples/Leukemia/leukemiatrain.txt
% 

%% first analyzed by Golub et al. (1999) 
% http://www.broadinstitute.org/mpr/publications/projects/Leukemia/Golub_et_al_1999.pdf
% original dataset contains 7129 genes and 38 training samples 
% http://www.broadinstitute.org/mpr/publications/projects/Leukemia/data_set_ALL_AML_train.tsv
% 

%% other users: Fan J. and Lv. J (2007)
% Sure Independence Screening for Ultra-High Dimensional Feature Space

%% loading and scaling
load Leukemia_data.mat data;   % size n by p+1, n: sample size; p: dimension
                            % last column is binary class info                             
textdata = textread('Leukemia_textdata.txt','%s');

fprintf('Leukemia data 3051 genes 72 samples\n');

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

[~, qvalues] = mafdr(pval, 'showplot', true);

idx = find(qvalues <= 0.05); % significance level

fprintf('sample size %d,\t dimension %d\n',size(data,1),length(idx));
new_data = data(:,idx); %size n by idx, n: sample size; idx: dimension
S = cov(new_data);
S = 0.5*(S+S');
new_text = textdata(idx);

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
textdata = new_text;
save Leukemia S data textdata;     