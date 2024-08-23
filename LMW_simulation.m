clear all; close all; clc;format compact;format shortG;
name = 'ryusiwon';
%name = 'hannu';
cd(['C:\Users\',name,'\Dropbox\1_Study\2. TA\2019-1 Topics\TA session material\simulation\LMW']);

global R ngrid alpha
alpha   = 0.05;
ngrid   = 40;


% Set design [mean of x1, sd of x1, mean of x2, sd of x2]
% H0 : x1 weakly sd x2 <=> F1 <= F2
% BDstat = sqrt(-)*sup(F1-F2)
design1 = [0 1 0 1]; %FSD, SSD
design2 = [-0.1 1 0 1] ; %FSD
design3 = [0.5 1 0 1]; % SSD
design4 = [0 0.5 0 1]; %None
design = [design1;design2;design3;design4];  

b1 = 10; b2 = 10;
%LMW test and choose subsample size
m = 1;
SDorder = 1;
sample1 = design(m,1)+design(m,2)*randn(100,1,1,1);
sample2 = design(m,3)+design(m,4)*randn(80,1,1,1);
grid = linspace(min(min([sample1;sample2])), max(max([sample1;sample2])),ngrid);
[lmw,pval] = lmwtest(sample1,sample2,grid,SDorder,b1,b2)
[lmw,pval] = lmwtest_select(sample1,sample2,grid,SDorder)


 
%%
% Simulation
tic; 
R = 1000; % Number of simulation
lmw_R = zeros(R,4,2); pval_R = zeros(R,4,2); 
for SDorder = [1]
    for m = [1]
        b1 = 10; b2 = 10; % Set subsample size
        sample1 = design(m,1)+design(m,2)*randn(100,1,1,R);
        sample2 = design(m,3)+design(m,4)*randn(80,1,1,R);
        grid = linspace(min(min([sample1;sample2])), max(max([sample1;sample2])),ngrid);
        [lmw_R(:,m,SDorder),pval_R(:,m,SDorder)] = lmwtest(sample1,sample2,grid,SDorder,b1,b2);
    end
end
if R > 1
    rej_prob = mean(pval_R<alpha,1);
    disp(['Rejection probabilty = '])
    disp(reshape(rej_prob,4,2,[]))
else
    disp('lmw = ') 
    disp(lmw_R)
    disp('pvalues = ')
    disp(pval_R)
end
save('result_LMW_new')
toc;
    

%%%%%% Function part %%%%%%

% Step 1: define function of lmwtest given subsample size b : lmwtest
function [lmw,pval] = lmwtest(sample1,sample2,grid,SDorder,b1,b2)
global R
n1 = size(sample1,1);
n2 = size(sample2,1);
n = n1*n2/(n1+n2);
lambda = n2/(n1+n2);
b = lambda * b1;
nsub = min(n1-b1+1, n2-b2+1);

operator = @(X,z) (X<=z).*(z-X).^(SDorder-1)/factorial(SDorder-1); %integral operator
ecdf    = @(X,z) mean(bsxfun(operator,X,z)); % n1 by R by grid_num
stat = @(n,sample1,sample2,grid) sqrt (n)*max (ecdf(sample1, grid)-ecdf(sample2, grid),[],2);
lmw = stat(n,sample1,sample2,grid);

% Block Subsampling
subindex1 = repmat(reshape(bsxfun(@plus, transpose(1:1:b1),(0:1:nsub-1)),b1,1,[]),[1,1,1,R]);
subindex2 = repmat(reshape(bsxfun(@plus, transpose(1:1:b2),(0:1:nsub-1)),b2,1,[]),[1,1,1,R]);

% Subsampling statistics
substat = @(grid) sqrt(b) * max (ecdf(sample1(subindex1), grid)-ecdf(sample2(subindex2), grid),[],2);
%substat = @(grid) sqrt(b) * max ((ecdf(sample1(subindex1), grid)-ecdf(sample1, grid)) -(ecdf(sample2(subindex2), grid)-ecdf(sample2, grid)),[],2);
pval = mean(substat(grid) > lmw);
end

% Step 2: define function of finding desired subsample size by repeating
% lmwtest in different subsample size b, and plot p-values over b
function [lmw,pvald] = lmwtest_select(sample1,sample2,grid,SDorder)
n1 = size(sample1,1);
n2 = size(sample2,1);
n = n1*n2/(n1+n2);
subsize_list = round(linspace(5,n,20));
pval = zeros(size(subsize_list,2),1);
for i = 1:size(subsize_list,2) %subsample size
    [lmw,pval_b] = lmwtest(sample1,sample2,grid,SDorder,subsize_list(i),subsize_list(i));
    pval(i,:) = pval_b;
end
pvald = [subsize_list',pval]; %p-values with subsample size
scatter(subsize_list,pval,'filled','black');
end