% Created by Wei Liu in 2020-02-23
% Updated date: 2020-3-24
clear;
warning('off');
%% Select some SNPs to GFM application based marginal screening.
s = tic;
data=importdata('NFBC/geneDat.txt');
stime = toc(s);stime
X = single(data); clear data
[n, p] = size(X);
phdat = load('NFBC/NFBC_filter_mph10.fam');
covriates = phdat(:,[6:8, 10:12, 14:15]);
bmi = single(phdat(1:n,13)); % extract the BMI variable.
pm_vec = sqrt(sum(X == 0) /n );
W = (X - repmat(pm_vec, n,1))./ repmat(sqrt(2*pm_vec.*(1-pm_vec) * p),n,1);
qf = 35;
hH1 = factorm(W, qf);
Pvec = zeros(1,p);

parfor j = 1:p
s = regstats(bmi, [hH1 covriates X(:,j)],'linear'); 
Pvec(j) = s.tstat.pval(end);
end
[Pvec_sort, order_Pvec] = sort(Pvec);
% Pvec_sort(1:10)
% save MargPvec_new Pvec_sort order_Pvec

%% Select number of factors
load usedSNP190510
% load NFBC/usedSNP190510
size(X)
% select q for GFM
r12 = [0, 0];
q_set = 1:10;
q_num = length(q_set);
allhH = cell(q_num,1);
allhB = cell(q_num,1);
for r = 1:q_num
    fprintf('r = %d \n ', r);
    [hH, hB] = nonlinfacest_gene(X, q_set(r), 1e-4, 5);
    allhH{r} = hH;
    allhB{r} = hB;
end
Vr = zeros(2, q_num);
for r = 1:q_num
    fprintf('r = %d \n ', r)
    Vr(:, r) = ICselect_gene(X, allhB{r}, allhH{r}, q_set(r), 1);
end
IC = sum(Vr);
[~, id] = min(IC);
r12(1) = q_set(id)

% select q for LFM.
q_set = 1:10;
q_num = length(q_set);
allhH1 = cell(q_num,1);
allhB1 = cell(q_num,1);
for r = 1:q_num
    fprintf('r = %d \n ', r);
    [hH, hB] = factorm(X, r, 0);
    allhH1{r} = hH;
    allhB1{r} = hB;
end
Vr1 = zeros(2, q_num);
for r = 1:q_num
    fprintf('r = %d \n ', r)
    Vr1(:, r) = ICselect_gene(X, allhB1{r}, allhH1{r}, q_set(r));
end
IC = sum(Vr1);

[~, id] = min(IC);
r12(2) = q_set(id)

q = 4;
s = tic;
[hH1, hB1] = factorm(X, q);
stime = toc(s); stime

%% Add  PCs of SNPs and estimate H and B
% load NFBC/usedSNP190510.mat
load usedSNP190510.mat
X_snp = X;
pc5 = load('NFBC/pc5_NFBC_filter_mph10.txt');
pc5 = load('pc5_NFBC_filter_mph10.txt');
pc5 = pc5(:,3:7); % the first five PCs.
[n, p1] = size(X_snp);

X_inter = single(pc5);
X = single([X_inter X_snp]);
p1 = size(X_inter,2);
%% Estimate the factors and loadings
q = 4;
maxIter = 10;
g1 = 1:p1; % identity, normal
g2 = (p1+1):size(X,2); % Bernoulli
p1 = length(g1);
p2 = length(g2);
clear history;
warning('off');
eps2=1e-6;
output = 1; eps1 = 1e-4; %
[n, p] = size(X);
% g1 = 1:p;
hB = 0;
dB = Inf; dH = Inf; dOmega = Inf;dc =Inf;
dOmega = max([dB, dH]);
omega = p^(-1/2);
H2 = zeros(n, q);
s = tic;
hH = factorm(X, q);
stime = toc(s); stime
tmpB = zeros(p,q);tmpH = hH; tmpc = 1e7;

w = ones(1, length(g1));
k = 1;
fprintf('start algorithm loop! \n')
% startmatlabpool;
while k <= maxIter && dOmega > eps1 && dc > eps2
    
    % given  H^(0), update B^(1)
   
    B3 = zeros(q, p);
    %startmatlabpool;
    parfor j = g1 % mj(u) = u
        B3(:,j) = glmfit(hH, X(:,j),'normal', 'constant', 'off');
    end
    parfor j = g2
        ntrail_j = length(unique(X(:,j)))-1;
        B3(:,j) = glmfit(hH, [X(:,j), ntrail_j*ones(n,1)], 'binomial', 'link', 'logit', 'constant', 'off');
    end
    hB = B3';
    % ensure indentifiability.
    [B0tmp, ~] = qr(hB, 0);
    B0= B0tmp * diag(sort(sqrt(eig(hB'*hB)), 'descend'));
    sB = sign(B0(1,:));
    hB = B0.*repmat(sB,p,1); % ensure B first nonzero is positive
    clear sB B0 B0tmp
    dB = norm(hB - tmpB, 'fro')/norm(hB, 'fro');
    tmpB = hB;
    
    fprintf('-------------------------------------------\n')
    fprintf('---------- B updation is finished!---------\n')
    
    % given B^(1), update H^(1)
    H3 = zeros(q,n);
    H1 = H3;
    parfor i = 1:n % 
        H1(:,i) = glmfit(hB(g1,:), X(i,g1)','normal', 'constant', 'off', 'weights', w);
        ntrail_i = length(unique(X(i,:)));
        H3(:,i) = glmfit(hB(g2,:), [X(i,g2)', ntrail_i*ones(length(g2),1)],'binomial', 'link', 'logit', 'constant', 'off'); 
    end
    % hH0 = (H1 + H3)'/2;
    hH0 = H3';
    s = tic;
    % hH0 = onestepupdateHfunreal1(X, hB, H4,g1,g2, w);
    stime = toc(s); stime
    
    [H0, ~] = qr(hH0, 0);
    hH1 = H0 * sqrt(n);
    sH0 = sign(hH0(1,:)).* sign(hH1(1,:));
    hH = hH1.* repmat(sH0,n,1);
    dH = norm(hH-tmpH, 'fro')/norm(hH, 'fro');
    clear hH1 H0 sH0
    tmpH = hH;
    
    fprintf('-------------------------------------------\n')
    fprintf('---------- H updation is finished!---------\n')
    fprintf('-------------------------------------------\n')
    dOmega = max([dB, dH]);
    % c = BHobjFun_gene(hH, hB, X, omega);
    
    numM = 12;
    M = p1/numM;
    Obj1 = zeros(numM,1);
    parfor j = 1: numM
        fprintf('j = %d \n', j);
        Obj1(j) = objfun_norm(hH, hB(((j-1)*M+1) : (j*M),:), X(:,((j-1)*M+1) : (j*M)), omega);
    end
    % parfor
    numM = 12;
    M = p2/numM;
    Obj2 = zeros(numM,1);
    parfor j = 1: numM
        fprintf('j = %d \n', j);
        Obj2(j) = BHobjFun_gene(hH, hB(p1+((j-1)*M+1) : (p1+j*M),:), X(:,(p1+(j-1)*M+1) : (p1+j*M)), omega);
    end
    c = sum(Obj1) + sum(Obj2);
    dc = abs(c - tmpc)/abs(tmpc);
    tmpc = c;
    if output
        fprintf('Iter %d \n', k);
        fprintf('dB= %4f, dH= %4f,dc= %4f, c=%4f \n', dB, dH,dc, c);
    end
    history.dB(k) = dB; history.dH(k) = dH; history.dc(k)=dc; history.c(k)=c;history.K = k;
    k = k+1;
end



%% Prediction
clear;
load chrom_allused190509.mat% not include pc5 or interaction
% load chrom_allused_pc1snps_190523
% load chrom_allused_pc52Asnps_190524.mat
load NFBC/usedSNP190510
load usedSNP190510
size(X)
q = 4;
s = tic;
[hH1, hB1] = factorm(X, q, 0);
stime = toc(s); stime
[n, p] = size(X);
pm_vec = sqrt(sum(X == 0) /n );
sum(pm_vec==0)
W = (X - repmat(pm_vec, n,1))./ repmat(sqrt(2*pm_vec.*(1-pm_vec) * p),n,1);
W(:,1:5) = X(:, 1:5);
% Ymat = load('NFBC/NFBC_filter_mph10.fam');
Ymat = load('NFBC_filter_mph10.fam');
bmi = Ymat(:,13);
n = length(bmi);
hHmat = hH; hH1mat = hH1(:,1);
measurefun(hH, hH1)
hH14mat = hH1(:,1:q);
addpath('glmnet', '-end');
N = 100;
test_errorMat = zeros(N, 5);
for i = 1:N
    rng(i);
    fprintf('i= %d\n',i);
    tr_ind = randsample(n, round(n/2));
    ts_ind = setdiff(1:n, tr_ind);
    % GFM_REG
    [hbeta, bint, r] = regress(bmi(tr_ind), hHmat(tr_ind,:));
    r_ts = bmi(ts_ind) - hHmat(ts_ind,:)*hbeta;
    nmse_gfm_ts = mean(r_ts.^2) / var(bmi(ts_ind));
    % LFM_REG
    [hbeta1, bint1, r2] = regress(bmi(tr_ind), hH1mat(tr_ind,:));
    r2_ts = bmi(ts_ind) - hH1mat(ts_ind,:)*hbeta1;
    nmse1_lfm_ts = mean(r2_ts.^2) / var(bmi(ts_ind));
    
    [hbeta14, bint1, r2] = regress(bmi(tr_ind), hH14mat(tr_ind,:));
    r2_ts = bmi(ts_ind) - hH14mat(ts_ind,:)*hbeta14;
    nmse14_lfm_ts = mean(r2_ts.^2) / var(bmi(ts_ind));
    % Ridge regression 
    lambda = 1e-3; n1 = length(tr_ind);
    hbeta_ridge = W(tr_ind,:)'*inv(W(tr_ind,:)*W(tr_ind,:)' + lambda*eye(n1))*bmi(tr_ind);
    r_ridge_ts = bmi(ts_ind) - W(ts_ind,:)*hbeta_ridge;
    nmse_rigde_ts = mean(r_ridge_ts.^2) / var(bmi(ts_ind))
    [nmse_gfm_ts nmse1_lfm_ts nmse_rigde_ts]
    % Lasso regression
    % slasso = cvglmnet(double(W(tr_ind,:)), bmi(tr_ind),'gaussian');
    opts.lambda = 0.05;
    options = glmnetSet(opts); 
    hbeta_lasso = glmnet(double(W(tr_ind,:)), bmi(tr_ind),'gaussian', options);
    r_lasso_ts = bmi(ts_ind) - W(ts_ind,:)*hbeta_lasso.beta;
    nmse_lasso_ts = mean(r_lasso_ts.^2) / var(bmi(ts_ind))
    
    test_errorMat(i,:) = [nmse_gfm_ts nmse1_lfm_ts nmse14_lfm_ts nmse_rigde_ts nmse_lasso_ts];
end
test_errorMat(1,:)
mean(test_errorMat)
save pc5-snppred_error20200225 test_errorMat N
load pc5-snppred_error20200225
% Run MRRR
load hH_mrrr0306.mat
size(hH_mrrr)
N2 = 100;
pred_err_mrrr = zeros(N2,1);

for i = 1:N2
    rng(i);
    fprintf('i= %d\n',i);
    tr_ind = randsample(n, round(n/2));
    ts_ind = setdiff(1:n, tr_ind);
    
    % MRRR-based liner regression
    [hbeta, bint, r] = regress(bmi(tr_ind), hH_mrrr(tr_ind,:));
    r_ts = bmi(ts_ind) - hH_mrrr(ts_ind,:)*hbeta;
    nmse_mrrr_ts = mean(r_ts.^2) / var(bmi(ts_ind));
    
    pred_err_mrrr(i) =  nmse_mrrr_ts;
end
pred_err_mrrr(1:2)
mean(pred_err_mrrr)
save pred_error_mrrr20200312 pred_err_mrrr N2


corX = corr(X);
% evaluate the likelihood values
g1 = 1:p1; % identity, normal
g2 = (p1+1):size(X,2); % Bernoulli
l_gfm = objfun_norm(hH, hB(g1,:), X(:,g1), n)+ BHobjFun_gene(hH, hB(g2,:), X(:,g2), n);
l_lfm = objfun_norm(hH1, hB1(g1,:), X(:,g1), n)+ BHobjFun_gene(hH1, hB1(g2,:), X(:,g2), n);
[l_gfm, l_lfm]

%% Evaluate the scaled hessian matrix and the eigvalues
[n, p] = size(X);
Bm = hB; Hm = hH;
q = size(Hm, 2);

type = {'normal', 'identity'};
g1 = 1: size(X,2);
gcell = {g1};
ng = size(type,1);
% mean matrix
mucell = cell(1,ng);
for j = 1:ng
    switch type{j,1}
    case 'normal'
        mucell{j} = Hm * Bm(gcell{j},:)';
    case 'poisson'
        mucell{j} = exp(Hm * Bm(gcell{j},:)');
    case 'binomial'
        Xj = X(:, gcell{j});
        ntrail_j = length(unique(Xj(:,1)))-1;
        mucell{j} =ntrail_j*1./(1 + exp(-Hm * Bm(gcell{j},:)'));
    end
end
% % Hessian matrix or information matrix
B = hB;
d2f = cell(n,1);
for i = 1:n
    Bng = zeros(q,q);
   for j = 1:ng
        switch type{j,1}
    case 'normal'
        %W = diag(1./ (std(X(:,gcell{j})).^2));
        Bng = Bng + B(gcell{j},:)'*B(gcell{j},:);
    case 'poisson'
        Bng = Bng + B(gcell{j},:)'* diag(mucell{j}(i,:)) * B(gcell{j},:);
        %Bng = Bng + (repmat(mucell{j}(i,:), q, 1)'.* B(gcell{j},:))'* B(gcell{j},:);
    case 'binomial'
        %Bng = Bng + (repmat(mucell{j}(i,:), q,1)' .* B(gcell{j},:))' *(repmat(1-mucell{j}(i,:), q,1)' .* B(gcell{j},:));
        Bng = Bng + (B(gcell{j},:))' * diag(mucell{j}(i,:).*(1-mucell{j}(i,:))) * B(gcell{j},:);

        end
   end
   d2f{i} = Bng;
end
d2f

% Hessian for hB
ng = size(type,1);

% mean matrix
d2fB = cell(p,1);
for j = 1:ng
    switch type{j,1}
    case 'normal'
        p11 = length(g1);
        for jj = 1:p11
            mutypej = Hm * Bm(gcell{j},:)'; % n*p1
            scorej = Hm' * (X(:, gcell{j}) -mutypej); % (q+1) * p1
            d2fB{jj} = - Hm'* Hm;
        end       
    case 'binomial'
        jvec = gcell{j};
        Xj = X(:, jvec);
        ntrail_j = length(unique(Xj(:,1)))-1;
        mutypej =ntrail_j*1./(1 + exp(-Hm * Bm(jvec,:)'));
        scorej = Hm' * (Xj -mutypej); % (q+1) * p1
        p1 = length(jvec);
        
        for jl = 1:p1
            fprintf('j = %d\n',j1);
            Hestypejjl = - Hm'* diag(mutypej(:,jl).*(1-mutypej(:,jl)))  * Hm;
            d2fB{p11+jl} =  Hestypejjl;
        end
    end
end
eigHvec = zeros(n,1);
for i = 1:n
    eigHvec(i) = min(eig(d2f{i}/n));
end
min(eigHvec)
eigBvec = zeros(p,1);
for j = 1:p
    eigBvec(j) = min(eig(-d2fB{j}/p));
end
min(eigBvec)