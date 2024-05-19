function [meavalue ] = GFM_gedtwoIC_norm_pois(i, q, n_seq, p_seq)
% The function is to conduct for different n and p as for normal + poisson
% It is used to validate the effeciency gain by one-step update.
% Created by Wei Liu on 2020/01/30.
% Updated date: 2020-01-30
% Copyright ? 2020 Wei Liu. All rights reserved.
%

if(~exist('q', 'var') || isempty(q))
    q = 5;
end

if(~exist('n_seq', 'var') || isempty(n_seq))
    n_seq = [ 100];
end
if(~exist('p_seq', 'var') || isempty(p_seq))
    p_seq = [100];
end

ex = 'norm_pois';
n_num = length(n_seq); p_num = length(p_seq);
nBH = 2; nIC = 2;
meavalue = zeros(n_num, p_num, nBH, nIC);
for jj = 1:n_num
    fprintf('n= %d\n', n_seq(jj));
    for kk = 1:p_num
        [X, H, B, mu] = gendata(i, n_seq(jj), p_seq(kk), ex, q);
        Bm = [B, mu'];
        p = p_seq(kk);
        group = [ones(1,floor(p/2)), 2*ones(1, p-floor(p/2))];
        type = {'normal', 'identity'; 'poisson', 'log'}; 
        X = single(X);
        eps2 = 1e-4; maxIter=10; output = 0;  
        dropout = 0;
        ICmethod = 'SVD';
        [hH_svd, hB_svd, hmu_svd] = gfm_eval2step2(X, group, type, q, dropout, eps2, maxIter, output,ICmethod);
        ICmethod = 'sep';
        [hH_sep, hB_sep, hmu_sep] = gfm_eval2step2(X, group, type, q, dropout, eps2, maxIter, output,ICmethod);
        try
        hBm_sep = [hB_sep, hmu_sep]; hBm_svd = [hB_svd, hmu_svd];
        meavalue(jj, kk, :,:) = [measurefun(H, hH_sep), measurefun(H, hH_svd) ; ...
            measurefun(Bm, hBm_sep), measurefun(Bm, hBm_svd), ];
        catch
            meavalue(jj, kk, :,:) = inf;
        end
     end
end