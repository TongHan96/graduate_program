function [meavluae] = GFM_gedeffi_norm_pois(i, q, n_seq, p_seq)
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
    n_seq = [50];
end
if(~exist('p_seq', 'var') || isempty(p_seq))
    p_seq = [50];
end

ex = 'norm_pois';
n_num = length(n_seq); p_num = length(p_seq);
nBH = 2; n_init= 3; n_os = 2;
meavluae = zeros(n_num, p_num, nBH, n_init, n_os);
for jj = 1:n_num
    fprintf('n= %d\n', n_seq(jj));
    for kk = 1:p_num
        [X, H, B, mu] = gendata(i, n_seq(jj), p_seq(kk), ex, q);
        p = p_seq(kk);
        group = [ones(1,floor(p/2)), 2*ones(1, p-floor(p/2))];
        type = {'normal', 'identity'; 'poisson', 'log'};
        X = single(X);
        eps2 = 1e-4; maxIter=10; output = 0;
        % only use normal to update H
        dropout = 2;
        [hH_norm, hB_norm, hmu_norm, history] = gfm_eval_intercept_init(X, group, type, q, dropout, eps2, maxIter,output);
        % only use poisson to update H
        dropout = 1;
        [hH_pois, hB_pois, hmu_pois, history] = gfm_eval_intercept_init(X, group, type, q, dropout, eps2, maxIter,output);
        
        % use two types of variables to update H
        dropout = 0;
        [hH_a, hB_a, hmu_a, history] = gfm_eval_intercept_init(X, group, type, q, dropout, eps2, maxIter,output);
        Bm = [B mu']; hBm_norm = [hB_norm, hmu_norm]; hBm_pois = [hB_pois, hmu_pois];
        hBm_a = [hB_a, hmu_a];
        meavluae(jj,kk, :,:, 1) = [measurefun(H, hH_norm),measurefun(H, hH_pois), measurefun(H, hH_a); ...
            measurefun(Bm, hBm_norm), measurefun(Bm, hBm_pois), measurefun(Bm, hBm_a)];
        try
        hH_init = hH_norm; hB_init= hB_norm; hmu_init= hmu_norm;
        [hH_final_norm, hB_final_norm, hmu_final_norm, objvalue] = gfm_eval_intercept_osfinal2(X, hH_init, hB_init,hmu_init, group, type);
        hH_init = hH_pois; hB_init= hB_pois; hmu_init= hmu_pois;
        [hH_final_pois, hB_final_pois, hmu_final_pois, objvalue] = gfm_eval_intercept_osfinal2(X, hH_init, hB_init,hmu_init, group, type);
        
        hH_init = hH_a; hB_init= hB_a; hmu_init= hmu_a;
        [hH_final_a, hB_final_a, hmu_final_a, objvalue] = gfm_eval_intercept_osfinal2(X, hH_init, hB_init,hmu_init, group, type);
        hBm_final_norm = [hB_final_norm, hmu_final_norm]; 
        hBm_final_pois = [hB_final_pois, hmu_final_pois];
        hBm_final_a = [hB_final_a, hmu_final_a];
        meavluae(jj,kk, :,:, 2) = [measurefun(H, hH_final_norm),measurefun(H, hH_final_pois), measurefun(H, hH_final_a); ...
            measurefun(Bm, hBm_final_norm), measurefun(Bm, hBm_final_pois), measurefun(Bm, hBm_final_a)];
        catch
            meavluae(jj,kk, :,:, 2) = inf;
        end
     end
end