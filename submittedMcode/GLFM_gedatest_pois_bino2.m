function [meavluae] = GLFM_gedatest_pois_bino2(i, q, n_seq, p_seq)

if(~exist('q', 'var') || isempty(q))
    q = 4;
end

if(~exist('n_seq', 'var') || isempty(n_seq))
    n_seq = [100, 200, 300];
end
if(~exist('p_seq', 'var') || isempty(p_seq))
    p_seq = [100, 200, 300, 400]; 
end

ex = 'pois_bino2';
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
meavluae = zeros(n_num, p_num, n_mea);
for jj = 1:n_num
    fprintf('n= %d\n', n_seq(jj));
    for kk = 1:p_num
        [X, H, B, mu] = gendata(i, n_seq(jj), p_seq(kk), ex, q);
        p =  p_seq(kk);
        group = [ones(1,floor(p/2)), 2*ones(1, p-floor(p/2))];
        type = cell(2,2);
        type{1,1} = 'poisson'; type{1,2} = 'log';
        type{2,1} = 'binomial';  type{2,2} = 'logit';
        X = single(X);
        dropout = 2; eps2 = 1e-4; maxIter=10; output = 0;  
        % estimate B,mu and H matrix.
        try 
            [hH1, hB1, hmu1] = factorm(X, q, 0);
            meavluae(jj, kk, 1) = measurefun(hH1, H);
            meavluae(jj, kk, 2) = measurefun([hB1, hmu1], [B, mu']);
        catch 
            meavluae(jj, kk, 1:2) = inf;
        end
        try
            [hH, hB, hmu] = gfm_eval2step(X, group, type, q, dropout, eps2, maxIter, output);
            meavluae(jj, kk, 3) = measurefun(hH, H);
            meavluae(jj, kk, 4) = measurefun([hB, hmu], [B, mu']);
        catch
            meavluae(jj, kk, 3:4) = inf;
        end
    end
end