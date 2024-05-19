function [meavalue] = GFM_gedatest_bino(i, q, n_seq, p_seq)

if(~exist('q', 'var') || isempty(q))
    q = 4;
end

if(~exist('n_seq', 'var') || isempty(n_seq))
    n_seq = [100, 150];
end
if(~exist('p_seq', 'var') || isempty(p_seq))
    p_seq = [100,150];
end

ex = 'bino';
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 2;
meavalue = zeros(n_num, p_num, n_mea);
for jj = 1:n_num
    for kk = 1:p_num
        [X, H, B, mu] = gendata(i, n_seq(jj), p_seq(kk), ex, q);
        group = [ones(1,p_seq(kk))];
        type = cell(1,2);
        type{1,1} = 'binomial';  
        type{1,2} = 'logit'; 
        eps2 = 1e-4; maxIter=10; output = 0; dropout = 0;
        X = single(X);
        
        % estimate B,mu and H matrix.
        [hH, hB, hmu] = gfm_eval2step(X, group, type, q, dropout, eps2, maxIter, output);
        meavalue(jj, kk, 1) = measurefun(hH, H);
        meavalue(jj, kk, 2) = measurefun([hB, hmu], [B, mu']);
    end
end