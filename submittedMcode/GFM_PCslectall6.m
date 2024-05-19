function [allIC] = GFM_PCslectall6(i, q, n_seq, p_seq, q_set)


if(~exist('q', 'var'))
    q = 6;
end

if(~exist('n_seq', 'var'))
    n_seq = [200, 300, 400];
end
if(~exist('p_seq', 'var'))
    p_seq = [200, 300, 400, 500]; 
end
if(~exist('q_set', 'var'))
    q_set = 1:7;
end
ex = 'npb2';
n_num = length(n_seq); p_num = length(p_seq);
q_num = length(q_set);
allIC = zeros(n_num, p_num);
for jj = 1:n_num
    % jj = 1;
    fprintf('n = %d \n ',  n_seq(jj))
    for kk = 1:p_num
        % kk = 1;
        X = gendata(i, n_seq(jj), p_seq(kk), ex, q);
        p = p_seq(kk);
        type = cell(3,2);
        type{1,1} = 'normal'; type{1,2} = 'identity';
        type{2,1} = 'poisson'; type{2,2} = 'log';
        type{3,1} = 'binomial';  type{3,2} = 'logit';
        group = [ones(1,floor(p/3)), 2*ones(1, floor(2*p/3)-floor(p/3)), 3*ones(1, p-floor(2*p/3))];
        dropout=[3];eps2 = 1e-4; maxIter=10; output = 0;  
        X = single(X);
        IC = zeros(1, q_num);
        try
            for r = 1:q_num
                fprintf('r = %d \n ', q_set(r))
                [hH, hB, hmu] = gfm_eval2step(X, group, type, q_set(r), dropout, eps2, maxIter, output);
                hHm = [ones(n_seq(jj),1), hH]; hBm = [hmu, hB];
                Vr = ICriteria(X, hBm, hHm, q_set(r), group, type);
                IC(r) = sum(Vr);
                fprintf('IC(r) = %f \n', IC(r))
            end
            [~, rid] = min(IC);
            allIC(jj, kk) = q_set(rid);
        catch
            allIC(jj, kk) = inf;
        end
    end
end