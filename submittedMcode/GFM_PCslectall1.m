function [allIC] = GFM_PCslectall1(i, q, n_seq, p_seq, q_set)

if(~exist('q', 'var'))
    q = 1;
end

if(~exist('n_seq', 'var'))
    n_seq = [30,50,100];
end
if(~exist('p_seq', 'var'))
    p_seq = [30,50,100,150];
end
if(~exist('q_set', 'var'))
    q_set = 1:7;
end
ex = 'homonorm2';
n_num = length(n_seq); p_num = length(p_seq);
q_num = length(q_set);
allIC = zeros(n_num, p_num);
for jj = 1:n_num
    fprintf('n = %d \n ',  n_seq(jj))
    for kk = 1:p_num
        X = gendata(i, n_seq(jj), p_seq(kk), ex, q);
        group = [ones(1,p_seq(kk))];
        type = cell(1,2);
        type{1,1} = 'normal'; type{1,2} = 'identity';
        X = single(X);
        dropout = 0; eps2 = 1e-4; maxIter=10;  output = 0;  
        IC = zeros(1, q_num);
        for r = 1:q_num
            fprintf('r = %d \n ', r)
            [hH, hB, hmu] = gfm_eval2step(X, group, type, q_set(r), dropout, eps2, maxIter, output);
            hHm = [ones(n_seq(jj),1), hH]; hBm = [hmu, hB];
            % Vr = ICselect_homonorm(X, hBm, hHm, q_set(r));
            Vr = ICriteria(X, hBm, hHm, q_set(r), group, type);
            IC(r) = sum(Vr);
            fprintf('IC(r) = %f \n', IC(r))
        end
        [~, rid] = min(IC);
        allIC(jj, kk) = q_set(rid);
    end
end