function [allIC] = GFM_PCslectall(i, ex, q, n_seq, p_seq, q_set)
if(~exist('ex', 'var'))
    ex = 'homonorm';
end
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
n_num = length(n_seq); p_num = length(p_seq);
q_num = length(q_set);
allIC = zeros(n_num, p_num);
for jj = 1:n_num
    for kk = 1:p_num
        X = gendata(i, n_seq(jj), p_seq(kk), ex, q);
        group = [ones(1,p)];
        type = cell(1,2);
        type{1,1} = 'normal'; type{1,2} = 'identity';
        X = single(X);
        dropout = 0; eps2 = 1e-4; maxIter=10;  output = 1;  
        IC = zeros(1, q_num);
        for r = 1:q_num
            %fprintf('r = %d \n ', r)
            [hH, hB] = nonlinfacest_homonorm(X, H, q_set(r));
            Vr = ICselect_homonorm(X, hB, hH, q_set(r));
            IC(r) = sum(Vr);
        end
        [~, rid] = min(IC);
        allIC(jj, kk) = q_set(rid);
    end
end