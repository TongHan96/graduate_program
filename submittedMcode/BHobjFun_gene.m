function obj = BHobjFun_gene(H, B, X, omega)
% normal + poisson + binary.
[n, p] = size(X);
eps1 = 1e-20;
% g1 = 1:p;
Bh = H * B';
clear B H;
% C = num2cell(X, 1);          
%p_ntrails = cellfun(@(x) length(unique(x)), C); 
% clear C;
% N = repmat(p_ntrails, n, 1);
N = 3;
% me3 = 1 ./(1+exp(-Bh(:,g1)));
me3 = 1 ./(1+exp(-Bh));
Q = -log(binopdf(X,N,me3) +eps1);
clear me3 N;

obj = 1/n*omega*sum(sum(Q));
