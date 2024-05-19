%% Arrhythmia dataset
clear;
warning('off');
load Arrhythmia0312.mat
X = final_clean_mat;
[n, p] = size(X);
q = 8;
hmu1 = mean(X);
X_center = X - repmat(hmu1, n, 1);
[hH1, hB1] = factorm(X_center, q);
X(X(:,166)==32,166) = 1;

% Iterative
clear history;
warning('off');
output = 1; eps1 = 1e-4;
eps2 = 1e-6; 
[n, p] = size(X);
g1 = 1:165; % identity, normal
g2 = 166:226; % Bernoulli
X = [X_center(:, g1), X(:,g2)];
hB = 0;
dB = Inf; dH = Inf; dOmega = Inf;dc =Inf;
dOmega = max([dB, dH]);
omega = p^(-1);
H2 = zeros(n, q);

hH = factorm(X, q);
tmpB = zeros(p,q);tmpH = hH; tmpc = 1e7;
maxIter = 100;
k = 1;
%w = 1./ (std(X(:, g1)).^2);
w = ones(1, length(g1));
tic;
while k <= maxIter && dOmega > eps1 && dc > eps2
    
    % given  H^(0), update B^(1)
    B1 = [];
    for j = g1 % mj(u) = u
        b1 = glmfit(hH, X(:,j),'normal', 'constant', 'off');
        B1 = [B1, b1];
    end
    B3 = [];
    for j =g2
        b3 = glmfit(hH, X(:,j), 'binomial', 'link', 'logit', 'constant', 'off');
        B3 = [B3, b3];
    end
    hB = [B1, B3]';
    % ensure indentifiability.
    [B0tmp, ~] = qr(hB, 0);
    B0= B0tmp * diag(sort(sqrt(eig(hB'*hB)), 'descend'));
    sB = sign(B0(1,:));
    hB = B0.*repmat(sB,p,1); % ensure B first nonzero is positive
    %hB(1:4,:), B(1:4,:)
    dB = norm(hB - tmpB, 'fro')/norm(hB, 'fro');
    tmpB = hB;
    fprintf('-------------------------------------------\n')
    fprintf('---------- B updation is finished!---------\n')
    
    % given B^(1), update H^(1)
    H1 = [];
    H3 = [];
    for i = 1:n % mj(u) = u
        h1 = glmfit(hB(g1,:), X(i,g1)','normal', 'constant', 'off', 'weights', w);
        H1 = [H1, h1];
        h3 = glmfit(hB(g2,:), X(i,g2)','binomial', 'link', 'logit', 'constant', 'off');
        H3 = [H3, h3];
    end
    %H4 = (H1 + H2+H3)'/3;
    H4 = (H1)';
    % H4 = H1';
    H5 = onestepupdateHfunreal1(X, hB, H4, g1,g2, w);
    hH0 = H5;
    % H1(:,1:5)', H2(:,1:5)', H3(1:5,:), H(1:5,:), hH(1:5,:)
    [H0, ~] = qr(hH0, 0);
    hH1 = H0 * sqrt(n);
    sH0 = sign(hH0(1,:)).* sign(hH1(1,:));
    hH = hH1.* repmat(sH0,n,1);
    dH = norm(hH-tmpH, 'fro')/norm(hH, 'fro');
    tmpH = hH;
    fprintf('-------------------------------------------\n')
    fprintf('---------- H updation is finished!---------\n')
    fprintf('-------------------------------------------\n')
    dOmega = max([dB, dH]);
    c = BHobjFunreal1(hH, hB, X, omega, g1,g2, w);
    dc = abs(c - tmpc)/abs(tmpc);
    tmpc = c;
    if output
        fprintf('Iter %d \n', k);
        fprintf('dB= %4f, dH= %4f,dc= %4f, c=%4f \n', dB, dH,dc, c);
    end
    history.dB(k) = dB; history.dH(k) = dH; history.dc(k)=dc; history.c(k)=c;history.K = k;
    k = k+1;
end
plot(history.c)
xlabel('iteration')
ylabel('obj. func. value')

save arrythmia_dat20200313 hH1 hB1 hH hB

%% Caculate the Hessian matrix
type = {'normal', 'identity'; 'binomial', 'logit'};
gcell = {g1, g2};

[n, p] = size(X);
Bm = hB; Hm = hH;
q = size(Hm, 2);
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
