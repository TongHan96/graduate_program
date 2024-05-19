function Vr = ICriteria(X, hB, hH, r, group, type, criteria)
warning('off');
% PC or IC creteria to choose number of factors.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.
if(~exist('criteria', 'var'))
    criteria = 'IC';
end

[n, p] = size(X);
% w = min(n^(1/7),p^(1/5));
w = log10(p);
% w = 1;
omega = 1/p;
ind_set = unique(group);
ng = length(ind_set);
gcell = cell(1, ng);
for j = 1:ng
    gcell{j} = find(group==j);
end
c = objfunc(hH, hB, X, omega, gcell, type);
fprintf('p = %d \n', p)
fprintf('omega = %f \n', omega)
if strcmp(criteria, 'IC')
   % fprintf('n = %d \n',n)
   % fprintf('p = %d \n',p)
   % Vr = [log(c) ,  r/min(sqrt(n), sqrt(p))^2*log(min(sqrt(n), sqrt(p))^2)];
   Vr = [c, r/(sqrt(n)/w^3 + sqrt(p)/w^2)^2 *log(min(sqrt(p)/w^2, sqrt(n)/w^3)^2)];
   fprintf('c = %f \n', Vr(1))
   fprintf('rg = %f \n', Vr(2))
elseif strcmp(criteria, 'PC')
   % Vr = [c , r * (n+p)/(n*p)*log(n*p/(n+p))];
   % Vr = [c,  r / min(p/w^4,n/w^6) * log(min(p/w^4,n/w^6))];
   % Vr = [c, r / min(p/w^4,n/w^6) * log(min(p/w^4,n/w^6))];
   % Vr = [c , 0.55 * r / (p/w^4+n/w^6) * log(min(n/w^6, p/w^4))];  
   % Vr = [c, 0.5 * r/(p/w^4) * log(n/w^6)];
   % Vr = [c, r/(p/w^4) * log(p/w^4)];
   fprintf('c = %f \n', Vr(1))
   fprintf('rg = %f \n', Vr(2))
% if strcmp(criteria, 'IC')
%    Vr = [c , r*(n+p)/(n*p)*log(n*p/(n+p))];    
%    % Vr = [log(c) , r*(n+p)/(n*p)*log(min(sqrt(n), sqrt(p))^2)];
% elseif strcmp(criteria, 'PC')
%    % Vr = [c , min(n^(1/7),p^(1/5))^6 * r*(n+p)/(n*p)*log(n*p/(n+p))];
%    Vr = [log(c), min(n^(1/7),p^(1/5))^6 * r*(n+p)/(n*p)*log(min(sqrt(n), sqrt(p))^2)];
end
% Vr = [log(c) , r*(n+p)/(n*p)*log(min(sqrt(n), sqrt(p))^2)];
% Vr = [log(c) , r*(n+p)/(n*p)*log(n*p/(+p))]n;
% Vr = [c , r*0.6*(n+p)/(n*p)*log(n*p/(n+p))];