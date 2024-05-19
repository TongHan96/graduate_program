% Created by Han Tong in 2024-05-04
clear;
pctRunOnAll warning('off', 'all');
%% %% %% Section 1: Simulation part
%% %% Part 1: Estimate B and Hm for the six examples.
%%  Example 1: for all j, X_ij is normal with homoskedasticity.
% Repeat N times 
s = tic;
[meavluae] = GLFM_gedatest_homonorm(1);
stime = toc(s);

N = 100;
% simu setting
q = 1; 
n_seq = [200,250,300];p_seq = [100,200,300,400];
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue = zeros(N, n_num, p_num, n_mea);
s = tic;
parfor i = 1:N  
    disp(i)
    allmeavalue(i,:,:,:) = GLFM_gedatest_homonorm(i, q, n_seq, p_seq);
end
size(allmeavalue)
stime = toc(s);

%%  Example 2: for all j, X_ij is normal with heteroskedasticity.
s = tic;
[meavluae] = GLFM_gedatest_heternorm(1);
stime = toc(s);

N = 1000
% simu setting
q = 2; 
q = 2; n_seq = [150, 200,250];p_seq = [50,75,100];
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue2 = zeros(N, n_num, p_num, n_mea);
s = tic;
parfor i = 1:N  
    disp(i)
    allmeavalue2(i,:,:,:) = GLFM_gedatest_heternorm(i, q, n_seq, p_seq);
end
size(allmeavalue2)
stime = toc(s);

%%  Example 2': for all j, X_ij is normal with heteroskedasticity.
s = tic;
[meavluae] = GLFM_gedatest_heternorm2(1);
stime = toc(s);

N = 1000
% simu setting
q = 2; n_seq = [150, 200,250];p_seq = [50,75,100];;
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue2 = zeros(N, n_num, p_num, n_mea);
s = tic;
parfor i = 1:N  
    disp(i)
    allmeavalue2(i,:,:,:) = GLFM_gedatest_heternorm2(i, q, n_seq, p_seq);
end
size(allmeavalue2)
stime = toc(s);


%%  Example 3:  Poisson variables
% Repeat N times 
s = tic;
[meavluae] = GLFM_gedatest_pois(1);
stime = toc(s);

N = 1000
% simu setting
q = 3; 
n_seq = [50, 100, 150];
p_seq = [100,150, 200, 250]; 
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue = zeros(N, n_num, p_num, n_mea);
s = tic;
parfor i = 1:N  
    allmeavalue(i,:,:,:) = GLFM_gedatest_pois(i, q, n_seq, p_seq);
end
size(allmeavalue)
stime = toc(s);
squeeze(allmeavalue(1,:,:,1)), squeeze(allmeavalue(1,:,:,3))
squeeze(allmeavalue(1,:,:,2)), squeeze(allmeavalue(1,:,:,4))

%%  Example 3':  Poisson variables
% Repeat N times 
s = tic;
[meavluae] = GLFM_gedatest_pois2(1);
stime = toc(s);

N = 1000
% simu setting
q = 3; 
n_seq = [50, 100, 150];
p_seq = [100,150, 200, 250]; 
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue = zeros(N, n_num, p_num, n_mea);
s = tic;
parfor i = 1:N  
    allmeavalue(i,:,:,:) = GLFM_gedatest_pois2(i, q, n_seq, p_seq);
end
size(allmeavalue)
stime = toc(s);
squeeze(allmeavalue(1,:,:,1)), squeeze(allmeavalue(1,:,:,3))
squeeze(allmeavalue(1,:,:,2)), squeeze(allmeavalue(1,:,:,4))


%% Example 4: mixture of poisson and binary variables
s = tic;
[meavluae] = GLFM_gedatest_pois_bino(1);
stime = toc(s);

N = 1000
% simu setting
q = 4; 
n_seq = [100, 200, 300];p_seq = [100, 200, 300, 400]; 
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue = zeros(N, n_num, p_num, n_mea);
s = tic;
parfor i = 1:N  
    allmeavalue(i,:,:,:) = GLFM_gedatest_pois_bino(i, q, n_seq, p_seq);
end
size(allmeavalue)
stime = toc(s);
squeeze(allmeavalue(1,:,:,1)), squeeze(allmeavalue(1,:,:,3))
squeeze(allmeavalue(1,:,:,2)), squeeze(allmeavalue(1,:,:,4))

%% Example 4': mixture of poisson and binary variables
s = tic;
[meavluae] = GLFM_gedatest_pois_bino2(1);
stime = toc(s);

N = 1000
% simu setting
q = 4; 
n_seq = [100, 200, 300];p_seq = [100, 200, 300, 400]; 
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue = zeros(N, n_num, p_num, n_mea);
s = tic;
parfor i = 1:N  
    allmeavalue(i,:,:,:) = GLFM_gedatest_pois_bino2(i, q, n_seq, p_seq);
end
size(allmeavalue)
stime = toc(s);
squeeze(allmeavalue(1,:,:,1)), squeeze(allmeavalue(1,:,:,3))
squeeze(allmeavalue(1,:,:,2)), squeeze(allmeavalue(1,:,:,4))

%% Example 5: Mixture of normal and poisson
s = tic;
[meavluae] = GLFM_gedatest_norm_pois(1)
stime = toc(s);

N = 1000
% simu setting
q = 5; 
n_seq = [100, 200, 300];p_seq = [100, 200, 300, 400]; 
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue = zeros(N, n_num, p_num, n_mea);
s = tic;
parfor i = 1:N  
    allmeavalue(i,:,:,:) = GLFM_gedatest_norm_pois(i, q, n_seq, p_seq);
end
size(allmeavalue)
stime = toc(s);
squeeze(allmeavalue(1,:,:,1)), squeeze(allmeavalue(1,:,:,3))
squeeze(allmeavalue(1,:,:,2)), squeeze(allmeavalue(1,:,:,4))
save norm_pois_res20200129
fprintf('elapsed time is %4f !\n', stime);

%% Example 5': Mixture of normal and poisson
s = tic;
[meavluae] = GLFM_gedatest_norm_pois2(1)
stime = toc(s);

N = 1000
% simu setting
q = 5; 
n_seq = [100, 200, 300];p_seq = [100, 200, 300, 400]; 
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue = zeros(N, n_num, p_num, n_mea);
s = tic;
parfor i = 1:N  
    allmeavalue(i,:,:,:) = GLFM_gedatest_norm_pois2(i, q, n_seq, p_seq);
end
size(allmeavalue)
stime = toc(s);
squeeze(allmeavalue(1,:,:,1)), squeeze(allmeavalue(1,:,:,3))
squeeze(allmeavalue(1,:,:,2)), squeeze(allmeavalue(1,:,:,4))
save norm_pois_res20200129
fprintf('elapsed time is %4f !\n', stime);
%% Example 6: Mixture of normal, poisson and binary.
s = tic;
[meavluae] = GLFM_gedatest_npb(1)
stime = toc(s);

N = 1000
% simu setting
q = 6;ex = 'npb';
n_seq = [200, 300, 400];p_seq = [200, 300, 400, 500]; 
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue = zeros(N, n_num, p_num, n_mea);
s = tic;
parfor i = 1:N  
    allmeavalue(i,:,:,:) = GLFM_gedatest_npb(i, q, n_seq, p_seq);
end
size(allmeavalue)
stime = toc(s)
squeeze(allmeavalue(1,:,:,1)), squeeze(allmeavalue(1,:,:,3))
squeeze(allmeavalue(1,:,:,2)), squeeze(allmeavalue(1,:,:,4))
% save npb_res20200129
fprintf('elapsed time is %4f !\n', stime);

%% Example 6': Mixture of normal, poisson and binary.
s = tic;
[meavluae] = GLFM_gedatest_npb2(1);
stime = toc(s);

N = 1000
% simu setting
q = 6;ex = 'npb';
n_seq = [200, 300, 400];p_seq = [200, 300, 400, 500]; 
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue = zeros(N, n_num, p_num, n_mea);
s = tic;
parfor i = 1:N  
    allmeavalue(i,:,:,:) = GLFM_gedatest_npb2(i, q, n_seq, p_seq);
end
size(allmeavalue)
stime = toc(s)
squeeze(allmeavalue(1,:,:,1)), squeeze(allmeavalue(1,:,:,3))
squeeze(allmeavalue(1,:,:,2)), squeeze(allmeavalue(1,:,:,4))
% save npb_res20200129
fprintf('elapsed time is %4f !\n', stime);


%% %% Part 2:  Model selection: select the number of factors q.
fprintf('Model selection now!\n')
% Example 1: q= 1. 
% Vr = [c,  r / min(p/w^4,n/w^6) * log(min(p/w^4,n/w^6))];
i=7; ex = 'homonorm2'; q = 1; n_seq = [30,50,100];p_seq = [30,50,100,150]; 
s1 = tic;
[allIC1] = GFM_PCslectall1(i, q, n_seq, p_seq, 1:3);
s2 = toc(s1)
allIC1
    % 2.8500    1.1900    1.0000    1.0000
    % 1.1400    1.7100    1.0000    1.0000
    % 1.0000    1.0000    1.0000    1.0000
%% BEFORE
% 2     1     1     1
% 1     1     1     1
% 1     1     1     1
%% OURS
% 1     1     1     1
% 1     1     1     1
% 1     1     1     1
% GFM  and LFM to select q
N = 100;
ex = 'homonorm2'; q = 1; n_seq = [30,50,100];p_seq = [30,50,100,150]; 
n_num = 3; p_num = 4; method = 2;
allIC= zeros(N, n_num, p_num,method);
s = tic;
parfor i = 1:N  
    i
    allIC(i,:,:, 1) = GFM_PCslectall1(i, q, n_seq, p_seq, 1:3);
end

N = 100;
ex = 'homonorm2'; q = 1; n_seq = [30,50,100];p_seq = [30,50,100,150]; 
n_num = 3; p_num = 4; method = 2;
allIC= zeros(N, n_num, p_num,method);
s = tic;
parfor i = 1:N  
    allIC(i,:,:, 1) = GFM_PCslectall1(i, q, n_seq, p_seq, 1:3,'PC');
end

% save example_homonorm_selq20200129
stime = toc(s);
stime

% Exmaple 2: q = 2.
fprintf('Model selection now!\n')
i=7; ex = 'heternorm2'; 
q = 2; n_seq = [100, 150, 200];p_seq = [100,150, 200, 250]; 
s1 = tic;
warning('off');
[allIC1] = GFM_PCslectall2(i, q, n_seq, p_seq, 1:3);
warning('off');
s2 = toc(s1)
allIC1
%% BEFORE
% 3     3     3     3
% 3     3     3     3
% 3     3     3     3

% GFM  and LFM to select q
N = 100
ex = 'heternorm2'; q = 2; n_seq = [100, 150, 200];p_seq = [100,150, 200, 250];
n_num = 3; p_num = 4; method = 2;
allIC= zeros(N, n_num, p_num, method);
s = tic;
parfor i = 1:N  
    allIC(i,:,:, 1) = GFM_PCslectall2(i, q, n_seq, p_seq, 1:3);
end
    % 1.8200    1.0100    1.0000    1.0000
    % 1.0000    1.1100    1.0000    1.0000
    % 1.0000    1.0000    1.0100    1.0000
% save example_heternorm_selq20200129
stime = toc(s);
stime

% Example 3: q=3
% Vr = [c , 0.6 * r / (w^6/n + w^4/p) * log(min(n/w^6, p/w^4))];  
i=7; ex = 'pois2'; q = 3; n_seq = [50, 100, 150];p_seq = [100,150, 200, 250]; 
s1 = tic;
[allIC1] = GFM_PCslectall3(i, q, n_seq, p_seq, 1:5);
s2 = toc(s1)
allIC1

%% BEFORE
% 5     5     5     5
% 5     5     5     5
% 5     5     5     5
%% OURS 
% 3     5     5     5
% 3     3     3     3
% 3     3     3     3

% 2.9600    1.0400    1.0000    1.0000
% 3.0000    1.5600    1.0100    1.0000
% 1.5300    2.1300    1.0200    1.0000
% GFM  and LFM to select q
N = 100
ex = 'pois2'; q = 3; n_seq = [50, 100, 150];p_seq = [100,150, 200, 250]; 
n_num = 3; p_num = 4; method = 2;
allIC= zeros(N, n_num, p_num, method);
s = tic;
parfor i = 1:N  
    allIC(i,:,:, 1) = GFM_PCslectall3(i, q, n_seq, p_seq, 1:5);
end
parfor i = 1:N  
    allIC(i,:,:, 2) = LFM_PCslectall(i, ex, q, n_seq, p_seq, 1:10);
end
size(allIC)
stime = toc(s);
stime

% Example 4:  q = 4.
% Vr = [c , 0.5 / r * (w^6/n + w^4/p) * log(min(n/w^6, p/w^4))];  
i=7; q = 4; ex =  'pois_bino2';n_seq = [100, 200, 300];p_seq = [100, 200, 300, 400]; 
s1 = tic;
[allIC1] = GFM_PCslectall4(i, q, n_seq, p_seq, 1:5);
s2 = toc(s1)
allIC1

% 3.4100    2.6600    1.0800    1.0100
% 2.0200    3.6100    1.6200    1.0000
% 1.9500    2.5200    2.4300    1.0000
%% BEFORE
% 5     5     5     5
% 5     5     5     5
% 5     5     5     5
%% OURS
% 4     4     4     4
% 2     2     3     4
% 1     1     1     1

% GFM  and LFM to select q
N = 100
q = 4; ex =  'pois_bino2';n_seq = [100, 200, 300];p_seq = [100, 200, 300, 400];  
n_num = 3; p_num = 4; method = 2;
allIC= zeros(N, n_num, p_num, method);
s = tic;
parfor i = 1:N  
    allIC(i,:,:, 1) = GFM_PCslectall4(i, q, n_seq, p_seq, 1:5);
end
parfor i = 1:N  
    allIC(i,:,:, 2) = LFM_PCslectall(i, ex, q, n_seq, p_seq, 1:10);
end
size(allIC)
% save example_pois_bino_selq20200129
stime = toc(s);
stime

% Example 5: q = 5
% Vr = [c, 0.4 * r / (p/w^4+n/w^6) * log(min(p/w^4,n/w^6))];
i=7; q = 5; ex =  'nrom_pois2';n_seq = [100, 200, 300];p_seq = [100, 200, 300, 400]; 
s1 = tic;
[allIC1] = GFM_PCslectall5(i, q, n_seq, p_seq, 3:6);
s2 = toc(s1)
allIC1

%% BEFORE
% 6     6     6     6
% 6     6     6     6
% 6     6     6     6
%% NOW
% 5     6     6     6
% 5     5     5     5
% 5     5     5     5

% GFM  and LFM to select q
N = 100
q = 5; ex =  'norm_pois';n_seq = [100, 200, 300];p_seq = [100, 200, 300, 400];  
n_num = 3; p_num = 4; method = 2;
allIC= zeros(N, n_num, p_num, method);
s = tic;
parfor i = 1:N  
    allIC(i,:,:, 1) = GFM_PCslectall5(i, q, n_seq, p_seq, 3:6);
end
parfor i = 1:N  
    allIC(i,:,:, 2) = LFM_PCslectall(i, ex, q, n_seq, p_seq, 1:10);
end
size(allIC)
save norm_pois_selq20200129
stime = toc(s);
stime

% EXample 6: q= 6.
% Vr = [c,  0.6*0.5*r / (p/w^4+n/w^6) * log(min(p/w^4,n/w^6))];
i=7; q = 6;ex = 'npb2';n_seq = [200, 300, 400];p_seq = [200, 300, 400, 500]; 
s1 = tic;
[allIC1] = GFM_PCslectall6(i, q, n_seq, p_seq, 4:8);
s2 = toc(s1)
allIC1
    % 4.0300    4.0300    4.0000    4.0100
    % 4.0000    4.0000    4.0100    4.0000
    % 4.0100    4.0300    4.0000    4.0100
    %% OURS
    % 4     6     6     6    
    % 4     4     5     6    
    % 4     4     4     5    
    %% BEFORE
    % 8     8     8     8
    % 8     8     8     8
    % 8     7     8     8

% GFM  and LFM to select q
N = 100
q = 6;ex = 'npb2';n_seq = [200, 300, 400];p_seq = [200, 300, 400, 500]; 
n_num = 3; p_num = 4; method = 2;
allIC= zeros(N, n_num, p_num, method);
s = tic;
parfor i = 1:N  
    allIC(i,:,:, 1) = GFM_PCslectall6(i, q, n_seq, p_seq, 4:8);
end
parfor i = 1:N  
    allIC(i,:,:, 2) = LFM_PCslectall(i, ex, q, n_seq, p_seq, 1:10);
end
size(allIC)
% save npb_selq20200129
stime = toc(s);
stime

%% %% Part 3: Running time
% Test the computation efficiency.
% The case for n = 500 and p=10000
clear;
N = 10;
cctMat = zeros(N, 3);
n = 500; p=10000; i = 1;
for i = 1:N
    fprintf('i = %d \n', i);
    [X, H, B, mu] = gendata(i, n, p, 'npb');
    G = [mu', B];
    % unified functions test
    type = cell(3,2);
    type{1,1} = 'normal'; type{1,2} = 'identity';
    type{2,1} = 'poisson'; type{2,2} = 'log';
    type{3,1} = 'binomial';  type{3,2} = 'logit';
    group = [ones(1,floor(p/3)), 2*ones(1, floor(2*p/3)-floor(p/3)), 3*ones(1, p-floor(2*p/3))];
    q= 6; dropout=[3]; output = 1; 
    X = single(X);
    s = tic;
    [hH, hB, hmu, history] = gfm_eval_intercept_init(X, group, type, q, dropout, 1e-4, 2,output);
    stime = toc(s);
    history
    hG = [hmu, hB];
    % Compare the estimation performance by GFM and LFM.
    cctMat(i,:) = [measurefun(H, hH),measurefun(G, hG), stime];
    
end

% The case for n = 5000 and p=100000.
clear;
N = 10;
cctMat = zeros(N, 3);
n = 5000; p=100000; i = 1;
for i = 1:N
    fprintf('i = %d \n', i);
    [X, H, B, mu] = gendata(i, n, p, 'npb');
    G = [mu', B];
    % unified functions test
    type = cell(3,2);
    type{1,1} = 'normal'; type{1,2} = 'identity';
    type{2,1} = 'poisson'; type{2,2} = 'log';
    type{3,1} = 'binomial';  type{3,2} = 'logit';
    group = [ones(1,floor(p/3)), 2*ones(1, floor(2*p/3)-floor(p/3)), 3*ones(1, p-floor(2*p/3))];
    q= 6; dropout=[3]; output = 1; 
    X = single(X);
    s = tic;
    [hH, hB, hmu, history] = gfm_eval_intercept_init(X, group, type, q, dropout, 1e-4, 2,output);
    stime = toc(s);
    history
    hG = [hmu, hB];
    % Compare the estimation performance by GFM and LFM.
    cctMat(i,:) = [measurefun(H, hH),measurefun(G, hG), stime];
end


%% %% Part 4: It is used to validate the effeciency gain by one-step update.
% Repeat N times
s = tic;
effmeasureVal = GFM_gedeffi_norm_pois(1);
stime = toc(s);

N = 1000
% simu setting
q = 5; n_seq = [50, 100]; p_seq = [50, 100];
n_num = length(n_seq); p_num = length(p_seq);
nBH = 2; n_init= 3; n_os = 2;
alleffmeasureVal = zeros(N, n_num, p_num, nBH, n_init, n_os);
s = tic;
parfor i = 1:N  
    alleffmeasureVal(i,:,:,:,:,:) = GFM_gedeffi_norm_pois(i, q, n_seq, p_seq);
end
size(alleffmeasureVal)
stime = toc(s);
squeeze(alleffmeasureVal(1,1,1,:,:, 1))
squeeze(alleffmeasureVal(1,1,1,:,:, 2))

%% %% Part 5: As a reviewer suggests, we will compare the two methods to
% exert the identifiability conditions  by Example 1 and  Example 5!
% Repeat N times
% -------------Example 1
s = tic;
[meavalue ] = GFM_gedtwoIC_homonorm(1);
squeeze(meavalue(1,1,:,:))
stime = toc(s);

N = 1000
% simu setting
q = 1; n_seq = [50, 100]; p_seq = [50, 100];
n_num = length(n_seq); p_num = length(p_seq);
nBH = 2; nIC = 2;
alltwoICval = zeros(N, n_num, p_num, nBH, nIC);
s = tic;
parfor i = 1:N  
    alltwoICval(i,:,:,:,:) = GFM_gedtwoIC_homonorm(i, q, n_seq, p_seq);
end
size(alltwoICval)
stime = toc(s);
squeeze(alltwoICval(1,1,1,:, :))
save norm_twoICmethod20200130

%-------------Example 5
s = tic;
[meavalue ] = GFM_gedtwoIC_norm_pois(1);
squeeze(meavalue(1,1,:,:))
stime = toc(s);

N = 1000;
% simu setting
q = 5; n_seq = [50, 100]; p_seq = [50, 100];
n_num = length(n_seq); p_num = length(p_seq);
nBH = 2; nIC = 2;
alltwoICval = zeros(N, n_num, p_num, nBH, nIC);
s = tic;
parfor i = 1:N  
    alltwoICval(i,:,:,:,:) = GFM_gedtwoIC_norm_pois(i, q, n_seq, p_seq);
end
size(alltwoICval)
stime = toc(s);
squeeze(alltwoICval(1:50,1,1,:, :))

%% %% Part 6: Generate data used for low rank model in R programming language.

ex = 'heternorm'; q = 2; n_seq = [100, 200];p_seq = [100, 200]
n_num = length(n_seq); p_num = length(p_seq); method = 2;
N = 1000;
for jj = 1:n_num
    fprintf('n= %d\n', n_seq(jj));
    for kk = 1:p_num
        tmp_cell = cell(2,1);
        tmp_X_Theta_array = zeros(N, n_seq(jj), p_seq(kk), 2);
        % tmp_Theta_array = tmp_X_array;
        for i = 1:N
            [X, H, B, mu] = gendata(i, n_seq(jj), p_seq(kk), ex, q);
            tmp_X_Theta_array(i,:,:,1) = X;
            tmp_X_Theta_array(i,:,:,2) = H* B';
        end
        eval(['save RmatlabData/heternorm_N100_n', num2str(n_seq(jj)), 'p', num2str( p_seq(kk)), ...
           ' tmp_X_Theta_array' ])
    end
end

% Generate data used for R 
ex = 'pois'; q = 3; n_seq = [100, 200];p_seq = [100, 200]
n_num = length(n_seq); p_num = length(p_seq); method = 2;
N = 1000;
for jj = 1:n_num
    fprintf('n= %d\n', n_seq(jj));
    for kk = 1:p_num
        tmp_cell = cell(2,1);
        tmp_X_Theta_array = zeros(N, n_seq(jj), p_seq(kk), 2);
        % tmp_Theta_array = tmp_X_array;
        for i = 1:N
            [X, H, B, mu] = gendata(i, n_seq(jj), p_seq(kk), ex, q);
            tmp_X_Theta_array(i,:,:,1) = X;
            tmp_X_Theta_array(i,:,:,2) = H* B';
        end
        eval(['save RmatlabData/pois_N100_n', num2str(n_seq(jj)), 'p', num2str( p_seq(kk)), ...
           ' tmp_X_Theta_array' ])
    end
end

% Generate data used for R 
ex = 'norm_pois'; q = 5; n_seq = [100, 200];p_seq = [100, 200]
n_num = length(n_seq); p_num = length(p_seq); method = 2;
N = 1000;
for jj = 1:n_num
    fprintf('n= %d\n', n_seq(jj));
    for kk = 1:p_num
        tmp_cell = cell(2,1);
        tmp_X_Theta_array = zeros(N, n_seq(jj), p_seq(kk), 2);
        % tmp_Theta_array = tmp_X_array;
        for i = 1:N
            [X, H, B, mu] = gendata(i, n_seq(jj), p_seq(kk), ex, q);
            tmp_X_Theta_array(i,:,:,1) = X;
            tmp_X_Theta_array(i,:,:,2) = H* B';
        end
        eval(['save RmatlabData/norm_pois_N100_n', num2str(n_seq(jj)), 'p', num2str( p_seq(kk)), ...
           ' tmp_X_Theta_array' ])
    end
end