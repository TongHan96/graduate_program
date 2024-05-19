clear;
pctRunOnAll warning('off', 'all');
%% %% Part 1: Estimate B and Hm for the six examples.
%%  Example 1: for all j, X_ij is normal with homoskedasticity.
% Repeat N times 
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

%% BEFORE （1.mat）
% results_H =
% 
%     0.8889    0.9111    0.9207    0.9250
%     0.8974    0.9075    0.9151    0.9205
%     0.8906    0.9006    0.9084    0.9183
% 
% 
% results_B =
% 
%     0.9295    0.9103    0.8891    0.8721
%     0.9506    0.9167    0.9003    0.8872
%     0.9533    0.9250    0.9083    0.8974

% results_H =
% 
%     0.0121    0.0097    0.0085    0.0073
%     0.0104    0.0115    0.0099    0.0080
%     0.0109    0.0099    0.0078    0.0077
% 
% 
% results_B =
% 
%     0.0114    0.0103    0.0097    0.0092
%     0.0073    0.0091    0.0094    0.0088
%     0.0077    0.0091    0.0078    0.0077

%% OURS （EX1_new.mat）
% results_H =
% 
%     0.7438    0.7767    0.7934    0.7974
%     0.7636    0.7652    0.7711    0.7788
%     0.7455    0.7430    0.7528    0.7694
% 
% 
% results_B =
% 
%     0.8089    0.7763    0.7420    0.7130
%     0.8587    0.7784    0.7450    0.7232
%     0.8572    0.7840    0.7535    0.7343

% results_H =
% 
%     0.0371    0.0302    0.0286    0.0265
%     0.0271    0.0363    0.0382    0.0306
%     0.0316    0.0314    0.0283    0.0288
% 
% 
% results_B =
% 
%     0.0382    0.0294    0.0288    0.0260
%     0.0278    0.0342    0.0348    0.0308
%     0.0303    0.0327    0.0270    0.0259

%%  Example 2: for all j, X_ij is normal with heterskedasticity.
N = 1000
% simu setting
q = 2; n_seq = [150, 200,250];p_seq = [40,60,80,100];
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


%% BEFORE (2.mat)
% results_H =
% 
%     0.8778    0.8910    0.9081    0.9059
%     0.8810    0.8970    0.8972    0.9023
%     0.8795    0.8889    0.9041    0.9127
% 
% 
% results_B =
% 
%     0.9395    0.9244    0.9169    0.9029
%     0.9059    0.9289    0.9152    0.9163
%     0.9365    0.9435    0.9448    0.9371

% results_H =
% 
%     0.0226    0.0184    0.0150    0.0149
%     0.0211    0.0163    0.0164    0.0151
%     0.0209    0.0180    0.0139    0.0121
% 
% 
% results_B =
% 
%     0.0153    0.0147    0.0137    0.0144
%     0.0255    0.0148    0.0150    0.0124
%     0.0175    0.0116    0.0094    0.0094

%% OURS （EX2_new.mat）
% results_H =
% 
%     0.7478    0.7655    0.7979    0.7853
%     0.7528    0.7781    0.7702    0.7694
%     0.7501    0.7565    0.7838    0.7979
% 
% 
% results_B =
% 
%     0.8333    0.8042    0.7993    0.7654
%     0.7639    0.8156    0.7816    0.7761
%     0.8259    0.8347    0.8452    0.8295

% results_H =
% 
%     0.0775    0.0632    0.0630    0.0648
%     0.0663    0.0576    0.0546    0.0823
%     0.0520    0.0601    0.0468    0.0510
% 
% 
% results_B =
% 
%     0.0861    0.0749    0.0671    0.0735
%     0.0843    0.0638    0.0621    0.0887
%     0.0631    0.0668    0.0518    0.0545

%%  Example 3:  Poisson variables
% Repeat N times 
N = 1000;
% simu setting
q = 3; 
n_seq = [50, 100, 150];
p_seq = [100,150, 200, 250]; 
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue3 = zeros(N, n_num, p_num, n_mea);
s = tic;
parfor i = 1:N  
    allmeavalue3(i,:,:,:) = GLFM_gedatest_pois(i, q, n_seq, p_seq);
end

allmeavalue3_ = zeros(N, n_num, p_num, n_mea);
s = tic;
parfor i = 1:N  
    allmeavalue3_(i,:,:,:) = GLFM_gedatest_pois2(i, q, n_seq, p_seq);
end

%% BEFORE (3.mat)
% results_H =
% 
%     0.9616    0.9751    0.9813    0.9788
%     0.9765    0.9789    0.9827    0.9843
%     0.9742    0.9754    0.9777    0.9819
% 
% 
% results_B =
% 
%     0.9046    0.9117    0.9271    0.9128
%     0.9591    0.9578    0.9567    0.9532
%     0.9704    0.9639    0.9644    0.9635
%
% results_H =
% 
%     0.0751    0.0484    0.0420    0.0527
%     0.0194    0.0148    0.0090    0.0018
%     0.0043    0.0083    0.0228    0.0137
% 
% 
% results_B =
% 
%     0.0607    0.0342    0.0291    0.0376
%     0.0212    0.0074    0.0074    0.0041
%     0.0059    0.0056    0.0098    0.0057

%% OURS
% results_H =
% 
%     0.9666    0.9700    0.9749    0.9781
%     0.9522    0.9558    0.9623    0.9641
%     0.9377    0.9426    0.9489    0.9566
% 
% 
% results_B =
% 
%     0.9008    0.8995    0.9028    0.8964
%     0.9402    0.9258    0.9171    0.9101
%     0.9451    0.9319    0.9215    0.9189

% results_H =
% 
%     0.0164    0.0051    0.0042    0.0123
%     0.0055    0.0055    0.0043    0.0040
%     0.0060    0.0058    0.0052    0.0041
% 
% 
% results_B =
% 
%     0.0175    0.0149    0.0097    0.0117
%     0.0077    0.0079    0.0076    0.0074
%     0.0071    0.0066    0.0076    0.0066

%% Example 4: mixture of poisson and binary variables
N = 1000;
% simu setting
q = 4; 
n_seq = [100, 200, 300];p_seq = [100, 200, 300, 400]; 
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue4 = zeros(N, n_num, p_num, n_mea);
parfor i = 1:N  
    allmeavalue4(i,:,:,:) = GLFM_gedatest_pois_bino(i, q, n_seq, p_seq);
end

allmeavalue4_ = zeros(N, n_num, p_num, n_mea);
parfor i = 1:N  
    allmeavalue4_(i,:,:,:) = GLFM_gedatest_pois_bino2(i, q, n_seq, p_seq);
end

%% BEFORE

%% OURS


%% Example 5: Mixture of normal and poisson
N = 1000;
% simu setting
q = 5; 
n_seq = [100, 200, 300]; p_seq = [100, 200, 300, 400]; 
n_num = length(n_seq); p_num = length(p_seq);
n_mea = 4;
allmeavalue5 = zeros(N, n_num, p_num, n_mea);
parfor i = 1:N  
    allmeavalue5(i,:,:,:) = GLFM_gedatest_norm_pois(i, q, n_seq, p_seq);
end
allmeavalue5_ = zeros(N, n_num, p_num, n_mea);
parfor i = 1:N  
    allmeavalue5_(i,:,:,:) = GLFM_gedatest_norm_pois2(i, q, n_seq, p_seq);
end

%% BEFORE

%% OURS

%% Example 6: Mixture of normal and poisson and binomial
N = 1000;
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
allmeavalue_ = zeros(N, n_num, p_num, n_mea);
parfor i = 1:N
allmeavalue_(i,:,:,:) = GLFM_gedatest_npb2(i, q, n_seq, p_seq);
end


%% BEFORE

%% OURS

%% %% Part 2:  Model selection: select the number of factors q.
pctRunOnAll warning('off', 'all');
fprintf('Model selection now!\n')
%% Example 1: q= 1. 
% Vr = [c, r/(sqrt(n)/w^3 + sqrt(p)/w^2)^2 *log(min(sqrt(p)/w^2, sqrt(n)/w^3)^2)];
% OURS (ex1_new.mat)
% 1.0000    1.0100    3.1600    5.0000
% 1.0000    1.0000    1.0000    1.0800
% 1.0000    1.0000    1.0000    1.0000
% 
%  0    0.1000    0.6150         0
%  0         0         0    0.2727
%  0         0         0         0

% BEFORE (ex1.mat)
% 5     5     5     5
% 5     5     5     5
% 5     5     5     5
% 
% 0     0     0     0
% 0     0     0     0
% 0     0     0     0

%
i=7; ex = 'homonorm2'; q = 1; n_seq = [200,250,300];p_seq = [50,100,200,300]; 
[allIC1] = GFM_PCslectall1(i, q, n_seq, p_seq, 1:3);

%
N = 100;
ex = 'homonorm2'; q = 1; n_seq = [200,250,300];p_seq = [100,200,300,400]; 
n_num = 3; p_num = 4;
allIC= zeros(N, n_num, p_num);
parfor i = 1:N  
    allIC(i,:,:) = GFM_PCslectall1(i, q, n_seq, p_seq, 1:5);
end


%% Example 2: q= 2. 

% Vr = [c, r/(sqrt(n)/w^3 + sqrt(p)/w^2)^2 *log(min(sqrt(p)/w^2, sqrt(n)/w^3)^2)];
% OURS (ex2_new.mat)
% 2.9500    2.9100    2.9400    2.9700
% 2.8300    2.6100    2.4400    2.3500
% 3.0000    2.0500    1.6800    1.5200
%
% 0.2611    0.3509    0.2778    0.2227
% 0.4507    0.6340    0.7696    0.8087
%      0    0.8333    0.7898    0.7314

% BEFORE (ex2.mat)
% 3.0000    2.9900    2.9900    2.9900
% 2.9900    3.0000    2.9900    3.0000
% 3.0000    2.9900    3.0000    3.0000
%
%      0    0.1000    0.1000    0.1000
% 0.1000         0    0.1000         0
%      0    0.1000         0         0
i=7; ex = 'heternorm2'; q = 2; n_seq = [150, 200,250];p_seq = [40,60,80,100];
[allIC1] = GFM_PCslectall2(i, q, n_seq, p_seq, 1:3);
%
N = 100;
ex = 'heternorm2'; q = 2; n_seq = [150, 200,250];p_seq = [40,60,80,100];
n_num = 3; p_num = 4;
allIC= zeros(N, n_num, p_num);
parfor i = 1:N  
    allIC(i,:,:) = GFM_PCslectall2(i, q, n_seq, p_seq, 1:3);
end

%% Example 3: q=3.
% ATTENTION: when p is big but n keeps same relative small value, the first term tend to be smaller and the second term keeps the same, then the penalty is smaller so the estimators of the factor numbers tend to be bigger. just like when n = 100 and p = 200,300,400 
% Vr = [c, r/(sqrt(n)/w^3 + sqrt(p)/w^2)^2 *log(min(sqrt(p)/w^2, sqrt(n)/w^3)^2)];
% OURS (ex3_new.mat)
% 3.0100    3.2700    4.6400    4.9400
% 3.0100    3.0000    3.0300    3.0000
% 3.0000    3.0000    3.0000    3.0000
%
% 0.1000    0.5096    0.5599    0.2778
% 0.1000         0    0.2227         0
%      0         0         0         0

% BEFORE (ex3.mat)
% 4.9600    4.9400    4.9700    4.9800
% 4.9600    4.9900    4.9600    4.9900
% 4.9500    4.9900    4.9900    4.9900
%
% 0.1969    0.2387    0.1714    0.1407
% 0.1969    0.1000    0.2429    0.1000
% 0.2611    0.1000    0.1000    0.1000
i=7; ex = 'pois2'; q = 3; n_seq = [50, 100, 150];p_seq = [100,150, 200, 250]; 
[allIC1] = GFM_PCslectall3(i, q, n_seq, p_seq, 1:5);

%
N = 100;
ex = 'pois2'; q = 3; n_seq = [50, 100, 150];p_seq = [100,150, 200, 250]; 
n_num = 3; p_num = 4;
allIC= zeros(N, n_num, p_num);
parfor i = 1:N  
    allIC(i,:,:) = GFM_PCslectall3(i, q, n_seq, p_seq, 1:5);
end

%% Example 4: q=4.
% Vr = [c, r/(sqrt(n)/w^3 + sqrt(p)/w^2)^2 *log(min(sqrt(p)/w^2, sqrt(n)/w^3)^2)];
% OURS (ex4_new.mat)
% 4.3900    5.0000    4.9800    4.9900
% 3.1100    4.9800    4.8200    5.0000
% 2.9900    2.9600    3.0700    4.0000
%
% 0.8864         0    0.1407    0.1000
% 0.3145    0.1407    0.3861         0
% 0.1000    0.1969    0.2564         0

% BEFORE (ex4.mat)
% 4.2200    4.9800    4.9700    4.9800
% 4.9900    4.9800    5.0000    5.0000
% 5.0000    4.9900    5.0000    5.0000
%
% 0.8236    0.2000    0.1714    0.1407
% 0.1000    0.1407         0         0
%      0    0.1000         0         0

i=7; q = 4; ex =  'pois_bino2';n_seq = [100, 200, 300];p_seq = [100, 200, 300, 400]; 
[allIC1] = GFM_PCslectall4(i, q, n_seq, p_seq, 1:5);

%
N = 100; q = 4;
ex =  'pois_bino2';n_seq = [100, 200, 300];p_seq = [100, 200, 300, 400]; 
n_num = 3; p_num = 4;
allIC4 = zeros(N, n_num, p_num);
parfor i = 1:N  
    allIC4(i,:,:) = GFM_PCslectall4(i, q, n_seq, p_seq, 1:5);
end

%% Example 5: q=5.
% ATTENTION: when p is big but n keeps same relative small value, the first term tend to be smaller and the second term keeps the same, then the penalty is smaller so the estimators of the factor numbers tend to be bigger. just like when n = 100 and p = 200,300,400 
% Vr = [c, r/(sqrt(n)/w^3 + sqrt(p)/w^2)^2 *log(min(sqrt(p)/w^2, sqrt(n)/w^3)^2)];
% OURS (ex5_new.mat)
% 5.0300    5.4400    5.9700    5.9400
% 4.9800    5.0000    5.0000    5.0100
% 5.0100    5.0000    5.0000    5.0000
%
% 0.1714    0.4989    0.1714    0.2778
% 0.2000         0         0    0.1000
% 0.1000         0         0    0.1421

% BEFORE (ex5.mat)
% 5.9500    5.9300    5.9700    5.9400
% 5.9600    6.0000    5.9900    5.9900
% 5.9900    6.0000    6.0000    5.9700
%
% 0.2190    0.2564    0.1714    0.2778
% 0.3153         0    0.1000    0.1000
% 0.1000         0         0    0.2227

%
i=7; q = 5; ex =  'nrom_pois2';n_seq = [100, 200, 300];p_seq = [100, 200, 300, 400]; 
[allIC1] = GFM_PCslectall5(i, q, n_seq, p_seq, 3:6);

%
N = 100; q=5;
ex =  'nrom_pois2';n_seq = [100, 200, 300];p_seq = [100, 200, 300, 400]; 
n_num = 3; p_num = 4;
allIC5 = zeros(N, n_num, p_num);
parfor i = 1:N  
    allIC5 (i,:,:) = GFM_PCslectall5(i, q, n_seq, p_seq, 3:6);
end

%% Example 6: q = 6.
% Vr = [c, r/(sqrt(n)/w^3 + sqrt(p)/w^2)^2 *log(min(sqrt(p)/w^2, sqrt(n)/w^3)^2)];
% OURS (ex6_new.mat)
% 5.9800    6.0100    5.9900    6.0000
% 5.9700    5.9800    6.0000    6.0100
% 4.7900    5.9800    5.9900    6.0100
%
% 0.2000    0.1000    0.1000         0
% 0.1714    0.1407         0    0.1000
% 0.5183    0.1407    0.1000    0.1000

% BEFORE (ex6.mat)
% 7.9200     7.9200     7.9200      7.9500 
% 7.9400     7.8300     7.9200      7.8800 
% 7.9700     7.9000     7.9500      7.9400 
%
% 0.3075    0.3075    0.3674    0.2611
% 0.2778    0.5515    0.3387    0.4330
% 0.1714    0.3333    0.2190    0.2387


i=7; q = 6;ex = 'npb2';n_seq = [200, 300, 400];p_seq = [200, 300, 400, 500]; 
[allIC1] = GFM_PCslectall6(i, q, n_seq, p_seq, 4:8);
%
N = 100;
ex =  'npb2'; n_seq = [200, 300, 400];p_seq = [200, 300, 400, 500]; 
q=6; n_num = 3; p_num = 4;
allIC6 = zeros(N, n_num, p_num);
parfor i = 1:N  
    allIC6(i,:,:) = GFM_PCslectall6(i, q, n_seq, p_seq, 4:8);
end



