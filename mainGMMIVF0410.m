% Part 5 (mainGMMIVF0410.m) is designed to examine an alternative classification scheme, dividing wealth and age into seven and five groups, respectively.
clearvars;
clc;

filename = 'C:\Users\USER\Desktop\moment condition_77.xlsx';
data = readtable(filename);
t   = data.t; % time period
n20 = data.n20; % no of 20 twins
n02 = data.n02; % no of 02 twins
n11 = data.n11; % no of 11 twins
n3 = data.n3; % no of multifetal
n_total = n20 + n02 + n11+ n3;
ivf_awt = data.ivf_awt; 
d20 = data.d20; 
d02 = data.d02; 
d11 = data.d11; 
d3 = data.d3;

delta_20 = n20 ./ n_total;
delta_02 = n02 ./ n_total;
delta_11 = n11 ./ n_total;
delta_3  = n3  ./ n_total;

% --- Upper & Lower Bounds ---
gamma1_min = 0.85;  gamma1_max = 0.9999;         % Range of gamma_1 
gamma2_min = 0.80;  gamma2_max = 0.95;         % Range of gamma_2 
gamma3_min = 0.63;  gamma3_max = 0.80;         % Range of gamma_3 
gamma4_min = 0.63;  gamma4_max = 0.75;         % Range of gamma_4 
theta_20_min = 0.22; theta_20_max = 0.26; % Range of  theta_20
theta_02_min = 0.19; theta_02_max = 0.26; % Range of  theta_02
theta_11_min = 0.50; theta_11_max = 0.60; % Range of  theta_11

lb = [theta_20_min; theta_02_min; theta_11_min; gamma1_min; gamma2_min; gamma3_min; gamma4_min]; % 下界
ub = [theta_20_max; theta_02_max; theta_11_max; gamma1_max; gamma2_max; gamma3_max; gamma4_max]; % 上界

num_restarts = 100;
best_solution = [];
best_obj = Inf;  

data_mean = mean(table2array(data));  
seed = sum(data_mean);
rng(seed);  % Set up the random seed

% Select randomized starting value
for i = 1:num_restarts
    
    % Set up the initial values for parameters
    theta_20 = rand(1) * (theta_20_max - theta_20_min) + theta_20_min;
    theta_02 = rand(1) * (theta_02_max - theta_02_min) + theta_02_min;
    theta_11 = rand(1) * (theta_11_max - theta_11_min) + theta_11_min;
    gamma1_initial = rand(1, 1) * (gamma1_max   - gamma1_min)   + gamma1_min;
    gamma2_initial = rand(1, 1) * (gamma2_max   - gamma2_min)   + gamma2_min;
    gamma3_initial = rand(1, 1) * (gamma3_max - gamma3_min) + gamma3_min;
    gamma4_initial = rand(1, 1) * (gamma4_max - gamma4_min) + gamma4_min;
    if i<=50
        x0_tilde=[0.2436;0.2180;0.5332;0.9418;0.9111;0.6967;0.6965]+rand(7,1)*0.05;
    else
        x0_tilde = [theta_20; theta_02; theta_11; gamma1_initial; gamma2_initial;gamma3_initial; gamma4_initial];
    end

    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'TolFun', 1e-6, 'TolX', 1e-6);
    objFun = @(vars_tilde) mygmmF1(vars_tilde, t, ivf_awt, d20, d02, d11, d3, ...
                              delta_20, delta_02, delta_11, delta_3, ...
                              n20, n02, n11, n3);
                                   
    [solution_tilde, obj_val] = fmincon(objFun, x0_tilde, [], [], [], [], lb, ub, [], options);
    solution_tilde
    obj_val
    checklb = sum(solution_tilde==lb);
    checkub = sum(solution_tilde==ub);

    % Updated the optimal solution 
    if and(obj_val<best_obj, and(checklb<1, checkub<1))
    % If obj_val<best_obj
        best_obj = obj_val;
        best_solution = solution_tilde;
    end    
end

% Show the optimal soultion
disp('Optimized Variables (theta, gamma):');
disp(best_solution);
% Show the seed used for randomization
disp(['Using seed: ', num2str(seed)]);

[SS, nn]=mygmmS(best_solution, t, ivf_awt, d20, d02, d11, d3, ...
                              delta_20, delta_02, delta_11, delta_3, ...
                              n20, n02, n11,n3);
Ga = mygmmSS(best_solution, t, ivf_awt, d20, d02, d11, d3, ...
                              delta_20, delta_02, delta_11, delta_3, ...
                              n20, n02, n11,n3);

AA=inv(Ga'*Ga);
BB=Ga'* (SS)*Ga;
varcov=  AA*BB*AA;

seA=sqrt(diag((varcov)/diag(nn)));
% se=sqrt(diag(varcov));
disp('Optimized (theta, gamma) and SE:');
resultX=[best_solution, seA]
disp('sum of theta 20 02 11, theta_other:');
[sum(best_solution(1:3)), 1-sum(best_solution(1:3))]
