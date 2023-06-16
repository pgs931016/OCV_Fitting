clc; close all; clear

% define the weighting 
% w = zeros(size(Cap')); % should have the same length with the data
% w(:)=1; % uniform weighting
    % prep: multi-start (tolerances)
%ms = MultiStart('UseParallel',true,'FunctionTolerance',1e-15,'XTolerance',1e-15);   
    % prep: optimset option (tolerances)
% options = optimoptions(@fmincon,'MaxIterations',5000,'StepTolerance',1e-15,'ConstraintTolerance', 1e-15, 'OptimalityTolerance', 1e-15);
    % initial guess and lower/upper bounds

% OCV_stoichiometry_model([0,0.0043,0.9,0.0043],OCP_n,OCP_p,OCV)
  
load('OCV_fit.mat')
%OCV_stoichiometry_model([0,0.0043,0.9,0.0043],OCP_n,OCP_p,OCV)
%import '/Users/g.park/Documents/gspark/MATLAB/Code/OCV_stoichiometry_model_06.m'
%Q_cell = 0.004265771728163;
%x_guess = [0,Q_cell,1,Q_cell];
x_guess = [0.01,Q_cell*1.2,0.9,Q_cell];
x_lb = [0,Q_cell*0.5,0,Q_cell*0.5];
x_ub = [1,Q_cell*2,1,Q_cell*2]; 
% Do something



%% Initial Guess
[~,OCV_guess] = OCV_stoichiometry_model_06(x_guess,OCP_n,OCP_p,OCV);

% compare initial guess and the data
figure(1) 
plot(OCV(:,1),OCV(:,2)); hold on
plot(OCV(:,1),OCV_guess)


    % optimization problem
        % obj = min(sum(OCV-OCV_model)^2 ; @stoichiometry_fit_vX.m
% problem = createOptimProblem('fmincon','x0',x_guess,'objective',@(x) OCV_stoichiometry_model(x, OCP_n, OCP_p, w),'lb',x_lb,'ub',x_ub,'options','MaxIter',100);
% run(problem);


% fmincon을 사용하여 최적화 수행
    options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 5000);
    % format: optimoptions(function_name, option_name, option_value, ...)
%     fhandle_cost = @(x)OCV_stoichiometry_model_06(x, OCP_n, OCP_p, OCV_golden);
   % function [x_id, fval] = fmincon_OCV_stoichiometry_model(x, OCP_n, OCP_p, OCV)
%    OCV_stoichiometry_model = @(x, OCP_n, OCP_p, OCV)OCV_stoichiometry_model_06(x, OCP_n, OCP_p, OCV);

%   OCV_stoichiometry_model = @(x, OCP_n, OCP_p, OCV) ...
%   OCV_stoichiometry_model_06(x, OCP_n, OCP_p, OCV);

%   [x_id, fval] = fmincon(OCV_stoichiometry_model, ...
%   x_guess, [], [], [], [], x_lb, x_ub, [],options);



    fhandle_cost = @(x)OCV_stoichiometry_model_06(x, OCP_n, OCP_p, OCV);
    [x_id, fval, exitflag, output] = fmincon(fhandle_cost, ...
        x_guess, [], [], [], [], x_lb, x_ub, [],options);



[cost_hat, OCV_hat] = OCV_stoichiometry_model_06(x_id,OCP_n,OCP_p,OCV);


figure(1)
plot(OCV(:,1),OCV_hat);
legend('data','guess','fit')
xlim([0,Q_cell])
ylim([2.6 4.2])
   
%end

% dataList.QN(k_list)=QN;
% dataList.QP(k_list)=QP;
% dataList.QLi(k_list)=TotalLi;
% dataList.x0(k_list)=x_0;
% dataList.x100(k_list)=x_100;
% dataList.y0(k_list)=y_0;
% dataList.y100(k_list)=y_100;
% dataList.OCVerr(k_list)=norm(OCV_sim - OCV);


    %% TO Do
    %1) OCV_data, OCV(x0,....), OCV(x_hat) overlay plot
