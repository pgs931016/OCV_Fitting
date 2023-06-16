clc;clear;close all

%input 지정(ocp,ocv데이터),function model,function cost,cost function
%minimization
%load the data
data1 = load('/Users/g.park/Library/CloudStorage/GoogleDrive-gspark@kentech.ac.kr/공유 드라이브/Battery Software Lab/Processed_Data/Hyundai_dataset/OCV/FCC_(5)_OCV_C100');
data2 = load('/Users/g.park/Library/CloudStorage/GoogleDrive-gspark@kentech.ac.kr/공유 드라이브/Battery Software Lab/Processed_Data/Hyundai_dataset/OCV/AHC_(5)_OCV_C100');
data3 = load('/Users/g.park/Library/CloudStorage/GoogleDrive-gspark@kentech.ac.kr/공유 드라이브/Battery Software Lab/Processed_Data/Hyundai_dataset/OCV/CHC_(5)_OCV_C100');

%half_cell_ocp(averaging the three) = golden ocp
 OCV   = data1.OCV_golden.OCVchg;
 OCP_n = data2.OCV_golden.OCVchg;
 OCP_p = data3.OCV_golden.OCVchg;
 
% %Variable name(dataList)
% %datalist = dataList;
 
%plot actual OCP data (parameter)
%for i=1:size(OCP_n,2) %MATLAB에서 배열의 열을 반복하는 데 사용
     figure(1); hold on; box on
     plot(OCP_n(:,1),OCP_n(:,2))
     figure(2); hold on; box on
     plot(OCP_p(:,1),OCP_p(:,2)) 
%end
figure(1)
title('Anode OCP')
xlabel('x in LixC6')
ylabel('OCP [V]')
figure(2)
title('Cathode OCP')
xlabel('y in LixMO2')
ylabel('OCP [V]')

% 
% % plot actual OCV data(Variable)
% %for i=1:size(datalist,1) %MATLAB에서 배열의 행을 반복하는 데 사용
    Cap = data1.OCV_golden.OCVchg(:,1); % [Ah] Discharged capacity
    Cap_end = Cap(end);
    Q_cell = abs(Cap_end);
   % Q_cell = Cap_end;
    OCV_chg = data1.OCV_golden.OCVchg(:,2); % [V] Cell Voltage     
    %measurement = [Cap,OCV]; % OCV measurement matrix [q [Ah], v[V]]
    figure(3); hold on; box on
    plot(Cap,OCV_chg)  
%end
figure(3)
title('Full Cell OCV')
xlabel('Cap[mAh]')
ylabel('OCV [V]')

save('OCV_fit.mat','Q_cell','OCP_n','OCP_p','OCV','Cap',"data1","data2","data3");

% % parameter 정의
% x0 = 0;
% y0 = 1;
% Qn = 0.0052;
% Qp = 0.0057;




% 
% % allocate the results
% x_0 = x_id(1);
% QN = x_id(2);
% y_0 = x_id(3);
% QP = x_id(4);
% 
% measurement = [Cap,OCV]; % OCV measurement matrix [q [Ah], v[V]]
%  % prep: multi-start (tolerances)
% ms = MultiStart('UseParallel',true,'FunctionTolerance',1e-15,'XTolerance',1e-15);
%    % prep: optimset option (tolerances)
% options = optimoptions(@fmincon,'MaxIterations',5000,'StepTolerance',1e-15,'ConstraintTolerance', 1e-15, 'OptimalityTolerance', 1e-15);
%     
%     % initial guess and lower/upper bounds
% x_guess = [0,Q_cell,1,Q_cell];
% x_lb = [0,Q_cell*0.5,0,Q_cell*0.5];
% x_ub = [1,Q_cell*2,1,Q_cell*2]; 
% 
%     % optimization problem
%         % obj = min(sum(OCV-OCV_model)^2 ; @stoichiometry_fit_vX.m
% problem = createOptimProblem('fmincon','x0',x_guess,'objective',@(x) ocv_fit(x, OCP_n, OCP_p, measurement, w),'lb',x_lb,'ub',x_ub,'options',options);
% [x_id,fval_ms,flag,outpt,allmins] = run(ms,problem,100);


%OCV_model = model_func(Cap,)



% % %fitting for a single measurement
% % %add another for loop for differetn samples
% %%initial guess and lower/upper bounds
% x_guess = [0,Q_cell,1,Q_cell];
% x_lb = [0,Q_cell*0.5,0,Q_cell*0.5];
% x_ub = [1,Q_cell*2,1,Q_cell*2];
% function [x, fval] = my_function (x) 
% 
% x = [x0, y0, Qn, Qp];
% 
% x0 = 0;
% y0 = 1;
% Qn = 0.0052;
% Qp = 0.0057;
% 
% %손실함수정의
% f = @(x0,y0,Qn,Qp) sum()
% 
% 
% % %select a cell
% for k_list = 1:size(data1.OCV_golden.OCVchg,1)
% 
% % %reduce the data
% %measurement = [Cap,OCV];
% % %end
% 
% %
% x0 =  
% 
% 
% % % define the weighting
% % w = zeros(size(SOC')); % should have the same length with the data
% % w(:)=1; % uniform weighting "(:)"는 배열의 모든 요소
% 
% 
% %prep: multi-start (tolerances)
% ms = MultiStart('UseParallel',true,'FunctionTolerance',1e-15,'XTolerance',1e-15);
% 
% 
% %  prep: optimset option (tolerances)
% options = optimoptions(@fmincon,'MaxIterations',5000,'StepTolerance',1e-15,'ConstraintTolerance', 1e-15, 'OptimalityTolerance', 1e-15);
% 
% 
% 
% % optimization problem\
%         % obj = min(sum(OCV-OCV_model)^2 ; @stoichiometry_fit_vX.m
% problem = createOptimProblem('fmincon','x0',x_guess,'objective',@(x) stoichiometry_fit_v3(x, OCP_n, OCP_p, measurement, w),'lb',x_lb,'ub',x_ub,'options',options);
% [x_id,fval_ms,flag,outpt,allmins] = run(ms,problem,100);
% 
% 
% % allocate the results
% x_0 = x_id(1);
% QN = x_id(2);
% y_0 = x_id(3);
% QP = x_id(4);
% 
% 
%  Cap = OCV_golden.OCVchg(:,1);
%     if (OCV_golden.OCVchg (end,1)<OCV_golden.OCVchg (1,1)) % Discharge OCV
%         x_sto = -(Cap - Cap(1))/QN + x_0;
%         y_sto = (Cap - Cap(1))/QP + y_0;
%     else  % Charge OCV
%         x_sto = (Cap - Cap(1))/QN + x_0;
%         y_sto = -(Cap - Cap(1))/QP + y_0;
%     end
% 
% 
% x_100 = x_sto(end);
% y_100 = y_sto(end);
% 
% OCP_n_sim = interp1(OCP_n(:,1), OCP_n(:,2), x_sto, 'linear','extrap');
% OCP_p_sim = interp1(OCP_p(:,1), OCP_p(:,2), y_sto, 'linear','extrap');
% 
% OCV_sim = OCP_p_sim - OCP_n_sim;
% TotalLi = QN*x_0 + QP*y_0;
% 
% dataList.QN(k_list)=QN;
% dataList.QP(k_list)=QP;
% dataList.QLi(k_list)=TotalLi;
% dataList.x0(k_list)=x_0;
% dataList.x100(k_list)=x_100;
% dataList.y0(k_list)=y_0;
% dataList.y100(k_list)=y_100;
% dataList.OCVerr(k_list)=norm(OCV_sim - OCV);
% 
% end
% 
% 
% %% FIGURE SET 2
%     CO = lines(8);
%     figure(103);
%     clf
%     plot(Cap, OCV, 'o');
%     hold on;
%     plot(Cap, OCV_sim,'linewidth',2);
%     grid on; box on;
%     xlabel('Capacity (Ah)');ylabel('OCV (V)');
%     legend('Experiment', 'Model');
%     set(gca,'fontsize',16);
%         title([strrep(dataList.testName{k_list},'_','\_'),  ' ', num2str(round(dataList.days(k_list))), ' days', ' C/150 OCP fitting'])
%         
%     figure(105);
%     clf
%     yyaxis left;
%     plot(Cap, OCP_n_sim, 'linewidth',2);
%     ylabel('Anode OCP (V)');
%     hold on;
%     yyaxis right;
%     plot(Cap, OCP_p_sim,'linewidth',2);
%     grid on; box on;
%     xlabel('Capacity (Ah)');ylabel('Cathode OCP (V)');
%     set(gca,'fontsize',20);
% 
%     figure(104);
%     clf
%     plot(Cap, OCV_sim - OCV);
%     grid on; box on;
%     xlabel('Capacity (Ah)');ylabel('OCV error (V)');
%     legend('Simulation - Experiment');
%     set(gca,'fontsize',20);
% 
% 
% 
% 
% 
% 
% % %function model
% % x_0 = x_id(1);
% % QN = x_id(2);
% % y_0 = x_id(3);
% % QP = x_id(4);
% % 
% % function [y] = func_model(  x_0,QN,y_0,QP)
% % 
% %   Cap = OCV_golden.OCVchg(:,1);
% % 
% %     if (OCV_golden.OCVchg (end)< OCV_golden.OCVchg (1) % Discharge OCV
% %         x_sto = -(Cap - Cap(1))/QN + x_0;
% %         y_sto = (Cap - Cap(1))/QP + y_0;
% %     else  % Charge OCV
% %         x_sto = (Cap - Cap(1))/QN + x_0;
% %         y_sto = -(Cap - Cap(1))/QP + y_0;
% %     end
% % 
% % x_100 = x_sto(end);
% % y_100 = y_sto(end);
% % 
% % OCP_n_sim = interp1(OCP_n(:,1), OCP_n(:,2), x_sto, 'linear','extrap');
% % OCP_p_sim = interp1(OCP_p(:,1), OCP_p(:,2), y_sto, 'linear','extrap');
% % 
% % OCV_sim = OCP_p_sim - OCP_n_sim;
% % TotalLi = QN*x_0 + QP*y_0;
% % 
% % dataList.QN(k_list)=QN;
% % dataList.QP(k_list)=QP;
% % dataList.QLi(k_list)=TotalLi;
% % dataList.x0(k_list)=x_0;
% % dataList.x100(k_list)=x_100;
% % dataList.y0(k_list)=y_0;
% % dataList.y100(k_list)=y_100;
% % dataList.OCVerr(k_list)=norm(OCV_sim - OCV);
% % 
% % end
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % %fmincon을 사용하여 최적화 
% % options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 100);
% % [opt_params, rms] = fmincon(@(params) cost_function(params, time_exp, deltaV_exp), ...
% % x_guess, [], [], [], [], [0, 0, 0], [], [], options);
% % 
% % % optimization function cost
% % % obj = min(sum(OCV-OCV_model)^2 ; @stoichiometry_fit_vX.m
% % problem = createOptimProblem('fmincon','x0',x_guess,'objective',@(x) stoichiometry_fit_v3(x, OCP_n, OCP_p, measurement, w),'lb',x_lb,'ub',x_ub,'options',options);
% % [x_id,fval_ms,flag,outpt,allmins] = run(ms,problem,100);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % % initial auess and lower/upperbound
% % x_guess = [0,Q_cell,1,Q_cell];
% % x_lb = [0,Q_cell*0.5,0,Q_cell*0.5];
% % x_ub = [1,Q_cell*2,1,Q_cell*2]; 
% % 
% %  % optimization problem
% %         % obj = min(sum(OCV_model-OCV_data)^2 ;
% % problem = createOptimProblem('fmincon','x0',x_guess,'objective',@(x) stoichiometry_fit_v3(x, OCP_n, OCP_p, measurement, w),'lb',x_lb,'ub',x_ub,'options',options);
% % [x_id,fval_ms,flag,outpt,allmins] = run(ms,problem,100);
% % 
% %     % allocate the results
% % x_0 = x_id(1);
% % QN = x_id(2);
% % y_0 = x_id(3);
% % QP = x_id(4);
% % 
% % Cap = OCV_golden.OCVchg(:,2);
% % 
% %     if (OCV_golden.OCVchg (end,2)<OCV(1,2)) % Discharge OCV
% %         x_sto = -(Cap - Cap(1))/QN + x_0;
% %         y_sto = (Cap - Cap(1))/QP + y_0;
% %     else  % Charge OCV
% %         x_sto = (Cap - Cap(1))/QN + x_0;
% %         y_sto = -(Cap - Cap(1))/QP + y_0;
% %     end
% % 
% % 
% % x_100 = x_sto(end);
% % y_100 = y_sto(end);
% % 
% % OCP_n_sim = interp1(OCP_n(:,1), OCP_n(:,2), x_sto, 'linear','extrap');
% % OCP_p_sim = interp1(OCP_p(:,1), OCP_p(:,2), y_sto, 'linear','extrap');
% % 
% % OCV_sim = OCP_p_sim - OCP_n_sim;
% % TotalLi = QN*x_0 + QP*y_0;
% % 
% % dataList.QN(k_list)=QN;
% % dataList.QP(k_list)=QP;
% % dataList.QLi(k_list)=TotalLi;
% % dataList.x0(k_list)=x_0;
% % dataList.x100(k_list)=x_100;
% % dataList.y0(k_list)=y_0;
% % dataList.y100(k_list)=y_100;
% % dataList.OCVerr(k_list)=norm(OCV_sim - OCV);
% % 
% % end
