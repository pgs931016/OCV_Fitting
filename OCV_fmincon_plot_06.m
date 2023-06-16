
%     % define the weighting 
% w = zeros(size(Cap')); % should have the same length with the data
% w(:)=1; % uniform weighting
% 
%     % prep: multi-start (tolerances)
% %ms = MultiStart('UseParallel',true,'FunctionTolerance',1e-15,'XTolerance',1e-15);
%     
%     % prep: optimset option (tolerances)
% % options = optimoptions(@fmincon,'MaxIterations',5000,'StepTolerance',1e-15,'ConstraintTolerance', 1e-15, 'OptimalityTolerance', 1e-15);
%     
%     % initial guess and lower/upper bounds
% x_guess = [0,Q_cell,1,Q_cell];
% x_lb = [0,Q_cell*0.5,0,Q_cell*0.5];
% x_ub = [1,Q_cell*2,1,Q_cell*2]; 
% 
%     % optimization problem
%         % obj = min(sum(OCV-OCV_model)^2 ; @stoichiometry_fit_vX.m
% % problem = createOptimProblem('fmincon','x0',x_guess,'objective',@(x) OCV_stoichiometry_model(x, OCP_n, OCP_p, w),'lb',x_lb,'ub',x_ub,'options','MaxIter',100);
% % run(problem);
% 
% 
% % fmincon을 사용하여 최적화 수행
%     options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 5000);
%     [x_id] = fmincon(@(x) OCV_stoichiometry_model_06(x, OCP_n, OCP_p, w), ...
%         x_guess, [], [], [], [], [0 0], [1 1],'lb',x_lb,'ub',x_ub, options);
% 
% 
% allocate the results
x_0 = x_id(1);
QN = x_id(2);
y_0 = x_id(3);
QP = x_id(4);

if (OCV(end)<OCV(1)) % Discharge OCV
    x_sto = -(Cap - Cap(1))/QN + x_0;
    y_sto = (Cap - Cap(1))/QP + y_0;
else  % Charge OCV
    x_sto = (Cap - Cap(1))/QN + x_0;
    y_sto = -(Cap - Cap(1))/QP + y_0;
end


x_100 = x_sto(end);
y_100 = y_sto(end);







%    
%     OCV_sim = OCP_p_sim - OCP_n_sim;
%     TotalLi = QN*x_0 + QP*y_0;
% 
% 
% 
% % dataList.QN(k_list)=QN;
% % dataList.QP(k_list)=QP;
% % dataList.QLi(k_list)=TotalLi;
% % dataList.x0(k_list)=x_0;
% % dataList.x100(k_list)=x_100;
% % dataList.y0(k_list)=y_0;
% % dataList.y100(k_list)=y_100;
% % dataList.OCVerr(k_list)=norm(OCV_sim - OCV);
% 
% %% FIGURE SET 
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
