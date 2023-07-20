clc;clear;close all

  
load ('OCV_fit.mat')

x_guess = [0.01,1*1.2,0.9,1];
x_lb = [0,1*0.5,0,1*0.5];
x_ub = [1,1*2,1,1*2]; 




%% Initial Guess
[~,OCV_guess] = OCV_stoichiometry_model_06(x_guess,OCP_n,OCP_p,OCV);


% fmincon을 사용하여 최적화 수행
  
options = optimoptions(@fmincon,'MaxIterations',5000,'StepTolerance',1e-15,'ConstraintTolerance', 1e-15, 'OptimalityTolerance', 1e-15);
   
% problem = createOptimProblem('fmincon', 'objective', @(x) OCV1_stoichiometry_model_06(x_id,OCP_n,OCP_p,OCV), ...
%             'x0', x_guess, 'lb', [0,1*0.5,0,1*0.5], 'ub', [1,1*2,1,1*2] , 'options', options);
%         ms = MultiStart('Display', 'iter');
%     
%         [x_id, fval, exitflag, output] = run(ms, problem, 20); 
 
fhandle_cost = @(x)OCV_stoichiometry_model_06(x, OCP_n, OCP_p, OCV);
    [x_id, fval, exitflag, output] = fmincon(fhandle_cost, ...
        x_guess, [], [], [], [], x_lb, x_ub, [],options);




[cost_hat, OCV_hat] = OCV_stoichiometry_model_06(x_id,OCP_n,OCP_p,OCV);

% plot
figure('position', [0 0 1000 400])
width = 6;     % Width in inches
height = 6;    % Height in inches
alw = 2;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 2;      % LineWidth
msz = 16;       % MarkerSize
 
plot(OCV(:,1),OCV(:,2),'b-','LineWidth',lw,'MarkerSize',msz); hold on
plot(OCV(:,1),OCV_hat,'r-','LineWidth',lw,'MarkerSize',msz);



pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

legend('FCC data','FCC fit')
xlabel('SOC');
ylabel('OCV (V)');
title('OCV1 (0.01C)');
print('OCV fig1','-dpng','-r300');



figure(3) %data dv/dq

window_size = 30;

x = OCV (:,1);
y = OCV (:,2);


x_values = [];
for i = 1:(length(x)-1)
    dvdq1(i) = (y(i + 1) - y(i)) / (x(i + 1) - x(i));
    x_values = [x_values; x(i)];
end
%dvdq1(end+1) = dvdq1(end);
   



x = OCV (:,1);
y = OCV_hat (:,1);


x_values2 = [];
for i = 1:(length(x) - 1)
    dvdq2(i) = (y(i + 1) - y(i)) / (x(i + 1) - x(i)); 
    x_values2 = [x_values2; x(i)];
end
% dvdq2(end+1) = dvdq2(end);


% plot(x(1:end),dvdq1(1:end),'b-','LineWidth',lw,'MarkerSize',msz); hold on
% plot(x(1:end),dvdq2(1:end),'r-','LineWidth',lw,'MarkerSize',msz);

dvdq1_moving_avg = movmean(dvdq1(1:end), window_size);
x_values_moving_avg = movmean(x_values, window_size);

% dvdq6에 이동 평균 적용
dvdq2_moving_avg = movmean(dvdq2(1:end), window_size);
x_values2_moving_avg = movmean(x_values2, window_size);

% 플롯 그리기
plot(x_values_moving_avg, dvdq1_moving_avg, 'b-', 'LineWidth', lw, 'MarkerSize', msz); hold on
plot(x_values2_moving_avg, dvdq2_moving_avg, 'r-', 'LineWidth', lw, 'MarkerSize', msz);



% width = 6;     % Width in inches
% height = 6;    % Height in inches
% alw = 2;    % AxesLineWidth
% fsz = 11;      % Fontsize
% lw = 2;      % LineWidth
% msz = 16;       % MarkerSize
% 
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
% set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
% 
% legend('FCC data','FCC fit')
% xlabel('SOC');
% ylabel('dV/dQ /  V (mAh)^-1');
% title('SOC vs. dV/dQ');
% ylim([-1.5 10])
% print('OCV fig3','-dpng','-r300');
% 
% 
% figure(40)
% w = ones(size(OCV(:,1)));
% greater_than_1_indices = find(OCV(:,1) > 0.1 & OCV(:,1) < 0.9);
% greater_than_2_indices = find(OCV(:,1) > 0.65 & OCV(:,1) < 0.8);
% 
% greater_than_1_values = OCV(greater_than_1_indices ,1);
% greater_than_2_values = OCV(greater_than_2_indices ,1);
% 
% 
% start_index = greater_than_1_indices(1,1); 
% end_index = greater_than_1_indices(end,1);
% start_index2 = greater_than_2_indices(1,1);
% end_index2 = greater_than_2_indices(end,1);
% 
% w(start_index:end_index) = dvdq1_moving_avg(start_index:end_index); 
% w(start_index2:end_index2) = dvdq1_moving_avg(start_index2:end_index2) +2
% 
% 
% % greater_than_1_indices = find(OCV(:,1) > 0.1 & OCV(:,1) < 0.9);
% % greater_than_1_values = OCV(greater_than_1_indices ,1);
% % w = ones(size(OCV(:,1)));
% % start_index = greater_than_1_indices(1,1); 
% % end_index = greater_than_1_indices(end,1);
% % w(start_index:end_index) = dvdq1_moving_avg(start_index:end_index); 
% 
% 
% 
% plot(OCV(:,1),OCV(:,2).*w,'b-','LineWidth',lw,'MarkerSize',msz); hold on
% 
% 
% save('ocv1w.mat','w');



figure('position', [0 0 500 400] );

width = 6;     % Width in inches
height = 6;    % Height in inches
alw = 2;    % AxesLineWidth
fsz = 20;      % Fontsize
lw = 2;      % LineWidth
msz = 16;       % MarkerSize


pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties



% 첫 번째 그래프
subplot(2, 1, 1);

plot(OCV(:,1),OCV(:,2),'b-','LineWidth',lw,'MarkerSize',msz);
hold on;
plot(OCV(:,1),OCV_hat,'r-','LineWidth',lw,'MarkerSize',msz);
% xlabel('SOC');
ylabel('OCV (V)');
title('OCV1 (0.01C)');
% yyaxis right;
% ax = gca;  % 현재 축 객체 가져오기
% ax.YColor = 'k';  % 검정색으로 설정
% ylabel('Weight')
% plot(OCV(1:end,1),w(1:end),'-g','LineWidth',lw,'MarkerSize',msz);
% ylim([0 7])
% 
legend('FCC data','FCC fit','Location', 'none', 'Position', [0.2 0.85 0.1 0.05],'FontSize', 6);

% 두 번째 그래프
subplot(2, 1, 2);

% subtightplot(2, 1, 2, [0.1 0.05], [0.05 0.1], [0.1 0.05]);
plot(x_values_moving_avg, dvdq1_moving_avg, 'b-', 'LineWidth', lw, 'MarkerSize', msz);
hold on;
plot(x_values2_moving_avg, dvdq2_moving_avg, 'r-', 'LineWidth', lw, 'MarkerSize', msz);
xlabel('SOC');
ylabel('dV/dQ /  V (mAh)^-1');
% title('SOC vs. dV/dQ');
ylim([0 2.5]);
print('OCV fig66','-dpng','-r300');


% % % 그래프 간격 조정
% % linkaxes('x');














