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




[cost_hat, OCV_hat] = OCV1_stoichiometry_model_06(x_id,OCP_n,OCP_p,OCV);

% plot
figure(1)
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
title('SOC vs. OCV (0.01C)');
print('OCV fig1','-dpng','-r300');


% figure(3) %data dv/dq

x = OCV (:,1);
y = OCV (:,2);

for i = 1:(length(x)-1)
    dvdq1(i) = (y(i + 1) - y(i)) / (x(i + 1) - x(i));
end
dvdq1(end+1) = dvdq1(end);
    


figure(3) %data dv/dq

x = OCV (:,1);
y = OCV_hat (:,1);

for i = 1:(length(x) - 1)
    dvdq2(i) = (y(i + 1) - y(i)) / (x(i + 1) - x(i));   
end
dvdq2(end+1) = dvdq2(end);

plot(x,dvdq1,'b-','LineWidth',lw,'MarkerSize',msz); hold on
plot(x,dvdq2,'r-','LineWidth',lw,'MarkerSize',msz);

width = 6;     % Width in inches
height = 6;    % Height in inches
alw = 2;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 2;      % LineWidth
msz = 16;       % MarkerSize

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

legend('FCC data','FCC fit')
xlabel('SOC');
ylabel('dV/dQ /  V (mAh)^-1');
title('SOC vs. dV/dQ');
print('OCV fig3','-dpng','-r300');


figure(40)
greater_than_1_indices = find(OCV(:,1) > 0.1 & OCV(:,1) < 0.9);
greater_than_1_values = OCV(greater_than_1_indices ,1);
w = ones(size(OCV(:,1)));
start_index = greater_than_1_indices(1,1); 
end_index = greater_than_1_indices(end,1);
w(start_index:end_index) = dvdq2(start_index:end_index); 



plot(OCV(:,1),OCV(:,2).*w,'b-','LineWidth',lw,'MarkerSize',msz); hold on


save('ocv1w.mat','w');













