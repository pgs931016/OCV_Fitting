clc; close all; clear

  
load ('OCV_fit.mat')

x_guess = [0.01,1*1.2,0.9,1];
x_lb = [0,1*0.5,0,1*0.5];
x_ub = [1,1*2,1,1*2]; 






%% Initial Guess
[~,OCV_guess] = OCV_stoichiometry_model_06(x_guess,OCP_n,OCP_p,OCV);


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

% plot
figure(1)
width = 6;     % Width in inches
height = 6;    % Height in inches
alw = 2;    % AxesLineWidth
fsz = 20;      % Fontsize
lw = 3.5;      % LineWidth
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


%0.05C NE figure(2)

OCV_chg2 = OCV_chg2(~isnan(OCV_chg2(:,1)),:);
OCV_chg3 = OCV_chg2;

x_guess2 = [0.01,1*1.2,0.9,1];
x_lb = [0,1*0.75,0,1*0.75];
x_ub = [1,1*1.25,1,1*1.25];  
options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 5000);


[~,OCV_guess2] = NE_OCV_stoichiometry_model_06(x_guess2,OCP_n,OCP_p,OCV_chg3);


    fhandle_cost2 = @(x)NE_OCV_stoichiometry_model_06(x, OCP_n, OCP_p,OCV_chg3);
    [x_id2, fval2, exitflag2, output2] = fmincon(fhandle_cost2, ...
        x_guess2, [], [], [], [], x_lb, x_ub, [],options);


[cost2_hat, OCV_chg3_hat] = NE_OCV_stoichiometry_model_06(x_id2,OCP_n,OCP_p,OCV_chg3);

% plot
figure(2)
width = 6;     % Width in inches
height = 6;    % Height in inches
alw = 2;    % AxesLineWidth
fsz = 20;      % Fontsize
lw = 3.5;      % LineWidth
msz = 16;       % MarkerSize


plot(OCV_chg3(:,1),OCV_chg3(:,2),'b-','LineWidth',lw,'MarkerSize',msz); hold on;
plot(OCV_chg3(:,1),OCV_chg3_hat,'r-','LineWidth',lw,'MarkerSize',msz);


pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

legend('NE data','NE fit')
xlabel('SOC');
ylabel('OCV (V)');
title('SOC vs. OCV (0.05C)');
print('OCV fig2','-dpng','-r300');



% figure(3) %data dv/dq

x = OCV (1:10:end,1);
y = OCV (1:10:end,2);

for i = 1:(length(x)-1)
    dvdq1(i) = (y(i + 1) - y(i)) / (x(i + 1) - x(i));
end
    dvdq1(end+1) = dvdq1(end);
    




figure(3) %data dv/dq

x = OCV (1:10:end,1);
y = OCV_hat (1:10:end,1);

for i = 1:(length(x) - 1)
    dvdq2(i) = (y(i + 1) - y(i)) / (x(i + 1) - x(i));   
end
 dvdq2(end+1) = dvdq2(end);

plot(x,dvdq1,'b-','LineWidth',lw,'MarkerSize',msz); hold on
plot(x,dvdq2,'r-','LineWidth',lw,'MarkerSize',msz);

width = 6;     % Width in inches
height = 6;    % Height in inches
alw = 2;    % AxesLineWidth
fsz = 20;      % Fontsize
lw = 3.5;      % LineWidth
msz = 16;       % MarkerSize

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

legend('FCC data','FCC fit')
xlabel('SOC');
ylabel('dV/dQ /  V (mAh)^-1');
title('SOC vs. dV/dQ');
print('OCV fig3','-dpng','-r300');



figure(4)  %NE data dv/dq

x = OCV_chg3 (:,1);
y = OCV_chg3 (:,2);

for i = 1:(length(x) - 1)
    dvdq3(i) = (y(i + 1) - y(i)) / (x(i + 1) - x(i));
    

end
dvdq3(end+1) = dvdq3(end);


figure(4)
x = OCV_chg3 (:,1);
y = OCV_chg3_hat (:,1);

for i = 1:(length(x) - 1)
    dvdq4(i) = (y(i + 1) - y(i)) / (x(i + 1) - x(i));   

end
dvdq4(end+1) = dvdq4(end);



plot(OCV_chg3(1:end,1),dvdq3,'b-','LineWidth',lw,'MarkerSize',msz); hold on
plot(OCV_chg3(1:end,1),dvdq4,'r-','LineWidth',lw,'MarkerSize',msz);




width = 6;     % Width in inches
height = 6;    % Height in inches
alw = 2;    % AxesLineWidth
fsz = 20;      % Fontsize
lw = 3.5;      % LineWidth
msz = 16;       % MarkerSize

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

legend('NE data','NE fit')
xlabel('SOC');
ylabel('dV/dQ /  V (mAh)^-1');
title('SOC vs. dV/dQ');

print('OCV fig4','-dpng','-r300');






figure(5) %OCP_n,OCP_p dvdq
x = OCP_n (1:10:end,1);
y = OCP_n (1:10:end,2);

for i = 1:(length(x) - 1)
    dvdq5(i) = (y(i + 1) - y(i)) / (x(i + 1) - x(i));   


end
dvdq5(end+1) = dvdq5(end);
plot(x,dvdq5,'b-','LineWidth',lw,'MarkerSize',msz);hold on


figure(5)
x = OCP_p (1:10:end,1);
y = OCP_p (1:10:end,2);

for i = 1:(length(x) - 1)
    dvdq6(i) = (y(i + 1) - y(i)) / (x(i + 1) - x(i));   
end
 dvdq6(end+1) = dvdq6(end);

width = 6;     % Width in inches
height = 6;    % Height in inches
alw = 2;    % AxesLineWidth
fsz = 20;      % Fontsize
lw = 3.5;      % LineWidth
msz = 16;       % MarkerSize

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties



plot(x,dvdq6,'r-','LineWidth',lw,'MarkerSize',msz);


legend('OCPn','OCPp')
xlabel('SOC');
ylabel('dV/dQ /  V (mAh)^-1');
title('SOC vs. dV/dQ');





print('OCV fig5','-dpng','-r300');




% dataList.QN(k_list)=QN;
% dataList.QP(k_list)=QP;
% dataList.QLi(k_list)=TotalLi;
% dataList.x0(k_list)=x_0;
% dataList.x100(k_list)=x_100;
% dataList.y0(k_list)=y_0;
% dataList.y100(k_list)=y_100;
% dataList.OCVerr(k_list)=norm(OCV_sim - OCV);


  
