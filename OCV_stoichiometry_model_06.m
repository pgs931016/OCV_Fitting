function [cost,OCV_sim] = OCV_stoichiometry_model_06(x, OCP_n, OCP_p, OCV,w)
%함수의 cost
%함수의 OCV_sim 값


  x_0 = x(1);
  QN = x(2);
  y_0 = x(3);
  QP = x(4);

Cap = OCV(:,1);
    if (OCV(end,2)<OCV (1,2)) % Discharge OCV
        x_sto =-(Cap - Cap(1))/QN + x_0;
        y_sto = (Cap - Cap(1))/QP + y_0;
    else  % Charge OCV
        x_sto = (Cap - Cap(1))/QN + x_0;
        y_sto =-(Cap - Cap(1))/QP + y_0;
    end


% Cap = OCV(:,1);
%     if (OCV_golden.OCVchg (end,2)<OCV_golden.OCVchg (1,2)) % Discharge OCV
%         x_sto = -(Cap - Cap(1))/QN + x_0;
%         y_sto = (Cap - Cap(1))/QP + y_0;
%     else  % Charge OCV
%         x_sto = (Cap - Cap(1))/QN + x_0;
%         y_sto = -(Cap - Cap(1))/QP + y_0;
%     end

    %model
    OCP_n_sim = interp1(OCP_n(:,1), OCP_n(:,2), x_sto, 'linear','extrap');
    OCP_p_sim = interp1(OCP_p(:,1), OCP_p(:,2), y_sto, 'linear','extrap');
    
    OCV_sim = OCP_p_sim - OCP_n_sim;

cost = sqrt(sum((OCV_sim - OCV(:,2)).^2.*w));
%cost는 매개변수의 비용


end


















    %OCV fullcell data
%     cost = sum((OCV_sim - OCV(:,2)).^2);





% %Fit QP, QN and x0, y0 (0% DOD stoichiometry) with the weight matrix w
% clc;clear;close all
% %function cost = stoichiometry_fit_v3(x, OCP_n, OCP_p, OCV, w)   
% %func_cost(x, OCP_n, OCP_p, OCV, w)는 계산된 오차의 제곱의 합을 반환합니다


%% Fit QP, QN and x0, y0 (0% DOD stoichiometry) with the weight matrix w



