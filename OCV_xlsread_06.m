%save_path
%폴더 경로 입력
clc;clear;close all
old_path = "/Users/g.park/Library/CloudStorage/GoogleDrive-gspark@kentech.ac.kr/공유 드라이브/Battery Software Lab/Data/NE_cell/";
new_path = regexprep(old_path, 'Data', 'Processed_data');


disp(new_path)

% 지정된 경로에 directory 생성
if ~exist(new_path, 'dir')
mkdir(new_path)
end

%data load
clc;clear;close all
n_hd = 6;
% excel_file = '/Users/g.park/Library/CloudStorage/GoogleDrive-gspark@kentech.ac.kr/공유 드라이브/Battery Software Lab/Data/Hyundai_dataset/NE_Cell_Characterization_performance.xlsx';
filename = '/Users/g.park/Library/CloudStorage/GoogleDrive-gspark@kentech.ac.kr/공유 드라이브/Battery Software Lab/Data/Hyundai_dataset/NE_Cell_Characterization_performance.xlsx';
sheet = 1;
data_NE= readtable(filename,...
                    'NumHeaderLines',n_hd,'readVariableNames',0);


%Engine
slash = filesep;

dataNE.SOC1 = data_NE.Var1;
dataNE.V1 = data_NE.Var2;
dataNE.V2 = data_NE.Var3;
dataNE.V3 = data_NE.Var4;

dataNE.SOC2 = data_NE.Var6;
dataNE.V4 = data_NE.Var7;
dataNE.V5 = data_NE.Var8;
dataNE.V6 = data_NE.Var9;

%충전,방전,평균 그래프
        figure(1)
        title ('SOC VS. OCV') 
        hold on; box on;
        plot(dataNE.SOC1,dataNE.V1,'-')
        xlabel('SOC (%)')
        ylabel('OCV (V)')
        legend ('0.1')

        figure(2)
        title ('SOC VS. OCV') 
        hold on; box on;
        plot(dataNE.SOC2,dataNE.V4,'-')
        xlabel('SOC (%)')
        ylabel('OCV (V)')
        legend ('0.05')
        % time_lo_bound = 7228.72 / 3600;
        % time_up_bound = 8128.72 / 3600;
        % xlim([time_lo_bound, time_up_bound])
%     end

    % make struct (output format)
 







% save_fullpath = [path slash files(i).name(1:end-4) '.mat'];
% save(save_fullpath,'data')
