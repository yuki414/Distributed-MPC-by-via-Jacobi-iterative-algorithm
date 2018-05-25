%%
clear all
close all
%% plot setting
set(0, 'defaultAxesFontSize', 24); % 軸のフォントサイズ
set(0, 'defaultTextFontSize', 24); % タイトル、注釈などのフォントサイズ
set(0, 'DefaultLineLineWidth', 1); % 軸の線の太さ
set(0, 'DefaultAxesLineWidth', 1);
%%
fig_plot = 
% 1:CMPC 2:DMPC(p=2) 3:すべて
%%
switch fig_plot
case 1
%% load data file
load data_cmpc
data = data_cx; % data of state
%%
figure(1)
set(1,'Position', [50 100 1000 500])
% hold on
for j = 1:RUNSTEP
    for i = 1:M
        plot(i,data(j,i),'o','Markersize',10)
        hold on
    end
    axis([0, M+1, -2, 2])
    drawnow
    pause(Ts)
    hold off
end
case 2
%% 自由応答
%% load data file
load ないよ
data = data_cx; % data of state
%% plot
figure(2)
set(2,'Position', [50 100 1000 500])
hold on
for j = 1:RUNSTEP
    for i = 1:M
        plot(i,data(j,i),'o','Markersize',10)
        hold on
    end
    axis([0, M+1, -2, 2])
    drawnow
    pause(Ts)
    hold off
end
case 3
%% 自由応答, CMPC, DMPC_2, DMPC_20 
%% load file
load data_cmpc
fig_1 = data_cx; % data of state
load data_free
fig_2 = data_cx;
load data_cmpc
fig_3 = data_cx; % data of state
load data_free
fig_4 = data_cx;
%% plot
figure(10)
set(10,'Position',[1 41 1920 969])
for j = 1:RUNSTEP
    for l = 1:4
        subplot(2,2,l)
        for i = 1:M
            if l == 1
                plot(i,fig_1(j,i),'o','Markersize',10)
            end
            if l == 2
                plot(i,fig_2(j,i),'o','Markersize',10)
            end
            if l == 3
                plot(i,fig_3(j,i),'o','Markersize',10)
            end
            if l == 4
                plot(i,fig_4(j,i),'o','Markersize',10)
            end
            hold on
        end
        axis([0, M+1, -2, 2])
        drawnow
%         pause(Ts)
        hold off
    end
end