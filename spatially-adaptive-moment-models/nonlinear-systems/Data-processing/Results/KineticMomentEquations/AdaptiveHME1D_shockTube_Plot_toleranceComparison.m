clear all
close all
clc
format long

prediction_solver = 'PRICE';
interface_solver = 'PRICE';
relaxation_time = '0p05';

smooth_par = 100;

toldown_s = "0p0001"; tolup_s = "0p00015";
toldown_m = "0p0003"; tolup_m = "0p00045";
toldown_l = "0p001"; tolup_l = "0p0015";

interpolated = 'interp';

name_s = strcat('pred',prediction_solver,'_interface',interface_solver,...
    '_smoothpar',string(smooth_par),'_toldown',toldown_s,'_tolup',tolup_s,...
    '_Kn',relaxation_time,'_',interpolated);type_s='hmeAdaptive';
name_m = strcat('pred',prediction_solver,'_interface',interface_solver,...
    '_smoothpar',string(smooth_par),'_toldown',toldown_m,'_tolup',tolup_m,...
    '_Kn',relaxation_time,'_',interpolated);type_m='hmeAdaptive';
name_l = strcat('pred',prediction_solver,'_interface',interface_solver,...
    '_smoothpar',string(smooth_par),'_toldown',toldown_l,'_tolup',tolup_l,...
    '_Kn',relaxation_time,'_',interpolated);type_l='hmeAdaptive';

[x1,rho1,u1,T1,f31,f41,f51,f61,f71,f81,f91,f101,f111,f121,moments1] = readDataAdaptiveHME1D(name_s,type_s);
[x2,rho2,u2,T2,f32,f42,f52,f62,f72,f82,f92,f102,f112,f122,moments2] = readDataAdaptiveHME1D(name_m,type_m);
[x3,rho3,u3,T3,f33,f43,f53,f63,f73,f83,f93,f103,f113,f123,moments3] = readDataAdaptiveHME1D(name_l,type_l);

blue = [0, 0.4470, 0.7410];
green = [0.4660, 0.6740, 0.1880];
yellow = [0.9290, 0.6940, 0.1250];
red = [0.8500, 0.3250, 0.0980];
brown = [171, 104, 87]./255;

fig = figure;

plot1 = plot(x1,rho1,x2,rho2,x3,rho3,x1,u1,x2,u2,x3,u3,x1,T1,x2,T2,x3,T3)

axis([-0.4,0.2,-0.25,7.25]);
xlabel('$x$','FontSize', 20,'Interpreter','latex')
ylabel('$\rho,u,\theta$','FontSize', 20,'Interpreter','latex')

text(-0.37,6.87,'$\rho$','FontSize',20,'Interpreter','latex')
text(-0.37,1.1,'$\theta$','FontSize',20,'Interpreter','latex')
text(-0.37,0,'$u$','FontSize',20,'Interpreter','latex')

linewidth1 = 3;
set(plot1(1:9),'LineWidth',linewidth1);
set(plot1(1),'LineStyle','-.');set(plot1(4),'LineStyle','-.');set(plot1(7),'LineStyle','-.');
set(plot1(2),'LineStyle','--');set(plot1(5),'LineStyle','--');set(plot1(8),'LineStyle','--');
set(plot1(3),'LineStyle',':');set(plot1(6),'LineStyle',':');set(plot1(9),'LineStyle',':');

set(plot1(1),'Color',red);set(plot1(4),'Color',red);set(plot1(7),'Color',red);
set(plot1(2),'Color',blue);set(plot1(5),'Color',blue);set(plot1(8),'Color',blue);
set(plot1(3),'Color','k');set(plot1(6),'Color','k');set(plot1(9),'Color','k');

leg = legend('$(\epsilon_\mathrm{C},\epsilon_\mathrm{R})_\mathrm{s}$',...
    '$(\epsilon_\mathrm{C},\epsilon_\mathrm{R})_\mathrm{m}$',...
    '$(\epsilon_\mathrm{C},\epsilon_\mathrm{R})_\mathrm{l}$','Location','northeast','Interpreter','latex');
set(leg,'FontSize',12); 

% export_name = strcat('Paper\shockTube','_smoothPar',string(smooth_par),'_Kn',relaxation_time,'toleranceComparison.pdf');
% 
% addpath('C:\Users\rikve\Github\PhD-RUG\spatially-adaptive-moment-models\nonlinear-systems\Data-processing\Results\export_fig\', '-end');
% export_fig(export_name, '-pdf','-transparent');
% 
% close(fig)
