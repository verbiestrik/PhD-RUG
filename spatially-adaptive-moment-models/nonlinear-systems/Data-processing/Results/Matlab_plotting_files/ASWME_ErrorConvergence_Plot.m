%% Clear statements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% These are examples of 2D ASWME error plot setups %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
format long

%% Import data

%h_error_ASWME = readmatrix('SmoothExpit_h_Error_Auto.csv');
%vr_error_ASWME = readmatrix('SmoothExpit_vr_Error_Auto.csv');
%vtheta_error_ASWME = readmatrix('SmoothExpit_vtheta_Error_Auto.csv');
%h_error_HASWME = readmatrix('SmoothExpit_h_Error_HASWME_Auto.csv');
%vr_error_HASWME = readmatrix('SmoothExpit_vr_Error_HASWME_Auto.csv');
%vtheta_error_HASWME = readmatrix('SmoothExpit_vtheta_Error_HASWME_Auto.csv');

%h_error_ASWME = readmatrix('SmoothExpit_h_Error_BoundaryExcluded.csv');
%vr_error_ASWME = readmatrix('SmoothExpit_vr_Error_BoundaryExcluded.csv');
%vtheta_error_ASWME = readmatrix('SmoothExpit_vtheta_Error_BoundaryExcluded.csv');
%h_error_HASWME = readmatrix('SmoothExpit_h_Error_HASWME_BoundaryExcluded.csv');
%vr_error_HASWME = readmatrix('SmoothExpit_vr_Error_HASWME_BoundaryExcluded.csv');
%vtheta_error_HASWME = readmatrix('SmoothExpit_vtheta_Error_HASWME_BoundaryExcluded.csv');

h_error_ASWME = readmatrix('ASWMEPaperRevised_SmoothExpit_h_Error_ASWME_BoundaryExcluded_Order4_4000x4_400x200.csv');
vr_error_ASWME = readmatrix('ASWMEPaperRevised_SmoothExpit_vr_Error_ASWME_BoundaryExcluded_Order4_4000x4_400x200.csv');
vtheta_error_ASWME = readmatrix('ASWMEPaperRevised_SmoothExpit_vtheta_Error_ASWME_BoundaryExcluded_Order4_4000x4_400x200.csv');
h_error_HASWME = readmatrix('ASWMEPaperRevised_SmoothExpit_h_Error_HASWME_BoundaryExcluded_Order4_4000x4_400x200.csv');
vr_error_HASWME = readmatrix('ASWMEPaperRevised_SmoothExpit_vr_Error_HASWME_BoundaryExcluded_Order4_4000x4_400x200.csv');
vtheta_error_HASWME = readmatrix('ASWMEPaperRevised_SmoothExpit_vtheta_Error_HASWME_BoundaryExcluded_Order4_4000x4_400x200.csv');

%h_error_ASWME = readmatrix('ASWMEPaperRevised_SmoothExpit_h_Error_SWME_BoundaryExcluded_Order4_2000x4.csv');
%vr_error_ASWME = readmatrix('ASWMEPaperRevised_SmoothExpit_vr_Error_SWME_BoundaryExcluded_Order4_2000x4.csv');
%vtheta_error_ASWME = readmatrix('ASWMEPaperRevised_SmoothExpit_vtheta_Error_SWME_BoundaryExcluded_Order4_2000x4.csv');
%h_error_HASWME = readmatrix('ASWMEPaperRevised_SmoothExpit_h_Error_HASWME_BoundaryExcluded_Order4_2000x4.csv');
%vr_error_HASWME = readmatrix('ASWMEPaperRevised_SmoothExpit_vr_Error_HASWME_BoundaryExcluded_Order4_2000x4.csv');
%vtheta_error_HASWME = readmatrix('ASWMEPaperRevised_SmoothExpit_vtheta_Error_HASWME_BoundaryExcluded_Order4_2000x4.csv');

%% Plot

x_axis = [0,1,2,3,4];
plot1 = semilogy(x_axis,h_error_ASWME,'-x',x_axis,h_error_HASWME,'-.x', ...
    x_axis,vr_error_ASWME,'-x',x_axis,vr_error_HASWME,'-.x', ...
    x_axis,vtheta_error_ASWME,'-x',x_axis,vtheta_error_HASWME,'-.x')

xticks([0 1 2 3 4])

set(plot1(1),'Color','blue')
set(plot1(1),'LineWidth',2)
set(plot1(2),'Color','red')
set(plot1(2),'LineWidth',2)
set(plot1(3),'Color','blue')
set(plot1(3),'LineWidth',2)
set(plot1(4),'Color','red')
set(plot1(4),'LineWidth',2)
set(plot1(5),'Color','blue')
set(plot1(5),'LineWidth',2)
set(plot1(6),'Color','red')
set(plot1(6),'LineWidth',2)

axis([-0.5,4.5,0.0005,1]);

xlabel('$N$','Interpreter','latex','FontSize',18)
ylabel('Relative error','FontSize',18)

%leg = legend('ASWME h','HASWME h', 'ASWME vr_m', 'HASWME vr_m','ASWME v\theta_m','HASWME v\theta_m');
leg = legend('ASWME','HASWME');
set(leg,'FontSize',14); 

text(4.1,0.009,'v_r','FontSize',14);
text(4.1,0.0038,'v_{\theta}','FontSize',14);
text(4.1,0.001,'h','FontSize',14);

%% Export figure
 %addpath('C:\Users\rikve\Prime\results\ASWME\export_fig\', '-end')
 %cd export_fig
 %export_fig('models_unstable_tend0p05_a1.pdf', '-pdf','-transparent');
 %export_fig('ASWMEPaperRevised_ErrorConvergence_SmoothExpit_Time1.0_BoundaryExcluded_Order4_RedBlue_4000x4_400x200_Publication.pdf', '-pdf','-transparent');

