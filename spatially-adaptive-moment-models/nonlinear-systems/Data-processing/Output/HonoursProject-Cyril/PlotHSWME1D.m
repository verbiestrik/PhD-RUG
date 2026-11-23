clear all
close all
clc
format long

basefolder = fullfile(fileparts(mfilename('fullpath')),'error_trainingData');

% initialHeight = 'hDamLarge';
% initialVelocity = 'uConst';
% velocityMagnitude = 'uSlow';
% lambda = 0.5;
% viscosity = 0.5;
% scheme = 'Roe';

initialHeight = 'hSmoothMid';
initialVelocity = 'uVarZnotvarX';
velocityMagnitude = 'uFast';
lambda = 0.05;
viscosity = 0.05;
scheme = 'Roe';

folder = strcat(initialHeight,'_',initialVelocity,'_',velocityMagnitude);

file0 = strcat('lambda',num2str(lambda),'_','viscosity',num2str(viscosity),'_','order0','_','FVM',scheme);
file1 = strcat('lambda',num2str(lambda),'_','viscosity',num2str(viscosity),'_','order1','_','FVM',scheme);
file2 = strcat('lambda',num2str(lambda),'_','viscosity',num2str(viscosity),'_','order2','_','FVM',scheme);
file3 = strcat('lambda',num2str(lambda),'_','viscosity',num2str(viscosity),'_','order3','_','FVM',scheme);
file4 = strcat('lambda',num2str(lambda),'_','viscosity',num2str(viscosity),'_','order4','_','FVM',scheme);
file5 = strcat('lambda',num2str(lambda),'_','viscosity',num2str(viscosity),'_','order5','_','FVM',scheme);
file6 = strcat('lambda',num2str(lambda),'_','viscosity',num2str(viscosity),'_','order6','_','FVM',scheme);
fileRef = strcat('lambda',num2str(lambda),'_','viscosity',num2str(viscosity),'_','ref.csv');

% Load data
name0 = fullfile(basefolder,folder,file0); type0 = 'swe';
name1 = fullfile(basefolder,folder,file1); type1 = 'swme';
name2 = fullfile(basefolder,folder,file2); type2 = 'swme';
name3 = fullfile(basefolder,folder,file3); type3 = 'swme';
name4 = fullfile(basefolder,folder,file4); type4 = 'swme';
name5 = fullfile(basefolder,folder,file5); type5 = 'swme';
name6 = fullfile(basefolder,folder,file6); type6 = 'swme';
nameRef = fullfile(basefolder,folder,fileRef);

ref = load(nameRef);

x_ref = ref(:,1); h_ref = ref(:,2); u_ref = ref(:,3); alpha1_ref = 3.*ref(:,4); 
alpha2_ref = 5.*ref(:,5); alpha3_ref = 7.*ref(:,6); alpha4_ref = 9.*ref(:,7); 
alpha5_ref = 11.*ref(:,8); alpha6_ref = 13.*ref(:,9);

[x0,h0,u0,alpha10,alpha20,alpha30,alpha40,alpha50,alpha60] = readDataHSWME(name0,type0);
[x1,h1,u1,alpha11,alpha21,alpha31,alpha41,alpha51,alpha61] = readDataHSWME(name1,type1);
[x2,h2,u2,alpha12,alpha22,alpha32,alpha42,alpha52,alpha62] = readDataHSWME(name2,type2);
[x3,h3,u3,alpha13,alpha23,alpha33,alpha43,alpha53,alpha63] = readDataHSWME(name3,type3);
[x4,h4,u4,alpha14,alpha24,alpha34,alpha44,alpha54,alpha64] = readDataHSWME(name4,type4);
[x5,h5,u5,alpha15,alpha25,alpha35,alpha45,alpha55,alpha65] = readDataHSWME(name5,type5);
[x6,h6,u6,alpha16,alpha26,alpha36,alpha46,alpha56,alpha66] = readDataHSWME(name6,type6);

plotting = 'h';
% plotting = 'u';
% plotting = 'alpha1';
% plotting = 'alpha2';
% plotting = 'alpha3';
% plotting = 'alpha4';
% plotting = 'alpha5';
% plotting = 'h';

brown = [171, 104, 87]./255;
blue = [0, 0.4470, 0.7410];
green = [0.4660, 0.6740, 0.1880];
yellow = [0.9290, 0.6940, 0.1250];
red = [0.8500, 0.3250, 0.0980];

if(strcmp(plotting,'h'))
    plot1 = plot(x0,h0,x1,h1,x2,h2,x3,h3,x4,h4,x5,h5,x6,h6,x_ref,h_ref,'*')

    axis([-1,1,0.8,3.2]);

    xlabel('x')
    ylabel('h')
    
    linewidth1 = 2;
    set(plot1(1:7),'LineWidth',linewidth1);
    set(plot1(1:7),'LineStyle','-');
    
    set(plot1(1),'Color',yellow);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color',green);
    set(plot1(4),'Color',red);
    set(plot1(5),'Color','b');
    set(plot1(6),'Color','cyan');
    set(plot1(7),'Color','magenta');
    
    set(plot1(8),'MarkerSize',3)

    leg = legend('M=0','M=1','M=2','M=3','M=4','M=5','M=6','ref','Location','southwest');
    set(leg,'FontSize',7); 
end
if(strcmp(plotting,'u'))
    plot1 = plot(x0,u0,x1,u1,x2,u2,x3,u3,x4,u4,x5,u5,x6,u6,x_ref,u_ref,'*')

    axis([-1,1,0,2.0]);

    xlabel('x')
    ylabel('u')
    
    linewidth1 = 2;
    set(plot1(1:7),'LineWidth',linewidth1);
    set(plot1(1:7),'LineStyle','-');
    
    set(plot1(1),'Color',yellow);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color',green);
    set(plot1(4),'Color',red);
    set(plot1(5),'Color','b');
    set(plot1(6),'Color','cyan');
    set(plot1(7),'Color','magenta');
    
    set(plot1(8),'MarkerSize',3)

    leg = legend('M=0','M=1','M=2','M=3','M=4','M=5','M=6','ref','Location','northwest');
    set(leg,'FontSize',7); 
end
if(strcmp(plotting,'alpha1'))
    % plot1 = plot(x1,alpha11,x2,alpha12,x3,alpha13,x4,alpha14,x_ref2000x200,alpha1_ref2000x200,'*')
    plot1 = plot(x1,alpha11,x2,alpha12,x3,alpha13,x4,alpha14)
    
    axis([-20,20,-0.335,0.2]);
    % axis([-3.75,3.25,-0.135,0]);
    xlabel('x')
    ylabel('\alpha_1')
    
    linewidth1 = 2;
    set(plot1(1:4),'LineWidth',linewidth1);

    grey = [0.4,0.4,0.4];
    set(plot1(1),'LineStyle','-');
    set(plot1(2),'LineStyle','-.');
    set(plot1(3),'LineStyle','--');
    set(plot1(4),'LineStyle',':');
    set(plot1(1),'Color','r');
    set(plot1(2),'Color','blue');
    set(plot1(3),'Color','black');
    set(plot1(4),'Color','cyan');
    
    % set(plot1(5),'Color',brown);
    % set(plot1(5),'MarkerSize',5);

    leg = legend('SWME_2','SWME_5','A-SWME-PBC','A-SWME-CIF','Location','northwest');
    % leg = legend('SWME_2','SWME_5','A-SWME-PBC','A-SWME-CIF','Reference','Location','southeast');
    set(leg,'FontSize',12); 
end
if(strcmp(plotting,'alpha2'))
    % plot1 = plot(x1,alpha21,x2,alpha22,x3,alpha23,x4,alpha24,x_ref2000x200,alpha2_ref2000x200,'*')
    plot1 = plot(x1,alpha21,x2,alpha22,x3,alpha23,x4,alpha24)
    
    axis([-20.,20,-0.09,0.015]);
    % axis([-3.0,3.0,-0.075,0.015]);
    xlabel('x')
    ylabel('\alpha_2')
    
    linewidth1 = 2;
    set(plot1(1:4),'LineWidth',linewidth1);

    grey = [0.4,0.4,0.4];
    set(plot1(1),'LineStyle','-');
    set(plot1(2),'LineStyle','-.');
    set(plot1(3),'LineStyle','--');
    set(plot1(4),'LineStyle',':');
    set(plot1(1),'Color','r');
    set(plot1(2),'Color','blue');
    set(plot1(3),'Color','black');
    set(plot1(4),'Color','cyan');
    
    % set(plot1(5),'Color',brown);
    % set(plot1(5),'MarkerSize',3);

    leg = legend('SWME_2','SWME_5','A-SWME-PBC','A-SWME-CIF','Location','southwest');
    % leg = legend('SWME_2','SWME_5','A-SWME-PBC','A-SWME-CIF','Reference','Location','southwest');
    set(leg,'FontSize',12); 
end
if(strcmp(plotting,'alpha3'))
    plot1 = plot(x1,alpha31,x2,alpha32,x3,alpha33,x4,alpha34)

    axis([-20.,20,-0.075,0.075]);
    xlabel('x')
    ylabel('\alpha_3')
    
    linewidth1 = 2;
    set(plot1(1:4),'LineWidth',linewidth1);

    grey = [0.4,0.4,0.4];
    set(plot1(1:4),'LineStyle','-.');
    set(plot1(1),'Color','r');
    set(plot1(2),'Color','blue');
    set(plot1(3),'Color','black');
    set(plot1(4),'Color','green');
    
    leg = legend('Low order','High order','Adaptive: noncons','Adaptive: cons','Location','northwest');
    set(leg,'FontSize',12); 
end
if(strcmp(plotting,'alpha4'))
    plot1 = plot(x1,alpha41,x2,alpha42,x3,alpha43,x4,alpha44)
%     legend(name1,name2,name3,name4)

    axis([-0.4,0.4,-0.8,0.0]);
    xlabel('x')
    ylabel('\alpha_4')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(2:4),'LineStyle','.');
    set(plot1(2:4),'MarkerSize',linewidth1+12);
    set(plot1(4),'LineStyle','.');
    set(plot1(4),'MarkerSize',linewidth1+12);

    grey = [0.4,0.4,0.4];
%     set(plot1(1),'Color','b');
    set(plot1(2),'Color','k');
%     set(plot1(3),'Color','b');
    set(plot1(3),'Color','k');
    set(plot1(4),'Color','k');
    
%     leg = legend('SWME','HSWME','\betaHSWME','reference','Location','southwest');
    leg = legend('SWME','reference','Location','southwest');
    set(leg,'FontSize',12); 
end
if(strcmp(plotting,'alpha5'))
    plot1 = plot(x1,alpha51,x2,alpha52,x3,alpha53,x4,alpha54)
%     legend(name1,name2,name3,name4)

    axis([-0.4,0.4,-0.8,0.0]);
    xlabel('x')
    ylabel('\alpha_5')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(2:4),'LineStyle','.');
    set(plot1(2:4),'MarkerSize',linewidth1+12);
    set(plot1(4),'LineStyle','.');
    set(plot1(4),'MarkerSize',linewidth1+12);

    grey = [0.4,0.4,0.4];
%     set(plot1(1),'Color','b');
    set(plot1(2),'Color','k');
%     set(plot1(3),'Color','b');
    set(plot1(3),'Color','k');
    set(plot1(4),'Color','k');
    
%     leg = legend('SWME','HSWME','\betaHSWME','reference','Location','southwest');
    leg = legend('SWME','reference','Location','southwest');
    set(leg,'FontSize',12); 
end


% addpath('C:\Users\rikve\Github\PhD-RUG\spatially-adaptive-moment-models\nonlinear-systems\Data-processing\Results\export_fig\', '-end');
% export_fig('AdaptiveSWME_Paper\damBreak-and-smooth_linear_lambda0.1_nu1.0_t5_orders_new.pdf', '-pdf','-transparent');