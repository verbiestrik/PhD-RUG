clear all
close all
clc
format long

prediction_solver = 'PRICE';
interface_solver = 'PRICE';
relaxation_time = '0p05';
smooth_par = 100;
toldown = '0p0003';
tolup = '0p00045';
order_class_low = 4;
order_class_high = 12;
interpolated = 'interp';

name1 = strcat('shockTube','_order',string(order_class_low),'_Kn',string(relaxation_time));type1='hme';
name2 = strcat('shockTube','_order',string(order_class_high),'_Kn',string(relaxation_time));type2='hme';
name3 = strcat('pred',prediction_solver,'_interface',interface_solver,...
    '_smoothpar',string(smooth_par),'_toldown',toldown,'_tolup',tolup,...
    '_Kn',relaxation_time,'_',interpolated);type3='hmeAdaptive';

[x1,rho1,u1,T1,f31,f41,f51,f61,f71,f81,f91,f101,f111,f121,moments1] = readDataAdaptiveHME1D(name1,type1);
[x2,rho2,u2,T2,f32,f42,f52,f62,f72,f82,f92,f102,f112,f122,moments2] = readDataAdaptiveHME1D(name2,type2);
[x3,rho3,u3,T3,f33,f43,f53,f63,f73,f83,f93,f103,f113,f123,moments3] = readDataAdaptiveHME1D(name3,type3);

plotting = 'rho';
% plotting = 'u';
% plotting = 'T';
% plotting = 'f3';
% plotting = 'f4';
% plotting = 'f5';
% plotting = 'f6';
% plotting = 'f7';
% plotting = 'f8';
% plotting = 'f9';
% plotting = 'f10';
% plotting = 'f11';
% plotting = 'f12';

blue = [0, 0.4470, 0.7410];
green = [0.4660, 0.6740, 0.1880];
yellow = [0.9290, 0.6940, 0.1250];
red = [0.8500, 0.3250, 0.0980];
brown = [171, 104, 87]./255;

fig = figure;

if(strcmp(plotting,'rho'))
    plot1 = plot(x1,rho1,x2,rho2,x3,rho3)
    
    axis([-1.0,1.25,0.75,7.25]);
    xlabel('$x$','FontSize', 20,'Interpreter','latex')
    ylabel('$\rho$','FontSize', 20,'Interpreter','latex')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(1:3),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color','k');
    
    leg = legend(strcat('$\mathrm{HME}_{', num2str(order_class_low), '}$'), ...
                    strcat('$\mathrm{HME}_{', num2str(order_class_high), '}$'), ...
                    'A-HME','Location','northeast', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'u'))
    plot1 = plot(x1,u1,x2,u2,x3,u3)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.0,1.25,-0.024,0.75]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.0,1.25,-0.024,0.75]);
    end

    xlabel('$x$','FontSize', 20,'Interpreter','latex')
    ylabel('$u$','FontSize', 20,'Interpreter','latex')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(1:3),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color','k');
    
    leg = legend(strcat('$\mathrm{HME}_{', num2str(order_class_low), '}$'), ...
                    strcat('$\mathrm{HME}_{', num2str(order_class_high), '}$'), ...
                    'A-HME','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'T'))
    plot1 = plot(x1,T1,x2,T2,x3,T3)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.0,1.25,0.45,1.85]);
    end
    
    if strcmp(relaxation_time,'0p05')
        axis([-1.0,1.25,0.45,1.85]);
    end

    xlabel('$x$','FontSize', 20,'Interpreter','latex')
    ylabel('$u$','FontSize', 20,'Interpreter','latex')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(1:3),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color','k');

    leg = legend(strcat('$\mathrm{HME}_{', num2str(order_class_low), '}$'), ...
                    strcat('$\mathrm{HME}_{', num2str(order_class_high), '}$'), ...
                    'A-HME','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f3'))
    plot1 = plot(x1,f31,x2,f32,x3,f33)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.3,0.2]);
    end

    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.14,0.11]);
    end

    xlabel('x')
    ylabel('f_3')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(1:3),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color','k');
    
    leg = legend(strcat('$\mathrm{HME}_{', num2str(order_class_low), '}$'), ...
                    strcat('$\mathrm{HME}_{', num2str(order_class_high), '}$'), ...
                    'A-HME','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f4'))
    plot1 = plot(x1,f41,x2,f42,x3,f43)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.15,0.15]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.045,0.03]);
    end

    xlabel('x')
    ylabel('f_4')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(1:3),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color','k');
    
    leg = legend(strcat('$\mathrm{HME}_{', num2str(order_class_low), '}$'), ...
                    strcat('$\mathrm{HME}_{', num2str(order_class_high), '}$'), ...
                    'A-HME','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f5'))
    plot1 = plot(x1,f51,x2,f52,x3,f53)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.05,0.06]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.016,0.013]);
    end
    xlabel('x')
    ylabel('f_5')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(1:3),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color','k');
    
    leg = legend(strcat('$\mathrm{HME}_{', num2str(order_class_low), '}$'), ...
                    strcat('$\mathrm{HME}_{', num2str(order_class_high), '}$'), ...
                    'A-HME','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f6'))
    plot1 = plot(x1,f61,x2,f62,x3,f63)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.01,0.02]);
    end

    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.002,0.005]);
    end

    xlabel('x')
    ylabel('f_6')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(1:3),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color','k');
    
    leg = legend(strcat('$\mathrm{HME}_{', num2str(order_class_low), '}$'), ...
                    strcat('$\mathrm{HME}_{', num2str(order_class_high), '}$'), ...
                    'A-HME','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f7'))
    plot1 = plot(x1,f71,x2,f72,x3,f73)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.009,0.009]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.0013,0.0013]);
    end
    xlabel('x')
    ylabel('f_7')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(1:3),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color','k');
    
    leg = legend(strcat('$\mathrm{HME}_{', num2str(order_class_low), '}$'), ...
                    strcat('$\mathrm{HME}_{', num2str(order_class_high), '}$'), ...
                    'A-HME','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f8'))
    plot1 = plot(x1,f81,x2,f82,x3,f83)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.004,0.003]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.0005,0.0005]);
    end
    xlabel('x')
    ylabel('f_8')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(1:3),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color','k');
    
    leg = legend(strcat('$\mathrm{HME}_{', num2str(order_class_low), '}$'), ...
                    strcat('$\mathrm{HME}_{', num2str(order_class_high), '}$'), ...
                    'A-HME','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f9'))
    plot1 = plot(x1,f91,x2,f92,x3,f93)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.0012,0.0013]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.0001,0.00011]);
    end
    xlabel('x')
    ylabel('f_9')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(1:3),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color','k');
    
    leg = legend(strcat('$\mathrm{HME}_{', num2str(order_class_low), '}$'), ...
                    strcat('$\mathrm{HME}_{', num2str(order_class_high), '}$'), ...
                    'A-HME','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f10'))
    plot1 = plot(x1,f101,x2,f102,x3,f103)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.0005,0.00065]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.00002,0.00004]);
    end

    xlabel('x')
    ylabel('f_{10}')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(1:3),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color','k');
    
    leg = legend(strcat('$\mathrm{HME}_{', num2str(order_class_low), '}$'), ...
                    strcat('$\mathrm{HME}_{', num2str(order_class_high), '}$'), ...
                    'A-HME','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f11'))
    plot1 = plot(x1,f111,x2,f112,x3,f113)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.0001,0.000075]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.000015,0.00001]);
    end
    xlabel('x')
    ylabel('f_{11}')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(1:3),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color','k');
    
    leg = legend(strcat('$\mathrm{HME}_{', num2str(order_class_low), '}$'), ...
                    strcat('$\mathrm{HME}_{', num2str(order_class_high), '}$'), ...
                    'A-HME','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f12'))
    plot1 = plot(x1,f121,x2,f122,x3,f123)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.00006,0.00005]);
    end 
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.000004,0.000002]);
    end 
    xlabel('x')
    ylabel('f_{12}')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(1:3),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color','k');
    
    leg = legend(strcat('$\mathrm{HME}_{', num2str(order_class_low), '}$'), ...
                    strcat('$\mathrm{HME}_{', num2str(order_class_high), '}$'), ...
                    'A-HME','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

% export_name = strcat('Paper\shockTube',plotting,'_smoothPar',string(smooth_par),...
%     '_toldown',toldown,'_tolup',tolup,'_Kn',relaxation_time,'_combined.pdf');
% 
% addpath('C:\Users\rikve\Github\PhD-RUG\spatially-adaptive-moment-models\nonlinear-systems\Data-processing\Results\export_fig\', '-end');
% export_fig(export_name, '-pdf','-transparent');
% 
% close(fig)
