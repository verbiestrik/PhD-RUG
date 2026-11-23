clear all
close all
clc
format long

relaxation_time = "0p05";

name1 = 'predPRICE_interfacePRICE_smoothpar100_toldown0p0003_tolup0p00045_Kn0p05_interp';type1='hmeAdaptive';
name2 = 'predPRICE_interfacePRICE_smoothpar200_toldown0p0003_tolup0p00045_Kn0p05_time1p0';type2='hmeAdaptive';
name3 = 'predPRICE_interfacePRICE_smoothpar200_toldown0p0003_tolup0p00045_Kn0p05_time3p0';type3='hmeAdaptive';

[x1,rho1,u1,T1,f31,f41,f51,f61,f71,f81,f91,f101,f111,f121,moments1] = readDataAdaptiveHME1D(name1,type1);
[x2,rho2,u2,T2,f32,f42,f52,f62,f72,f82,f92,f102,f112,f122,moments2] = readDataAdaptiveHME1D(name2,type2);
[x3,rho3,u3,T3,f33,f43,f53,f63,f73,f83,f93,f103,f113,f123,moments3] = readDataAdaptiveHME1D(name3,type3);

% plotting = 'rho';
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
plotting = 'rho_and_f3';

blue = [0, 0.4470, 0.7410];
green = [0.4660, 0.6740, 0.1880];
yellow = [0.9290, 0.6940, 0.1250];
red = [0.8500, 0.3250, 0.0980];
brown = [171, 104, 87]./255;

fig = figure;

if(strcmp(plotting,'rho'))
    plot1 = plot(x2,rho2,x3,rho3)
    
    axis([-3.0,3.0,0.8,7.2]);
    xlabel('$x$','FontSize', 20,'Interpreter','latex')
    ylabel('$\rho$','FontSize', 20,'Interpreter','latex')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1),'LineStyle','-');
    set(plot1(2),'LineStyle','-.');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);

    leg = legend('$t = 1.0$','$t = 3.0$','Location','northeast', 'Interpreter', 'latex');
    set(leg,'FontSize',14); 
end

if(strcmp(plotting,'u'))
    plot1 = plot(x2,u2,x3,u3)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.0,1.25,-0.024,0.75]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-3.0,3.0,-0.024,0.75]);
    end

    xlabel('$x$','FontSize', 20,'Interpreter','latex')
    ylabel('$u$','FontSize', 20,'Interpreter','latex')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('$t = 1.0$','$t = 3.0$','Location','northeast', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'T'))
    plot1 = plot(x2,T2,x3,T3)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.0,1.25,0.45,1.85]);
    end
    
    if strcmp(relaxation_time,'0p05')
        axis([-3.0,3.0,0.45,1.85]);
    end

    xlabel('$x$','FontSize', 20,'Interpreter','latex')
    ylabel('$u$','FontSize', 20,'Interpreter','latex')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('$t = 1.0$','$t = 3.0$','Location','northeast', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f3'))
    plot1 = plot(x2,f32,x3,f33)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.3,0.2]);
    end

    if strcmp(relaxation_time,'0p05')
        axis([0,1.5,-0.085,0.03]);
    end

    xlabel('x')
    ylabel('f_3')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('$t = 1.0$','$t = 3.0$','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f4'))
    plot1 = plot(x1,f41,x2,f42)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.15,0.15]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.045,0.03]);
    end

    xlabel('x')
    ylabel('f_4')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('interface PRICE','interface Osher','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f5'))
    plot1 = plot(x1,f51,x2,f52)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.05,0.06]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.016,0.013]);
    end
    xlabel('x')
    ylabel('f_5')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('interface PRICE','interface Osher','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f6'))
    plot1 = plot(x1,f61,x2,f62)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.01,0.02]);
    end

    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.002,0.005]);
    end

    xlabel('x')
    ylabel('f_6')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('interface PRICE','interface Osher','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f7'))
    plot1 = plot(x1,f71,x2,f72)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.009,0.009]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.0013,0.0013]);
    end
    xlabel('x')
    ylabel('f_7')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('interface PRICE','interface Osher','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f8'))
    plot1 = plot(x1,f81,x2,f82)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.004,0.003]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.0005,0.0005]);
    end
    xlabel('x')
    ylabel('f_8')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('interface PRICE','interface Osher','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f9'))
    plot1 = plot(x1,f91,x2,f92)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.0012,0.0013]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.0001,0.00011]);
    end
    xlabel('x')
    ylabel('f_9')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('interface PRICE','interface Osher','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f10'))
    plot1 = plot(x1,f101,x2,f102)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.0005,0.00065]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.00002,0.00004]);
    end

    xlabel('x')
    ylabel('f_{10}')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('interface PRICE','interface Osher','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f11'))
    plot1 = plot(x1,f111,x2,f112)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.0001,0.000075]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.000015,0.00001]);
    end
    xlabel('x')
    ylabel('f_{11}')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('interface PRICE','interface Osher','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f12'))
    plot1 = plot(x2,f122,x3,f123)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.00006,0.00005]);
    end 
    if strcmp(relaxation_time,'0p05')
        axis([-3.0,3.0,-0.000004,0.000002]);
    end 
    xlabel('x')
    ylabel('f_{12}')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('interface PRICE','interface Osher','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'rho_and_f3'))
    figure('Color','w');
    ax = axes;
    hold(ax,'on');
    
    % Colors
    cF3  = [0.85 0.33 0.10];   % red/orange
    cRho = [0.00 0.45 0.74];   % blue
    
    % --- Left axis: f3 ---
    yyaxis left
    hF3_1 = plot(x2,f32,'-','LineWidth',2,'Color',cF3);
    hF3_2 = plot(x3,f33,'--','LineWidth',2,'Color',cF3);
    ylabel('$f_3$','Interpreter','latex','FontSize',18)
    if strcmp(relaxation_time,'0p5')
        ylim([-0.3 0.2])
    elseif strcmp(relaxation_time,'0p05')
        ylim([-0.085 0.03])
    end
    ax.YAxis(1).Color = cF3;
    
    % --- Right axis: rho ---
    yyaxis right
    hRho_1 = plot(x2,rho2,'-','LineWidth',2,'Color',cRho);
    hRho_2 = plot(x3,rho3,'--','LineWidth',2,'Color',cRho);
    ylabel('$\rho$','Interpreter','latex','FontSize',18)
    ylim([0.8 7.2])
    ax.YAxis(2).Color = cRho;
    
    % Shared x-axis and formatting
    xlabel('$x$','Interpreter','latex','FontSize',18)
    xlim([-3 3])
    set(ax,'FontSize',14,'LineWidth',1,'TickLabelInterpreter','latex')
    grid on
    box on
    
    legend([hF3_1 hF3_2 hRho_1 hRho_2], ...
        {'$f_3,\ t=1.0$','$f_3,\ t=3.0$','$\rho,\ t=1.0$','$\rho,\ t=3.0$'}, ...
        'Interpreter','latex','FontSize',12,'Location','northwest');
end

% export_name = 'interface_viscosity_comparison_shockTube_rho.pdf';
% 
% addpath('C:\Users\rikve\Github\PhD-RUG\spatially-adaptive-moment-models\nonlinear-systems\Data-processing\Results\export_fig\', '-end');
% export_fig(export_name, '-pdf','-transparent');
% 
% close(fig)
