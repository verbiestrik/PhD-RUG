clear all
close all
clc
format long

relaxation_time = "0p05";

% name1 = 'predPRICE_interfacePRICE_smoothpar200_toldown0p0003_tolup0p00045_Kn0p5';type1='hmeAdaptive';
% name2 = 'predPRICE_interfaceOsher_smoothpar200_toldown0p0003_tolup0p00045_Kn0p5';type2='hmeAdaptive';
% name1 = 'predPRICE_interfacePRICE_smoothpar200_toldown0p0003_tolup0p00045_Kn0p5';type1='hmeAdaptive';
% name2 = 'predPRICE_interfaceOsher_smoothpar200_toldown0p0003_tolup0p00045_Kn0p5';type2='hmeAdaptive';

% name1 = 'predPRICE_interfacePRICE_smoothpar200_toldown0p0003_tolup0p00045_Kn0p05';type1='hmeAdaptive';
% name2 = 'predPRICE_interfaceOsher_smoothpar200_toldown0p0003_tolup0p00045_Kn0p05';type2='hmeAdaptive';
% name1 = 'predPRICE_interfacePRICE_smoothpar200_toldown0p001_tolup0p0015_Kn0p05';type1='hmeAdaptive';
% name2 = 'predPRICE_interfaceOsher_smoothpar200_toldown0p001_tolup0p0015_Kn0p05';type2='hmeAdaptive';

name1 = 'predPRICE_interfacePRICE_smoothpar200_toldown0p0003_tolup0p00045_Kn0p05_time3p0';type1='hmeAdaptive';
name2 = 'predPRICE_interfaceOsher_smoothpar200_toldown0p0003_tolup0p00045_Kn0p05_time3p0';type2='hmeAdaptive';

[x1,rho1,u1,T1,f31,f41,f51,f61,f71,f81,f91,f101,f111,f121,moments1] = readDataAdaptiveHME1D(name1,type1);
[x2,rho2,u2,T2,f32,f42,f52,f62,f72,f82,f92,f102,f112,f122,moments2] = readDataAdaptiveHME1D(name2,type2);

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
    plot1 = plot(x1,rho1,x2,rho2)
    
    axis([-0.15,-0.01,4.7,5.025]);
    xlabel('$x$','FontSize', 20,'Interpreter','latex')
    ylabel('$\rho$','FontSize', 20,'Interpreter','latex')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1),'LineStyle','-');
    set(plot1(2),'LineStyle','-.');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('PRICE','Osher-Solomon','Location','northeast', 'Interpreter', 'latex');
    set(leg,'FontSize',14); 
end

if(strcmp(plotting,'u'))
    plot1 = plot(x1,u1,x2,u2)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.0,1.25,-0.024,0.75]);
    end
    if strcmp(relaxation_time,'0p05')
        axis([-1.0,1.25,-0.024,0.75]);
    end

    xlabel('$x$','FontSize', 20,'Interpreter','latex')
    ylabel('$u$','FontSize', 20,'Interpreter','latex')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('interface PRICE','interface Osher','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'T'))
    plot1 = plot(x1,T1,x2,T2)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.0,1.25,0.45,1.85]);
    end
    
    if strcmp(relaxation_time,'0p05')
        axis([-1.0,1.25,0.45,1.85]);
    end

    xlabel('$x$','FontSize', 20,'Interpreter','latex')
    ylabel('$u$','FontSize', 20,'Interpreter','latex')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('interface PRICE','interface Osher','Location','northwest', 'Interpreter', 'latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'f3'))
    plot1 = plot(x1,f31,x2,f32)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.3,0.2]);
    end

    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.14,0.11]);
    end

    xlabel('x')
    ylabel('f_3')
    
    linewidth1 = 2;
    set(plot1(1:2),'LineWidth',linewidth1);
    set(plot1(1:2),'LineStyle','-');

    set(plot1(1),'Color',red);
    set(plot1(2),'Color',blue);
    
    leg = legend('interface PRICE','interface Osher','Location','northwest', 'Interpreter', 'latex');
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
    plot1 = plot(x1,f121,x2,f122)
    
    if strcmp(relaxation_time,'0p5')
        axis([-1.5,1.5,-0.00006,0.00005]);
    end 
    if strcmp(relaxation_time,'0p05')
        axis([-1.5,1.5,-0.000004,0.000002]);
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
    % figure('Color','w');
    % ax = axes;
    % hold(ax,'on');
    % set(ax,'FontSize',14,'LineWidth',1,'TickLabelInterpreter','latex');
    % 
    % % Colors for the two methods
    % cPRICE = red;
    % cOSHER = blue;
    % 
    % % Left axis: rho
    % yyaxis left
    % hRho1 = plot(x1,rho1,'-','Color',cPRICE,'LineWidth',2);
    % hRho2 = plot(x2,rho2,'-.','Color',cOSHER,'LineWidth',2);
    % ylabel('$\rho$','Interpreter','latex','FontSize',18)
    % ylim([4.7 5.025])
    % text(-0.02,4.785,'$\rho$','Interpreter','latex','FontSize',16, ...
    %     'FontWeight','bold','Color','k','VerticalAlignment','top');
    % 
    % % Right axis: f3
    % yyaxis right
    % hF31 = plot(x1,f31,'-','Color',cPRICE,'LineWidth',2);
    % hF32 = plot(x2,f32,'-.','Color',cOSHER,'LineWidth',2);
    % ylabel('$f_3$','Interpreter','latex','FontSize',18)
    % 
    % if strcmp(relaxation_time,'0p5')
    %     ylim([-0.3 0.2]);
    % elseif strcmp(relaxation_time,'0p05')
    %     ylim([-0.03 0.015]);
    %     text(-0.02,0.0025,'$f_3$','Interpreter','latex','FontSize',16, ...
    %     'FontWeight','bold','Color','k','VerticalAlignment','top');
    % end
    % 
    % xlabel('$x$','Interpreter','latex','FontSize',18)
    % 
    % % Use the common x-range here
    % xlim([-0.15 -0.01])
    % 
    % box on
    % grid on
    % 
    % % Keep axes neutral so the colors belong to the curves
    % ax.YAxis(1).Color = 'k';
    % ax.YAxis(2).Color = 'k';
    % 
    % % Compact legend for the methods
    % legend([hRho1 hRho2],{'PRICE','Osher-Solomon'}, ...
    %     'Location','southwest','Interpreter','latex','FontSize',14);
    % 
    % % drawnow;                          % flush pending rendering
    % % fig = gcf;
    % % % Vector PDF (good for Illustrator / high quality)
    % % exportgraphics(fig, 'interface_comparison_rho_and_f3_t0p3.pdf', 'ContentType', 'vector');




    figure('Color','w');
    ax = axes;
    hold(ax,'on');
    set(ax,'FontSize',14,'LineWidth',1,'TickLabelInterpreter','latex');

    % Colors for the two methods
    cPRICE = red;
    cOSHER = blue;

    % Left axis: rho
    yyaxis left
    hRho1 = plot(x1,rho1,'-','Color',cPRICE,'LineWidth',2);
    hRho2 = plot(x2,rho2,'-.','Color',cOSHER,'LineWidth',2);
    ylabel('$\rho$','Interpreter','latex','FontSize',18)
    ylim([1.1 5.8])
    text(-2.5,5.615,'$\rho$','Interpreter','latex','FontSize',16, ...
        'FontWeight','bold','Color','k','VerticalAlignment','top');

    % Right axis: f3
    yyaxis right
    hF31 = plot(x1,f31,'-','Color',cPRICE,'LineWidth',2);
    hF32 = plot(x2,f32,'-.','Color',cOSHER,'LineWidth',2);
    ylabel('$f_3$','Interpreter','latex','FontSize',18)

    if strcmp(relaxation_time,'0p5')
        ylim([-0.3 0.2]);
    elseif strcmp(relaxation_time,'0p05')
        ylim([-0.05 0.025]);
        text(-2.5,-0.001,'$f_3$','Interpreter','latex','FontSize',16, ...
        'FontWeight','bold','Color','k','VerticalAlignment','top');
    end

    xlabel('$x$','Interpreter','latex','FontSize',18)

    % Use the common x-range here
    xlim([-3 3])

    box on
    grid on

    % Keep axes neutral so the colors belong to the curves
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';

    % Compact legend for the methods
    legend([hRho1 hRho2],{'PRICE','Osher-Solomon'}, ...
        'Location','southwest','Interpreter','latex','FontSize',14);

    drawnow;                          % flush pending rendering
    fig = gcf;
    % Vector PDF (good for Illustrator / high quality)
    exportgraphics(fig, 'interface_comparison_rho_and_f3_t3p0.pdf', 'ContentType', 'vector');
    
end

% export_name = 'interface_viscosity_comparison_shockTube_rho.pdf';
% 
% addpath('C:\Users\rikve\Github\PhD-RUG\SpatiallyAdaptiveMomentModels\Nonlinear-systems\Data-processing\Results\export_fig\', '-end');
% export_fig(export_name, '-pdf','-transparent');
% 
% close(fig)
