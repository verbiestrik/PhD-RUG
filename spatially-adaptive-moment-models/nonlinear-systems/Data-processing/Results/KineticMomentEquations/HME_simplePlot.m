clear all
close all
clc
format long

name = 'shockPlusSmooth_order12_t0_Kn0p5';type='hme';

[x,rho,u,T,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,moments] = readDataAdaptiveHME1D(name,type);

% plotting = 'rho';
plotting = 'rho_u_T';

blue = [0, 0.4470, 0.7410];
green = [0.4660, 0.6740, 0.1880];
yellow = [0.9290, 0.6940, 0.1250];
red = [0.8500, 0.3250, 0.0980];
brown = [171, 104, 87]./255;

fig = figure;

if(strcmp(plotting,'rho_u_T'))
    plot1 = plot(x,rho,x,u,x,T)
    set(plot1(1:3),'LineWidth',3)
    set(plot1(1),'LineStyle','-.')
    set(plot1(2),'LineStyle','--')    
    set(plot1(3),'LineStyle',':')

    set(plot1(1),'Color',red)
    set(plot1(2),'Color',blue)    
    set(plot1(3),'Color',brown)

    xlabel('$x$','FontSize',15,'Interpreter','latex')
    ylabel('\(\rho,u,\theta\)','FontSize',15,'interpreter','latex')
    ylim([0.25 2.25])
    xlim([-2.75,3.25])

    leg = legend('$\rho$','$u$','$\theta$','Location','best','interpreter','latex');
    set(leg,'FontSize',12); 
end

if(strcmp(plotting,'rho'))
    plot1 = plot(x,rho)
    
    axis([-2.75,3.25,0.75,2.25]);
    xlabel('$x$','FontSize', 20,'Interpreter','latex')
    ylabel('$\rho$','FontSize', 20,'Interpreter','latex')
    
    linewidth1 = 3;
    set(plot1(1),'LineWidth',linewidth1);

    set(plot1(1),'Color',red);

end

% export_name = 'shockPlusSmooth_initialConditions.pdf';
% 
% addpath('C:\Users\rikve\Github\PhD-RUG\spatially-adaptive-moment-models\nonlinear-systems\Data-processing\Results\export_fig\', '-end');
% export_fig(export_name, '-pdf','-transparent');