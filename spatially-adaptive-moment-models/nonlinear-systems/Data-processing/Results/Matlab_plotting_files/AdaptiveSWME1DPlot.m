clear all
close all
clc
format long

% Load data
name1 = 'smooth_lowOrder'; type1 = 'swe';
name2 = 'smooth_highOrder'; type2 = 'swme';
name3 = 'smooth_adaptive_conservative'; type3 = 'swme';
name4 = 'smooth_adaptive_nonConservative'; type4 = 'swme';

[x1,h1,u1,alpha11,alpha21,alpha31,alpha41,alpha51] = readDataSWME1D(name1,type1);
[x2,h2,u2,alpha12,alpha22,alpha32,alpha42,alpha52] = readDataSWME1D(name2,type2);
[x3,h3,u3,alpha13,alpha23,alpha33,alpha43,alpha53] = readDataSWME1D(name3,type3);
[x4,h4,u4,alpha14,alpha24,alpha34,alpha44,alpha54] = readDataSWME1D(name4,type4);

% plotting = 'all';
 plotting = 'h';
% plotting = 'u';
% plotting = 'alpha1';
% plotting = 'alpha2';
% plotting = 'alpha3';
% plotting = 'alpha4';
% plotting = 'alpha5';


if(strcmp(plotting,'all'))
    ax1 = subplot(2,2,1);
    ax2 = subplot(2,2,2);
    ax3 = subplot(2,2,3);
    ax4 = subplot(2,2,4);
    ax5 = subplot(2,2,5);
    ax6 = subplot(2,2,6);
    ax7 = subplot(2,2,7);

    plot(ax1,x1,h1,x2,h2,x3,h3,x4,h4)
    legend(ax1,name1,name2,name3,name4)
    plot(ax2,x1,u1,x2,u2,x3,u3,x4,u4)

    plot(ax3,x1,alpha11,x2,alpha12,x3,alpha13,x4,alpha14)

    plot(ax4,x1,alpha21,x2,alpha22,x3,alpha23,x4,alpha24)

    plot(ax5,x1,alpha31,x2,alpha32,x3,alpha33,x4,alpha34)

    plot(ax6,x1,alpha41,x2,alpha42,x3,alpha43,x4,alpha44)

    plot(ax7,x1,alpha51,x2,alpha52,x3,alpha53,x4,alpha54)
end

if(strcmp(plotting,'h'))
    plot1 = plot(x1,h1,x2,h2,x3,h3,x4,h4)
%     legend(name1,name2,name3,name4)

    axis([-20.,20,2.75,4.25]);
    xlabel('x')
    ylabel('h')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    %set(plot1(2:4),'LineStyle','.');
    %set(plot1(2:4),'MarkerSize',linewidth1+12);
    %set(plot1(4),'LineStyle','.');
    %set(plot1(4),'MarkerSize',linewidth1+12);

    grey = [0.4,0.4,0.4];
    set(plot1(1:2),'LineStyle','-.');
    set(plot1(3),'LineStyle','--');
    set(plot1(1),'Color','r');
    set(plot1(2),'Color','blue');
    set(plot1(3),'Color','black');
    
%     leg = legend('SWME','HSWME','\betaHSWME','reference','Location','southwest');
%     leg = legend('SWME','reference','Location','southwest');
    leg = legend('Low order','High order','Adaptive','Location','northwest');
    set(leg,'FontSize',12); 
end
if(strcmp(plotting,'u'))
    plot1 = plot(x1,u1,x2,u2,x3,u3,x4,u4)
%     legend(name1,name2,name3,name4)

    axis([-10.,10.,-0.15,0.15]);
    xlabel('x')
    ylabel('u_m')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    %set(plot1(2:4),'LineStyle','.');
    %set(plot1(2:4),'MarkerSize',linewidth1+12);
    %set(plot1(4),'LineStyle','.');
    %set(plot1(4),'MarkerSize',linewidth1+12);

    grey = [0.4,0.4,0.4];
    set(plot1(1:2),'LineStyle','-.');
    set(plot1(3),'LineStyle','--');
    set(plot1(1),'Color','r');
    set(plot1(2),'Color','blue');
    set(plot1(3),'Color','black');
    
%     leg = legend('SWME','HSWME','\betaHSWME','reference','Location','southwest');
%     leg = legend('SWME','reference','Location','southwest');
    leg = legend('Low order','High order','Adaptive','Location','northwest');
    set(leg,'FontSize',12); 
end
if(strcmp(plotting,'alpha1'))
    plot1 = plot(x1,alpha11,x2,alpha12,x3,alpha13,x4,alpha14)
%     legend(name1,name2,name3,name4)

    axis([-0.4,0.4,-0.8,0.0]);
    xlabel('x')
    ylabel('\alpha_1')
    
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
if(strcmp(plotting,'alpha2'))
    plot1 = plot(x1,alpha21,x2,alpha22,x3,alpha23,x4,alpha24)
%     legend(name1,name2,name3,name4)

    axis([-0.4,0.4,-0.8,0.0]);
    xlabel('x')
    ylabel('\alpha_2')
    
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
if(strcmp(plotting,'alpha3'))
    plot1 = plot(x1,alpha31,x2,alpha32,x3,alpha33,x4,alpha34)
%     legend(name1,name2,name3,name4)

    axis([-0.4,0.4,-0.8,0.0]);
    xlabel('x')
    ylabel('\alpha_3')
    
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

%addpath('C:\Users\rikve\Github\PhD-RUG\spatially-adaptive-moment-models\nonlinear-systems\Data-processing\Results\export_fig\', '-end')
%export_fig('smooth_h.pdf', '-pdf','-transparent');