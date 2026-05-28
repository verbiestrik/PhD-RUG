clear all
close all
clc
format long

% Load data
name1 = 'damBreak-and-smooth_linear_lambda1.0_nu0.1_order1_t5_10000'; type1 = 'swme';
name2 = 'damBreak-and-smooth_linear_lambda1.0_nu0.1_order5_t5_10000'; type2 = 'swme';
name3 = 'damBreak-and-smooth_linear_lambda1.0_nu0.1_adaptiveNonConservative_t5_10000'; type3 = 'swmeAdaptive';
name4 = 'damBreak-and-smooth_linear_lambda1.0_nu0.1_adaptiveConservative_t5_10000'; type4 = 'swmeAdaptive';
reference_data2000x200 = load("damBreak-and-smooth_linear_lambda1.0_nu0.1_t5_2000x200.csv");

% name1 = 'damBreak-and-smooth_linear_lambda0.1_nu0.1_order1_t5_10000'; type1 = 'swme';
% name2 = 'damBreak-and-smooth_linear_lambda0.1_nu0.1_order5_t5_10000'; type2 = 'swme';
% name3 = 'damBreak-and-smooth_linear_lambda0.1_nu0.1_adaptiveNonConservative_t5_10000'; type3 = 'swmeAdaptive';
% name4 = 'damBreak-and-smooth_linear_lambda0.1_nu0.1_adaptiveConservative_t5_10000'; type4 = 'swmeAdaptive';
% % name4 = 'test_newFluctuationsStorage2'; type4 = 'swmeAdaptive';
% reference_data2000x200 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_2000x200.csv");

% name0 = 'damBreak-and-smooth_linear_lambda0.1_nu1.0_order1_t5_10000'; type0 = 'swme';
% name1 = 'damBreak-and-smooth_linear_lambda0.1_nu1.0_order2_t5_10000'; type1 = 'swme';
% name2 = 'damBreak-and-smooth_linear_lambda0.1_nu1.0_order5_t5_10000'; type2 = 'swme';
% name3 = 'damBreak-and-smooth_linear_lambda0.1_nu1.0_adaptiveNonConservative_t5_10000'; type3 = 'swmeAdaptive';
% name4 = 'damBreak-and-smooth_linear_lambda0.1_nu1.0_adaptiveConservative_t5_10000'; type4 = 'swmeAdaptive';
% reference_data2000x200 = load("damBreak-and-smooth_linear_lambda0.1_nu1.0_t5_2000x200.csv");

t_m_evol = load("discreteTimes_couplingCIF.csv");

moments_at_Min15 = load("smoothPlusDam_timeEvolution_order_at_position_Min15_CIF.csv");
moments_at_Min10 = load("smoothPlusDam_timeEvolution_order_at_position_Min10_CIF.csv");
moments_at_Min5 = load("smoothPlusDam_timeEvolution_order_at_position_Min5_CIF.csv");
moments_at_0 = load("smoothPlusDam_timeEvolution_order_at_position_0_CIF.csv");
moments_at_5 = load("smoothPlusDam_timeEvolution_order_at_position_5_CIF.csv");
moments_at_10 = load("smoothPlusDam_timeEvolution_order_at_position_10_CIF.csv");
moments_at_15 = load("smoothPlusDam_timeEvolution_order_at_position_15_CIF.csv");

% name1 = 'damBreak-and-smooth_linear_lambda0.1_nu0.1_adaptiveConservative_t0.5_10000'; type1 = 'swmeAdaptive';
% name2 = 'damBreak-and-smooth_linear_lambda0.1_nu0.1_adaptiveConservative_t2.5_10000'; type2 = 'swmeAdaptive';
% name3 = 'damBreak-and-smooth_linear_lambda0.1_nu0.1_adaptiveConservative_t5_10000'; type3 = 'swmeAdaptive';
% name4 = 'damBreak-and-smooth_linear_lambda0.1_nu0.1_adaptiveConservative_t3.5_10000'; type4 = 'swmeAdaptive';
% reference_data2000x200 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_2000x200.csv");

x_ref2000x200 = reference_data2000x200(:,1); h_ref2000x200 = reference_data2000x200(:,2); u_ref2000x200 = reference_data2000x200(:,3);
alpha1_ref2000x200 = 3.*reference_data2000x200(:,4); alpha2_ref2000x200 = 5.*reference_data2000x200(:,5);

x_ref2000x200_short = reference_data2000x200(1:4:end,1); 
h_ref2000x200_short = reference_data2000x200(1:4:end,2); 
u_ref2000x200_short = reference_data2000x200(1:4:end,3);
alpha1_ref2000x200_short = 3.*reference_data2000x200(1:4:end,4); 
alpha2_ref2000x200_short = 5.*reference_data2000x200(1:4:end,5);

% [x0,h0,u0,alpha10,alpha20,alpha30,alpha40,alpha50,moments0] = readDataAdaptiveSWME1D(name0,type0);
[x1,h1,u1,alpha11,alpha21,alpha31,alpha41,alpha51,moments1] = readDataAdaptiveSWME1D(name1,type1);
[x2,h2,u2,alpha12,alpha22,alpha32,alpha42,alpha52,moments2] = readDataAdaptiveSWME1D(name2,type2);
[x3,h3,u3,alpha13,alpha23,alpha33,alpha43,alpha53,moments3] = readDataAdaptiveSWME1D(name3,type3);
[x4,h4,u4,alpha14,alpha24,alpha34,alpha44,alpha54,moments4] = readDataAdaptiveSWME1D(name4,type4);

% plotting = 'all';
% plotting = 'h';
% plotting = 'u';
plotting = 'alpha1';
% plotting = 'alpha2';
% plotting = 'alpha3';
% plotting = 'alpha4';
% plotting = 'alpha5';
% plotting = 'number_of_moments';
% plotting = 'number_of_moments-evolution';

brown = [171, 104, 87]./255;

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
    % plot1 = plot(x1,h1,x2,h2,x3,h3,x4,h4,x_ref2000x200,h_ref2000x200,'*')
    plot1 = plot(x1,h1,x2,h2,x3,h3,x4,h4)
    % plot1 = plot(x0,h0,x1,h1,x2,h2,x3,h3,x4,h4)

    axis([-20,20,2.95,4.1]);
    % axis([-3.2,3.2,2.95,3.7]);
    xlabel('$x$','Interpreter','latex','FontSize',18)
    ylabel('$h$','Interpreter','latex','FontSize',18)

    linewidth1 = 2;

    set(plot1(1),'LineWidth',2);
    set(plot1(2),'LineWidth',2);
    set(plot1(3),'LineWidth',2);
    set(plot1(4),'LineWidth',2);
    % set(plot1(5),'LineWidth',2);

    grey = [0.4,0.4,0.4];
    set(plot1(1),'LineStyle','-');
    % set(plot1(2),'LineStyle','-');
    set(plot1(2),'LineStyle','-.');
    set(plot1(3),'LineStyle','--');
    set(plot1(4),'LineStyle',':');
    set(plot1(1),'Color','r');
    % set(plot1(2),'Color','magenta');
    set(plot1(2),'Color','blue');
    set(plot1(3),'Color','black');
    set(plot1(4),'Color','cyan');

    % set(plot1(5),'Color',brown);
    % set(plot1(5),'MarkerSize',3);

    % leg = legend('SWME_1','SWME_2','SWME_5','A-SWME-PBC','A-SWME-PCIF','Location','best');
    leg = legend('SWME_1','SWME_5','A-SWME-PBC','A-SWME-PCIF','Location','best');
    % leg = legend('SWME_2','SWME_5','A-SWME-PBC','A-SWME-PCIF','Reference','Location','northeast');
    set(leg,'FontSize',12); 
end
if(strcmp(plotting,'u'))
    % plot1 = plot(x1,u1,x2,u2,x3,u3,x4,u4,x_ref2000x200,u_ref2000x200,'*')
    plot1 = plot(x1,u1,x2,u2,x3,u3,x4,u4)
    % plot1 = plot(x0,u0,x1,u1,x2,u2,x3,u3,x4,u4)

    axis([-20.,20,0.01,0.425]);
    % axis([-3.25,3.25,-0.025,0.2]);
    % plot1 = plot(x1,h1,x2,h2,x3,h3,x4,h4,x_ref2000x200,h_ref2000x200,'*')
    % plot1 = plot(x1,h1,x2,h2,x3,h3,x4,h4)

    % axis([-3.2,3.2,2.95,3.7]);
    xlabel('$x$','Interpreter','latex','FontSize',18)
    ylabel('$u_m$','Interpreter','latex','FontSize',18)
    
    linewidth1 = 2;

    set(plot1(1),'LineWidth',2);
    set(plot1(2),'LineWidth',2);
    set(plot1(3),'LineWidth',2);
    set(plot1(4),'LineWidth',2);
    % set(plot1(5),'LineWidth',2);

    grey = [0.4,0.4,0.4];
    set(plot1(1),'LineStyle','-');
    % set(plot1(2),'LineStyle','-');
    set(plot1(2),'LineStyle','-.');
    set(plot1(3),'LineStyle','--');
    set(plot1(4),'LineStyle',':');
    set(plot1(1),'Color','r');
    % set(plot1(2),'Color','magenta');
    set(plot1(2),'Color','blue');
    set(plot1(3),'Color','black');
    set(plot1(4),'Color','cyan');

    % set(plot1(5),'Color',brown);
    % set(plot1(5),'MarkerSize',3);

    % leg = legend('SWME_1','SWME_2','SWME_5','A-SWME-PBC','A-SWME-PCIF','Location','best');
    leg = legend('SWME_1','SWME_5','A-SWME-PBC','A-SWME-PCIF','Location','northeast');
    % leg = legend('SWME_2','SWME_5','A-SWME-PBC','A-SWME-PCIF','Reference','Location','northeast');
    set(leg,'FontSize',12); 
end
if(strcmp(plotting,'alpha1'))
    % plot1 = plot(x1,alpha11,x2,alpha12,x3,alpha13,x4,alpha14,x_ref2000x200,alpha1_ref2000x200,'*')
    plot1 = plot(x1,alpha11,x2,alpha12,x3,alpha13,x4,alpha14)
    % plot1 = plot(x0,alpha10,x1,alpha11,x2,alpha12,x3,alpha13,x4,alpha14)

    axis([-20,20,-0.1,0.001]);
    % axis([-3.75,3.25,-0.135,0]);

    xlabel('$x$','Interpreter','latex','FontSize',18)
    ylabel('$\alpha_1$','Interpreter','latex','FontSize',18)

    linewidth1 = 2;

    set(plot1(1),'LineWidth',2);
    set(plot1(2),'LineWidth',2);
    set(plot1(3),'LineWidth',2);
    set(plot1(4),'LineWidth',2);
    % set(plot1(5),'LineWidth',2);

    grey = [0.4,0.4,0.4];
    set(plot1(1),'LineStyle','-');
    % set(plot1(2),'LineStyle','-');
    set(plot1(2),'LineStyle','-.');
    set(plot1(3),'LineStyle','--');
    set(plot1(4),'LineStyle',':');
    set(plot1(1),'Color','r');
    % set(plot1(2),'Color','magenta');
    set(plot1(2),'Color','blue');
    set(plot1(3),'Color','black');
    set(plot1(4),'Color','cyan');

    % set(plot1(5),'Color',brown);
    % set(plot1(5),'MarkerSize',3);

    % leg = legend('SWME_1','SWME_2','SWME_5','A-SWME-PBC','A-SWME-PCIF','Location','best');
    leg = legend('SWME_1','SWME_5','A-SWME-PBC','A-SWME-PCIF','Location','best');
    % leg = legend('SWME_2','SWME_5','A-SWME-PBC','A-SWME-PCIF','Reference','Location','northeast');
    set(leg,'FontSize',12); 
end
if(strcmp(plotting,'alpha2'))
    % plot1 = plot(x1,alpha21,x2,alpha22,x3,alpha23,x4,alpha24,x_ref2000x200,alpha2_ref2000x200,'*')
    % plot1 = plot(x1,alpha21,x2,alpha22,x3,alpha23,x4,alpha24)
    plot1 = plot(x0,alpha20,x1,alpha21,x2,alpha22,x3,alpha23,x4,alpha24)

    axis([-20.,20,-0.1,0.015]);
    % axis([-3.0,3.0,-0.075,0.015]);
    xlabel('$x$','Interpreter','latex','FontSize',18)
    ylabel('$\alpha_2$','Interpreter','latex','FontSize',18)
    
    linewidth1 = 2;

    set(plot1(1),'LineWidth',2);
    set(plot1(2),'LineWidth',2);
    set(plot1(3),'LineWidth',2);
    set(plot1(4),'LineWidth',2);
    set(plot1(5),'LineWidth',2);

    grey = [0.4,0.4,0.4];
    set(plot1(1),'LineStyle','-');
    set(plot1(2),'LineStyle','-');
    set(plot1(3),'LineStyle','-.');
    set(plot1(4),'LineStyle','--');
    set(plot1(5),'LineStyle',':');
    set(plot1(1),'Color','r');
    set(plot1(2),'Color','magenta');
    set(plot1(3),'Color','blue');
    set(plot1(4),'Color','black');
    set(plot1(5),'Color','cyan');

    % set(plot1(5),'Color',brown);
    % set(plot1(5),'MarkerSize',3);

    leg = legend('SWME_1','SWME_2','SWME_5','A-SWME-PBC','A-SWME-PCIF','Location','best');
    % leg = legend('SWME_1','SWME_5','A-SWME-PBC','A-SWME-PCIF','Location','best');
    % leg = legend('SWME_2','SWME_5','A-SWME-PBC','A-SWME-PCIF','Reference','Location','northeast');
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
if(strcmp(plotting,'number_of_moments'))
    yyaxis right
    % plot1 = plot(x1,moments1,x2,moments2,x3,moments3,x4,moments4)
    plot1 = plot(x3,moments3,x4,moments4)
    ylabel('Order $M$','FontSize',18,'Interpreter','latex');
    ylim([-0.5 5.5])

    yticks([0 1 2 3 4 5])

    set(plot1(1:2),'LineStyle','None');
    set(plot1(1),'marker','x');
    set(plot1(2),'marker','+');
    set(plot1(1),'Color','red');
    set(plot1(2),'Color','blue');
    set(plot1(1), 'MarkerSize',8, 'LineWidth',2)
    set(plot1(2), 'MarkerSize',4.5, 'LineWidth',1.2)

    yyaxis left
    % plot(x1,u1,'LineWidth',2)
    plot(x_ref2000x200,h_ref2000x200,'LineWidth',3)
    
    % axis([-20,20,-0.5,5.5]);
    % axis([-3.5,0,3.4,3.9]);
    xlabel('$x$','Interpreter','latex','FontSize',18)
    ylabel('$h$','Interpreter','latex','FontSize',18)
    ylim([2.9 4.1])
    % ylim([-0.02 0.3])

    leg = legend('$h$ (reference)','Adaptive order: PBC','Adaptive order: PCIF','Location','southwest','Interpreter','latex');
    % leg = legend('h','Adaptive order: PBC','Adaptive order: PCIF-group','Location','southwest');
    set(leg,'FontSize',12); 
end
if(strcmp(plotting,'number_of_moments-evolution'))
    t=t_m_evol;
    n = 5293;
    
    idx = 1:10:n;

    % plot1 = plot( ...
    %     t(1:n), moments_at_Min15(1:n), ...
    %     t(1:n), moments_at_Min10(1:n), ...
    %     t(1:n), moments_at_Min5(1:n), ...
    %     t(1:n), moments_at_0(1:n), ...
    %     t(1:n), moments_at_5(1:n), ...
    %     t(1:n), moments_at_10(1:n), ...
    %     t(1:n), moments_at_15(1:n));

    plot1 = plot(t(idx), moments_at_Min5(idx),t(idx), moments_at_0(idx),t(idx), moments_at_5(idx));

    % axis([-1,1,0.85,5.15]);
    xlabel('$t$', 'Interpreter','latex','FontSize',18)
    ylabel('Order $M$', 'Interpreter','latex')
    ylim([-0.3 5.5])

    yticks([0 1 2 3 4 5])

    cols = [
        0.0000 0.4470 0.7410
        0.8500 0.3250 0.0980
        0.9290 0.6940 0.1250
        0.4940 0.1840 0.5560
        0.4660 0.6740 0.1880
        0.3010 0.7450 0.9330
        0.6350 0.0780 0.1840
    ];
    
    styles = {'-','--',':','-','--',':','-.'};
    markers = {'none','none','none','o','s','^','d'};
    
    for k = 1:3
        set(plot1(k), ...
            'LineStyle', styles{k}, ...
            'Color', cols(k,:), ...
            'LineWidth', 2.5, ...
            'Marker', markers{k}, ...
            'MarkerSize', 5);
    end
    % set(plot1(5:7),'LineWidth',0.75);

    yticks_points = [0 1 2 3 4 5];
    
    yticks(yticks_points)
    
    ax = gca;
    ax.TickLabelInterpreter = 'latex';

    % leg = legend({ ...
    %     '$x=-15$', '$x=-10$', '$x=-5$', '$x=0$', ...
    %     '$x=5$', '$x=10$', '$x=15$'}, ...
    %     'Location','best', 'Interpreter','latex');
    
    leg = legend({'$x=-5$', '$x=0$','$x=5$'}, ...
        'Location','best', 'Interpreter','latex');

    leg.FontSize = 13;
    % leg.IconColumnWidth = 10;      % narrower legend box
end

% addpath('C:\Users\rikve\Github\PhD-RUG\SpatiallyAdaptiveMomentModels\Nonlinear-systems\Data-processing\Results\export_fig\', '-end');
% export_fig('lambda1.0_nu0.1_t5_alpha1_full_newName.pdf', '-pdf','-transparent');