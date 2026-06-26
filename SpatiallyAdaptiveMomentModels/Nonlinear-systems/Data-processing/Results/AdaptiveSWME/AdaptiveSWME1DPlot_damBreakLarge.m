clear all
close all
clc
format long

% name1 = 'damBreak_classical_order4_time0p4'; type1 = 'swme';
% name2 = 'damBreak_classical_order5_time0p4'; type2 = 'swme';
% name3 = 'damBreak_adaptive_couplingPBC_time0p4'; type3 = 'swmeAdaptive';
% name4 = 'damBreak_adaptive_couplingCIF_time0p4'; type4 = 'swmeAdaptive';

name1 = 'damBreak_classical_order0_time1p0'; type1 = 'swe';
name2 = 'damBreak_classical_order5_time1p0'; type2 = 'swme';
name3 = 'damBreak_adaptive_couplingPBC_time1p0'; type3 = 'swmeAdaptive';
name4 = 'damBreak_adaptive_couplingCIF_time1p0'; type4 = 'swmeAdaptive';

t_m_evol = load("discreteTimes_couplingCIF.csv");

moments_at_0p1 = load("timeEvolution_order_at_position_0p1_CIF.csv");
moments_at_0p25 = load("timeEvolution_order_at_position_0p25_CIF.csv");
moments_at_0p5 = load("timeEvolution_order_at_position_0p5_CIF.csv");
moments_at_0p75 = load("timeEvolution_order_at_position_0p75_CIF.csv");
moments_at_Min0p1 = load("timeEvolution_order_at_position_Min0p1_CIF.csv");
moments_at_Min0p25 = load("timeEvolution_order_at_position_Min0p25_CIF.csv");
moments_at_Min0p5 = load("timeEvolution_order_at_position_Min0p5_CIF.csv");
moments_at_Min0p75 = load("timeEvolution_order_at_position_Min0p75_CIF.csv");

% ref = load("damBreakLarge_lambda0p1_nu1p0_t0p4_250x100.csv");
ref = load("damBreakLarge_lambda0p1_nu1p0_t1p0_800x200_allVars.csv");
x_ref = ref(:,1); h_ref = ref(:,2); u_ref = ref(:,3);alpha1_ref = 3.*ref(:,4); alpha2_ref = 5.*ref(:,5);

ref0p4 = load("damBreakLarge_lambda0p1_nu1p0_t0p4_250x100.csv");
ref0p35 = load("damBreakLarge_lambda0p1_nu1p0_t0p35_250x100.csv");
ref0p3 = load("damBreakLarge_lambda0p1_nu1p0_t0p3_250x100.csv");
ref0p25 = load("damBreakLarge_lambda0p1_nu1p0_t0p25_250x100.csv");
ref0p2 = load("damBreakLarge_lambda0p1_nu1p0_t0p2_250x100.csv");
ref0p15 = load("damBreakLarge_lambda0p1_nu1p0_t0p15_250x100.csv");
ref0p1 = load("damBreakLarge_lambda0p1_nu1p0_t0p1_250x100.csv");
ref0p05 = load("damBreakLarge_lambda0p1_nu1p0_t0p05_250x100.csv");
h_ref0p05=ref0p05(:,2);u_ref0p05=ref0p05(:,3);alpha1_ref0p05=3.*ref0p05(:,4);alpha20p05=5.*ref0p05(:,5);
h_ref0p1=ref0p1(:,2);u_ref0p1=ref0p1(:,3);alpha1_ref0p1=3.*ref0p1(:,4);alpha20p1=5.*ref0p1(:,5);
h_ref0p15=ref0p15(:,2);u_ref0p15=ref0p15(:,3);alpha1_ref0p15=3.*ref0p15(:,4);alpha20p15=5.*ref0p15(:,5);
h_ref0p2=ref0p2(:,2);u_ref0p2=ref0p2(:,3);alpha1_ref0p2=3.*ref0p2(:,4);alpha20p2=5.*ref0p2(:,5);
h_ref0p25=ref0p25(:,2);u_ref0p25=ref0p05(:,3);alpha1_ref0p25=3.*ref0p25(:,4);alpha20p25=5.*ref0p25(:,5);
h_ref0p3=ref0p3(:,2);u_ref0p3=ref0p3(:,3);alpha1_ref0p3=3.*ref0p3(:,4);alpha20p3=5.*ref0p3(:,5);
h_ref0p35=ref0p35(:,2);u_ref0p35=ref0p35(:,3);alpha1_ref0p35=3.*ref0p35(:,4);alpha20p35=5.*ref0p35(:,5);
h_ref0p4=ref0p4(:,2);u_ref0p4=ref0p4(:,3);alpha1_ref0p4=3.*ref0p4(:,4);alpha20p4=5.*ref0p4(:,5);

[x1,h1,u1,alpha11,alpha21,alpha31,alpha41,alpha51,moments1] = readDataAdaptiveSWME1D(name1,type1);
[x2,h2,u2,alpha12,alpha22,alpha32,alpha42,alpha52,moments2] = readDataAdaptiveSWME1D(name2,type2);
[x3,h3,u3,alpha13,alpha23,alpha33,alpha43,alpha53,moments3] = readDataAdaptiveSWME1D(name3,type3);
[x4,h4,u4,alpha14,alpha24,alpha34,alpha44,alpha54,moments4] = readDataAdaptiveSWME1D(name4,type4);

plotting_mode = 'model_comparison';
% plotting_mode = 'variables_time_evolution';
% plotting_mode = 'orders_time_evolution';

% plotting = 'h';
% plotting = 'u';
% plotting = 'alpha1';
plotting = 'alpha2';
% plotting = 'number_of_moments';
% plotting = 'number_of_moments-evolution';

brown = [171, 104, 87]./255;

if(strcmp(plotting,'h'))
    if(strcmp(plotting_mode,'model_comparison'))
        % plot1 = plot(x1,h1,x2,h2,x3,h3,x4,h4,x_ref,h_ref,'*')
        plot1 = plot(x3,h3,x4,h4,x_ref,h_ref)
    
        axis([-3.075,2.075,0.9,5.1]);
        xlabel('$x$','FontSize',18, 'Interpreter','latex')
        ylabel('$h$','FontSize',18, 'Interpreter','latex')

        linewidth1 = 2;
        set(plot1(1),'LineWidth',4);
        set(plot1(2),'LineWidth',3);
        set(plot1(3),'LineWidth',2);

        grey = [0.4,0.4,0.4];
        set(plot1(1),'LineStyle',':');
        set(plot1(2),'LineStyle','--');
        set(plot1(3),'LineStyle','-');
        % set(plot1(3),'LineStyle','--');
        % set(plot1(4),'LineStyle',':');
        set(plot1(1),'Color','r');
        set(plot1(2),'Color','blue');
        set(plot1(3),'Color','k');
        % set(plot1(3),'Color','black');
        % set(plot1(4),'Color','cyan');

        % set(plot1(3),'Color','k');
        % set(plot1(3),'MarkerSize',8);

        % leg = legend('SWME_0','SWME_5','A-SWME-PBC','A-SWME-PCIF','Reference','Location','best');
        leg = legend('A-SWME-PBC','A-SWME-PCIF','Reference','Location','best');
        set(leg,'FontSize',12);
    else
        plot1 = plot(x_ref,h_ref0p05,x_ref,h_ref0p1,x_ref,h_ref0p15,x_ref,h_ref0p2, ...
            x_ref,h_ref0p25,x_ref,h_ref0p3,x_ref,h_ref0p35,x_ref,h_ref0p4)
    
        axis([-1,1,0.85,5.15]);
        xlabel('$x$','FontSize',18, 'Interpreter','latex')
        ylabel('$h$','FontSize',18, 'Interpreter','latex')
        
        linewidth1 = 5;
        set(plot1(1:8),'LineWidth',linewidth1);
    
        cols = [
            0.0000 0.4470 0.7410
            0.8500 0.3250 0.0980
            0.9290 0.6940 0.1250
            0.4940 0.1840 0.5560
            0.4660 0.6740 0.1880
            0.3010 0.7450 0.9330
            0.6350 0.0780 0.1840
            0.2000 0.2000 0.2000
        ];
        
        styles = {'-','--',':','-','--',':','-.','-'};
        markers = {'none','none','none','none','o','s','^','d'};
        
        for k = 1:8
            set(plot1(k), ...
                'LineStyle', styles{k}, ...
                'Color', cols(k,:), ...
                'LineWidth', 2.5, ...
                'Marker', markers{k}, ...
                'MarkerSize', 5);
        end
        set(plot1(5:8),'LineWidth',0.75);

        leg = legend({ ...
            '$t=0.05$', '$t=0.1$', '$t=0.15$', '$t=0.2$', ...
            '$t=0.25$', '$t=0.3$', '$t=0.35$', '$t=0.4$' }, ...
            'Location','best', 'Interpreter','latex');
        
        xs = [-0.75 -0.5 -0.25 -0.1 0.1 0.25 0.5 0.75];
        xticks(xs)
        xtickformat('%.2f')
        
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
        
        hold on
        for k = 1:numel(xs)
            xline(xs(k), '--', 'HandleVisibility', 'off', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        end
        hold off

        leg.FontSize = 14;
        leg.IconColumnWidth = 10;      % narrower legend box
    end
end
if(strcmp(plotting,'u'))
    if(strcmp(plotting_mode,'model_comparison'))
        % plot1 = plot(x1,u1,x2,u2,x3,u3,x4,u4,x_ref,u_ref,'*')
        plot1 = plot(x3,u3,x4,u4,x_ref,u_ref)

        axis([-3.075,2.075,-0.2,1.2]);
        xlabel('$x$','FontSize',18, 'Interpreter','latex')
        ylabel('$u_m$','FontSize',18, 'Interpreter','latex')

        linewidth1 = 2;
        set(plot1(1),'LineWidth',4);
        set(plot1(2),'LineWidth',3);
        set(plot1(3),'LineWidth',2);

        grey = [0.4,0.4,0.4];
        set(plot1(1),'LineStyle',':');
        set(plot1(2),'LineStyle','--');
        set(plot1(3),'LineStyle','-');
        % set(plot1(3),'LineStyle','--');
        % set(plot1(4),'LineStyle',':');
        set(plot1(1),'Color','r');
        set(plot1(2),'Color','blue');
        set(plot1(3),'Color','k');
        % set(plot1(3),'Color','black');
        % set(plot1(4),'Color','cyan');

        % set(plot1(3),'Color','k');
        % set(plot1(3),'MarkerSize',8);

        % leg = legend('SWME_0','SWME_5','A-SWME-PBC','A-SWME-PCIF','Reference','Location','best');
        leg = legend('A-SWME-PBC','A-SWME-PCIF','Reference','Location','best');
        set(leg,'FontSize',12); 
    else
        plot1 = plot(x_ref,u_ref0p05,x_ref,u_ref0p1,x_ref,u_ref0p15,x_ref,u_ref0p2, ...
            x_ref,u_ref0p25,x_ref,u_ref0p3,x_ref,u_ref0p35,x_ref,u_ref0p4)
    
        axis([-1,1,-0.2,1.25]);
        xlabel('$x$','FontSize',18, 'Interpreter','latex')
        ylabel('$u_m$','FontSize',18, 'Interpreter','latex')
        
        linewidth1 = 2;
        set(plot1(1:8),'LineWidth',linewidth1);
    
        cols = [
            0.0000 0.4470 0.7410
            0.8500 0.3250 0.0980
            0.9290 0.6940 0.1250
            0.4940 0.1840 0.5560
            0.4660 0.6740 0.1880
            0.3010 0.7450 0.9330
            0.6350 0.0780 0.1840
            0.2000 0.2000 0.2000
        ];
        
        styles = {'-','--',':','-','--',':','-.','-'};
        markers = {'none','none','none','none','o','s','^','d'};
        
        for k = 1:8
            set(plot1(k), ...
                'LineStyle', styles{k}, ...
                'Color', cols(k,:), ...
                'LineWidth', 1.5, ...
                'Marker', markers{k}, ...
                'MarkerSize', 5);
        end
    
        leg = legend({ ...
            '$t=0.05$', '$t=0.1$', '$t=0.15$', '$t=0.2$', ...
            '$t=0.25$', '$t=0.3$', '$t=0.35$', '$t=0.4$' }, ...
            'Location','best', 'Interpreter','latex');
        
        leg.FontSize = 14;
        leg.IconColumnWidth = 10;      % narrower legend box
    end
end
if(strcmp(plotting,'alpha1'))
    if(strcmp(plotting_mode,'model_comparison'))
        plot1 = plot(x1,alpha11,x2,alpha12,x3,alpha13,x4,alpha14,x_ref,alpha1_ref,'*')
    
        axis([-3,2,-0.7,0.12]);
        xlabel('$x$','FontSize',18, 'Interpreter','latex')
        ylabel('$\alpha_1$','FontSize',18, 'Interpreter','latex')
        
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
        
        set(plot1(5),'Color',brown);
        set(plot1(5),'MarkerSize',3);
    
        leg = legend('SWME_0','SWME_5','A-SWME-PBC','A-SWME-PCIF','Reference','Location','best');
        set(leg,'FontSize',12); 
    end
end
if(strcmp(plotting,'alpha2'))
    if(strcmp(plotting_mode,'model_comparison'))
        plot1 = plot(x1,alpha21,x2,alpha22,x3,alpha23,x4,alpha24,x_ref,alpha2_ref,'*')
    
        axis([-3,2,-0.7,0.12]);
        xlabel('$x$','FontSize',18, 'Interpreter','latex')
        ylabel('$\alpha_2$','FontSize',18, 'Interpreter','latex')
        
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
        
        set(plot1(5),'Color',brown);
        set(plot1(5),'MarkerSize',3);
    
        leg = legend('SWME_0','SWME_5','A-SWME-PBC','A-SWME-PCIF','Reference','Location','best');
        set(leg,'FontSize',12); 
    end
end
if(strcmp(plotting,'number_of_moments'))
    yyaxis right
    % plot1 = plot(x1,moments1,x2,moments2,x3,moments3,x4,moments4)
    plot1 = plot(x3,moments3,x4,moments4)
    ylabel('Order')
    ylim([-0.3 5.5])

    set(plot1(1:2),'LineStyle','None');
    set(plot1(1),'marker','x');
    set(plot1(2),'marker','+');
    set(plot1(1),'Color','red');
    set(plot1(2),'Color','blue');
    set(plot1(1:2), 'MarkerSize',3)

    yyaxis left
    % plot(x1,u1,'LineWidth',2)
    plot(x2,h2,'LineWidth',2)
    
    % axis([-20,20,-0.5,5.5]);
    % axis([-3.5,0,3.4,3.9]);
    xlabel('x')
    % ylabel('u_m')
    ylabel('h')
    ylim([2.9 4.1])
    % ylim([-0.02 0.3])

    leg = legend('h (SWME_5)','Adaptive order: PBC','Adaptive order: PCIF','Location','southwest');
    % leg = legend('h','Adaptive order: PBC','Adaptive order: PCIF-group','Location','southwest');
    set(leg,'FontSize',13); 
end
if(strcmp(plotting,'number_of_moments-evolution'))
    t=t_m_evol;
    plot1 = plot(t,moments_at_0p1,t,moments_at_0p25,t,moments_at_0p5,t,moments_at_0p75, ...
        t,moments_at_Min0p1,t,moments_at_Min0p25,t,moments_at_Min0p5,t,moments_at_Min0p75)

    % axis([-1,1,0.85,5.15]);
    xlabel('$t$', 'Interpreter','latex')
    ylabel('Order $M$', 'Interpreter','latex')
    
    linewidth1 = 2;
    set(plot1(1:4),'LineWidth',linewidth1);

    cols = [
        0.0000 0.4470 0.7410
        0.8500 0.3250 0.0980
        0.9290 0.6940 0.1250
        0.4940 0.1840 0.5560
        0.4660 0.6740 0.1880
        0.3010 0.7450 0.9330
        0.6350 0.0780 0.1840
        0.2000 0.2000 0.2000
    ];
    
    styles = {'-','--',':','-','--',':','-.','-'};
    markers = {'none','none','none','none','o','s','^','d'};
    
    for k = 1:8
        set(plot1(k), ...
            'LineStyle', styles{k}, ...
            'Color', cols(k,:), ...
            'LineWidth', 2.5, ...
            'Marker', markers{k}, ...
            'MarkerSize', 5);
    end
    set(plot1(5:8),'LineWidth',0.75);

    yticks_points = [0 1 2 3 4 5];
    
    yticks(yticks_points)
    
    ax = gca;
    ax.TickLabelInterpreter = 'latex';

    leg = legend({ ...
        '$x=0.1$', '$x=0.25$', '$x=0.5$', '$x=0.75$', ...
        '$x=-0.1$', '$x=-0.25$', '$x=-0.5$', '$x=-0.75$'}, ...
        'Location','best', 'Interpreter','latex');
    
    leg.FontSize = 12;
    leg.IconColumnWidth = 10;      % narrower legend box
end

% addpath('C:\Users\rikve\Github\PhD-RUG\SpatiallyAdaptiveMomentModels\Nonlinear-systems\Data-processing\Results\export_fig\', '-end');
% export_fig('AdaptiveSWME_Paper\damBreakLarge_h_newName.pdf', '-pdf','-transparent');