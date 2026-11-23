clear all
close all
clc
format long

prediction_solver = 'PRICE';
interface_solver = 'PRICE';
relaxation_time = '0p05';

% smooth_par = 200;
% toldown = '0p001';
% tolup = '0p0015';
% interpolated = 'interp';

% smooth_pars = [50,100,200];
% toldowns = ["0p0001","0p0005","0p001"];
% tolups = ["0p00015","0p00075","0p0015"];
% interpolated = 'interp';
% plottings = ["rho","u","T","adaptiveOrders"];
order_class_low = 4;
order_class_high = 12;
smooth_pars = [100];
toldowns = ["0p0001"];
tolups = ["0p00015"];
% toldowns = ["0p001"];
% tolups = ["0p0015"];
interpolated = 'interp';
% plottings = ["rho","u","T","adaptiveOrders"];
% plottings = ["rho","u","T"];
plottings = ["comparison_all"];
for i = 1:length(smooth_pars)
    for j = 1:length(toldowns)
        for k = 1:length(plottings)
            plotting = plottings(k);
            smooth_par = smooth_pars(i);
            toldown = toldowns(j);
            tolup = tolups(j);
            name1 = strcat('shockTube','_order',string(order_class_low),'_Kn',string(relaxation_time));type1='hme';
            name2 = strcat('shockTube','_order',string(order_class_high),'_Kn',string(relaxation_time));type2='hme';
            name3 = strcat('pred',prediction_solver,'_interface',interface_solver,...
                '_smoothpar',string(smooth_par),'_toldown',toldown,'_tolup',tolup,...
                '_Kn',relaxation_time,'_',interpolated);type3='hmeAdaptive';
            
            [x1,rho1,u1,T1,f31,f41,f51,f61,f71,f81,f91,f101,f111,f121,moments1] = readDataAdaptiveHME1D(name1,type1);
            [x2,rho2,u2,T2,f32,f42,f52,f62,f72,f82,f92,f102,f112,f122,moments2] = readDataAdaptiveHME1D(name2,type2);
            [x3,rho3,u3,T3,f33,f43,f53,f63,f73,f83,f93,f103,f113,f123,moments3] = readDataAdaptiveHME1D(name3,type3);
            
            % rho = load(strcat('DVM_density','_Kn',relaxation_time,'.dat'))';
            % u = load(strcat('DVM_velocity','_Kn',relaxation_time,'.dat'))';
            % p = load(strcat('DVM_pressure','_Kn',relaxation_time,'.dat'))';
            
            xA = -2; xE = 2; deltaX = 0.002;
            x_dvm = xA+deltaX/2:deltaX:xE-deltaX/2;
            rho_dvm = load(strcat('rho','_Kn',relaxation_time,'.dat'))';
            u_dvm =  load(strcat('rho_u','_Kn',relaxation_time,'.dat'))'./load(strcat('rho','_Kn',relaxation_time,'.dat'))';
            
            E_dvm = load(strcat('E','_Kn',relaxation_time,'.dat'))';
            theta_dvm =  (E_dvm - 0.5 * rho_dvm .* u_dvm .* u_dvm)./(0.5* rho_dvm);
            
            f3 = 0*rho_dvm;
            f4 = 0*rho_dvm;
            Q = 0*rho_dvm;
            p = rho_dvm .* theta_dvm;
            
            % % plotting = 'rho';
            % % plotting = 'u';
            % % plotting = 'T';
            % % plotting = 'f4';
            % % plotting = 'adaptiveOrders';
            % plotting = 'comparison_all';
            
            blue = [0, 0.4470, 0.7410];
            green = [0.4660, 0.6740, 0.1880];
            yellow = [0.9290, 0.6940, 0.1250];
            red = [0.8500, 0.3250, 0.0980];
            brown = [171, 104, 87]./255;
            
            fig = figure;

            if(strcmp(plotting,'rho'))
                plot1 = plot(x1,rho1,x2,rho2,x3,rho3,x_dvm,rho_dvm)
                
                axis([-1.0,1.25,0.75,7.25]);
                xlabel('$x$','FontSize', 20,'Interpreter','latex')
                ylabel('$\rho$','FontSize', 20,'Interpreter','latex')
                
                linewidth1 = 3;
                set(plot1(1:3),'LineWidth',linewidth1);
                set(plot1(1),'LineStyle','-.');
                set(plot1(2),'LineStyle','--');
                set(plot1(3),'LineStyle',':');
                set(plot1(4),'LineStyle','-');
                set(plot1(4),'LineWidth',2)
            
                set(plot1(1),'Color',red);
                set(plot1(2),'Color',blue);
                set(plot1(3),'Color','k');
                set(plot1(4),'Color',brown);
                
                leg = legend('HME_{8}','HME_{12}','A-HME','DVM','Location','northeast');
                set(leg,'FontSize',12); 
            end
            
            if(strcmp(plotting,'u'))
                plot1 = plot(x1,u1,x2,u2,x3,u3,x_dvm,u_dvm)
                
                axis([-1.0,1.25,-0.024,0.75]);
                xlabel('$x$','FontSize', 20,'Interpreter','latex')
                ylabel('$u$','FontSize', 20,'Interpreter','latex')
                
                linewidth1 = 3;
                set(plot1(1:3),'LineWidth',linewidth1);
                set(plot1(1),'LineStyle','-.');
                set(plot1(2),'LineStyle','--');
                set(plot1(3),'LineStyle',':');
                set(plot1(4),'LineStyle','-');
                set(plot1(4),'LineWidth',2)
            
                set(plot1(1),'Color',red);
                set(plot1(2),'Color',blue);
                set(plot1(3),'Color','k');
                set(plot1(4),'Color',brown);
                
                leg = legend('HME_{8}','HME_{12}','A-HME','DVM','Location','northwest');
                set(leg,'FontSize',12); 
            end
            
            if(strcmp(plotting,'T'))
                plot1 = plot(x1,T1,x2,T2,x3,T3,x_dvm,theta_dvm)
                
                axis([-1.0,1.25,0.45,1.85]);
                xlabel('$x$','FontSize', 20,'Interpreter','latex')
                ylabel('$u$','FontSize', 20,'Interpreter','latex')
                
                linewidth1 = 3;
                set(plot1(1:3),'LineWidth',linewidth1);
                set(plot1(1),'LineStyle','-.');
                set(plot1(2),'LineStyle','--');
                set(plot1(3),'LineStyle',':');
                set(plot1(4),'LineStyle','-');
                set(plot1(4),'LineWidth',2)
            
                set(plot1(1),'Color',red);
                set(plot1(2),'Color',blue);
                set(plot1(3),'Color','k');
                set(plot1(4),'Color',brown);
            
                leg = legend('HME_{8}','HME_{12}','A-HME','DVM','Location','northwest');
                set(leg,'FontSize',12); 
            end
            
            if(strcmp(plotting,'f3'))
                plot1 = plot(x1,f31,x2,f32,x3,f33)
                
                axis([-1.5,1.5,-0.4,0.4]);
                xlabel('x')
                ylabel('f_3')
                
                linewidth1 = 2;
                set(plot1(1:3),'LineWidth',linewidth1);
                set(plot1(1:3),'LineStyle','-');
                
                set(plot1(1),'Color',yellow);
                set(plot1(2),'Color',blue);
                set(plot1(3),'Color',green);
                
                leg = legend('M=2','M=12','Adaptive','Location','northwest');
                set(leg,'FontSize',7); 
            end
            
            if(strcmp(plotting,'f4'))
                plot1 = plot(x1,f41,x2,f42,x3,f43)
                
                axis([-1.5,1.5,-0.1,0.1]);
                xlabel('x')
                ylabel('f_4')
                
                linewidth1 = 2;
                set(plot1(1:3),'LineWidth',linewidth1);
                set(plot1(1:3),'LineStyle','-');
                
                set(plot1(1),'Color',yellow);
                set(plot1(2),'Color',blue);
                set(plot1(3),'Color',green);
                
                leg = legend('M=2','M=12','Adaptive','Location','northwest');
                set(leg,'FontSize',7); 
            end
            
            if(strcmp(plotting,'adaptiveOrders'))
                yyaxis right
                plot1 = plot(x3,moments3)
                ylim([-0.3 12.5])
                
                ax = gca;           
                ax.YColor = 'k'; 
               
                ylabel('Order \it{M}', 'FontSize', 15, 'Color', 'k')   

                set(plot1(1),'LineStyle','None');
                set(plot1(1),'marker','.');
                set(plot1(1),'MarkerSize',20)
                set(plot1(1),'Color','k');
            
                yyaxis left
                % plot(x3,rho3,'LineWidth',2)
                plot2 = plot(x3,rho3,x3,u3,x3,T3)
                set(plot2(1:3),'LineWidth',3)
                set(plot2(1),'LineStyle','-.')
                set(plot2(2),'LineStyle','--')    
                set(plot2(3),'LineStyle',':')
            
                set(plot2(1),'Color',red)
                set(plot2(2),'Color',blue)    
                set(plot2(3),'Color',brown)
            
                % axis([-20,20,-0.5,5.5]);
                % axis([-3.5,0,3.4,3.9]);
                xlabel('$x$','FontSize',20,'Interpreter','latex')
                % ylabel('u_m')
                ylabel('\(\rho,u,\theta\)','FontSize',20,'interpreter','latex')
                ylim([-0.85 7.5])
                % ylim([-0.02 0.3])
            
                leg = legend('$\rho$','$u$','$\theta$','$M$','Location','west','interpreter','latex');
                set(leg,'FontSize',12); 
            end
            
            if(strcmp(plotting,'comparison_all'))
                plot1 = plot(...
                    x_dvm,rho_dvm,x1,rho1,x2,rho2,x3,rho3,...
                    x_dvm,u_dvm,x1,u1,x2,u2,x3,u3,...
                    x_dvm,theta_dvm,x1,T1,x2,T2,x3,T3)
                
                axis([-0.7,1.0,-0.15,7.15]);
                xlabel('$x$','FontSize', 20,'Interpreter','latex')
                ylabel('$\rho,u,\theta$','FontSize', 20,'Interpreter','latex')
                
                text(-0.65,6.65,'$\rho$','FontSize',20,'Interpreter','latex')
                text(-0.65,1.3,'$\theta$','FontSize',20,'Interpreter','latex')
                text(-0.65,0.3,'$u$','FontSize',20,'Interpreter','latex')

                linewidth1 = 3;
                set(plot1(1:12),'LineWidth',linewidth1);
                set(plot1(2),'LineStyle','-.');set(plot1(6),'LineStyle','-.');set(plot1(10),'LineStyle','-.');
                set(plot1(3),'LineStyle','--');set(plot1(7),'LineStyle','--');set(plot1(11),'LineStyle','--');
                set(plot1(4),'LineStyle',':');set(plot1(8),'LineStyle',':');set(plot1(12),'LineStyle',':');
                set(plot1(1),'LineStyle','-');set(plot1(5),'LineStyle','-');set(plot1(9),'LineStyle','-');
                % set(plot1(1),'LineWidth',3);set(plot1(5),'LineWidth',2);set(plot1(9),'LineWidth',2);

                set(plot1(2),'Color',red);set(plot1(6),'Color',red);set(plot1(10),'Color',red);
                set(plot1(3),'Color',blue);set(plot1(7),'Color',blue);set(plot1(11),'Color',blue);
                set(plot1(4),'Color','k');set(plot1(8),'Color','k');set(plot1(12),'Color','k');
                set(plot1(1),'Color',brown);set(plot1(5),'Color',brown);set(plot1(9),'Color',brown);

                leg = legend('DVM',strcat('HME_{',string(order_class_low),'}'),...
                    strcat('HME_{',string(order_class_high),'}'),'A-HME','Location','northeast');
                set(leg,'FontSize',12); 
            end

            if(strcmp(plotting,'adaptive_all'))
                plot1 = plot(x3,rho3,x3,u3,x3,T3)
                
                axis([-1.0,1.2,-0.15,7.15]);
                xlabel('$x$','FontSize', 20,'Interpreter','latex')
                ylabel('$\rho,u,\theta$','FontSize', 20,'Interpreter','latex')
                
                text(-0.95,6.65,'$\rho$','FontSize',20,'Interpreter','latex')
                text(-0.95,1.3,'$\theta$','FontSize',20,'Interpreter','latex')
                text(-0.95,0.3,'$u$','FontSize',20,'Interpreter','latex')

                linewidth1 = 2;
                set(plot1(1:3),'LineWidth',linewidth1);
                set(plot1(1:3),'LineStyle','-');


                set(plot1(1),'Color',red);
                set(plot1(2),'Color',blue);
                set(plot1(3),'Color','k');

                % leg = legend(strcat('HME_{',string(order_class_low),'}'),...
                %     strcat('HME_{',string(order_class_high),'}'),'A-HME','DVM','Location','northeast');
                % set(leg,'FontSize',12); 
            end
            
            % export_name = strcat('Paper\shockTube',plotting,'_smoothPar',string(smooth_par),...
            %     '_toldown',toldown,'_tolup',tolup,'_Kn',relaxation_time,'_combined.pdf');
            % 
            % drawnow;                          % flush pending rendering
            % fig = gcf;
            % % Vector PDF (good for Illustrator / high quality)
            % exportgraphics(fig, export_name, 'ContentType', 'vector');
            % 
            % 
            % addpath('C:\Users\rikve\Github\PhD-RUG\spatially-adaptive-moment-models\nonlinear-systems\Data-processing\Results\export_fig\', '-end');
            % export_fig(export_name, '-pdf','-transparent');
            % 
            % close(fig)
        end
    end
end