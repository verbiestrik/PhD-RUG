clear all
close all
clc
format long

% name1 = 'smoothAndShockTube_order2_relaxation2.0_time1.0';
% name2 = 'smoothAndShockTube_order3_relaxation2.0_time1.0';
% name3 = 'smoothAndShockTube_order4_relaxation2.0_time1.0';
% name4 = 'smoothAndShockTube_order5_relaxation2.0_time1.0';
% name5 = 'smoothAndShockTube_order6_relaxation2.0_time1.0';
% name6 = 'smoothAndShockTube_order7_relaxation2.0_time1.0';
% name7 = 'smoothAndShockTube_order8_relaxation2.0_time1.0';
% name8 = 'smoothAndShockTube_order9_relaxation2.0_time1.0';
% name9 = 'smoothAndShockTube_order10_relaxation2.0_time1.0';

name1 = 'smoothAndShockTube_order8_relaxation0.1_time1.0';
name2 = 'smoothAndShockTube_order10_relaxation0.1_time1.0_2000';
name3 = 'smoothAndShockTube_order10_relaxation0.1_time1.0_3000';
name4 = 'smoothAndShockTube_order10_relaxation0.1_time1.0_4000';
name5 = 'smoothAndShockTube_order10_relaxation0.1_time1.0_5000';
name6 = 'smoothAndShockTube_order10_relaxation0.1_time1.0_6000';
name7 = 'smoothAndShockTube_order10_relaxation0.1_time1.0_7000';
name8 = 'smoothAndShockTube_order10_relaxation0.1_time1.0_8000';
name9 = 'smoothAndShockTube_order10_relaxation0.1_time1.0_9000';
name10 = 'smoothAndShockTube_order10_relaxation0.1_time1.0_10000';

% name1 = 'smoothAndShockTube_order2_relaxation1.0_time1.0_PRICE';
% name2 = 'smoothAndShockTube_order3_relaxation1.0_time1.0_PRICE';
% name3 = 'smoothAndShockTube_order4_relaxation1.0_time1.0_PRICE';
% name4 = 'smoothAndShockTube_order5_relaxation1.0_time1.0_PRICE';
% name5 = 'smoothAndShockTube_order6_relaxation1.0_time1.0_PRICE';
% name6 = 'smoothAndShockTube_order7_relaxation1.0_time1.0_PRICE';
% name7 = 'smoothAndShockTube_order8_relaxation1.0_time1.0_PRICE';
% name8 = 'smoothAndShockTube_order9_relaxation1.0_time1.0_PRICE';
% name9 = 'smoothAndShockTube_order10_relaxation1.0_time1.0_PRICE';

[x1,rho1,u1,T1,f31,f41,f51,f61,f71,f81,f91,f101] = readDataHME(name2);
[x2,rho2,u2,T2,f32,f42,f52,f62,f72,f82,f92,f102] = readDataHME(name2);
[x3,rho3,u3,T3,f33,f43,f53,f63,f73,f83,f93,f103] = readDataHME(name3);
[x4,rho4,u4,T4,f34,f44,f54,f64,f74,f84,f94,f104] = readDataHME(name4);
[x5,rho5,u5,T5,f35,f45,f55,f65,f75,f85,f95,f105] = readDataHME(name5);
[x6,rho6,u6,T6,f36,f46,f56,f66,f76,f86,f96,f106] = readDataHME(name6);
[x7,rho7,u7,T7,f37,f47,f57,f67,f77,f87,f97,f107] = readDataHME(name7);
[x8,rho8,u8,T8,f38,f48,f58,f68,f78,f88,f98,f108] = readDataHME(name8);
[x9,rho9,u9,T9,f39,f49,f59,f69,f79,f89,f99,f109] = readDataHME(name9);
[x10,rho10,u10,T10,f310,f410,f510,f610,f710,f810,f910,f1010] = readDataHME(name10);

% ax3 = subplot(4,2,1);
% ax4 = subplot(4,2,2);

% plotting = 'rho';
% plotting = 'u';
plotting = 'T';


 blue = [0, 0.4470, 0.7410];
 green = [0.4660, 0.6740, 0.1880];
 yellow = [0.9290, 0.6940, 0.1250];
 red = [0.8500, 0.3250, 0.0980];
 brown = [171, 104, 87]./255;

if(strcmp(plotting,'rho'))
    % plot1 = plot(x5,rho6,x7,rho7,x8,rho8,x9,rho9)
    plot1 = plot(x2,rho2,x3,rho3,x4,rho4,x5,rho5,x6,rho6,x7,rho7,x8,rho8,x9,rho9,x10,rho10)
%     legend(name1,name2,name3,name4)
    
    axis([-2.0,2.0,0.9,2.4]);
    xlabel('x')
    ylabel('rho')
    
    linewidth1 = 2;
    set(plot1(1:9),'LineWidth',linewidth1);
    set(plot1(1:9),'LineStyle','-');
    
    set(plot1(1),'Color',yellow);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color',green);
    set(plot1(4),'Color',red);
    set(plot1(5),'Color','b');
    set(plot1(6),'Color','cyan');
    set(plot1(7),'Color','magenta');
    set(plot1(8),'Color','green');
    set(plot1(9),'Color','black');
    
%     leg = legend('SWME','HSWME','\betaHSWME','reference','Location','southwest');
%     leg = legend('SWME','reference','Location','southwest');
    % leg = legend('M=7','M=8','M=9','M=10','Location','northwest');
    leg = legend('M=8','2000','3000','4000','5000','6000','7000','8000','9000','Location','northwest');
    set(leg,'FontSize',7); 
end

if(strcmp(plotting,'u'))
    % plot1 = plot(x5,rho6,x7,rho7,x8,rho8,x9,rho9)
    plot1 = plot(x2,u2,x3,u3,x4,u4,x5,u5,x6,u6,x7,u7,x8,u8,x9,u9,x10,u10)
%     legend(name1,name2,name3,name4)
    
    axis([-2.0,2.0,-0.2,0.2]);
    xlabel('x')
    ylabel('u')
    
    linewidth1 = 2;
    set(plot1(1:9),'LineWidth',linewidth1);
    set(plot1(1:9),'LineStyle','-');
    
    set(plot1(1),'Color',yellow);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color',green);
    set(plot1(4),'Color',red);
    set(plot1(5),'Color','b');
    set(plot1(6),'Color','cyan');
    set(plot1(7),'Color','magenta');
    set(plot1(8),'Color','green');
    set(plot1(9),'Color','black');
    
%     leg = legend('SWME','HSWME','\betaHSWME','reference','Location','southwest');
%     leg = legend('SWME','reference','Location','southwest');
    % leg = legend('M=7','M=8','M=9','M=10','Location','northwest');
    leg = legend('M=8','2000','3000','4000','5000','6000','7000','8000','9000','Location','northwest');
    set(leg,'FontSize',7); 
end

if(strcmp(plotting,'T'))
    % plot1 = plot(x5,rho6,x7,rho7,x8,rho8,x9,rho9)
    plot1 = plot(x2,T2,x3,T3,x4,T4,x5,T5,x6,T6,x7,T7,x8,T8,x9,T9,x10,T10,x10,T10)
%     legend(name1,name2,name3,name4)
    
    axis([-2.0,2.0,0.5,1.5]);
    xlabel('x')
    ylabel('T')
    
    linewidth1 = 2;
    set(plot1(1:10),'LineWidth',linewidth1);
    set(plot1(1:10),'LineStyle','-');
    
    set(plot1(1),'Color',yellow);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color',green);
    set(plot1(4),'Color',red);
    set(plot1(5),'Color','b');
    set(plot1(6),'Color','cyan');
    set(plot1(7),'Color','magenta');
    set(plot1(8),'Color','green');
    set(plot1(9),'Color','black');
    set(plot1(10),'Color','yellow');
    
%     leg = legend('SWME','HSWME','\betaHSWME','reference','Location','southwest');
%     leg = legend('SWME','reference','Location','southwest');
    % leg = legend('M=7','M=8','M=9','M=10','Location','northwest');
    leg = legend('M=8','2000','3000','4000','5000','6000','7000','8000','9000','10000','Location','southwest');
    set(leg,'FontSize',7); 
end

% % export_fig('models_unstable_tend0p05_a1.pdf', '-pdf','-transparent');
% export_fig('SWME_unstable_tend0p05_h.pdf', '-pdf','-transparent');