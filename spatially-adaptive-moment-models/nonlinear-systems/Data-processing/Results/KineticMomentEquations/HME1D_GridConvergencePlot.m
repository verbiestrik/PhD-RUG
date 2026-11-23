clear all
close all
clc
format long

% name1 = 'resolution100numerical_scheme_Osherrelaxation_time0.05';
% name2 = 'resolution200numerical_scheme_Osherrelaxation_time0.05';
% name3 = 'resolution400numerical_scheme_Osherrelaxation_time0.05';
% name4 = 'resolution800numerical_scheme_Osherrelaxation_time0.05';
% name5 = 'resolution1600numerical_scheme_Osherrelaxation_time0.05';
% name6 = 'resolution3200numerical_scheme_Osherrelaxation_time0.05';
% name7 = 'resolution6400numerical_scheme_Osherrelaxation_time0.05';

name1 = 'resolution100numerical_scheme_Osherrelaxation_time0.5';
name2 = 'resolution200numerical_scheme_Osherrelaxation_time0.5';
name3 = 'resolution400numerical_scheme_Osherrelaxation_time0.5';
name4 = 'resolution800numerical_scheme_Osherrelaxation_time0.5';
name5 = 'resolution1600numerical_scheme_Osherrelaxation_time0.5';
name6 = 'resolution3200numerical_scheme_Osherrelaxation_time0.5';
name7 = 'resolution6400numerical_scheme_Osherrelaxation_time0.5';

% name1 = 'resolution100numerical_scheme_PRICErelaxation_time0.05';
% name2 = 'resolution200numerical_scheme_PRICErelaxation_time0.05';
% name3 = 'resolution400numerical_scheme_PRICErelaxation_time0.05';
% name4 = 'resolution800numerical_scheme_PRICErelaxation_time0.05';
% name5 = 'resolution1600numerical_scheme_PRICErelaxation_time0.05';
% name6 = 'resolution3200numerical_scheme_PRICErelaxation_time0.05';
% name7 = 'resolution6400numerical_scheme_PRICErelaxation_time0.05';

% name1 = 'resolution100numerical_scheme_PRICErelaxation_time0.5';
% name2 = 'resolution200numerical_scheme_PRICErelaxation_time0.5';
% name3 = 'resolution400numerical_scheme_PRICErelaxation_time0.5';
% name4 = 'resolution800numerical_scheme_PRICErelaxation_time0.5';
% name5 = 'resolution1600numerical_scheme_PRICErelaxation_time0.5';
% name6 = 'resolution3200numerical_scheme_PRICErelaxation_time0.5';
% name7 = 'resolution6400numerical_scheme_PRICErelaxation_time0.5';

[x1,rho1,u1,T1,f31,f41,f51,f61,f71,f81,f91,f101] = readDataHME(name1);
[x2,rho2,u2,T2,f32,f42,f52,f62,f72,f82,f92,f102] = readDataHME(name2);
[x3,rho3,u3,T3,f33,f43,f53,f63,f73,f83,f93,f103] = readDataHME(name3);
[x4,rho4,u4,T4,f34,f44,f54,f64,f74,f84,f94,f104] = readDataHME(name4);
[x5,rho5,u5,T5,f35,f45,f55,f65,f75,f85,f95,f105] = readDataHME(name5);
[x6,rho6,u6,T6,f36,f46,f56,f66,f76,f86,f96,f106] = readDataHME(name6);
[x7,rho7,u7,T7,f37,f47,f57,f67,f77,f87,f97,f107] = readDataHME(name7);

% plotting = 'rho';
% plotting = 'u';
plotting = 'T';

 blue = [0, 0.4470, 0.7410];
 green = [0.4660, 0.6740, 0.1880];
 yellow = [0.9290, 0.6940, 0.1250];
 red = [0.8500, 0.3250, 0.0980];
 brown = [171, 104, 87]./255;

if(strcmp(plotting,'rho'))
    plot1 = plot(x1,rho1,x2,rho2,x3,rho3,x4,rho4,x5,rho5,x6,rho6,x7,rho7)
    
    axis([-1.5,1.5,0.9,7.4]);
    xlabel('x')
    ylabel('rho')
    
    linewidth1 = 2;
    set(plot1(1:7),'LineWidth',linewidth1);
    set(plot1(1:7),'LineStyle','-');
    
    set(plot1(1),'Color',yellow);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color',green);
    set(plot1(4),'Color',red);
    set(plot1(5),'Color','b');
    set(plot1(6),'Color','cyan');
    set(plot1(7),'Color',brown);
    
    leg = legend('100','200','400','800','1600','3200','6400','Location','southwest');
    set(leg,'FontSize',7); 
end

if(strcmp(plotting,'u'))
    plot1 = plot(x1,u1,x2,u2,x3,u3,x4,u4,x5,u5,x6,u6,x7,u7)
    
    axis([-1.5,1.5,-0.2,0.8]);
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
    set(plot1(7),'Color',brown);
    
    leg = legend('100','200','400','800','1600','3200','6400','Location','southwest');
    set(leg,'FontSize',7); 
end

if(strcmp(plotting,'T'))
    plot1 = plot(x1,T1,x2,T2,x3,T3,x4,T4,x5,T5,x6,T6,x7,T7)
    
    axis([-1.5,1.5,0.3,1.9]);
    xlabel('x')
    ylabel('T')
    
    linewidth1 = 2;
    set(plot1(1:7),'LineWidth',linewidth1);
    set(plot1(1:7),'LineStyle','-');
    
    set(plot1(1),'Color',yellow);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color',green);
    set(plot1(4),'Color',red);
    set(plot1(5),'Color','b');
    set(plot1(6),'Color','cyan');
    set(plot1(7),'Color',brown);
    
    leg = legend('100','200','400','800','1600','3200','6400','Location','southwest');
    set(leg,'FontSize',7); 
end

% % export_fig('models_unstable_tend0p05_a1.pdf', '-pdf','-transparent');
% export_fig('SWME_unstable_tend0p05_h.pdf', '-pdf','-transparent');