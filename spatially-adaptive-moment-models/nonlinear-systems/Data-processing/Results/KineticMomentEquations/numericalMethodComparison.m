clear all
close all
clc
format long

resolutionPRICE = 1600;
resolutionOsher = 800;

% name1 = 'resolution'+string(resolutionPRICE)+'numerical_scheme_PRICErelaxation_time0.05';
% name2 = 'resolution'+string(resolutionOsher)+'numerical_scheme_Osherrelaxation_time0.05';
% name3 = 'resolution3200numerical_scheme_PRICErelaxation_time0.05';
% name4 = 'resolution3200numerical_scheme_Osherrelaxation_time0.05';

name1 = 'resolution'+string(resolutionPRICE)+'numerical_scheme_PRICErelaxation_time0.5';
name2 = 'resolution'+string(resolutionOsher)+'numerical_scheme_Osherrelaxation_time0.5';
name3 = 'resolution6400numerical_scheme_PRICErelaxation_time0.5';
name4 = 'resolution6400numerical_scheme_Osherrelaxation_time0.5';


[x1,rho1,u1,T1,f31,f41,f51,f61,f71,f81,f91,f101] = readDataHME(name1);
[x2,rho2,u2,T2,f32,f42,f52,f62,f72,f82,f92,f102] = readDataHME(name2);
[x3,rho3,u3,T3,f33,f43,f53,f63,f73,f83,f93,f103] = readDataHME(name3);
[x4,rho4,u4,T4,f34,f44,f54,f64,f74,f84,f94,f104] = readDataHME(name4);

% plotting = 'rho';
plotting = 'u';
% plotting = 'T';

 blue = [0, 0.4470, 0.7410];
 green = [0.4660, 0.6740, 0.1880];
 yellow = [0.9290, 0.6940, 0.1250];
 red = [0.8500, 0.3250, 0.0980];
 brown = [171, 104, 87]./255;

if(strcmp(plotting,'rho'))
    plot1 = plot(x1,rho1,x2,rho2,x3,rho3,x4,rho4)
    
    axis([-1.5,1.5,0.9,7.4]);
    xlabel('x')
    ylabel('rho')
    
    linewidth1 = 2;
    set(plot1(1:4),'LineWidth',linewidth1);
    set(plot1(1:4),'LineStyle','-');
    
    set(plot1(1),'Color',yellow);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color',green);
    set(plot1(4),'Color',red);

    leg = legend('PRICE, nx: '+string(resolutionPRICE),'Osher, nx: '+string(resolutionOsher),'PRICE, nx: 6400','Osher, nx: 6400','Location','southwest');
    set(leg,'FontSize',7); 
end

if(strcmp(plotting,'u'))
    plot1 = plot(x1,u1,x2,u2,x3,u3,x4,u4)
    
    axis([-1.5,1.5,-0.2,0.8]);
    xlabel('x')
    ylabel('u')
    
    linewidth1 = 2;
    set(plot1(1:4),'LineWidth',linewidth1);
    set(plot1(1:4),'LineStyle','-');
    
    set(plot1(1),'Color',brown);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color',green);
    set(plot1(4),'Color',red);
    
    leg = legend('PRICE, nx: '+string(resolutionPRICE),'Osher, nx: '+string(resolutionOsher),'PRICE, nx: 6400','Osher, nx: 6400','Location','southwest');
    set(leg,'FontSize',7); 
end

if(strcmp(plotting,'T'))
    plot1 = plot(x1,T1,x2,T2,x3,T3,x4,T4)
    
    axis([-1.5,1.5,0.5,1.7]);
    xlabel('x')
    ylabel('T')
    
    linewidth1 = 2;
    set(plot1(1:4),'LineWidth',linewidth1);
    set(plot1(1:4),'LineStyle','-');
    
    set(plot1(1),'Color',brown);
    set(plot1(2),'Color',blue);
    set(plot1(3),'Color',green);
    set(plot1(4),'Color',red);
    
    leg = legend('PRICE, nx: '+string(resolutionPRICE),'Osher, nx: '+string(resolutionOsher),'PRICE, nx: 6400','Osher, nx: 6400','Location','southwest');
    set(leg,'FontSize',7); 
end

% % export_fig('models_unstable_tend0p05_a1.pdf', '-pdf','-transparent');
% export_fig('SWME_unstable_tend0p05_h.pdf', '-pdf','-transparent');