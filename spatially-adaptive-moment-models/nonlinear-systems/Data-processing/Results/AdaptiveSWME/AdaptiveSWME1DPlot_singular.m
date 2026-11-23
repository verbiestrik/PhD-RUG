clear all
close all
clc
format long

name = 'damBreak-and-smooth_linear_init'; type = 'swme';
[x,h,u,alpha1,alpha2,alpha3,alpha4,alpha5,moments] = readDataAdaptiveSWME1D(name,type);

plotting = 'h';

brown = [171, 104, 87]./255;

if(strcmp(plotting,'h'))
    plot1 = plot(x,h)
    
    axis([-20,20,2.8,4.2]);
    xlabel('x')
    ylabel('h')
    
    linewidth1 = 2;
    set(plot1(1),'LineWidth',linewidth1);
end

% addpath('C:\Users\rikve\Github\PhD-RUG\spatially-adaptive-moment-models\nonlinear-systems\Data-processing\Results\export_fig\', '-end')
% export_fig('AdaptiveSWME_Paper\damBreak-and-smooth_linear_init.pdf', '-pdf','-transparent');