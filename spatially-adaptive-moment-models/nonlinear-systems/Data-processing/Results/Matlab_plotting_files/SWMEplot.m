clear all
close all
clc
format long

% name1 = 'Friction_R0p01_xi0p01_L'; type1 = 'ref';
% name2 = 'Friction_R0p1_xiINF_L'; type2 = 'ref';
% name3 = 'Friction_R0p1_xi0p1_L'; type3 = 'ref';
% name4 = 'Friction_R0p01_xiINF_L'; type4 = 'ref';

% name1 = 'Friction_R0p1_xi0p1_L'; type1 = 'ref';
% name2 = 'R0p1smoothSWMEN1dt0p001'; type2 = 'swme';
% name3 = 'R0p1smoothSWMEN2dt0p001'; type3 = 'swme';
% name4 = 'R0p1smoothSWMEN3dt0p001'; type4 = 'swme';

% name1 = 'Friction_R0p1_xi0p1_L'; type1 = 'ref';
% name2 = 'R0p1smoothHSWMEN3dt0p001'; type2 = 'swme';
% name3 = 'R0p1smoothHSWMEN4dt0p001'; type3 = 'swme';
% name4 = 'R0p1smoothHSWMEN5dt0p001'; type4 = 'swme';

% name1 = 'Friction_R0p1_xi0p1_L'; type1 = 'ref';
% name2 = 'R0p1smoothbetaHSWMEN3dt0p001'; type2 = 'swme';
% name3 = 'R0p1smoothbetaHSWMEN4dt0p001'; type3 = 'swme';
% name4 = 'R0p1smoothbetaHSWMEN5dt0p001'; type4 = 'swme';

% name1 = 'Friction_R0p1_xi0p1_L'; type1 = 'ref';
% name2 = 'R0p1smoothSWMEN2dt0p001'; type2 = 'swme';
% name3 = 'R0p1smoothSWMEN3dt0p001'; type3 = 'swme';
% name4 = 'R0p1smoothSWMEN4dt0p001'; type4 = 'swme';

% name1 = 'Friction_R0p1_xi0p1_L'; type1 = 'ref';
% name2 = 'R0p1smoothSWMEN3dt0p001'; type2 = 'swme';
% name3 = 'R0p1smoothHSWMEN3dt0p001'; type3 = 'swme';
% name4 = 'R0p1smoothbetaHSWMEN3dt0p001'; type4 = 'swme';

% name1 = 'Kinematic_L'; type1 = 'ref';
% name2 = 'kinematic_smoothSWMEN1dt0p001'; type2 = 'swme';
% name3 = 'kinematic_smoothSWMEN2dt0p001'; type3 = 'swme';
% name4 = 'kinematic_smoothSWMEN3dt0p001'; type4 = 'swme';

% SWME_unstable_tend0p2_a1
name1 = 'DamLargeSWME_dh5tend0p05'; type1 = 'swme';
name2 = 'Friction_R0p001_xi0p01tend0p05_DamNonHypnx160nv40'; type2 = 'ref';
name3 = 'Friction_R0p001_xi0p01tend0p05_DamNonHypnx160nv40'; type3 = 'ref';
name4 = 'Friction_R0p001_xi0p01tend0p05_DamNonHypnx160nv40'; type4 = 'ref';

% models_unstable_tend0p05_h
% name1 = 'DamLargeSWME_dh5tend0p05'; type1 = 'swme';
% name2 = 'DamLargeSWME_dh5tend0p05H'; type2 = 'swme';
% name3 = 'DamLargeSWME_dh5tend0p05betaH'; type3 = 'swme';
% name4 = 'Friction_R0p001_xi0p01tend0p05_DamNonHypnx160nv40'; type4 = 'ref';

% name1 = 'Friction_R0p001_xi0p01tend0p05_DamNonHypnx80nv80'; type1 = 'ref';
% name2 = 'Friction_R0p001_xi0p01tend0p05_DamNonHypnx160nv40'; type2 = 'ref';
% name3 = 'Friction_R0p001_xi0p01tend0p05_DamNonHypnx160nv80'; type3 = 'ref';
% name4 = 'DamLargeSWME_dh5tend0p05betaH'; type4 = 'swme';

% name1 = 'DamLargeR0p1SWME4_dh5tend0p05'; type1 = 'swme';
% name2 = 'DamLargeR0p1SWME4_dh5tend0p05H'; type2 = 'swme';
% name3 = 'DamLargeR0p1SWME4_dh5tend0p05betaH'; type3 = 'swme';
% name4 = 'Friction_R0p1_xi0p1tend0p2_Damnx160nv40'; type4 = 'swme';

% name1 = 'DamLargeR0p1SWME1_dh5tend0p05'; type1 = 'swme';
% name2 = 'DamLargeR0p1SWME2_dh5tend0p05'; type2 = 'swme';
% name3 = 'DamLargeR0p1SWME3_dh5tend0p05'; type3 = 'swme';
% name4 = 'DamLargeR0p1SWME4_dh5tend0p05'; type4 = 'swme';

% name1 = 'DamLargeR0p1SWME2_dh5tend0p05H'; type1 = 'swme';
% name2 = 'DamLargeR0p1SWME3_dh5tend0p05H'; type2 = 'swme';
% name3 = 'DamLargeR0p1SWME4_dh5tend0p05H'; type3 = 'swme';
% name4 = 'DamLargeR0p1SWME5_dh5tend0p05H'; type4 = 'swme';

% name1 = 'DamLargeR0p1SWME2_dh5tend0p05betaH'; type1 = 'swme';
% name2 = 'DamLargeR0p1SWME3_dh5tend0p05betaH'; type2 = 'swme';
% name3 = 'DamLargeR0p1SWME4_dh5tend0p05betaH'; type3 = 'swme';
% name4 = 'DamLargeR0p1SWME5_dh5tend0p05betaH'; type4 = 'swme';

% name1 = 'DamSmallSWME1_dh1p5tend0p2'; type1 = 'swme';
% name2 = 'DamSmallSWME2_dh1p5tend0p2H'; type2 = 'swme';
% name3 = 'DamSmallSWME4_dh1p5tend0p2H'; type3 = 'swme';
% name4 = 'DamSmallSWME5_dh1p5tend0p2H'; type4 = 'swme';

% name1 = 'DamSmallSWME2_dh1p5tend0p2H'; type1 = 'swme';
% name2 = 'DamSmallSWME3_dh1p5tend0p2betaH'; type2 = 'swme';
% name3 = 'DamSmallSWME4_dh1p5tend0p2betaH'; type3 = 'swme';
% name4 = 'DamSmallSWME5_dh1p5tend0p2betaH'; type4 = 'swme';

% name1 = 'DamSmallSWME4_dh1p5tend0p2'; type1 = 'swme';
% name2 = 'DamSmallSWME4_dh1p5tend0p2H'; type2 = 'swme';
% name3 = 'DamSmallSWME4_dh1p5tend0p2betaH'; type3 = 'swme';
% name4 = 'DamSmallSWME4_dh1p5tend0p2betaH'; type4 = 'swme';

[x1,h1,u1,alpha11,alpha21] = readDataSWME(name1,type1);
[x2,h2,u2,alpha12,alpha22] = readDataSWME(name2,type2);
[x3,h3,u3,alpha13,alpha23] = readDataSWME(name3,type3);
[x4,h4,u4,alpha14,alpha24] = readDataSWME(name4,type4);

% alpha13 = alpha13 ./ 3;
% alpha23 = alpha23 ./ 5;
% alpha14 = alpha14 ./ 3;
% alpha24 = alpha24 ./ 5;

% ax3 = subplot(4,2,1);
% ax4 = subplot(4,2,2);

% plotting = 'all';
plotting = 'h';
% plotting = 'u';
% plotting = 'alpha1';


if(strcmp(plotting,'all'))
    ax1 = subplot(2,2,1);
    ax2 = subplot(2,2,2);
    ax3 = subplot(2,2,3);
    ax4 = subplot(2,2,4);

    plot(ax1,x1,h1,x2,h2,x3,h3,x4,h4)
    legend(ax1,name1,name2,name3,name4)
    plot(ax2,x1,u1,x2,u2,x3,u3,x4,u4)

    plot(ax3,x1,alpha11,x2,alpha12,x3,alpha13,x4,alpha14)

    plot(ax4,x1,alpha21,x2,alpha22,x3,alpha23,x4,alpha24)
end
if(strcmp(plotting,'h'))
    plot1 = plot(x1,h1,x2,h2,x3,h3,x4,h4)
%     legend(name1,name2,name3,name4)

    axis([-0.4,0.4,0.9,5.1]);
    xlabel('x')
    ylabel('h')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(2:4),'LineStyle','.');
    set(plot1(2:4),'MarkerSize',linewidth1+12);
    set(plot1(4),'LineStyle','.');
    set(plot1(4),'MarkerSize',linewidth1+12);

    grey = [0.4,0.4,0.4];
%     set(plot1(1),'Color','b');
    set(plot1(2),'Color','k');
    set(plot1(3),'Color','k');
    set(plot1(4),'Color','k');
    
%     leg = legend('SWME','HSWME','\betaHSWME','reference','Location','southwest');
%     leg = legend('SWME','reference','Location','southwest');
    leg = legend('SWME','reference','Location','southwest');
    set(leg,'FontSize',12); 
end
if(strcmp(plotting,'u'))
    plot1 = plot(x1,u1,x2,u2,x3,u3,x4,u4)
%     legend(name1,name2,name3,name4)

    axis([-0.4,0.4,0.2,1.7]);
    xlabel('x')
    ylabel('u_m')
    
    linewidth1 = 2;
    set(plot1(1:3),'LineWidth',linewidth1);
    set(plot1(2:4),'LineStyle','.');
    set(plot1(2:4),'MarkerSize',linewidth1+12);
    set(plot1(4),'LineStyle','.');
    set(plot1(4),'MarkerSize',linewidth1+12);
    
    grey = [0.4,0.4,0.4];
%     set(plot1(1),'Color','b');
    set(plot1(2),'Color','k');
    set(plot1(3),'Color','k');
    set(plot1(4),'Color','k');
    
%     leg = legend('SWME','HSWME','\betaHSWME','reference','Location','northwest');
%     leg = legend('SWME','reference','Location','northwest');
    leg = legend('SWME','reference','Location','northwest');
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

% export_fig('models_unstable_tend0p05_a1.pdf', '-pdf','-transparent');
export_fig('SWME_unstable_tend0p05_h.pdf', '-pdf','-transparent');