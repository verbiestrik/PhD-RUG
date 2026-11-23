%% Clear statements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% These are examples of 2D ASWME plot setups %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
format long

%% Simulation files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Specify the used simulation files %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%name1 = 'RadialDamBreak_NoVelocity_Time0.025_ThirdOrder'; type1 = 'swme';
%name2 = 'RadialDamBreak_NoVelocity_Time0.5_ThirdOrder'; type2 = 'swme';
%name3 = 'RadialDamBreak_NoVelocity_Time0_ThirdOrder'; type3 = 'swme';
%name4 = 'RadialDamBreak_ConstantVelocity_Time0_ThirdOrder'; type4 = 'swme';

%name1 = 'SmoothExpit_ASWMEPaper_Time0.5_FourthOrder'; type1 = 'swme';
%name2 = 'SmoothExpit_ASWMEPaper_Time0.5_FourthOrder_HSWME'; type2 = 'swme';
%name3 = 'SmoothExpit_ASWMEPaper_Time0_ASWE'; type3 = 'swe';
%name4 = 'SmoothExpit_ASWMEPaper_Time0_ASWE'; type4 = 'swe';

%name1 = 'SmoothExpit_ASWMEPaper_Time0.5_ASWE_5000x4'; type1 = 'swe';
%name2 = 'SmoothExpit_ASWMEPaper_Time0.5_FirstOrder_5000x4'; type2 = 'swme';
%name3 = 'SmoothExpit_ASWMEPaper_Time0.5_SecondOrder_5000x4'; type3 = 'swme';
%name4 = 'SmoothExpit_ASWMEPaper_Time0.5_ThirdOrder_5000x4'; type4 = 'swme';

%name1 = 'SmoothExpit_ASWMEPaper_Time0.5_ASWE'; type1 = 'swe';
%name2 = 'SmoothExpit_ASWMEPaper_Time0.5_FirstOrder'; type2 = 'swme';
%name3 = 'SmoothExpit_ASWMEPaper_Time0.5_SecondOrder'; type3 = 'swme';
%name4 = 'SmoothExpit_ASWMEPaper_Time0.5_ThirdOrder'; type4 = 'swme';

%name1 = 'SmoothExpit_ASWMEPaper_Time0.5_ASWE'; type1 = 'swe';
%name2 = 'SmoothExpit_ASWMEPaper_Time0.5_FirstOrder_Auto'; type2 = 'swme';
%name3 = 'SmoothExpit_ASWMEPaper_Time0.5_SecondOrder_Auto'; type3 = 'swme';
%name4 = 'SmoothExpit_ASWMEPaper_Time0.5_ThirdOrder_Auto'; type4 = 'swme';

%name1 = 'RadialDamBreak_NoVelocity_Time0.5_ThirdOrder'; type1 = 'swme';
%name2 = 'RadialDamBreak_NoVelocity_Time0.5_ThirdOrder_HSWME'; type2 = 'swme';
%name3 = 'RadialDamBreak_NoVelocity_Time0_ThirdOrder'; type3 = 'swme';
%name4 = 'RadialDamBreak_NoVelocity_Time0_ThirdOrder'; type4 = 'swme';

%name1 = 'RadialDamBreak_ASWMEPaper_Time0.15_ThirdOrder'; type1 = 'swme';
%name2 = 'RadialDamBreak_ASWMEPaper_Time0.2_ThirdOrder'; type2 = 'swme';
%name3 = 'RadialDamBreak_ASWMEPaper_Time0.25_ThirdOrder'; type3 = 'swme';
%name4 = 'RadialDamBreak_ASWMEPaper_Time0.3_ThirdOrder'; type4 = 'swme';

%name1 = 'RadialDamBreak_AngularVelocity_ZeroethOrder_Time0.5_1'; type1 = 'swe';
%name2 = 'RadialDamBreak_AngularVelocity_FirstOrder_Time0.5_1'; type2 = 'swme';
%name3 = 'RadialDamBreak_AngularVelocity_SecondOrder_Time0.5_1'; type3 = 'swme';
%name4 = 'RadialDamBreak_AngularVelocity_ThirdOrder_Time0.5_1'; type4 = 'swme';

%name1 = 'ASWMEPaperRevised_RadialDamBreak_Unstable_AngularVelocity_ThirdOrder_SWME_Time0.3'; type1 = 'swme';
%name2 = 'ASWMEPaperRevised_RadialDamBreak_Unstable_AngularVelocity_ThirdOrder_HSWME_Time0.5'; type2 = 'swme';
%name2 = 'test_1'; type2 = 'swme';
%name3 = 'ASWMEPaperRevised_RadialDamBreak_Unstable_AngularVelocity_ThirdOrder_SWME_Time0.025'; type3 = 'swme';
%name4 = 'ASWMEPaperRevised_RadialDamBreak_Unstable_AngularVelocity_ThirdOrder_HSWME_Time0.025'; type4 = 'swme';

%name1 = 'SmoothExpit_InitialAngularVelocity_ThirdOrder_Time0.5_1'; type1 = 'swme';
%name2 = 'SmoothExpit_InitialAngularVelocity_FirstOrder_Time0.5_1'; type2 = 'swme';
%name3 = 'ASWMEPaperRevised_RadialDamBreak_Unstable_AngularVelocity_ThirdOrder_SWME_Time0.025'; type3 = 'swme';
%name4 = 'ASWMEPaperRevised_RadialDamBreak_Unstable_AngularVelocity_ThirdOrder_HSWME_Time0.025'; type4 = 'swme';






name1 = 'ASWMEPaperRevised_RadialDamBreak_Unstable_AngularVelocity26_ThirdOrder_SWME_Time0.3'; type1 = 'swme';
name2 = 'ASWMEPaperRevised_RadialDamBreak_Unstable_AngularVelocity26_ThirdOrder_HSWME_Time0.3'; type2 = 'swme';
name3 = 'ASWMEPaperRevised_RadialDamBreak_Unstable_AngularVelocity26_ThirdOrder_HSWME_Time0.3'; type3 = 'swme';
name4 = 'ASWMEPaperRevised_RadialDamBreak_Unstable_AngularVelocity26_ThirdOrder_HSWME_Time0.3'; type4 = 'swme';

%name1 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_ZeroethOrder_Time1.0'; type1 = 'swe';
%name2 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_FirstOrder_Time1.0'; type2 = 'swme';
%name3 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_SecondOrder_SWME_Time1.0'; type3 = 'swme';
%name4 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_ThirdOrder_SWME_Time1.0'; type4 = 'swme';

%name1 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_ZeroethOrder_Time1.0_4000'; type1 = 'swe';
%name2 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_FirstOrder_Time1.0_4000'; type2 = 'swme';
%name3 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_SecondOrder_HSWME_Time1.0_4000'; type3 = 'swme';
%name4 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_ThirdOrder_HSWME_Time1.0_4000'; type4 = 'swme';

% Moment model data

[r1,theta1,h1,vr1,alpha11,alpha21,alpha31,vtheta1,gamma11,gamma21,gamma31] = readDataASWME(name1,type1);
[r2,theta2,h2,vr2,alpha12,alpha22,alpha32,vtheta2,gamma12,gamma22,gamma32] = readDataASWME(name2,type2);
[r3,theta3,h3,vr3,alpha13,alpha23,alpha33,vtheta3,gamma13,gamma23,gamma33] = readDataASWME(name3,type3);
[r4,theta4,h4,vr4,alpha14,alpha24,alpha34,vtheta4,gamma14,gamma24,gamma34] = readDataASWME(name4,type4);

% Reference data

%reference_data = load("RadialDamBreakUnstable_ASWMEPaper_Time0.1_Reference_400x100_CFL0.5.csv");
%reference_data = load("SmoothExpit_ASWMEPaper_Time0.5_Reference.csv");
%reference_data = load("RadialDamBreakUnstable_ReferenceTest2.csv");
%reference_data = load("ASWMEPaperRevised_RadialDamBreakUnstable_Time0.3_Reference_400x100_CFL0.5_New.csv");
%reference_data=load("referenceTestRevised.csv");
%reference_data=load("ASWMEPaperRevised_SmoothExpit[1,11]_Time0.5_Reference_400x100_CFL0.5.csv");
%reference_data=load("ASWMEPaperRevised_SmoothExpit[1,8]_Time1.0_Reference_200x100_CFL0.5.csv");
%reference_data=load("ASWMEPaperRevised_SmoothExpit[1,8]_Time1.0_Reference_400x200_CFL0.5_New.csv");
reference_data=load("ASWMEPaperRevised_RadialDamBreakUnstable[2,6]_Time0.3_Reference_200x100_CFL0.5.csv");



r_ref = reference_data(:,1);
h_ref = reference_data(:,2);
vr_ref = reference_data(:,3);
alpha1_ref = reference_data(:,4);
alpha1_ref = 3.*alpha1_ref;
alpha2_ref = reference_data(:,5);
alpha2_ref = 5.*alpha2_ref;
vtheta_ref = reference_data(:,6);

%% Grid information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Operations on the specified grid %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grid = '100x4Axisymmetric';
% grid = '100x4Axisymmetric2';
% grid = '100x10Axisymmetric';
% grid = '100x10Axisymmetric2';
% grid = '10000x4Axisymmetric';
% grid = '5000x4Axisymmetric';
% grid = '1000x4Axisymmetric';
 grid = '2000x4Axisymmetric';
% grid = '4000x4Axisymmetric';

if(strcmp(grid,'100x4Axisymmetric'))
    sameAngle = 100;
end

if(strcmp(grid,'100x4Axisymmetric2'))
    sameAngle = 100;
end

if(strcmp(grid,'100x10Axisymmetric'))
    sameAngle = 100;
end

if(strcmp(grid,'100x10Axisymmetric2'))
    sameAngle = 100;
end

if(strcmp(grid,'10000x4Axisymmetric'))
    sameAngle = 10000;
end

if(strcmp(grid,'1000x4Axisymmetric'))
    sameAngle = 1000;
end

if(strcmp(grid,'5000x4Axisymmetric'))
    sameAngle = 5000;
end

if(strcmp(grid,'2000x4Axisymmetric'))
    sameAngle = 2000;
end

if(strcmp(grid,'4000x4Axisymmetric'))
    sameAngle = 4000;
end

r1_sameAngle = r1(1:sameAngle);
r2_sameAngle = r2(1:sameAngle);
r3_sameAngle = r3(1:sameAngle);
r4_sameAngle = r4(1:sameAngle);

theta1_sameAngle = theta1(1:sameAngle);
theta2_sameAngle = theta2(1:sameAngle);
theta3_sameAngle = theta3(1:sameAngle);
theta4_sameAngle = theta4(1:sameAngle);

h1_sameAngle = h1(1:sameAngle);
h2_sameAngle = h2(1:sameAngle);
h3_sameAngle = h3(1:sameAngle);
h4_sameAngle = h4(1:sameAngle);

vr1_sameAngle = vr1(1:sameAngle);vr1_sameAngle(1)=vr1_sameAngle(2);
vr2_sameAngle = vr2(1:sameAngle);vr2_sameAngle(1)=vr2_sameAngle(2);
vr3_sameAngle = vr3(1:sameAngle);vr3_sameAngle(1)=vr3_sameAngle(2);
vr4_sameAngle = vr4(1:sameAngle);vr4_sameAngle(1)=vr4_sameAngle(2);

alpha11_sameAngle = alpha11(1:sameAngle);
alpha12_sameAngle = alpha12(1:sameAngle);
alpha13_sameAngle = alpha13(1:sameAngle);
alpha14_sameAngle = alpha14(1:sameAngle);

alpha21_sameAngle = alpha21(1:sameAngle);
alpha22_sameAngle = alpha22(1:sameAngle);
alpha23_sameAngle = alpha23(1:sameAngle);
alpha24_sameAngle = alpha24(1:sameAngle);

alpha31_sameAngle = alpha31(1:sameAngle);
alpha32_sameAngle = alpha32(1:sameAngle);
alpha33_sameAngle = alpha33(1:sameAngle);
alpha34_sameAngle = alpha34(1:sameAngle);

vtheta1_sameAngle = vtheta1(1:sameAngle);vtheta1_sameAngle(1)=vtheta1_sameAngle(2);
vtheta2_sameAngle = vtheta2(1:sameAngle);vtheta2_sameAngle(1)=vtheta2_sameAngle(2);
vtheta3_sameAngle = vtheta3(1:sameAngle);vtheta3_sameAngle(1)=vtheta3_sameAngle(2);
vtheta4_sameAngle = vtheta4(1:sameAngle);vtheta4_sameAngle(1)=vtheta4_sameAngle(2);

gamma11_sameAngle = gamma11(1:sameAngle);
gamma12_sameAngle = gamma12(1:sameAngle);
gamma13_sameAngle = gamma13(1:sameAngle);
gamma14_sameAngle = gamma14(1:sameAngle);

gamma21_sameAngle = gamma21(1:sameAngle);
gamma22_sameAngle = gamma22(1:sameAngle);
gamma23_sameAngle = gamma23(1:sameAngle);
gamma24_sameAngle = gamma24(1:sameAngle);

gamma31_sameAngle = gamma31(1:sameAngle);
gamma32_sameAngle = gamma32(1:sameAngle);
gamma33_sameAngle = gamma33(1:sameAngle);
gamma34_sameAngle = gamma34(1:sameAngle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alpha13 = alpha13 ./ 3;
% alpha23 = alpha23 ./ 5;
% alpha14 = alpha14 ./ 3;
% alpha24 = alpha24 ./ 5;

% ax3 = subplot(4,2,1);
% ax4 = subplot(4,2,2);

%% Plotting
 blue = [0, 0.4470, 0.7410];
 green = [0.4660, 0.6740, 0.1880];
 yellow = [0.9290, 0.6940, 0.1250];
 red = [0.8500, 0.3250, 0.0980];
 brown = [171, 104, 87]./255;

% plotting = 'all';
 plotting = 'h';
% plotting = 'vr';
% plotting = 'alpha1';
% plotting = 'alpha2';
% plotting = 'alpha3';
% plotting = 'vtheta';
% plotting = 'gamma1';
% plotting = 'gamma2';
% plotting = 'gamma3';

% plotting_mode = 'time_evolution';
 plotting_mode = 'model_comparison';
% plotting_mode = 'order_comparison';

if(strcmp(plotting,'all'))
    ax1 = subplot(3,3,1);
    ax2 = subplot(3,3,2);
    ax3 = subplot(3,3,3);
    ax4 = subplot(3,3,4);
    ax5 = subplot(3,3,5);
    ax6 = subplot(3,3,6);
    ax7 = subplot(3,3,7);
    ax8 = subplot(3,3,8);
    ax9 = subplot(3,3,9);

    if(strcmp(plotting_mode,'time_evolution'))

        subplot1 = plot(ax1,r1_sameAngle,h1_sameAngle,r2_sameAngle,h2_sameAngle, ...
            r3_sameAngle,h3_sameAngle,r4_sameAngle,h4_sameAngle)
        legend(ax1,'t=0','t=0.025','t=0.1','t=0.25') % Add legend to othr subplots too if necessary

        subplot2 = plot(ax2,r1_sameAngle,vr1_sameAngle,r2_sameAngle,vr2_sameAngle, ...
            r3_sameAngle,vr3_sameAngle,r4_sameAngle,vr4_sameAngle)

        subplot3 = plot(ax3,r1_sameAngle,alpha11_sameAngle,r2_sameAngle,alpha12_sameAngle, ...
            r3_sameAngle,alpha13_sameAngle,r4_sameAngle,alpha14_sameAngle)

        subplot4 = plot(ax4,r1_sameAngle,alpha21_sameAngle,r2_sameAngle,alpha22_sameAngle, ...
            r3_sameAngle,alpha23_sameAngle,r4_sameAngle,alpha24_sameAngle)

        subplot5 = plot(ax5,r1_sameAngle,alpha31_sameAngle,r2_sameAngle,alpha32_sameAngle, ...
            r3_sameAngle,alpha33_sameAngle,r4_sameAngle,alpha34_sameAngle)

        subplot6 = plot(ax6,r1_sameAngle,vtheta1_sameAngle,r2_sameAngle,vtheta2_sameAngle, ...
            r3_sameAngle,vtheta3_sameAngle,r4_sameAngle,vtheta4_sameAngle)

        subplot7 = plot(ax7,r1_sameAngle,gamma11_sameAngle,r2_sameAngle,gamma12_sameAngle, ...
            r3_sameAngle,gamma13_sameAngle,r4_sameAngle,gamma14_sameAngle)

        subplot8 = plot(ax8,r1_sameAngle,gamma21_sameAngle,r2_sameAngle,gamma22_sameAngle, ...
            r3_sameAngle,gamma23_sameAngle,r4_sameAngle,gamma24_sameAngle)

        subplot9 = plot(ax9,r1_sameAngle,gamma31_sameAngle,r2_sameAngle,gamma32_sameAngle, ...
            r3_sameAngle,gamma33_sameAngle,r4_sameAngle,gamma34_sameAngle)

        set(subplot1(1),'Color',blue);
        set(subplot1(2),'Color',green);
        set(subplot1(3),'Color',yellow);
        set(subplot1(4),'Color',red);

        set(subplot2(1),'Color',blue);
        set(subplot2(2),'Color',green);
        set(subplot2(3),'Color',yellow);
        set(subplot2(4),'Color',red);

        set(subplot3(1),'Color',blue);
        set(subplot3(2),'Color',green);
        set(subplot3(3),'Color',yellow);
        set(subplot3(4),'Color',red);

        set(subplot4(1),'Color',blue);
        set(subplot4(2),'Color',green);
        set(subplot4(3),'Color',yellow);
        set(subplot4(4),'Color',red);

        set(subplot5(1),'Color',blue);
        set(subplot5(2),'Color',green);
        set(subplot5(3),'Color',yellow);
        set(subplot5(4),'Color',red);

        set(subplot6(1),'Color',blue);
        set(subplot6(2),'Color',green);
        set(subplot6(3),'Color',yellow);
        set(subplot6(4),'Color',red);

        set(subplot7(1),'Color',blue);
        set(subplot7(2),'Color',green);
        set(subplot7(3),'Color',yellow);
        set(subplot7(4),'Color',red);

        set(subplot8(1),'Color',blue);
        set(subplot8(2),'Color',green);
        set(subplot8(3),'Color',yellow);
        set(subplot8(4),'Color',red);

        set(subplot9(1),'Color',blue);
        set(subplot9(2),'Color',green);
        set(subplot9(3),'Color',yellow);
        set(subplot9(4),'Color',red);
    end
    if(strcmp(plotting_mode,'model_comparison'))
        subplot1 = plot(ax1,r1_sameAngle,h1_sameAngle,r2_sameAngle,h2_sameAngle)
        legend(ax1,'ASWME','HASWME') % Add legend to other subplots too if necessary
        subplot2 = plot(ax2,r1_sameAngle,vr1_sameAngle,r2_sameAngle,vr2_sameAngle)

        subplot3 = plot(ax3,r1_sameAngle,alpha11_sameAngle,r2_sameAngle,alpha12_sameAngle)

        subplot4 = plot(ax4,r1_sameAngle,alpha21_sameAngle,r2_sameAngle,alpha22_sameAngle)

        subplot5 = plot(ax5,r1_sameAngle,alpha31_sameAngle,r2_sameAngle,alpha32_sameAngle)

        subplot6 = plot(ax6,r1_sameAngle,vtheta1_sameAngle,r2_sameAngle,vtheta2_sameAngle)

        subplot7 = plot(ax7,r1_sameAngle,gamma11_sameAngle,r2_sameAngle,gamma12_sameAngle)

        subplot8 = plot(ax8,r1_sameAngle,gamma21_sameAngle,r2_sameAngle,gamma22_sameAngle)

        subplot9 = plot(ax9,r1_sameAngle,gamma31_sameAngle,r2_sameAngle,gamma32_sameAngle)

        set(subplot1(1),'Color',blue);
        set(subplot1(2),'Color',red);

        set(subplot2(1),'Color',blue);
        set(subplot2(2),'Color',red);

        set(subplot3(1),'Color',blue);
        set(subplot3(2),'Color',red);

        set(subplot4(1),'Color',blue);
        set(subplot4(2),'Color',red);

        set(subplot5(1),'Color',blue);
        set(subplot5(2),'Color',red);

        set(subplot6(1),'Color',blue);
        set(subplot6(2),'Color',red);

        set(subplot7(1),'Color',blue);
        set(subplot7(2),'Color',red);

        set(subplot8(1),'Color',blue);
        set(subplot8(2),'Color',red);

        set(subplot9(1),'Color',blue);
        set(subplot9(2),'Color',red);
    end
    if(strcmp(plotting_mode,'order_comparison'))
        subplot1 = plot(ax1,r1_sameAngle,h1_sameAngle,r2_sameAngle,h2_sameAngle, ...
            r3_sameAngle,h3_sameAngle,r4_sameAngle,h4_sameAngle)
        legend(ax1,'N=0','N=1','N=2','N=3') % Add legend to othr subplots too if necessary

        subplot2 = plot(ax2,r1_sameAngle,vr1_sameAngle,r2_sameAngle,vr2_sameAngle, ...
            r3_sameAngle,vr3_sameAngle,r4_sameAngle,vr4_sameAngle)

        subplot3 = plot(ax3,r1_sameAngle,alpha11_sameAngle,r2_sameAngle,alpha12_sameAngle, ...
            r3_sameAngle,alpha13_sameAngle,r4_sameAngle,alpha14_sameAngle)

        subplot4 = plot(ax4,r1_sameAngle,alpha21_sameAngle,r2_sameAngle,alpha22_sameAngle, ...
            r3_sameAngle,alpha23_sameAngle,r4_sameAngle,alpha24_sameAngle)

        subplot5 = plot(ax5,r1_sameAngle,alpha31_sameAngle,r2_sameAngle,alpha32_sameAngle, ...
            r3_sameAngle,alpha33_sameAngle,r4_sameAngle,alpha34_sameAngle)

        subplot6 = plot(ax6,r1_sameAngle,vtheta1_sameAngle,r2_sameAngle,vtheta2_sameAngle, ...
            r3_sameAngle,vtheta3_sameAngle,r4_sameAngle,vtheta4_sameAngle)

        subplot7 = plot(ax7,r1_sameAngle,gamma11_sameAngle,r2_sameAngle,gamma12_sameAngle, ...
            r3_sameAngle,gamma13_sameAngle,r4_sameAngle,gamma14_sameAngle)

        subplot8 = plot(ax8,r1_sameAngle,gamma21_sameAngle,r2_sameAngle,gamma22_sameAngle, ...
            r3_sameAngle,gamma23_sameAngle,r4_sameAngle,gamma24_sameAngle)

        subplot9 = plot(ax9,r1_sameAngle,gamma31_sameAngle,r2_sameAngle,gamma32_sameAngle, ...
            r3_sameAngle,gamma33_sameAngle,r4_sameAngle,gamma34_sameAngle)

        set(subplot1(1),'Color',blue);
        set(subplot1(2),'Color',green);
        set(subplot1(3),'Color',yellow);
        set(subplot1(4),'Color',red);

        set(subplot2(1),'Color',blue);
        set(subplot2(2),'Color',green);
        set(subplot2(3),'Color',yellow);
        set(subplot2(4),'Color',red);

        set(subplot3(1),'Color',blue);
        set(subplot3(2),'Color',green);
        set(subplot3(3),'Color',yellow);
        set(subplot3(4),'Color',red);

        set(subplot4(1),'Color',blue);
        set(subplot4(2),'Color',green);
        set(subplot4(3),'Color',yellow);
        set(subplot4(4),'Color',red);

        set(subplot5(1),'Color',blue);
        set(subplot5(2),'Color',green);
        set(subplot5(3),'Color',yellow);
        set(subplot5(4),'Color',red);

        set(subplot6(1),'Color',blue);
        set(subplot6(2),'Color',green);
        set(subplot6(3),'Color',yellow);
        set(subplot6(4),'Color',red);

        set(subplot7(1),'Color',blue);
        set(subplot7(2),'Color',green);
        set(subplot7(3),'Color',yellow);
        set(subplot7(4),'Color',red);

        set(subplot8(1),'Color',blue);
        set(subplot8(2),'Color',green);
        set(subplot8(3),'Color',yellow);
        set(subplot8(4),'Color',red);

        set(subplot9(1),'Color',blue);
        set(subplot9(2),'Color',green);
        set(subplot9(3),'Color',yellow);
        set(subplot9(4),'Color',red);
    end
    xlabel(ax1,'r')
    ylabel(ax1,'h')

    xlabel(ax2,'r')
    ylabel(ax2,'vr')

    xlabel(ax3,'r')
    ylabel(ax3,'\alpha_1')

    xlabel(ax4,'r')
    ylabel(ax4,'\alpha_2')

    xlabel(ax5,'r')
    ylabel(ax5,'\alpha_3')

    xlabel(ax6,'r')
    ylabel(ax6,'v\theta')

    xlabel(ax7,'r')
    ylabel(ax7,'\gamma_1')

    xlabel(ax8,'r')
    ylabel(ax8,'\gamma_2')

    xlabel(ax9,'r')
    ylabel(ax9,'\gamma_3')

end
if(strcmp(plotting,'h'))
    if(strcmp(plotting_mode,'time_evolution'))
        plot1 = plot(r1_sameAngle,h1_sameAngle,r2_sameAngle,h2_sameAngle, ...
            r1_sameAngle,h3_sameAngle,r2_sameAngle,h4_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',green);
        set(plot1(3),'Color',yellow);
        set(plot1(4),'Color',red);
    
        leg = legend('t=0.15','t=0.2', 't=0.25', 't=0.3');
        set(leg,'FontSize',12); 
    end
    if(strcmp(plotting_mode,'model_comparison'))
        plot1 = plot(r1_sameAngle,h1_sameAngle,r2_sameAngle,h2_sameAngle, r_ref,h_ref,'*')
        %plot1 = plot(r1_sameAngle,h1_sameAngle,r2_sameAngle,h2_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(1),'LineWidth',2);

        set(plot1(2),'Color',red);
        set(plot1(2),'LineWidth',2);

        set(plot1(3),'Color',brown);
        set(plot1(3),'MarkerSize',4);
    
        leg = legend('ASWME','HASWME','Reference');
        %leg = legend('t=0','t=t_{end}');
        set(leg,'FontSize',14); 
    end
    if(strcmp(plotting_mode,'order_comparison'))
        plot1 = plot(r1_sameAngle,h1_sameAngle,r2_sameAngle,h2_sameAngle,'-.', ...
            r1_sameAngle,h3_sameAngle,'--',r2_sameAngle,h4_sameAngle,':',r_ref,h_ref,'*')

        set(plot1(1),'Color',blue);
        set(plot1(1),'LineWidth',2);

        set(plot1(2),'Color',green);
        set(plot1(2),'LineWidth',2);

        set(plot1(3),'Color',yellow);
        set(plot1(3),'LineWidth',2);

        set(plot1(4),'Color',red);
        set(plot1(4),'LineWidth',2);

        set(plot1(5),'Color',brown);
        %set(plot1(5),'LineWidth',1);
        set(plot1(5),'MarkerSize',5);
     
        leg = legend('N=0','N=1', 'N=2', 'N=3','Reference','Location','Northwest');
        set(leg,'FontSize',14); 
    end
    axis([3,5,0.8,5.2]); % Change if necessary
    yticks([1.0 2.0 3.0 4.0 5.0])
    xlabel('$r$','Interpreter','latex','FontSize',18)
    ylabel('$h$','Interpreter','latex','FontSize',18) 
end
if(strcmp(plotting,'vr'))
    if(strcmp(plotting_mode,'time_evolution'))
        plot1 = plot(r1_sameAngle,vr1_sameAngle,r2_sameAngle,vr2_sameAngle, ...
            r1_sameAngle,vr3_sameAngle,r2_sameAngle,vr4_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',green);
        set(plot1(3),'Color',yellow);
        set(plot1(4),'Color',red);
    
        leg = legend('t=0.025','t=0.1', 't=0.25', 't=0.275');
        set(leg,'FontSize',12); 
    end
    if(strcmp(plotting_mode,'model_comparison'))
        plot1 = plot(r1_sameAngle,vr1_sameAngle,r2_sameAngle,vr2_sameAngle, r_ref,vr_ref,'*')

        set(plot1(1),'Color',blue);
        set(plot1(1),'LineWidth',2);

        set(plot1(2),'Color',red);
        set(plot1(2),'LineWidth',2);

        set(plot1(3),'Color',brown);
        %set(plot1(5),'LineWidth',1);
        set(plot1(3),'MarkerSize',4);

        leg = legend('ASWME','HASWME','Reference','Location','Northwest');
        set(leg,'FontSize',14); 
    end
    if(strcmp(plotting_mode,'order_comparison'))
        plot1 = plot(r1_sameAngle,vr1_sameAngle,r2_sameAngle,vr2_sameAngle,'-.', ...
            r1_sameAngle,vr3_sameAngle,'--',r2_sameAngle,vr4_sameAngle,':',r_ref,vr_ref,'*')

        set(plot1(1),'Color',blue);
        set(plot1(1),'LineWidth',2);

        set(plot1(2),'Color',green);
        set(plot1(2),'LineWidth',2);

        set(plot1(3),'Color',yellow);
        set(plot1(3),'LineWidth',2);

        set(plot1(4),'Color',red);
        set(plot1(4),'LineWidth',2);

        set(plot1(5),'Color',brown);
        %set(plot1(5),'LineWidth',1);
        set(plot1(5),'MarkerSize',5);

        leg = legend('N=0','N=1', 'N=2', 'N=3','Reference','Location','southeast');
        set(leg,'FontSize',12); 
    end
    axis([3,5,0,1.8]); % Change if necessary
    yticks([0 0.5 1.0 1.5])
    xlabel('$r$','Interpreter','latex','FontSize',18)
    ylabel('$v_{r,m}$','Interpreter','latex','FontSize',18) 
end
if(strcmp(plotting,'alpha1'))
    if(strcmp(plotting_mode,'time_evolution'))
        plot1 = plot(r1_sameAngle,alpha11_sameAngle,r2_sameAngle,alpha12_sameAngle, ...
            r1_sameAngle,alpha13_sameAngle,r2_sameAngle,alpha14_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',green);
        set(plot1(3),'Color',yellow);
        set(plot1(4),'Color',red);
    
        leg = legend('t=0.025','t=0.1', 't=0.25', 't=0.275');
        set(leg,'FontSize',12); 
    end
    if(strcmp(plotting_mode,'model_comparison'))
        plot1 = plot(r1_sameAngle,alpha11_sameAngle,'-.',r2_sameAngle,alpha12_sameAngle,'--' ...
            ,r_ref,alpha1_ref,'*')

        set(plot1(1),'Color',blue);
        set(plot1(1),'LineWidth',2);

        set(plot1(2),'Color',red);
        set(plot1(2),'LineWidth',2);

        set(plot1(3),'Color',brown);
%         set(plot1(5),'LineWidth',1);
        set(plot1(3),'MarkerSize',4);
    
        %plot2 = plot(r_ref,alpha1_ref,'.');
        %set(plot2,'Color',green)

        leg = legend('ASWME','HASWME','reference');
        set(leg,'FontSize',14); 
    end
    axis([3.5,5,-0.8,0.05]); % Change if necessary
    %xticks([0 5 10])
    xlabel('$r$','Interpreter','latex','FontSize',18)
    ylabel('$\alpha_1$','Interpreter','latex','FontSize',18) 
end
if(strcmp(plotting,'alpha2'))
    if(strcmp(plotting_mode,'time_evolution'))
        plot1 = plot(r1_sameAngle,alpha21_sameAngle,r2_sameAngle,alpha22_sameAngle, ...
            r1_sameAngle,alpha23_sameAngle,r2_sameAngle,alpha24_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',green);
        set(plot1(3),'Color',yellow);
        set(plot1(4),'Color',red);
    
        leg = legend('t=0.15','t=0.2', 't=0.25', 't=0.3');
        set(leg,'FontSize',12); 
    end
    if(strcmp(plotting_mode,'model_comparison'))
        plot1 = plot(r1_sameAngle,alpha21_sameAngle,r2_sameAngle,alpha22_sameAngle)

        hold on;

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',red);
    
        %plot2 = plot(r_ref,alpha2_ref,'.');
        s%et(plot2,'Color',green)

        leg = legend('ASWME','HASWME');
        set(leg,'FontSize',12); 
    end
    axis([1,11,-0.6,0.2]); % Change if necessary
    xlabel('$r$','Interpreter','latex','FontSize',18)
    ylabel('$\alpha_2$','Interpreter','latex','FontSize',18) 
end
if(strcmp(plotting,'alpha3'))
    if(strcmp(plotting_mode,'time_evolution'))
        plot1 = plot(r1_sameAngle,alpha31_sameAngle,r2_sameAngle,alpha32_sameAngle, ...
            r1_sameAngle,alpha33_sameAngle,r2_sameAngle,alpha34_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',green);
        set(plot1(3),'Color',yellow);
        set(plot1(4),'Color',red);
    
        leg = legend('t=0.025','t=0.1', 't=0.25', 't=0.275');
        set(leg,'FontSize',12); 
    end
    if(strcmp(plotting_mode,'model_comparison'))
        plot1 = plot(r1_sameAngle,alpha31_sameAngle,r2_sameAngle,alpha32_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',red);
    
        leg = legend('ASWME','HASWME');
        set(leg,'FontSize',12); 
    end
    axis([1,11,-0.6,0.6]); % Change if necessary
    xlabel('$r$','Interpreter','latex','FontSize',18)
    ylabel('$\alpha_3$','Interpreter','latex','FontSize',18)  
end
if(strcmp(plotting,'vtheta'))
    if(strcmp(plotting_mode,'time_evolution'))
        plot1 = plot(r1_sameAngle,vtheta1_sameAngle,r2_sameAngle,vtheta2_sameAngle, ...
            r1_sameAngle,vtheta3_sameAngle,r2_sameAngle,vtheta4_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',green);
        set(plot1(3),'Color',yellow);
        set(plot1(4),'Color',red);

        leg = legend('t=0.025','t=0.1', 't=0.25', 't=0.275');
        set(leg,'FontSize',12); 
    end
    if(strcmp(plotting_mode,'model_comparison'))
        plot1 = plot(r1_sameAngle,vtheta1_sameAngle,r2_sameAngle,vtheta2_sameAngle, ...
            r_ref,vtheta_ref,'*')

        set(plot1(1),'Color',blue);
        set(plot1(1),'LineWidth',1.2);

        set(plot1(2),'Color',red);
        set(plot1(2),'LineWidth',1.2);

        set(plot1(3),'Color',brown);
        %set(plot1(5),'LineWidth',1);
        set(plot1(3),'MarkerSize',4);

        leg = legend('ASWME','HASWME','Reference','Location','northwest');
        set(leg,'FontSize',12); 
    end
    if(strcmp(plotting_mode,'order_comparison'))
        plot1 = plot(r1_sameAngle,vtheta1_sameAngle,r2_sameAngle,vtheta2_sameAngle,'-.', ...
            r1_sameAngle,vtheta3_sameAngle,'--',r2_sameAngle,vtheta4_sameAngle,':',r_ref,vtheta_ref,'*')

        set(plot1(1),'Color',blue);
        set(plot1(1),'LineWidth',2);

        set(plot1(2),'Color',green);
        set(plot1(2),'LineWidth',2);

        set(plot1(3),'Color',yellow);
        set(plot1(3),'LineWidth',2);

        set(plot1(4),'Color',red);
        set(plot1(4),'LineWidth',2);

        set(plot1(5),'Color',brown);
        %set(plot1(5),'LineWidth',1);
        set(plot1(5),'MarkerSize',5);

        leg = legend('N=0','N=1', 'N=2', 'N=3', 'Reference','Location','northwest');
        set(leg,'FontSize',14);
        %set(leg,'Location','northwest');
    end
    axis([1.05,8,-0.1,0.55]); % Change if necessary
    %yline(0);
    %yticklabels({'','0.02', '','0.03','','0.04','','0.05','','0.06',''})
    xlabel('$r$','Interpreter','latex','FontSize',18)
    ylabel('$v_{\theta,m}$','Interpreter','latex','FontSize',18)  
end
if(strcmp(plotting,'gamma1'))
    if(strcmp(plotting_mode,'time_evolution'))
        plot1 = plot(r1_sameAngle,gamma11_sameAngle,r2_sameAngle,gamma12_sameAngle, ...
            r1_sameAngle,gamma13_sameAngle,r2_sameAngle,gamma14_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',green);
        set(plot1(3),'Color',yellow);
        set(plot1(4),'Color',red);
    
        leg = legend('t=0.025','t=0.1', 't=0.25', 't=0.275');
        set(leg,'FontSize',12); 
    end
    if(strcmp(plotting_mode,'model_comparison'))
        plot1 = plot(r1_sameAngle,gamma11_sameAngle,r2_sameAngle,gamma12_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',red);
    
        leg = legend('ASWME','HASWME');
        set(leg,'FontSize',12); 
    end
    axis([1,11,-0.5,0.5]); % Change if necessary
    xlabel('$r$','Interpreter','latex','FontSize',18)
    ylabel('$\gamma_1$','Interpreter','latex','FontSize',18) 
end
if(strcmp(plotting,'gamma2'))
    if(strcmp(plotting_mode,'time_evolution'))
        plot1 = plot(r1_sameAngle,gamma21_sameAngle,r2_sameAngle,gamma22_sameAngle, ...
            r1_sameAngle,gamma23_sameAngle,r2_sameAngle,gamma24_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',green);
        set(plot1(3),'Color',yellow);
        set(plot1(4),'Color',red);
    
        leg = legend('t=0.025','t=0.1', 't=0.25', 't=0.275');
        set(leg,'FontSize',12); 
    end
    if(strcmp(plotting_mode,'model_comparison'))
        plot1 = plot(r1_sameAngle,gamma21_sameAngle,r2_sameAngle,gamma22_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',red);
    
        leg = legend('ASWME','HASWME');
        set(leg,'FontSize',12); 
    end
    axis([1,11,-0.5,0.5]); % Change if necessary
    xlabel('$r$','Interpreter','latex','FontSize',18)
    ylabel('$\gamma_2$','Interpreter','latex','FontSize',18) 
end
if(strcmp(plotting,'gamma3'))
    if(strcmp(plotting_mode,'time_evolution'))
        plot1 = plot(r1_sameAngle,gamma31_sameAngle,r2_sameAngle,gamma32_sameAngle, ...
            r1_sameAngle,gamma33_sameAngle,r2_sameAngle,gamma34_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',green);
        set(plot1(3),'Color',yellow);
        set(plot1(4),'Color',red);
    
        leg = legend('t=0.025','t=0.1', 't=0.25', 't=0.275');
        set(leg,'FontSize',12); 
    end
    if(strcmp(plotting_mode,'model_comparison'))
        plot1 = plot(r1_sameAngle,gamma31_sameAngle,r2_sameAngle,gamma32_sameAngle)

        set(plot1(1),'Color',blue);
        set(plot1(2),'Color',red);
    
        leg = legend('ASWME','HASWME');
        set(leg,'FontSize',12); 
    end
    axis([1,11,-0.5,0.5]); % Change if necessary
    xlabel('$r$','Interpreter','latex','FontSize',18)
    ylabel('$\gamma_3$','Interpreter','latex','FontSize',18) 
end


%% Export figure
 %addpath('C:\Users\rikve\Prime\results\ASWME\export_fig\', '-end')
 %cd export_fig
 %export_fig('ASWMEPaperRevised_RadialDamBreakUnstable_h_Time0.3_SubDomain_Publication.pdf', '-pdf','-transparent');
