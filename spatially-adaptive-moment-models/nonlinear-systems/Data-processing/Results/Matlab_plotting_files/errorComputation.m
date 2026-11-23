%% Clear statements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Computation of the model errors %%%%%%%%%%%%%%%%
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

%name1 = 'SmoothExpit_ASWMEPaper_Time0.5_FirstOrder'; type1 = 'swme';
%name2 = 'SmoothExpit_ASWMEPaper_Time0.5_ThirdOrder'; type2 = 'swme';
%name3 = 'SmoothExpit_ASWMEPaper_Time0_ASWE'; type3 = 'swe';
%name4 = 'SmoothExpit_ASWMEPaper_Time0_ASWE'; type4 = 'swe';

%name1 = 'SmoothExpit_ASWMEPaper_Time0.5_ASWE'; type1 = 'swe';
%name2 = 'SmoothExpit_ASWMEPaper_Time0.5_FirstOrder_Auto'; type2 = 'swme';
%name3 = 'SmoothExpit_ASWMEPaper_Time0.5_SecondOrder_Auto'; type3 = 'swme';
%name4 = 'SmoothExpit_ASWMEPaper_Time0.5_ThirdOrder_Auto'; type4 = 'swme';
%name5 = 'SmoothExpit_ASWMEPaper_Time0.5_FourthOrder'; type5 = 'swme';
%name6 = 'SmoothExpit_ASWMEPaper_Time0.5_FifthOrder'; type6 = 'swme';

name1 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_ZeroethOrder_Time1.0_4000'; type1 = 'swe';
name2 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_FirstOrder_Time1.0_4000'; type2 = 'swme';
name3 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_SecondOrder_SWME_Time1.0_4000'; type3 = 'swme';
name4 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_ThirdOrder_SWME_Time1.0_4000'; type4 = 'swme';
name5 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_FourthOrder_SWME_Time1.0_4000'; type5 = 'swme';
%name6 = 'ASWMEPaperRevised_SmoothExpit18_InitialAngularVelocity_FourthOrder_SWME_Time1.0'; type6 = 'swme';

%name1 = 'SmoothExpit_ASWMEPaper_Time0.5_ASWE_5000x4'; type1 = 'swe';
%name2 = 'SmoothExpit_ASWMEPaper_Time0.5_FirstOrder_5000x4'; type2 = 'swme';
%name3 = 'SmoothExpit_ASWMEPaper_Time0.5_SecondOrder_5000x4'; type3 = 'swme';
%name4 = 'SmoothExpit_ASWMEPaper_Time0.5_ThirdOrder_5000x4'; type4 = 'swme';
%name5 = 'SmoothExpit_ASWMEPaper_Time0.5_FourthOrder_5000x4'; type5 = 'swme';

%name1 = 'RadialDamBreak_NoVelocity_Time0.5_ThirdOrder'; type1 = 'swme';
%name2 = 'RadialDamBreak_NoVelocity_Time0.5_ThirdOrder_HSWME'; type2 = 'swme';
%name3 = 'RadialDamBreak_NoVelocity_Time0_ThirdOrder'; type3 = 'swme';
%name4 = 'RadialDamBreak_NoVelocity_Time0_ThirdOrder'; type4 = 'swme';

%name1 = 'RadialDamBreak_ASWMEPaper_Time0.15_ThirdOrder'; type1 = 'swme';
%name2 = 'RadialDamBreak_ASWMEPaper_Time0.2_ThirdOrder'; type2 = 'swme';
%name3 = 'RadialDamBreak_ASWMEPaper_Time0.25_ThirdOrder'; type3 = 'swme';
%name4 = 'RadialDamBreak_ASWMEPaper_Time0.3_ThirdOrder'; type4 = 'swme';

%name1 = 'RadialDamBreakUnstable_ASWMEPaper_Time0.15_ThirdOrder'; type1 = 'swme';
%name2 = 'RadialDamBreakUnstable_ASWMEPaper_Time0.15_ThirdOrder_HSWME'; type2 = 'swme';
%name3 = 'RadialDamBreakUnstable_ASWMEPaper_Time0_ThirdOrder'; type3 = 'swme';
%name4 = 'RadialDamBreakUnstable_ASWMEPaper_Time0_ThirdOrder_HSWME'; type4 = 'swme';

% Moment model data

[r1,theta1,h1,vr1,alpha11,alpha21,alpha31,vtheta1,gamma11,gamma21,gamma31] = readDataASWME(name1,type1);
[r2,theta2,h2,vr2,alpha12,alpha22,alpha32,vtheta2,gamma12,gamma22,gamma32] = readDataASWME(name2,type2);
[r3,theta3,h3,vr3,alpha13,alpha23,alpha33,vtheta3,gamma13,gamma23,gamma33] = readDataASWME(name3,type3);
[r4,theta4,h4,vr4,alpha14,alpha24,alpha34,vtheta4,gamma14,gamma24,gamma34] = readDataASWME(name4,type4);
[r5,theta5,h5,vr5,alpha15,alpha25,alpha35,vtheta5,gamma15,gamma25,gamma35] = readDataASWME(name5,type5);
%[r6,theta6,h6,vr6,alpha16,alpha26,alpha36,vtheta6,gamma16,gamma26,gamma36] = readDataASWME(name6,type6);

% Reference data

%reference_data = load("RadialDamBreakUnstable_ASWMEPaper_Time0.15_Reference.csv");
%reference_data = load("ASWMEPaperRevised_SmoothExpit[1,8]_Time1.0_Reference_200x100_CFL0.5.csv");
reference_data = load("ASWMEPaperRevised_SmoothExpit[1,8]_Time1.0_Reference_500x200_CFL0.5.csv");


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
% grid = '2000x4Axisymmetric';
 grid = '4000x4Axisymmetric';

if(strcmp(grid,'100x4Axisymmetric'))
    sameAngle = 100;
    numberOfCells_MomentModel = 100;
end

if(strcmp(grid,'100x4Axisymmetric2'))
    sameAngle = 100;
    numberOfCells_MomentModel = 100;
end

if(strcmp(grid,'100x10Axisymmetric'))
    sameAngle = 100;
    numberOfCells_MomentModel = 100;
end

if(strcmp(grid,'100x10Axisymmetric2'))
    sameAngle = 100;
    numberOfCells_MomentModel = 100;
end

if(strcmp(grid,'10000x4Axisymmetric'))
    sameAngle = 10000;
    numberOfCells_MomentModel = 10000;
end

if(strcmp(grid,'5000x4Axisymmetric'))
    sameAngle = 5000;
    numberOfCells_MomentModel = 5000;
end

if(strcmp(grid,'2000x4Axisymmetric'))
    sameAngle = 2000;
    numberOfCells_MomentModel = 2000;
end

if(strcmp(grid,'4000x4Axisymmetric'))
    sameAngle = 4000;
    numberOfCells_MomentModel = 4000;
end

% Change this if necessary
numberOfCells_Reference = 500;
r_1 = 1.0;
r_2 = 8.0;
delta_MomentModel = (r_2-r_1)/numberOfCells_MomentModel;
delta_Reference = (r_2-r_1)/numberOfCells_Reference;

delta_ratio = delta_Reference/delta_MomentModel;

r1 = r1(1:sameAngle);
r2 = r2(1:sameAngle);
r3 = r3(1:sameAngle);
r4 = r4(1:sameAngle);
r5 = r5(1:sameAngle);
%r6 = r6(1:sameAngle);

theta1 = theta1(1:sameAngle);
theta2 = theta2(1:sameAngle);
theta3 = theta3(1:sameAngle);
theta4 = theta4(1:sameAngle);
theta5 = theta5(1:sameAngle);
%theta6 = theta6(1:sameAngle);

h1 = h1(1:sameAngle);
h2 = h2(1:sameAngle);
h3 = h3(1:sameAngle);
h4 = h4(1:sameAngle);
h5 = h5(1:sameAngle);
%h6 = h6(1:sameAngle);

vr1 = vr1(1:sameAngle);
vr2 = vr2(1:sameAngle);
vr3 = vr3(1:sameAngle);
vr4 = vr4(1:sameAngle);
vr5 = vr5(1:sameAngle);
%vr6 = vr6(1:sameAngle);

%% Error convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Error convergence analysis %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_interpol1 = zeros(length(r_ref),1);
h_interpol2 = zeros(length(r_ref),1);
h_interpol3 = zeros(length(r_ref),1);
h_interpol4 = zeros(length(r_ref),1);
h_interpol5 = zeros(length(r_ref),1);
%h_interpol6 = zeros(length(r_ref),1);

vr_interpol1 = zeros(length(r_ref),1);
vr_interpol2 = zeros(length(r_ref),1);
vr_interpol3 = zeros(length(r_ref),1);
vr_interpol4 = zeros(length(r_ref),1);
vr_interpol5 = zeros(length(r_ref),1);
%vr_interpol6 = zeros(length(r_ref),1);

vtheta_interpol1 = zeros(length(r_ref),1);
vtheta_interpol2 = zeros(length(r_ref),1);
vtheta_interpol3 = zeros(length(r_ref),1);
vtheta_interpol4 = zeros(length(r_ref),1);
vtheta_interpol5 = zeros(length(r_ref),1);
%vtheta_interpol6 = zeros(length(r_ref),1);

start=4;

parity=(-1)^delta_ratio;

if(parity==1)
    for i = start:length(r_ref)
        h_interpol1(i) = (h1((i-1/2)*delta_ratio)+h1((i-1/2)*delta_ratio+1))/2;
        h_interpol2(i) = (h2((i-1/2)*delta_ratio)+h2((i-1/2)*delta_ratio+1))/2;
        h_interpol3(i) = (h3((i-1/2)*delta_ratio)+h3((i-1/2)*delta_ratio+1))/2;
        h_interpol4(i) = (h4((i-1/2)*delta_ratio)+h4((i-1/2)*delta_ratio+1))/2;
        h_interpol5(i) = (h5((i-1/2)*delta_ratio)+h5((i-1/2)*delta_ratio+1))/2;
        %h_interpol6(i) = (h6((i-1/2)*delta_ratio)+h6((i-1/2)*delta_ratio+1))/2;
    
        vr_interpol1(i) = (vr1((i-1/2)*delta_ratio)+vr1((i-1/2)*delta_ratio+1))/2;
        vr_interpol2(i) = (vr2((i-1/2)*delta_ratio)+vr2((i-1/2)*delta_ratio+1))/2;
        vr_interpol3(i) = (vr3((i-1/2)*delta_ratio)+vr3((i-1/2)*delta_ratio+1))/2;
        vr_interpol4(i) = (vr4((i-1/2)*delta_ratio)+vr4((i-1/2)*delta_ratio+1))/2;
        vr_interpol5(i) = (vr5((i-1/2)*delta_ratio)+vr5((i-1/2)*delta_ratio+1))/2;
        %vr_interpol6(i) = (vr6((i-1/2)*delta_ratio)+vr6((i-1/2)*delta_ratio+1))/2;
    
        vtheta_interpol1(i) = (vtheta1((i-1/2)*delta_ratio)+vtheta1((i-1/2)*delta_ratio+1))/2;
        vtheta_interpol2(i) = (vtheta2((i-1/2)*delta_ratio)+vtheta2((i-1/2)*delta_ratio+1))/2;
        vtheta_interpol3(i) = (vtheta3((i-1/2)*delta_ratio)+vtheta3((i-1/2)*delta_ratio+1))/2;
        vtheta_interpol4(i) = (vtheta4((i-1/2)*delta_ratio)+vtheta4((i-1/2)*delta_ratio+1))/2;
        vtheta_interpol5(i) = (vtheta5((i-1/2)*delta_ratio)+vtheta5((i-1/2)*delta_ratio+1))/2;
        %vtheta_interpol6(i) = (vtheta6((i-1/2)*delta_ratio)+vtheta6((i-1/2)*delta_ratio+1))/2;
    end
end

if(parity==-1)
    for i = start:length(r_ref)
        h_interpol1(i) = h1((i-1/2)*delta_ratio+1/2);
        h_interpol2(i) = h2((i-1/2)*delta_ratio+1/2);
        h_interpol3(i) = h3((i-1/2)*delta_ratio+1/2);
        h_interpol4(i) = h4((i-1/2)*delta_ratio+1/2);
        h_interpol5(i) = h5((i-1/2)*delta_ratio+1/2);
        %h_interpol6(i) = (h6((i-1/2)*delta_ratio)+h6((i-1/2)*delta_ratio+1))/2;
    
        vr_interpol1(i) = vr1((i-1/2)*delta_ratio+1/2);
        vr_interpol2(i) = vr2((i-1/2)*delta_ratio+1/2);
        vr_interpol3(i) = vr3((i-1/2)*delta_ratio+1/2);
        vr_interpol4(i) = vr4((i-1/2)*delta_ratio+1/2);
        vr_interpol5(i) = vr5((i-1/2)*delta_ratio+1/2);
        %vr_interpol6(i) = (vr6((i-1/2)*delta_ratio)+vr6((i-1/2)*delta_ratio+1))/2;
    
        vtheta_interpol1(i) = vtheta1((i-1/2)*delta_ratio+1/2);
        vtheta_interpol2(i) = vtheta2((i-1/2)*delta_ratio+1/2);
        vtheta_interpol3(i) = vtheta3((i-1/2)*delta_ratio+1/2);
        vtheta_interpol4(i) = vtheta4((i-1/2)*delta_ratio+1/2);
        vtheta_interpol5(i) = vtheta5((i-1/2)*delta_ratio+1/2);
        %vtheta_interpol6(i) = (vtheta6((i-1/2)*delta_ratio)+vtheta6((i-1/2)*delta_ratio+1))/2;
    end
end

h_diff_norm1 = norm(h_interpol1(start:end) - h_ref(start:end))/norm(h_ref(start:end));
h_diff_norm2 = norm(h_interpol2(start:end) - h_ref(start:end))/norm(h_ref(start:end));
h_diff_norm3 = norm(h_interpol3(start:end) - h_ref(start:end))/norm(h_ref(start:end));
h_diff_norm4 = norm(h_interpol4(start:end) - h_ref(start:end))/norm(h_ref(start:end));
h_diff_norm5 = norm(h_interpol5(start:end) - h_ref(start:end))/norm(h_ref(start:end));
%h_diff_norm6 = norm(h_interpol6(2:end) - h_ref(2:end))/norm(h_ref(2:end));

vr_diff_norm1 = norm(vr_interpol1(start:end) - vr_ref(start:end))/norm(vr_ref(start:end));
vr_diff_norm2 = norm(vr_interpol2(start:end) - vr_ref(start:end))/norm(vr_ref(start:end));
vr_diff_norm3 = norm(vr_interpol3(start:end) - vr_ref(start:end))/norm(vr_ref(start:end));
vr_diff_norm4 = norm(vr_interpol4(start:end) - vr_ref(start:end))/norm(vr_ref(start:end));
vr_diff_norm5 = norm(vr_interpol5(start:end) - vr_ref(start:end))/norm(vr_ref(start:end));
%vr_diff_norm6 = norm(vr_interpol6(2:end) - vr_ref(2:end))/norm(vr_ref(2:end));

vtheta_diff_norm1 = norm(vtheta_interpol1(start:end) - vtheta_ref(start:end))/norm(vtheta_ref(start:end));
vtheta_diff_norm2 = norm(vtheta_interpol2(start:end) - vtheta_ref(start:end))/norm(vtheta_ref(start:end));
vtheta_diff_norm3 = norm(vtheta_interpol3(start:end) - vtheta_ref(start:end))/norm(vtheta_ref(start:end));
vtheta_diff_norm4 = norm(vtheta_interpol4(start:end) - vtheta_ref(start:end))/norm(vtheta_ref(start:end));
vtheta_diff_norm5 = norm(vtheta_interpol5(start:end) - vtheta_ref(start:end))/norm(vtheta_ref(start:end));
%vtheta_diff_norm6 = norm(vtheta_interpol6(2:end) - vtheta_ref(2:end))/norm(vtheta_ref(2:end));

h_diff_norms = [h_diff_norm1,h_diff_norm2,h_diff_norm3,h_diff_norm4,h_diff_norm5];
vr_diff_norms = [vr_diff_norm1,vr_diff_norm2,vr_diff_norm3,vr_diff_norm4,vr_diff_norm5];
vtheta_diff_norms = [vtheta_diff_norm1,vtheta_diff_norm2,vtheta_diff_norm3,vtheta_diff_norm4,vtheta_diff_norm5];

%writematrix(h_diff_norms,'ASWMEPaperRevised_SmoothExpit_h_Error_ASWME_BoundaryExcluded_Order4_4000x4_500x200.csv')
%writematrix(vr_diff_norms,'ASWMEPaperRevised_SmoothExpit_vr_Error_ASWME_BoundaryExcluded_Order4_4000x4_500x200.csv')
%writematrix(vtheta_diff_norms,'ASWMEPaperRevised_SmoothExpit_vtheta_Error_ASWME_BoundaryExcluded_Order4_4000x4_500x200.csv')

x_axis = [0,1,2,3,4];

plot(x_axis,vtheta_diff_norms)