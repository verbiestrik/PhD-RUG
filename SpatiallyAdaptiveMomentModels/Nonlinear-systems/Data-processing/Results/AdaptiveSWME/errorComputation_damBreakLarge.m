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

% name0 = 'damBreak_classical_order0_time0p4'; type0 = 'swe';
% name1 = 'damBreak_classical_order1_time0p4'; type1 = 'swme';
% name2 = 'damBreak_classical_order2_time0p4'; type2 = 'swme';
% name3 = 'damBreak_classical_order3_time0p4'; type3 = 'swme';
% name4 = 'damBreak_classical_order4_time0p4'; type4 = 'swme';
% name5 = 'damBreak_classical_order5_time0p4'; type5 = 'swme';
% namePBC = 'damBreak_adaptive_couplingPBC_time0p4'; typePBC = 'swmeAdaptive';
% nameCIF = 'damBreak_adaptive_couplingCIF_time0p4'; typeCIF = 'swmeAdaptive';
% reference_data = load("damBreakLarge_lambda0p1_nu1p0_t0p4_250x100.csv");

name0 = 'damBreak_classical_order0_time1p0'; type0 = 'swe';
name1 = 'damBreak_classical_order1_time1p0'; type1 = 'swme';
name2 = 'damBreak_classical_order2_time1p0'; type2 = 'swme';
name3 = 'damBreak_classical_order3_time1p0'; type3 = 'swme';
name4 = 'damBreak_classical_order4_time1p0'; type4 = 'swme';
name5 = 'damBreak_classical_order5_time1p0'; type5 = 'swme';
namePBC = 'damBreak_adaptive_couplingPBC_time1p0'; typePBC = 'swmeAdaptive';
nameCIF = 'damBreak_adaptive_couplingCIF_time1p0'; typeCIF = 'swmeAdaptive';
reference_data = load("damBreakLarge_lambda0p1_nu1p0_t1p0_800x200_new.csv");

time0 = load(strcat("time",name0,".csv"));
time1 = load(strcat("time",name1,".csv"));
time2 = load(strcat("time",name2,".csv"));
time3 = load(strcat("time",name3,".csv"));
time4 = load(strcat("time",name4,".csv"));
time5 = load(strcat("time",name5,".csv"));
timePBC = load(strcat("time",namePBC,".csv"));
timeCIF = load(strcat("time",nameCIF,".csv"));

times_rel = [time0,time1,time2,time3,time4,time5,timePBC,timeCIF]/time5

% Moment model data

[x0,h0,u0,alpha10,alpha20,alpha30,alpha40,alpha50,moments0] = readDataAdaptiveSWME1D(name0,type0);
[x1,h1,u1,alpha11,alpha21,alpha31,alpha41,alpha51,moments1] = readDataAdaptiveSWME1D(name1,type1);
[x2,h2,u2,alpha12,alpha22,alpha32,alpha42,alpha52,moments2] = readDataAdaptiveSWME1D(name2,type2);
[x3,h3,u3,alpha13,alpha23,alpha33,alpha43,alpha53,moments3] = readDataAdaptiveSWME1D(name3,type3);
[x4,h4,u4,alpha14,alpha24,alpha34,alpha44,alpha54,moments4] = readDataAdaptiveSWME1D(name4,type4);
[x5,h5,u5,alpha15,alpha25,alpha35,alpha45,alpha55,moments5] = readDataAdaptiveSWME1D(name5,type5);
[xPBC,hPBC,uPBC,alpha1PBC,alpha2PBC,alpha3PBC,alpha4PBC,alpha5PBC,momentsPBC] = readDataAdaptiveSWME1D(namePBC,typePBC);
[xCIF,hCIF,uCIF,alpha1CIF,alpha2CIF,alpha3CIF,alpha4CIF,alpha5CIF,momentsCIF] = readDataAdaptiveSWME1D(nameCIF,typeCIF);

% Reference data

x_ref = reference_data(:,1); h_ref = reference_data(:,2); u_ref = reference_data(:,3);
alpha1_ref = 3.*reference_data(:,4); alpha2_ref = 5.*reference_data(:,5);
alpha3_ref = 7.*reference_data(:,6); alpha4_ref = 9.*reference_data(:,7);
alpha5_ref = 11.*reference_data(:,8);


%% Error convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Error convergence analysis %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a common grid over the interval
xmin = max(min(x1), min(x_ref));
xmax = min(max(x1), max(x_ref));

N = 10000;                         % Number of grid points
x_common = linspace(xmin, xmax, N);

% Interpolate both functions onto the common grid
h0_interp = interp1(x0, h0, x_common, 'linear');
h1_interp = interp1(x1, h1, x_common, 'linear');
h2_interp = interp1(x1, h2, x_common, 'linear');
h3_interp = interp1(x1, h3, x_common, 'linear');
h4_interp = interp1(x1, h4, x_common, 'linear');
h5_interp = interp1(x1, h5, x_common, 'linear');
hPBC_interp = interp1(x1, hPBC, x_common, 'linear');
hCIF_interp = interp1(x1, hCIF, x_common, 'linear');
h_ref_interp = interp1(x_ref, h_ref, x_common, 'linear');

u0_interp = interp1(x1, u0, x_common, 'linear');
u1_interp = interp1(x1, u1, x_common, 'linear');
u2_interp = interp1(x1, u2, x_common, 'linear');
u3_interp = interp1(x1, u3, x_common, 'linear');
u4_interp = interp1(x1, u4, x_common, 'linear');
u5_interp = interp1(x1, u5, x_common, 'linear');
uPBC_interp = interp1(x1, uPBC, x_common, 'linear');
uCIF_interp = interp1(x1, uCIF, x_common, 'linear');
u_ref_interp = interp1(x_ref, u_ref, x_common, 'linear');

alpha10_interp = interp1(x1, alpha10, x_common, 'linear');
alpha11_interp = interp1(x1, alpha11, x_common, 'linear');
alpha12_interp = interp1(x1, alpha12, x_common, 'linear');
alpha13_interp = interp1(x1, alpha13, x_common, 'linear');
alpha14_interp = interp1(x1, alpha14, x_common, 'linear');
alpha15_interp = interp1(x1, alpha15, x_common, 'linear');
alpha1PBC_interp = interp1(x1, alpha1PBC, x_common, 'linear');
alpha1CIF_interp = interp1(x1, alpha1CIF, x_common, 'linear');
alpha1_ref_interp = interp1(x_ref, alpha1_ref, x_common, 'linear');

alpha20_interp = interp1(x1, alpha20, x_common, 'linear');
alpha21_interp = interp1(x1, alpha21, x_common, 'linear');
alpha22_interp = interp1(x1, alpha22, x_common, 'linear');
alpha23_interp = interp1(x1, alpha23, x_common, 'linear');
alpha24_interp = interp1(x1, alpha24, x_common, 'linear');
alpha25_interp = interp1(x1, alpha25, x_common, 'linear');
alpha2PBC_interp = interp1(x1, alpha2PBC, x_common, 'linear');
alpha2CIF_interp = interp1(x1, alpha2CIF, x_common, 'linear');
alpha2_ref_interp = interp1(x_ref, alpha2_ref, x_common, 'linear');

alpha30_interp = interp1(x1, alpha30, x_common, 'linear');
alpha31_interp = interp1(x1, alpha31, x_common, 'linear');
alpha32_interp = interp1(x1, alpha32, x_common, 'linear');
alpha33_interp = interp1(x1, alpha33, x_common, 'linear');
alpha34_interp = interp1(x1, alpha34, x_common, 'linear');
alpha35_interp = interp1(x1, alpha35, x_common, 'linear');
alpha3PBC_interp = interp1(x1, alpha3PBC, x_common, 'linear');
alpha3CIF_interp = interp1(x1, alpha3CIF, x_common, 'linear');
alpha3_ref_interp = interp1(x_ref, alpha3_ref, x_common, 'linear');

alpha40_interp = interp1(x1, alpha40, x_common, 'linear');
alpha41_interp = interp1(x1, alpha41, x_common, 'linear');
alpha42_interp = interp1(x1, alpha42, x_common, 'linear');
alpha43_interp = interp1(x1, alpha43, x_common, 'linear');
alpha44_interp = interp1(x1, alpha44, x_common, 'linear');
alpha45_interp = interp1(x1, alpha45, x_common, 'linear');
alpha4PBC_interp = interp1(x1, alpha4PBC, x_common, 'linear');
alpha4CIF_interp = interp1(x1, alpha4CIF, x_common, 'linear');
alpha4_ref_interp = interp1(x_ref, alpha4_ref, x_common, 'linear');

alpha50_interp = interp1(x1, alpha50, x_common, 'linear');
alpha51_interp = interp1(x1, alpha51, x_common, 'linear');
alpha52_interp = interp1(x1, alpha52, x_common, 'linear');
alpha53_interp = interp1(x1, alpha53, x_common, 'linear');
alpha54_interp = interp1(x1, alpha54, x_common, 'linear');
alpha55_interp = interp1(x1, alpha55, x_common, 'linear');
alpha5PBC_interp = interp1(x1, alpha5PBC, x_common, 'linear');
alpha5CIF_interp = interp1(x1, alpha5CIF, x_common, 'linear');
alpha5_ref_interp = interp1(x_ref, alpha5_ref, x_common, 'linear');

% Compute the L2 norm over the interval
L2_norm_h0 = sqrt(trapz(x_common, (h0_interp - h_ref_interp).^2));
L2_norm_h1 = sqrt(trapz(x_common, (h1_interp - h_ref_interp).^2));
L2_norm_h2 = sqrt(trapz(x_common, (h2_interp - h_ref_interp).^2));
L2_norm_h3 = sqrt(trapz(x_common, (h3_interp - h_ref_interp).^2));
L2_norm_h4 = sqrt(trapz(x_common, (h4_interp - h_ref_interp).^2));
L2_norm_h5 = sqrt(trapz(x_common, (h5_interp - h_ref_interp).^2));
L2_norm_hPBC = sqrt(trapz(x_common, (hPBC_interp - h_ref_interp).^2));
L2_norm_hCIF = sqrt(trapz(x_common, (hCIF_interp - h_ref_interp).^2));

L2_norm_u0 = sqrt(trapz(x_common, (u0_interp - u_ref_interp).^2));
L2_norm_u1 = sqrt(trapz(x_common, (u1_interp - u_ref_interp).^2));
L2_norm_u2 = sqrt(trapz(x_common, (u2_interp - u_ref_interp).^2));
L2_norm_u3 = sqrt(trapz(x_common, (u3_interp - u_ref_interp).^2));
L2_norm_u4 = sqrt(trapz(x_common, (u4_interp - u_ref_interp).^2));
L2_norm_u5 = sqrt(trapz(x_common, (u5_interp - u_ref_interp).^2));
L2_norm_uPBC = sqrt(trapz(x_common, (uPBC_interp - u_ref_interp).^2));
L2_norm_uCIF = sqrt(trapz(x_common, (uCIF_interp - u_ref_interp).^2));

L2_norm_alpha10 = sqrt(trapz(x_common, (alpha10_interp - alpha1_ref_interp).^2));
L2_norm_alpha11 = sqrt(trapz(x_common, (alpha11_interp - alpha1_ref_interp).^2));
L2_norm_alpha12 = sqrt(trapz(x_common, (alpha12_interp - alpha1_ref_interp).^2));
L2_norm_alpha13 = sqrt(trapz(x_common, (alpha13_interp - alpha1_ref_interp).^2));
L2_norm_alpha14 = sqrt(trapz(x_common, (alpha14_interp - alpha1_ref_interp).^2));
L2_norm_alpha15 = sqrt(trapz(x_common, (alpha15_interp - alpha1_ref_interp).^2));
L2_norm_alpha1PBC = sqrt(trapz(x_common, (alpha1PBC_interp - alpha1_ref_interp).^2));
L2_norm_alpha1CIF = sqrt(trapz(x_common, (alpha1CIF_interp - alpha1_ref_interp).^2));

L2_norm_alpha20 = sqrt(trapz(x_common, (alpha20_interp - alpha2_ref_interp).^2));
L2_norm_alpha21 = sqrt(trapz(x_common, (alpha21_interp - alpha2_ref_interp).^2));
L2_norm_alpha22 = sqrt(trapz(x_common, (alpha22_interp - alpha2_ref_interp).^2));
L2_norm_alpha23 = sqrt(trapz(x_common, (alpha23_interp - alpha2_ref_interp).^2));
L2_norm_alpha24 = sqrt(trapz(x_common, (alpha24_interp - alpha2_ref_interp).^2));
L2_norm_alpha25 = sqrt(trapz(x_common, (alpha25_interp - alpha2_ref_interp).^2));
L2_norm_alpha2PBC = sqrt(trapz(x_common, (alpha2PBC_interp - alpha2_ref_interp).^2));
L2_norm_alpha2CIF = sqrt(trapz(x_common, (alpha2CIF_interp - alpha2_ref_interp).^2));

L2_norm_alpha30 = sqrt(trapz(x_common, (alpha30_interp - alpha3_ref_interp).^2));
L2_norm_alpha31 = sqrt(trapz(x_common, (alpha31_interp - alpha3_ref_interp).^2));
L2_norm_alpha32 = sqrt(trapz(x_common, (alpha32_interp - alpha3_ref_interp).^2));
L2_norm_alpha33 = sqrt(trapz(x_common, (alpha33_interp - alpha3_ref_interp).^2));
L2_norm_alpha34 = sqrt(trapz(x_common, (alpha34_interp - alpha3_ref_interp).^2));
L2_norm_alpha35 = sqrt(trapz(x_common, (alpha35_interp - alpha3_ref_interp).^2));
L2_norm_alpha3PBC = sqrt(trapz(x_common, (alpha3PBC_interp - alpha3_ref_interp).^2));
L2_norm_alpha3CIF = sqrt(trapz(x_common, (alpha3CIF_interp - alpha3_ref_interp).^2));

L2_norm_alpha40 = sqrt(trapz(x_common, (alpha40_interp - alpha4_ref_interp).^2));
L2_norm_alpha41 = sqrt(trapz(x_common, (alpha41_interp - alpha4_ref_interp).^2));
L2_norm_alpha42 = sqrt(trapz(x_common, (alpha42_interp - alpha4_ref_interp).^2));
L2_norm_alpha43 = sqrt(trapz(x_common, (alpha43_interp - alpha4_ref_interp).^2));
L2_norm_alpha44 = sqrt(trapz(x_common, (alpha44_interp - alpha4_ref_interp).^2));
L2_norm_alpha45 = sqrt(trapz(x_common, (alpha45_interp - alpha4_ref_interp).^2));
L2_norm_alpha4PBC = sqrt(trapz(x_common, (alpha4PBC_interp - alpha4_ref_interp).^2));
L2_norm_alpha4CIF = sqrt(trapz(x_common, (alpha4CIF_interp - alpha4_ref_interp).^2));

L2_norm_alpha50 = sqrt(trapz(x_common, (alpha50_interp - alpha5_ref_interp).^2));
L2_norm_alpha51 = sqrt(trapz(x_common, (alpha51_interp - alpha5_ref_interp).^2));
L2_norm_alpha52 = sqrt(trapz(x_common, (alpha52_interp - alpha5_ref_interp).^2));
L2_norm_alpha53 = sqrt(trapz(x_common, (alpha53_interp - alpha5_ref_interp).^2));
L2_norm_alpha54 = sqrt(trapz(x_common, (alpha54_interp - alpha5_ref_interp).^2));
L2_norm_alpha55 = sqrt(trapz(x_common, (alpha55_interp - alpha5_ref_interp).^2));
L2_norm_alpha5PBC = sqrt(trapz(x_common, (alpha5PBC_interp - alpha5_ref_interp).^2));
L2_norm_alpha5CIF = sqrt(trapz(x_common, (alpha5CIF_interp - alpha5_ref_interp).^2));

h_diff_norms = [L2_norm_h0,L2_norm_h1,L2_norm_h2,L2_norm_h3,L2_norm_h4,L2_norm_h5,L2_norm_hPBC,L2_norm_hCIF]/norm(h_ref);
u_diff_norms = [L2_norm_u0,L2_norm_u1,L2_norm_u2,L2_norm_u3,L2_norm_u4,L2_norm_u5,L2_norm_uPBC,L2_norm_uCIF]/norm(u_ref);
alpha1_diff_norms = [L2_norm_alpha10,L2_norm_alpha11,L2_norm_alpha12,L2_norm_alpha13,L2_norm_alpha14,L2_norm_alpha15,L2_norm_alpha1PBC,L2_norm_alpha1CIF]/norm(alpha1_ref);
alpha2_diff_norms = [L2_norm_alpha20,L2_norm_alpha21,L2_norm_alpha22,L2_norm_alpha23,L2_norm_alpha24,L2_norm_alpha25,L2_norm_alpha2PBC,L2_norm_alpha2CIF]/norm(alpha2_ref);
alpha3_diff_norms = [L2_norm_alpha30,L2_norm_alpha31,L2_norm_alpha32,L2_norm_alpha33,L2_norm_alpha34,L2_norm_alpha35,L2_norm_alpha3PBC,L2_norm_alpha3CIF]/norm(alpha3_ref);
alpha4_diff_norms = [L2_norm_alpha40,L2_norm_alpha41,L2_norm_alpha42,L2_norm_alpha43,L2_norm_alpha44,L2_norm_alpha45,L2_norm_alpha4PBC,L2_norm_alpha4CIF]/norm(alpha4_ref);
alpha5_diff_norms = [L2_norm_alpha50,L2_norm_alpha51,L2_norm_alpha52,L2_norm_alpha53,L2_norm_alpha54,L2_norm_alpha55,L2_norm_alpha5PBC,L2_norm_alpha5CIF]/norm(alpha5_ref);

% writematrix(h_diff_norms,'AdaptiveSWME_Paper\damBreak-and-smooth_linear_lambda1.0_nu0.1_error_h.csv')
% writematrix(u_diff_norms,'AdaptiveSWME_Paper\damBreak-and-smooth_linear_lambda1.0_nu0.1_error_u.csv')
% writematrix(alpha1_diff_norms,'AdaptiveSWME_Paper\damBreak-and-smooth_linear_lambda1.0_nu0.1_error_alpha1.csv')
% writematrix(alpha2_diff_norms,'AdaptiveSWME_Paper\damBreak-and-smooth_linear_lambda0.1_nu1.0_error_alpha2.csv')

x_axis = [0,1,2,3];

h_diff_norms

u_diff_norms

alpha1_diff_norms

alpha2_diff_norms

alpha3_diff_norms

alpha4_diff_norms

alpha5_diff_norms

times_rel