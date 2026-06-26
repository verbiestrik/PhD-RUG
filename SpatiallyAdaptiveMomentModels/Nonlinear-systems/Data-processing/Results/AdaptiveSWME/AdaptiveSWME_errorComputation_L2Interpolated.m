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

name1 = 'damBreak-and-smooth_linear_lambda0.1_nu1.0_order2_t5_10000'; type1 = 'swme';
name2 = 'damBreak-and-smooth_linear_lambda0.1_nu1.0_order5_t5_10000'; type2 = 'swme';
name3 = 'damBreak-and-smooth_linear_lambda0.1_nu1.0_adaptiveNonConservative_t5_10000'; type3 = 'swmeAdaptive';
name4 = 'damBreak-and-smooth_linear_lambda0.1_nu1.0_adaptiveConservative_t5_10000'; type4 = 'swmeAdaptive';
reference_data = load("damBreak-and-smooth_linear_lambda0.1_nu1.0_t5_2000x200.csv");

% name1 = 'damBreak-and-smooth_linear_lambda0.1_nu0.1_order1_t5_10000'; type1 = 'swme';
% name2 = 'damBreak-and-smooth_linear_lambda0.1_nu0.1_order5_t5_10000'; type2 = 'swme';
% name3 = 'damBreak-and-smooth_linear_lambda0.1_nu0.1_adaptiveNonConservative_t5_10000'; type3 = 'swmeAdaptive';
% name4 = 'damBreak-and-smooth_linear_lambda0.1_nu0.1_adaptiveConservative_t5_10000'; type4 = 'swmeAdaptive';
% reference_data = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_2000x200.csv");

% name1 = 'damBreak-and-smooth_linear_lambda1.0_nu0.1_order1_t5_10000'; type1 = 'swme';
% name2 = 'damBreak-and-smooth_linear_lambda1.0_nu0.1_order5_t5_10000'; type2 = 'swme';
% name3 = 'damBreak-and-smooth_linear_lambda1.0_nu0.1_adaptiveNonConservative_t5_10000'; type3 = 'swmeAdaptive';
% name4 = 'damBreak-and-smooth_linear_lambda1.0_nu0.1_adaptiveConservative_t5_10000'; type4 = 'swmeAdaptive';
% reference_data = load("damBreak-and-smooth_linear_lambda1.0_nu0.1_t5_2000x200.csv");


% Moment model data

[x1,h1,u1,alpha11,alpha21,alpha31,alpha41,alpha51,moments1] = readDataAdaptiveSWME1D(name1,type1);
[x2,h2,u2,alpha12,alpha22,alpha32,alpha42,alpha52,moments2] = readDataAdaptiveSWME1D(name2,type2);
[x3,h3,u3,alpha13,alpha23,alpha33,alpha43,alpha53,moments3] = readDataAdaptiveSWME1D(name3,type3);
[x4,h4,u4,alpha14,alpha24,alpha34,alpha44,alpha54,moments4] = readDataAdaptiveSWME1D(name4,type4);

% Reference data

x_ref = reference_data(:,1); h_ref = reference_data(:,2); u_ref = reference_data(:,3);
alpha1_ref = 3.*reference_data(:,4); alpha2_ref = 5.*reference_data(:,5);


%% Grid information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Operations on the specified grid %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numberOfCells_MomentModel = 10000;

% Change this if necessary
numberOfCells_Reference = 2000;

x_1 = -20.0;
x_2 = 20.0;
delta_MomentModel = (x_2-x_1)/numberOfCells_MomentModel;
delta_Reference = (x_2-x_1)/numberOfCells_Reference;

delta_ratio = delta_Reference/delta_MomentModel;

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
h1_interp = interp1(x1, h1, x_common, 'linear');
h2_interp = interp1(x1, h2, x_common, 'linear');
h3_interp = interp1(x1, h3, x_common, 'linear');
h4_interp = interp1(x1, h4, x_common, 'linear');
h_ref_interp = interp1(x_ref, h_ref, x_common, 'linear');

u1_interp = interp1(x1, u1, x_common, 'linear');
u2_interp = interp1(x1, u2, x_common, 'linear');
u3_interp = interp1(x1, u3, x_common, 'linear');
u4_interp = interp1(x1, u4, x_common, 'linear');
u_ref_interp = interp1(x_ref, u_ref, x_common, 'linear');

alpha11_interp = interp1(x1, alpha11, x_common, 'linear');
alpha12_interp = interp1(x1, alpha12, x_common, 'linear');
alpha13_interp = interp1(x1, alpha13, x_common, 'linear');
alpha14_interp = interp1(x1, alpha14, x_common, 'linear');
alpha1_ref_interp = interp1(x_ref, alpha1_ref, x_common, 'linear');

alpha21_interp = interp1(x1, alpha21, x_common, 'linear');
alpha22_interp = interp1(x1, alpha22, x_common, 'linear');
alpha23_interp = interp1(x1, alpha23, x_common, 'linear');
alpha24_interp = interp1(x1, alpha24, x_common, 'linear');
alpha2_ref_interp = interp1(x_ref, alpha2_ref, x_common, 'linear');

% Compute the L2 norm over the interval
L2_norm_h1 = sqrt(trapz(x_common, (h1_interp - h_ref_interp).^2));
L2_norm_h2 = sqrt(trapz(x_common, (h2_interp - h_ref_interp).^2));
L2_norm_h3 = sqrt(trapz(x_common, (h3_interp - h_ref_interp).^2));
L2_norm_h4 = sqrt(trapz(x_common, (h4_interp - h_ref_interp).^2));

L2_norm_u1 = sqrt(trapz(x_common, (u1_interp - u_ref_interp).^2));
L2_norm_u2 = sqrt(trapz(x_common, (u2_interp - u_ref_interp).^2));
L2_norm_u3 = sqrt(trapz(x_common, (u3_interp - u_ref_interp).^2));
L2_norm_u4 = sqrt(trapz(x_common, (u4_interp - u_ref_interp).^2));

L2_norm_alpha11 = sqrt(trapz(x_common, (alpha11_interp - alpha1_ref_interp).^2));
L2_norm_alpha12 = sqrt(trapz(x_common, (alpha12_interp - alpha1_ref_interp).^2));
L2_norm_alpha13 = sqrt(trapz(x_common, (alpha13_interp - alpha1_ref_interp).^2));
L2_norm_alpha14 = sqrt(trapz(x_common, (alpha14_interp - alpha1_ref_interp).^2));

L2_norm_alpha21 = sqrt(trapz(x_common, (alpha21_interp - alpha2_ref_interp).^2));
L2_norm_alpha22 = sqrt(trapz(x_common, (alpha22_interp - alpha2_ref_interp).^2));
L2_norm_alpha23 = sqrt(trapz(x_common, (alpha23_interp - alpha2_ref_interp).^2));
L2_norm_alpha24 = sqrt(trapz(x_common, (alpha24_interp - alpha2_ref_interp).^2));

norm_h_ref = sqrt(trapz(x_common, (h_ref_interp).^2));
norm_u_ref = sqrt(trapz(x_common, (u_ref_interp).^2));
norm_alpha1_ref = sqrt(trapz(x_common, (alpha1_ref_interp).^2));
norm_alpha2_ref = sqrt(trapz(x_common, (alpha2_ref_interp).^2));

h_diff_norms = [L2_norm_h1,L2_norm_h2,L2_norm_h3,L2_norm_h4]/norm_h_ref;
u_diff_norms = [L2_norm_u1,L2_norm_u2,L2_norm_u3,L2_norm_u4]/norm_u_ref;
alpha1_diff_norms = [L2_norm_alpha11,L2_norm_alpha12,L2_norm_alpha13,L2_norm_alpha14]/norm_alpha1_ref;
alpha2_diff_norms = [L2_norm_alpha21,L2_norm_alpha22,L2_norm_alpha23,L2_norm_alpha24]/norm_alpha2_ref;
 
% writematrix(h_diff_norms,'AdaptiveSWME_Paper\damBreak-and-smooth_linear_lambda1.0_nu0.1_error_h.csv')
% writematrix(u_diff_norms,'AdaptiveSWME_Paper\damBreak-and-smooth_linear_lambda1.0_nu0.1_error_u.csv')
% writematrix(alpha1_diff_norms,'AdaptiveSWME_Paper\damBreak-and-smooth_linear_lambda1.0_nu0.1_error_alpha1.csv')
% writematrix(alpha2_diff_norms,'AdaptiveSWME_Paper\damBreak-and-smooth_linear_lambda0.1_nu1.0_error_alpha2.csv')

x_axis = [0,1,2,3];

h_diff_norms

u_diff_norms

alpha1_diff_norms

alpha2_diff_norms