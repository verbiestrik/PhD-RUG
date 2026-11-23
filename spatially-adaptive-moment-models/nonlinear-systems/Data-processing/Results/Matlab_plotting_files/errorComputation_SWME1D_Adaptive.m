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

% smooth Gaussian
name1 = 'smooth_lowOrder'; type1 = 'swe';
name2 = 'smooth_highOrder'; type2 = 'swme';
name3 = 'smooth_adaptive'; type3 = 'swme';

% Moment model data

[x1,h1,u1,alpha11,alpha21,alpha31,alpha41,alpha51] = readDataSWME1D(name1,type1);
[x2,h2,u2,alpha12,alpha22,alpha32,alpha42,alpha52] = readDataSWME1D(name2,type2);
[x3,h3,u3,alpha13,alpha23,alpha33,alpha43,alpha53] = readDataSWME1D(name3,type3);

%% Error convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Error calculation %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_diff_norm1 = norm(h1 - h2)/norm(h2);
h_diff_norm2 = norm(h3 - h2)/norm(h2);

u_diff_norm1 = norm(u1 - u2)/norm(u2);
u_diff_norm2 = norm(u3 - u2)/norm(u2);

alpha1_diff_norm1 = norm(alpha11 - alpha12)/norm(alpha12);
alpha1_diff_norm2 = norm(alpha13 - alpha12)/norm(alpha12);

h_diff_norms = [h_diff_norm1,h_diff_norm2];
u_diff_norms = [u_diff_norm1,u_diff_norm2];
alpha1_diff_norms = [alpha1_diff_norm1,alpha1_diff_norm2];

%writematrix(h_diff_norms,'ASWMEPaperRevised_SmoothExpit_h_Error_ASWME_BoundaryExcluded_Order4_4000x4_500x200.csv')
%writematrix(vr_diff_norms,'ASWMEPaperRevised_SmoothExpit_vr_Error_ASWME_BoundaryExcluded_Order4_4000x4_500x200.csv')
%writematrix(vtheta_diff_norms,'ASWMEPaperRevised_SmoothExpit_vtheta_Error_ASWME_BoundaryExcluded_Order4_4000x4_500x200.csv')

h_diff_norms
u_diff_norms
alpha1_diff_norms

x_axis = [0,1];

plot(x_axis,h_diff_norms)