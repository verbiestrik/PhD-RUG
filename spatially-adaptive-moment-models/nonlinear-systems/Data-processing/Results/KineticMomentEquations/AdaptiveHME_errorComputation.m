%% Clear statements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Computation of the HME model errors %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
format long

%% Simulation files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Specify the used simulation files %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prediction_solver = 'PRICE';
interface_solver = 'PRICE';
smooth_par = 200;
toldown = '0p001';
tolup = '0p0015';
relaxation_time = '0p5';
interpolated = 'interp';
order_class_low = 10;
order_class_high = 12;

name1 = strcat('shockTube','_order',string(order_class_low),'_Kn',string(relaxation_time));type1='hme';
name2 = strcat('shockTube','_order',string(order_class_high),'_Kn',string(relaxation_time));type2='hme';
name3 = strcat('pred',prediction_solver,'_interface',interface_solver,...
    '_smoothpar',string(smooth_par),'_toldown',toldown,'_tolup',tolup,...
    '_Kn',relaxation_time,'_',interpolated);type3='hmeAdaptive';

[x1,rho1,u1,T1,f31,f41,f51,f61,f71,f81,f91,f101,f111,f121,moments1] = readDataAdaptiveHME1D(name1,type1);
[x2,rho2,u2,T2,f32,f42,f52,f62,f72,f82,f92,f102,f112,f122,moments2] = readDataAdaptiveHME1D(name2,type2);
[x3,rho3,u3,T3,f33,f43,f53,f63,f73,f83,f93,f103,f113,f123,moments3] = readDataAdaptiveHME1D(name3,type3);

xA = -2; xE = 2; deltaX = 0.002;
x_dvm = xA+deltaX/2:deltaX:xE-deltaX/2;
rho_dvm = load(strcat('rho','_Kn',relaxation_time,'.dat'))';
u_dvm =  load(strcat('rho_u','_Kn',relaxation_time,'.dat'))'./load(strcat('rho','_Kn',relaxation_time,'.dat'))';

E_dvm = load(strcat('E','_Kn',relaxation_time,'.dat'))';
theta_dvm =  (E_dvm - 0.5 * rho_dvm .* u_dvm .* u_dvm)./(0.5* rho_dvm);

f3 = 0*rho_dvm;
f4 = 0*rho_dvm;
Q = 0*rho_dvm;
p = rho_dvm .* theta_dvm;

%% Error convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Error convergence analysis %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho_dvm = rho_dvm(251:1750);
rho_dvm = rho_dvm(:);
u_dvm = u_dvm(251:1750);
u_dvm = u_dvm(:);
theta_dvm = theta_dvm(251:1750);
theta_dvm = theta_dvm(:);
x_dvm = x_dvm(251:1750);
x_dvm = x_dvm(:);

x1 = x1(:); x2 = x2(:); x3 = x3(:);
rho1 = rho1(:); rho2 = rho2(:); rho3 = rho3(:);
u1 = u1(:); u2 = u2(:); u3 = u3(:); 
T1 = T1(:); T2 = T2(:); T3 = T3(:);

% compute overlap interval
x_min = max(min(x1), min(x_dvm));
x_max = min(max(x1), max(x_dvm));

% choose a grid density (e.g. max length of the two, or a fixed high resolution)
Ngrid = 10000;
x_grid = linspace(x_min, x_max, Ngrid);

rho1_on_grid = interp1(x1, rho1, x_grid, 'spline');
rho2_on_grid = interp1(x2, rho2, x_grid, 'spline');
rho3_on_grid = interp1(x3, rho3, x_grid, 'spline');
rho_dvm_on_grid = interp1(x_dvm,  rho_dvm, x_grid, 'spline');

u1_on_grid = interp1(x1, u1, x_grid, 'spline');
u2_on_grid = interp1(x2, u2, x_grid, 'spline');
u3_on_grid = interp1(x3, u3, x_grid, 'spline');
u_dvm_on_grid = interp1(x_dvm,  u_dvm, x_grid, 'spline');

T1_on_grid = interp1(x1, T1, x_grid, 'spline');
T2_on_grid = interp1(x2, T2, x_grid, 'spline');
T3_on_grid = interp1(x3, T3, x_grid, 'spline');
theta_dvm_on_grid = interp1(x_dvm,  theta_dvm, x_grid, 'spline');

% % drop NaNs if any, then compute norms
% mask = ~isnan(v1_on_grid) & ~isnan(v2_on_grid);
% L2_grid = norm(v1_on_grid(mask) - v2_on_grid(mask));
% rel_L2_grid = L2_grid / norm(v1_on_grid(mask));

rho_diff_norm_highOrder = norm(rho2_on_grid-rho_dvm_on_grid);
rho_diff_norm_lowOrder = norm(rho1_on_grid-rho_dvm_on_grid);
rho_diff_norm_adaptive = norm(rho3_on_grid-rho_dvm_on_grid);

u_diff_norm_highOrder = norm(u2_on_grid-u_dvm_on_grid);
u_diff_norm_lowOrder = norm(u1_on_grid-u_dvm_on_grid);
u_diff_norm_adaptive = norm(u3_on_grid-u_dvm_on_grid);

theta_diff_norm_highOrder = norm(T2_on_grid-theta_dvm_on_grid);
theta_diff_norm_lowOrder = norm(T1_on_grid-theta_dvm_on_grid);
theta_diff_norm_adaptive = norm(T3_on_grid-theta_dvm_on_grid);

rho_diff_norms = [rho_diff_norm_lowOrder,rho_diff_norm_highOrder,rho_diff_norm_adaptive];
u_diff_norms = [u_diff_norm_lowOrder,u_diff_norm_highOrder,u_diff_norm_adaptive];
theta_diff_norms = [theta_diff_norm_lowOrder,theta_diff_norm_highOrder,theta_diff_norm_adaptive];

rho_diff_norms_rel = rho_diff_norms/norm(rho_dvm_on_grid);
u_diff_norms_rel = u_diff_norms/norm(u_dvm_on_grid);
theta_diff_norms_rel = theta_diff_norms/norm(theta_dvm_on_grid)

smooth_par
toldown
tolup

% diff_norms_rel_sum = sqrt(rho_diff_norms_rel.*rho_diff_norms_rel...
%     +u_diff_norms_rel.*u_diff_norms_rel+theta_diff_norms_rel.*theta_diff_norms_rel)
diff_norms_rel_sum = rho_diff_norms_rel+u_diff_norms_rel+theta_diff_norms_rel

% x_axis = [0,1,2];

% plot(x_axis,theta_diff_norms_rel)
% plot(x_grid,T1_on_grid)