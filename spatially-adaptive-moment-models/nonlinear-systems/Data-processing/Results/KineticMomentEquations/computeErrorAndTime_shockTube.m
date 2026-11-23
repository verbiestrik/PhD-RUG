clear all
close all
clc
format long

smooth_par = 50;
toldown = '0p001';
tolup = '0p0015';
relaxation_time = '0p05';

prediction_solver = 'PRICE';
interface_solver = 'PRICE';
interpolated = 'interp';
order_class_low = 10;
order_class_high = 12;

name2 = strcat('shockTube','_order',string(2),'_Kn',string(relaxation_time));type2='hme';
name4 = strcat('shockTube','_order',string(4),'_Kn',string(relaxation_time));type4='hme';
name6 = strcat('shockTube','_order',string(6),'_Kn',string(relaxation_time));type6='hme';
name8 = strcat('shockTube','_order',string(8),'_Kn',string(relaxation_time));type8='hme';
name10 = strcat('shockTube','_order',string(10),'_Kn',string(relaxation_time));type10='hme';
name12 = strcat('shockTube','_order',string(12),'_Kn',string(relaxation_time));type12='hme';
nameA = strcat('pred',prediction_solver,'_interface',interface_solver,...
    '_smoothpar',string(smooth_par),'_toldown',toldown,'_tolup',tolup,...
    '_Kn',relaxation_time,'_',interpolated);typeA='hmeAdaptive';

[x2,rho2,u2,T2,f32,f42,f52,f62,f72,f82,f92,f102,f112,f122,moments2] = readDataAdaptiveHME1D(name2,type2);
[x4,rho4,u4,T4,f34,f44,f54,f64,f74,f84,f94,f104,f114,f124,moments4] = readDataAdaptiveHME1D(name4,type4);
[x6,rho6,u6,T6,f36,f46,f56,f66,f76,f86,f96,f106,f116,f126,moments6] = readDataAdaptiveHME1D(name6,type6);
[x8,rho8,u8,T8,f38,f48,f58,f68,f78,f88,f98,f108,f118,f128,moments8] = readDataAdaptiveHME1D(name8,type8);
[x10,rho10,u10,T10,f310,f410,f510,f610,f710,f810,f910,f1010,f1110,f1210,moments10] = readDataAdaptiveHME1D(name10,type10);
[x12,rho12,u12,T12,f312,f412,f512,f612,f712,f812,f912,f1012,f1112,f1212,moments12] = readDataAdaptiveHME1D(name12,type12);
[xA,rhoA,uA,TA,f3A,f4A,f5A,f6A,f7A,f8A,f9A,f10A,f11A,f12A,momentsA] = readDataAdaptiveHME1D(nameA,typeA);

xS = -2; xE = 2; deltaX = 0.002;
x_dvm = xS+deltaX/2:deltaX:xE-deltaX/2;
rho_dvm = load(strcat('rho','_Kn',relaxation_time,'.dat'))';
u_dvm =  load(strcat('rho_u','_Kn',relaxation_time,'.dat'))'./load(strcat('rho','_Kn',relaxation_time,'.dat'))';

E_dvm = load(strcat('E','_Kn',relaxation_time,'.dat'))';
theta_dvm =  (E_dvm - 0.5 * rho_dvm .* u_dvm .* u_dvm)./(0.5* rho_dvm);

f3 = 0*rho_dvm;
f4 = 0*rho_dvm;
Q = 0*rho_dvm;
p = rho_dvm .* theta_dvm;

time2_arr = load(strcat('time',name2,'.csv'));
time4_arr = load(strcat('time',name4,'.csv'));
time6_arr = load(strcat('time',name6,'.csv'));
time8_arr = load(strcat('time',name8,'.csv'));
time10_arr = load(strcat('time',name10,'.csv'));
time12_arr = load(strcat('time',name12,'.csv'));
timeA_arr = load(strcat('time',nameA,'.csv'));

time2 = time2_arr;
time4 = time4_arr;
time6 = time6_arr;
time8 = time8_arr;
time10 = time10_arr;
time12 = time12_arr;
timeA = timeA_arr;

relTime2 = time2/time12;
relTime4 = time4/time12;
relTime6 = time6/time12;
relTime8 = time8/time12;
relTime10 = time10/time12;
relTime12 = time12/time12;
relTimeA = timeA/time12;

relTimes = [relTime2,relTime4,relTime6,relTime8,relTime10,relTime12,relTimeA];

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

x2=x2(:);x4=x4(:);x6=x6(:);x8=x8(:);x10=x10(:);x12=x12(:);
rho2=rho2(:);rho4=rho4(:);rho6=rho6(:);rho8=rho8(:);rho10=rho10(:);rho12=rho12(:);
u2=u2(:);u4=u4(:);u6=u6(:);u8=u8(:);u10=u10(:);u12=u12(:); 
T2=T2(:);T4=T4(:);T6=T6(:);T8=T8(:);T10=T10(:);T12=T12(:);

% compute overlap interval
x_min = max(min(x2), min(x_dvm));
x_max = min(max(x2), max(x_dvm));

% choose a grid density (e.g. max length of the two, or a fixed high resolution)
Ngrid = 10000;
x_grid = linspace(x_min, x_max, Ngrid);

rho2_on_grid = interp1(x2, rho2, x_grid, 'spline');
rho4_on_grid = interp1(x4, rho4, x_grid, 'spline');
rho6_on_grid = interp1(x6, rho6, x_grid, 'spline');
rho8_on_grid = interp1(x8, rho8, x_grid, 'spline');
rho10_on_grid = interp1(x10, rho10, x_grid, 'spline');
rho12_on_grid = interp1(x12, rho12, x_grid, 'spline');
rhoA_on_grid = interp1(xA, rhoA, x_grid, 'spline');
rho_dvm_on_grid = interp1(x_dvm,  rho_dvm, x_grid, 'spline');

u2_on_grid = interp1(x2, u2, x_grid, 'spline');
u4_on_grid = interp1(x4, u4, x_grid, 'spline');
u6_on_grid = interp1(x6, u6, x_grid, 'spline');
u8_on_grid = interp1(x8, u8, x_grid, 'spline');
u10_on_grid = interp1(x10, u10, x_grid, 'spline');
u12_on_grid = interp1(x12, u12, x_grid, 'spline');
uA_on_grid = interp1(xA, uA, x_grid, 'spline');
u_dvm_on_grid = interp1(x_dvm,  u_dvm, x_grid, 'spline');

T2_on_grid = interp1(x2, T2, x_grid, 'spline');
T4_on_grid = interp1(x4, T4, x_grid, 'spline');
T6_on_grid = interp1(x6, T6, x_grid, 'spline');
T8_on_grid = interp1(x8, T8, x_grid, 'spline');
T10_on_grid = interp1(x10, T10, x_grid, 'spline');
T12_on_grid = interp1(x12, T12, x_grid, 'spline');
TA_on_grid = interp1(xA, TA, x_grid, 'spline');
theta_dvm_on_grid = interp1(x_dvm,  theta_dvm, x_grid, 'spline');

% % drop NaNs if any, then compute norms
% mask = ~isnan(v1_on_grid) & ~isnan(v2_on_grid);
% L2_grid = norm(v1_on_grid(mask) - v2_on_grid(mask));
% rel_L2_grid = L2_grid / norm(v1_on_grid(mask));

rho_diff_norm_2 = norm(rho2_on_grid-rho_dvm_on_grid);
rho_diff_norm_4 = norm(rho4_on_grid-rho_dvm_on_grid);
rho_diff_norm_6 = norm(rho6_on_grid-rho_dvm_on_grid);
rho_diff_norm_8 = norm(rho8_on_grid-rho_dvm_on_grid);
rho_diff_norm_10 = norm(rho10_on_grid-rho_dvm_on_grid);
rho_diff_norm_12 = norm(rho12_on_grid-rho_dvm_on_grid);
rho_diff_norm_adaptive = norm(rhoA_on_grid-rho_dvm_on_grid);

u_diff_norm_2 = norm(u2_on_grid-u_dvm_on_grid);
u_diff_norm_4 = norm(u4_on_grid-u_dvm_on_grid);
u_diff_norm_6 = norm(u6_on_grid-u_dvm_on_grid);
u_diff_norm_8 = norm(u8_on_grid-u_dvm_on_grid);
u_diff_norm_10 = norm(u10_on_grid-u_dvm_on_grid);
u_diff_norm_12 = norm(u12_on_grid-u_dvm_on_grid);
u_diff_norm_adaptive = norm(uA_on_grid-u_dvm_on_grid);

T_diff_norm_2 = norm(T2_on_grid-theta_dvm_on_grid);
T_diff_norm_4 = norm(T4_on_grid-theta_dvm_on_grid);
T_diff_norm_6 = norm(T6_on_grid-theta_dvm_on_grid);
T_diff_norm_8 = norm(T8_on_grid-theta_dvm_on_grid);
T_diff_norm_10 = norm(T10_on_grid-theta_dvm_on_grid);
T_diff_norm_12 = norm(T12_on_grid-theta_dvm_on_grid);
T_diff_norm_adaptive = norm(TA_on_grid-theta_dvm_on_grid);

rho_diff_norms = [rho_diff_norm_2,rho_diff_norm_4,rho_diff_norm_6,rho_diff_norm_8,...
    rho_diff_norm_10,rho_diff_norm_12,rho_diff_norm_adaptive];
u_diff_norms = [u_diff_norm_2,u_diff_norm_4,u_diff_norm_6,u_diff_norm_8,...
    u_diff_norm_10,u_diff_norm_12,u_diff_norm_adaptive];
theta_diff_norms = [T_diff_norm_2,T_diff_norm_4,T_diff_norm_6,T_diff_norm_8,...
    T_diff_norm_10,T_diff_norm_12,T_diff_norm_adaptive];

rho_diff_norms_rel = rho_diff_norms/norm(rho_dvm_on_grid);
u_diff_norms_rel = u_diff_norms/norm(u_dvm_on_grid);
theta_diff_norms_rel = theta_diff_norms/norm(theta_dvm_on_grid);

smooth_par
toldown
tolup

% relTimes

% diff_norms_rel_sum = sqrt(rho_diff_norms_rel.*rho_diff_norms_rel...
%     +u_diff_norms_rel.*u_diff_norms_rel+theta_diff_norms_rel.*theta_diff_norms_rel)
diff_norms_rel_sum = rho_diff_norms_rel+u_diff_norms_rel+theta_diff_norms_rel

