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

name1 = 'smoothAndShockTube_order2_relaxation0.05_time1.0_Roe';
name2 = 'smoothAndShockTube_order3_relaxation0.05_time1.0_Roe';
name3 = 'smoothAndShockTube_order4_relaxation0.05_time1.0_Roe';
name4 = 'smoothAndShockTube_order5_relaxation0.05_time1.0_Roe';
name5 = 'smoothAndShockTube_order6_relaxation0.05_time1.0_Roe';
name6 = 'smoothAndShockTube_order7_relaxation0.05_time1.0_Roe';
name7 = 'smoothAndShockTube_order8_relaxation0.05_time1.0_Roe';
name8 = 'smoothAndShockTube_order9_relaxation0.05_time1.0_Roe';
nameRef = 'smoothAndShockTube_order10_relaxation0.05_time1.0_Roe';

% name1 = 'smoothAndShockTube_order2_relaxation0.05_time1.0_PRICE';
% name2 = 'smoothAndShockTube_order3_relaxation0.05_time1.0_PRICE';
% name3 = 'smoothAndShockTube_order4_relaxation0.05_time1.0_PRICE';
% name4 = 'smoothAndShockTube_order5_relaxation0.05_time1.0_PRICE';
% name5 = 'smoothAndShockTube_order6_relaxation0.05_time1.0_PRICE';
% name6 = 'smoothAndShockTube_order7_relaxation0.05_time1.0_PRICE';
% name7 = 'smoothAndShockTube_order8_relaxation0.05_time1.0_PRICE';
% name8 = 'smoothAndShockTube_order9_relaxation0.05_time1.0_PRICE';
% nameRef = 'smoothAndShockTube_order10_relaxation0.05_time1.0_PRICE';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% name1 = 'smoothAndShockTube_order2_relaxation0.1_time1.0_Roe';
% name2 = 'smoothAndShockTube_order3_relaxation0.1_time1.0_Roe';
% name3 = 'smoothAndShockTube_order4_relaxation0.1_time1.0_Roe';
% name4 = 'smoothAndShockTube_order5_relaxation0.1_time1.0_Roe';
% name5 = 'smoothAndShockTube_order6_relaxation0.1_time1.0_Roe';
% name6 = 'smoothAndShockTube_order7_relaxation0.1_time1.0_Roe';
% name7 = 'smoothAndShockTube_order8_relaxation0.1_time1.0_Roe';
% name8 = 'smoothAndShockTube_order9_relaxation0.1_time1.0_Roe';
% nameRef = 'smoothAndShockTube_order10_relaxation0.1_time1.0_Roe';

% name1 = 'smoothAndShockTube_order2_relaxation0.1_time1.0_PRICE';
% name2 = 'smoothAndShockTube_order3_relaxation0.1_time1.0_PRICE';
% name3 = 'smoothAndShockTube_order4_relaxation0.1_time1.0_PRICE';
% name4 = 'smoothAndShockTube_order5_relaxation0.1_time1.0_PRICE';
% name5 = 'smoothAndShockTube_order6_relaxation0.1_time1.0_PRICE';
% name6 = 'smoothAndShockTube_order7_relaxation0.1_time1.0_PRICE';
% name7 = 'smoothAndShockTube_order8_relaxation0.1_time1.0_PRICE';
% name8 = 'smoothAndShockTube_order9_relaxation0.1_time1.0_PRICE';
% nameRef = 'smoothAndShockTube_order10_relaxation0.1_time1.0_PRICE';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% name1 = 'smoothAndShockTube_order2_relaxation0.5_time1.0_Roe';
% name2 = 'smoothAndShockTube_order3_relaxation0.5_time1.0_Roe';
% name3 = 'smoothAndShockTube_order4_relaxation0.5_time1.0_Roe';
% name4 = 'smoothAndShockTube_order5_relaxation0.5_time1.0_Roe';
% name5 = 'smoothAndShockTube_order6_relaxation0.5_time1.0_Roe';
% name6 = 'smoothAndShockTube_order7_relaxation0.5_time1.0_Roe';
% name7 = 'smoothAndShockTube_order8_relaxation0.5_time1.0_Roe';
% name8 = 'smoothAndShockTube_order9_relaxation0.5_time1.0_Roe';
% nameRef = 'smoothAndShockTube_order10_relaxation0.5_time1.0_Roe';

% name1 = 'smoothAndShockTube_order2_relaxation0.5_time1.0_PRICE';
% name2 = 'smoothAndShockTube_order3_relaxation0.5_time1.0_PRICE';
% name3 = 'smoothAndShockTube_order4_relaxation0.5_time1.0_PRICE';
% name4 = 'smoothAndShockTube_order5_relaxation0.5_time1.0_PRICE';
% name5 = 'smoothAndShockTube_order6_relaxation0.5_time1.0_PRICE';
% name6 = 'smoothAndShockTube_order7_relaxation0.5_time1.0_PRICE';
% name7 = 'smoothAndShockTube_order8_relaxation0.5_time1.0_PRICE';
% name8 = 'smoothAndShockTube_order9_relaxation0.5_time1.0_PRICE';
% nameRef = 'smoothAndShockTube_order10_relaxation0.5_time1.0_PRICE';

% Moment model data

A = readDataHME_matrixForm(name1);
nr_of_cells = size(A,1);

data_matrix = zeros(9,nr_of_cells,12);

data_matrix(1,:,:) = readDataHME_matrixForm(name1);
data_matrix(2,:,:) = readDataHME_matrixForm(name2);
data_matrix(3,:,:) = readDataHME_matrixForm(name3);
data_matrix(4,:,:) = readDataHME_matrixForm(name4);
data_matrix(5,:,:) = readDataHME_matrixForm(name5);
data_matrix(6,:,:) = readDataHME_matrixForm(name6);
data_matrix(7,:,:) = readDataHME_matrixForm(name7);
data_matrix(8,:,:) = readDataHME_matrixForm(name8);
data_matrix(9,:,:) = readDataHME_matrixForm(nameRef);

% [x2,rho2,u2,T2,f32,f42,f52,f62,f72,f82,f92,f102] = readDataHME(name2);
% [x3,rho3,u3,T3,f33,f43,f53,f63,f73,f83,f93,f103] = readDataHME(name3);
% [x4,rho4,u4,T4,f34,f44,f54,f64,f74,f84,f94,f104] = readDataHME(name4);
% [x5,rho5,u5,T5,f35,f45,f55,f65,f75,f85,f95,f105] = readDataHME(name5);
% [x6,rho6,u6,T6,f36,f46,f56,f66,f76,f86,f96,f106] = readDataHME(name6);
% [x7,rho7,u7,T7,f37,f47,f57,f67,f77,f87,f97,f107] = readDataHME(name7);
% [x8,rho8,u8,T8,f38,f48,f58,f68,f78,f88,f98,f108] = readDataHME(name8);

% Reference data
% 
% [x_ref,rho_ref,u_ref,T_ref,f3_ref,f4_ref,f5_ref,f6_ref,f7_ref,f8_ref,f9_ref,f10_ref] = readDataHME_matrixForm(nameRef);


%% Error convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Error convergence analysis %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error_norms = zeros(8,10);

for i = 1:8
    for j = 1:10
        error_norms(i,j) = norm(data_matrix(i,:,j+1)-data_matrix(9,:,j+1))
    end
end

x_axis = [0,1,2,3,4,5,6,7];

plot(x_axis,error_norms(:,4))