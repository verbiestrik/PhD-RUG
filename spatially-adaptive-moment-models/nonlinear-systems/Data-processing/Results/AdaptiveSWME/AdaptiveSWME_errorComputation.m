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

h_interpol1 = zeros(length(x_ref),1);
h_interpol2 = zeros(length(x_ref),1);
h_interpol3 = zeros(length(x_ref),1);
h_interpol4 = zeros(length(x_ref),1);

u_interpol1 = zeros(length(x_ref),1);
u_interpol2 = zeros(length(x_ref),1);
u_interpol3 = zeros(length(x_ref),1);
u_interpol4 = zeros(length(x_ref),1);

alpha1_interpol1 = zeros(length(x_ref),1);
alpha1_interpol2 = zeros(length(x_ref),1);
alpha1_interpol3 = zeros(length(x_ref),1);
alpha1_interpol4 = zeros(length(x_ref),1);

alpha2_interpol1 = zeros(length(x_ref),1);
alpha2_interpol2 = zeros(length(x_ref),1);
alpha2_interpol3 = zeros(length(x_ref),1);
alpha2_interpol4 = zeros(length(x_ref),1);

start=1;

parity=(-1)^delta_ratio;

if(parity==1)
    for i = start:length(x_ref)
        h_interpol1(i) = (h1((i-1/2)*delta_ratio)+h1((i-1/2)*delta_ratio+1))/2;
        h_interpol2(i) = (h2((i-1/2)*delta_ratio)+h2((i-1/2)*delta_ratio+1))/2;
        h_interpol3(i) = (h3((i-1/2)*delta_ratio)+h3((i-1/2)*delta_ratio+1))/2;
        h_interpol4(i) = (h4((i-1/2)*delta_ratio)+h4((i-1/2)*delta_ratio+1))/2;
    
        u_interpol1(i) = (u1((i-1/2)*delta_ratio)+u1((i-1/2)*delta_ratio+1))/2;
        u_interpol2(i) = (u2((i-1/2)*delta_ratio)+u2((i-1/2)*delta_ratio+1))/2;
        u_interpol3(i) = (u3((i-1/2)*delta_ratio)+u3((i-1/2)*delta_ratio+1))/2;
        u_interpol4(i) = (u4((i-1/2)*delta_ratio)+u4((i-1/2)*delta_ratio+1))/2;

        alpha1_interpol1(i) = (alpha11((i-1/2)*delta_ratio)+alpha11((i-1/2)*delta_ratio+1))/2;
        alpha1_interpol2(i) = (alpha12((i-1/2)*delta_ratio)+alpha12((i-1/2)*delta_ratio+1))/2;
        alpha1_interpol3(i) = (alpha13((i-1/2)*delta_ratio)+alpha13((i-1/2)*delta_ratio+1))/2;
        alpha1_interpol4(i) = (alpha14((i-1/2)*delta_ratio)+alpha14((i-1/2)*delta_ratio+1))/2;

        alpha2_interpol1(i) = (alpha21((i-1/2)*delta_ratio)+alpha21((i-1/2)*delta_ratio+1))/2;
        alpha2_interpol2(i) = (alpha22((i-1/2)*delta_ratio)+alpha22((i-1/2)*delta_ratio+1))/2;
        alpha2_interpol3(i) = (alpha23((i-1/2)*delta_ratio)+alpha23((i-1/2)*delta_ratio+1))/2;
        alpha2_interpol4(i) = (alpha24((i-1/2)*delta_ratio)+alpha24((i-1/2)*delta_ratio+1))/2;

    end
end

if(parity==-1)
    for i = start:length(x_ref)
        h_interpol1(i) = h1((i-1/2)*delta_ratio+1/2);
        h_interpol2(i) = h2((i-1/2)*delta_ratio+1/2);
        h_interpol3(i) = h3((i-1/2)*delta_ratio+1/2);
        h_interpol4(i) = h4((i-1/2)*delta_ratio+1/2);
    
        u_interpol1(i) = u1((i-1/2)*delta_ratio+1/2);
        u_interpol2(i) = u2((i-1/2)*delta_ratio+1/2);
        u_interpol3(i) = u3((i-1/2)*delta_ratio+1/2);
        u_interpol4(i) = u4((i-1/2)*delta_ratio+1/2);
    
        alpha1_interpol1(i) = alpha11((i-1/2)*delta_ratio+1/2);
        alpha1_interpol2(i) = alpha12((i-1/2)*delta_ratio+1/2);
        alpha1_interpol3(i) = alpha13((i-1/2)*delta_ratio+1/2);
        alpha1_interpol4(i) = alpha14((i-1/2)*delta_ratio+1/2);

        alpha2_interpol1(i) = alpha21((i-1/2)*delta_ratio+1/2);
        alpha2_interpol2(i) = alpha22((i-1/2)*delta_ratio+1/2);
        alpha2_interpol3(i) = alpha23((i-1/2)*delta_ratio+1/2);
        alpha2_interpol4(i) = alpha24((i-1/2)*delta_ratio+1/2);
    end
end

h_diff_norm1 = norm(h_interpol1(start:end) - h_ref(start:end))/norm(h_ref(start:end));
h_diff_norm2 = norm(h_interpol2(start:end) - h_ref(start:end))/norm(h_ref(start:end));
h_diff_norm3 = norm(h_interpol3(start:end) - h_ref(start:end))/norm(h_ref(start:end));
h_diff_norm4 = norm(h_interpol4(start:end) - h_ref(start:end))/norm(h_ref(start:end));

u_diff_norm1 = norm(u_interpol1(start:end) - u_ref(start:end))/norm(u_ref(start:end));
u_diff_norm2 = norm(u_interpol2(start:end) - u_ref(start:end))/norm(u_ref(start:end));
u_diff_norm3 = norm(u_interpol3(start:end) - u_ref(start:end))/norm(u_ref(start:end));
u_diff_norm4 = norm(u_interpol4(start:end) - u_ref(start:end))/norm(u_ref(start:end));

alpha1_diff_norm1 = norm(alpha1_interpol1(start:end) - alpha1_ref(start:end))/norm(alpha1_ref(start:end));
alpha1_diff_norm2 = norm(alpha1_interpol2(start:end) - alpha1_ref(start:end))/norm(alpha1_ref(start:end));
alpha1_diff_norm3 = norm(alpha1_interpol3(start:end) - alpha1_ref(start:end))/norm(alpha1_ref(start:end));
alpha1_diff_norm4 = norm(alpha1_interpol4(start:end) - alpha1_ref(start:end))/norm(alpha1_ref(start:end));

alpha2_diff_norm1 = norm(alpha2_interpol1(start:end) - alpha2_ref(start:end))/norm(alpha2_ref(start:end));
alpha2_diff_norm2 = norm(alpha2_interpol2(start:end) - alpha2_ref(start:end))/norm(alpha2_ref(start:end));
alpha2_diff_norm3 = norm(alpha2_interpol3(start:end) - alpha2_ref(start:end))/norm(alpha2_ref(start:end));
alpha2_diff_norm4 = norm(alpha2_interpol4(start:end) - alpha2_ref(start:end))/norm(alpha2_ref(start:end));

h_diff_norms = [h_diff_norm1,h_diff_norm2,h_diff_norm3,h_diff_norm4];
u_diff_norms = [u_diff_norm1,u_diff_norm2,u_diff_norm3,u_diff_norm4];
alpha1_diff_norms = [alpha1_diff_norm1,alpha1_diff_norm2,alpha1_diff_norm3,alpha1_diff_norm4];
alpha2_diff_norms = [alpha2_diff_norm1,alpha2_diff_norm2,alpha2_diff_norm3,alpha2_diff_norm4];
% 
% writematrix(h_diff_norms,'AdaptiveSWME_Paper\damBreak-and-smooth_linear_lambda0.1_nu1.0_error_h.csv')
% writematrix(u_diff_norms,'AdaptiveSWME_Paper\damBreak-and-smooth_linear_lambda0.1_nu1.0_error_u.csv')
% writematrix(alpha1_diff_norms,'AdaptiveSWME_Paper\damBreak-and-smooth_linear_lambda0.1_nu1.0_error_alpha1.csv')
% writematrix(alpha2_diff_norms,'AdaptiveSWME_Paper\damBreak-and-smooth_linear_lambda0.1_nu1.0_error_alpha2.csv')

x_axis = [0,1,2,3];

plot(x_axis,alpha2_diff_norms)