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

% Moment model data

[x1,h1,u1,alpha11,alpha21,alpha31,alpha41,alpha51,moments1] = readDataAdaptiveSWME1D(name1,type1);
[x2,h2,u2,alpha12,alpha22,alpha32,alpha42,alpha52,moments2] = readDataAdaptiveSWME1D(name2,type2);
[x3,h3,u3,alpha13,alpha23,alpha33,alpha43,alpha53,moments3] = readDataAdaptiveSWME1D(name3,type3);
[x4,h4,u4,alpha14,alpha24,alpha34,alpha44,alpha54,moments4] = readDataAdaptiveSWME1D(name4,type4);

% Reference data

reference_data200x20 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_200x20.csv");
reference_data200x50 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_200x50.csv");
reference_data200x100 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_200x100.csv");
reference_data400x40 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_400x40.csv");
reference_data400x100 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_400x100.csv");
reference_data400x200 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_400x200.csv");
reference_data500x50 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_500x50.csv");
reference_data500x125 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_500x125.csv");
reference_data500x250 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_500x250.csv");
reference_data1000x100 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_1000x100.csv");
reference_data1000x250 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_1000x250.csv");
reference_data1000x500 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_1000x500.csv");
reference_data2000x200 = load("damBreak-and-smooth_linear_lambda0.1_nu0.1_t5_2000x200.csv");

x_ref200x20 = reference_data200x20(:,1); h_ref200x20 = reference_data200x20(:,2); u_ref200x20 = reference_data200x20(:,3);
alpha1_ref200x20 = 3.*reference_data200x20(:,4); alpha2_ref200x20 = 5.*reference_data200x20(:,5);

x_ref200x50 = reference_data200x50(:,1); h_ref200x50 = reference_data200x50(:,2); u_ref200x50 = reference_data200x50(:,3);
alpha1_ref200x50 = 3.*reference_data200x50(:,4); alpha2_ref200x50 = 5.*reference_data200x50(:,5);

x_ref200x100 = reference_data200x100(:,1); h_ref200x100 = reference_data200x100(:,2); u_ref200x100 = reference_data200x100(:,3);
alpha1_ref200x100 = 3.*reference_data200x100(:,4); alpha2_ref200x100 = 5.*reference_data200x100(:,5);

%%%%

x_ref400x40 = reference_data400x40(:,1); h_ref400x40 = reference_data400x40(:,2); u_ref400x40 = reference_data400x40(:,3);
alpha1_ref400x40 = 3.*reference_data400x40(:,4); alpha2_ref400x40 = 5.*reference_data400x40(:,5);

x_ref400x100 = reference_data400x100(:,1); h_ref400x100 = reference_data400x100(:,2); u_ref400x100 = reference_data400x100(:,3);
alpha1_ref400x100 = 3.*reference_data400x100(:,4); alpha2_ref400x100 = 5.*reference_data400x100(:,5);

x_ref400x200 = reference_data400x200(:,1); h_ref400x200 = reference_data400x200(:,2); u_ref400x200 = reference_data400x200(:,3);
alpha1_ref400x200 = 3.*reference_data400x200(:,4); alpha2_ref400x200 = 5.*reference_data400x200(:,5);

%%%%

x_ref500x50 = reference_data500x50(:,1); h_ref500x50 = reference_data500x50(:,2); u_ref500x50 = reference_data500x50(:,3);
alpha1_ref500x50 = 3.*reference_data500x50(:,4); alpha2_ref500x50 = 5.*reference_data500x50(:,5);

x_ref500x125 = reference_data500x125(:,1); h_ref500x125 = reference_data500x125(:,2); u_ref500x125 = reference_data500x125(:,3);
alpha1_ref500x125 = 3.*reference_data500x125(:,4); alpha2_ref500x125 = 5.*reference_data500x125(:,5);

x_ref500x250 = reference_data500x250(:,1); h_ref500x250 = reference_data500x250(:,2); u_ref500x250 = reference_data500x250(:,3);
alpha1_ref500x250 = 3.*reference_data500x250(:,4); alpha2_ref500x250 = 5.*reference_data500x250(:,5);

%%%%

x_ref1000x100 = reference_data1000x100(:,1); h_ref1000x100 = reference_data1000x100(:,2); u_ref1000x100 = reference_data1000x100(:,3);
alpha1_ref1000x100 = 3.*reference_data1000x100(:,4); alpha2_ref1000x100 = 5.*reference_data1000x100(:,5);

x_ref1000x250 = reference_data1000x250(:,1); h_ref1000x250 = reference_data1000x250(:,2); u_ref1000x250 = reference_data1000x250(:,3);
alpha1_ref1000x250 = 3.*reference_data1000x250(:,4); alpha2_ref1000x250 = 5.*reference_data1000x250(:,5);

x_ref1000x500 = reference_data1000x500(:,1); h_ref1000x500 = reference_data1000x500(:,2); u_ref1000x500 = reference_data1000x500(:,3);
alpha1_ref1000x500 = 3.*reference_data1000x500(:,4); alpha2_ref1000x500 = 5.*reference_data1000x500(:,5);

%%%%

x_ref2000x200 = reference_data2000x200(:,1); h_ref2000x200 = reference_data2000x200(:,2); u_ref2000x200 = reference_data2000x200(:,3);
alpha1_ref2000x200 = 3.*reference_data2000x200(:,4); alpha2_ref2000x200 = 5.*reference_data2000x200(:,5);

%% Plotting to check if the reference simulation has converged

plotting = 'h';
% plotting = 'u';
% plotting = 'alpha1';
% plotting = 'alpha2';

if(strcmp(plotting,'h'))
    % plot1 = plot(x_ref200x50,h_ref200x50, ...
    %     x_ref200x100,h_ref200x100, ...
    %     x_ref400x100,h_ref400x100, ...
    %     x_ref400x200,h_ref400x200, ...
    %     x_ref500x125,h_ref500x125, ...
    %     x_ref500x250,h_ref500x250, ...
    %     x_ref1000x250,h_ref1000x250, ...
    %     x_ref1000x500,h_ref1000x500)
    
    plot1 = plot(x_ref200x20,h_ref200x20, ...
        x_ref400x40,h_ref400x40, ...
        x_ref500x50,h_ref500x50, ...
        x_ref1000x100,h_ref1000x100, ...
        x_ref2000x200,h_ref2000x200)

    axis([-20,20,2.8,4.25]);
    xlabel('x')
    ylabel('h')
    
    linewidth1 = 1;
    set(plot1(1:5),'LineWidth',linewidth1);

    grey = [0.4,0.4,0.4];
    % set(plot1(1:2),'LineStyle','-');
    % set(plot1(3:4),'LineStyle','-.');
    % set(plot1(5:6),'LineStyle','--');
    % set(plot1(7:8),'LineStyle',':');
    % set(plot1(1),'Color','r');
    % set(plot1(3),'Color','r');
    % set(plot1(5),'Color','r');
    % set(plot1(7),'Color','r');
    % set(plot1(2),'Color','blue');
    % set(plot1(4),'Color','blue');
    % set(plot1(6),'Color','blue');
    % set(plot1(8),'Color','blue');
    
    set(plot1(1),'LineStyle','-');
    set(plot1(2),'LineStyle','-.');
    set(plot1(3),'LineStyle','--');
    set(plot1(4),'LineStyle',':');
    set(plot1(5),'LineStyle','-');
    set(plot1(1),'Color','r');
    set(plot1(2),'Color','g');
    set(plot1(3),'Color','b');
    set(plot1(4),'Color','c');
    set(plot1(5),'Color','m');


    % leg = legend('200x50','200x100','400x100','400x200','500x125','500x250','1000x250','1000x500','Location','northeast');
    leg = legend('200x20','400x40','500x50','1000x100','2000x200','Location','northeast');
    set(leg,'FontSize',8); 
end
if(strcmp(plotting,'u'))
    plot1 = plot(x_ref200x50,u_ref200x50, ...
        x_ref200x100,u_ref200x100, ...
        x_ref400x100,u_ref400x100, ...
        x_ref400x200,u_ref400x200, ...
        x_ref500x125,u_ref500x125, ...
        x_ref500x250,u_ref500x250, ...
        x_ref1000x250,u_ref1000x250, ...
        x_ref1000x500,u_ref1000x500)
    
    axis([-20,20,2.8,4.25]);
    xlabel('x')
    ylabel('u')
    
    linewidth1 = 2;
    set(plot1(1:8),'LineWidth',linewidth1);

    grey = [0.4,0.4,0.4];
    set(plot1(1:2),'LineStyle','-');
    set(plot1(3:4),'LineStyle','-.');
    set(plot1(5:6),'LineStyle','--');
    set(plot1(7:8),'LineStyle',':');
    set(plot1(1),'Color','r');
    set(plot1(3),'Color','r');
    set(plot1(5),'Color','r');
    set(plot1(7),'Color','r');
    set(plot1(2),'Color','blue');
    set(plot1(4),'Color','blue');
    set(plot1(6),'Color','blue');
    set(plot1(8),'Color','blue');

    leg = legend('200x50','200x100','400x100','400x200','500x125','500x250','1000x250','1000x500','Location','northeast');
    set(leg,'FontSize',8); 
end
if(strcmp(plotting,'alpha1'))
    plot1 = plot(x_ref200x50,alpha1_ref200x50, ...
        x_ref200x100,alpha1_ref200x100, ...
        x_ref400x100,alpha1_ref400x100, ...
        x_ref400x200,alpha1_ref400x200, ...
        x_ref500x125,alpha1_ref500x125, ...
        x_ref500x250,alpha1_ref500x250, ...
        x_ref1000x250,alpha1_ref1000x250, ...
        x_ref1000x500,alpha1_ref1000x500)
    
    axis([-20,20,2.8,4.25]);
    xlabel('x')
    ylabel('\alpha_1')
    
    linewidth1 = 2;
    set(plot1(1:8),'LineWidth',linewidth1);

    grey = [0.4,0.4,0.4];
    set(plot1(1:2),'LineStyle','-');
    set(plot1(3:4),'LineStyle','-.');
    set(plot1(5:6),'LineStyle','--');
    set(plot1(7:8),'LineStyle',':');
    set(plot1(1),'Color','r');
    set(plot1(3),'Color','r');
    set(plot1(5),'Color','r');
    set(plot1(7),'Color','r');
    set(plot1(2),'Color','blue');
    set(plot1(4),'Color','blue');
    set(plot1(6),'Color','blue');
    set(plot1(8),'Color','blue');
    
    leg = legend('200x50','200x100','400x100','400x200','500x125','500x250','1000x250','1000x500','Location','northeast');
    set(leg,'FontSize',8); 
end
if(strcmp(plotting,'alpha2'))
    plot1 = plot(x_ref200x50,alpha2_ref200x50, ...
        x_ref200x100,alpha2_ref200x100, ...
        x_ref400x100,alpha2_ref400x100, ...
        x_ref400x200,alpha2_ref400x200, ...
        x_ref500x125,alpha2_ref500x125, ...
        x_ref500x250,alpha2_ref500x250, ...
        x_ref1000x250,alpha2_ref1000x250, ...
        x_ref1000x500,alpha2_ref1000x500)
    
    axis([-20,20,2.8,4.25]);
    xlabel('x')
    ylabel('\alpha_2')
    
    linewidth1 = 2;
    set(plot1(1:8),'LineWidth',linewidth1);

    grey = [0.4,0.4,0.4];
    set(plot1(1:2),'LineStyle','-');
    set(plot1(3:4),'LineStyle','-.');
    set(plot1(5:6),'LineStyle','--');
    set(plot1(7:8),'LineStyle',':');
    set(plot1(1),'Color','r');
    set(plot1(3),'Color','r');
    set(plot1(5),'Color','r');
    set(plot1(7),'Color','r');
    set(plot1(2),'Color','blue');
    set(plot1(4),'Color','blue');
    set(plot1(6),'Color','blue');
    set(plot1(8),'Color','blue');
    
    leg = legend('200x50','200x100','400x100','400x200','500x125','500x250','1000x250','1000x500','Location','northeast');
    set(leg,'FontSize',8); 
end


%% Grid information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Operations on the specified grid %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numberOfCells_MomentModel = 2000;

% Change this if necessary
numberOfCells_Reference = 1000;

x_1 = -20.0;
x_2 = 20.0;
delta_MomentModel = (x_2-x_1)/numberOfCells_MomentModel;
delta_Reference = (x_2-x_1)/numberOfCells_Reference;

delta_ratio = delta_Reference/delta_MomentModel;

%% Error convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Error convergence analysis %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_ref = x_ref1000x500;u_ref = u_ref1000x500;alpha1_ref = alpha1_ref1000x500;alpha2_ref = alpha2_ref1000x500;

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
    for i = start:length(r_ref)
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
u_diff_norms = [vr_diff_norm1,vr_diff_norm2,vr_diff_norm3,vr_diff_norm4,vr_diff_norm5];
alpha1_diff_norms = [alpha1_diff_norm1,alpha1_diff_norm2,alpha1_diff_norm3,alpha1_diff_norm4];
alpha2_diff_norms = [alpha2_diff_norm1,alpha2_diff_norm2,alpha2_diff_norm3,alpha2_diff_norm4];

%writematrix(h_diff_norms,'ASWMEPaperRevised_SmoothExpit_h_Error_ASWME_BoundaryExcluded_Order4_4000x4_500x200.csv')
%writematrix(vr_diff_norms,'ASWMEPaperRevised_SmoothExpit_vr_Error_ASWME_BoundaryExcluded_Order4_4000x4_500x200.csv')
%writematrix(vtheta_diff_norms,'ASWMEPaperRevised_SmoothExpit_vtheta_Error_ASWME_BoundaryExcluded_Order4_4000x4_500x200.csv')

x_axis = [0,1,2,3];

plot(x_axis,u_diff_norms)