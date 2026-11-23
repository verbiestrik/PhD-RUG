reference_data1 = load("ASWMEPaperRevised_SmoothExpit[1,8]_Time1.0_Reference_400x200_CFL0.5_New.csv");
reference_data2 = load("ASWMEPaperRevised_SmoothExpit[1,8]_Time1.0_Reference_600x200_CFL0.5.csv");

r_ref1 = reference_data1(:,1);
r_ref2 = reference_data2(:,1);

h_ref1 = reference_data1(:,2);
vr_ref1 = reference_data1(:,3);
alpha1_ref1 = reference_data1(:,4);
alpha1_ref1 = 3.*alpha1_ref1;
alpha2_ref1 = reference_data1(:,5);
alpha2_ref1 = 5.*alpha2_ref1;
vtheta_ref1 = reference_data1(:,6);

h_ref2 = reference_data2(:,2);
vr_ref2 = reference_data2(:,3);
alpha1_ref2 = reference_data2(:,4);
alpha1_ref2 = 3.*alpha1_ref2;
alpha2_ref2 = reference_data2(:,5);
alpha2_ref2 = 5.*alpha2_ref2;
vtheta_ref2 = reference_data2(:,6);

blue = [0, 0.4470, 0.7410];
green = [0.4660, 0.6740, 0.1880];
yellow = [0.9290, 0.6940, 0.1250];
red = [0.8500, 0.3250, 0.0980];
brown = [171, 104, 87]./255;

plot1 = plot(r_ref1,vr_ref1,r_ref2,vr_ref2)

        set(plot1(1),'Color',blue);
        set(plot1(1),'LineWidth',1.2);

        set(plot1(2),'Color',red);
        set(plot1(2),'LineWidth',1.2);

        leg = legend('ReferenceCoarse','RerferenceFine','Location','northwest');
        set(leg,'FontSize',12); 
    axis([1,8,-0.5,0.1]); % Change if necessary