%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Spatially Adaptive Moment error computing %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear statements

clear all
close all
clc
format long

%% Simulation files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Specify the used simulation files %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name_Adaptive1 = "linear+shockTube_Relaxation10.0_Mesh100+200_spatiallyAdaptiveNeighbour_Moments46.csv";
name_Reference1 = "linear+shockTube_Relaxation10.0_Mesh100+200_Reference_Moments46.csv";
name_LowOrder1 = "linear+shockTube_Relaxation10.0_Mesh100+200_LowOrder_Moments46.csv";

name_Adaptive2 = "linear+shockTube_Relaxation1.0_Mesh100+200_spatiallyAdaptiveNeighbour_Moments46.csv";
name_Reference2 = "linear+shockTube_Relaxation1.0_Mesh100+200_Reference_Moments46.csv";
name_LowOrder2 = "linear+shockTube_Relaxation1.0_Mesh100+200_LowOrder_Moments46.csv";

name_Adaptive3 = "linear+shockTube_Relaxation0.1_Mesh100+200_spatiallyAdaptiveNeighbour_Moments46.csv";
name_Reference3 = "linear+shockTube_Relaxation0.1_Mesh100+200_Reference_Moments46.csv";
name_LowOrder3 = "linear+shockTube_Relaxation0.1_Mesh100+200_LowOrder_Moments46.csv";

data_Adaptive1 = load(name_Adaptive1);
data_Reference1 = load(name_Reference1);
data_LowOrder1 = load(name_LowOrder1);

data_Adaptive2 = load(name_Adaptive2);
data_Reference2 = load(name_Reference2);
data_LowOrder2 = load(name_LowOrder2);

data_Adaptive3 = load(name_Adaptive3);
data_Reference3 = load(name_Reference3);
data_LowOrder3 = load(name_LowOrder3);

x_Adaptive1 = data_Adaptive1(:,1);
rho_Adaptive1 = data_Adaptive1(:,2);
u_Adaptive1 = data_Adaptive1(:,3);
T_Adaptive1 = data_Adaptive1(:,4);

x_Reference1 = data_Reference1(:,1);
rho_Reference1 = data_Reference1(:,2);
u_Reference1 = data_Reference1(:,3);
T_Reference1 = data_Reference1(:,4);

x_LowOrder1 = data_LowOrder1(:,1);
rho_LowOrder1 = data_LowOrder1(:,2);
u_LowOrder1 = data_LowOrder1(:,3);
T_LowOrder1 = data_LowOrder1(:,4);

x_Adaptive2 = data_Adaptive2(:,1);
rho_Adaptive2 = data_Adaptive2(:,2);
u_Adaptive2 = data_Adaptive2(:,3);
T_Adaptive2 = data_Adaptive2(:,4);

x_Reference2 = data_Reference2(:,1);
rho_Reference2 = data_Reference2(:,2);
u_Reference2 = data_Reference2(:,3);
T_Reference2 = data_Reference2(:,4);

x_LowOrder2 = data_LowOrder2(:,1);
rho_LowOrder2 = data_LowOrder2(:,2);
u_LowOrder2 = data_LowOrder2(:,3);
T_LowOrder2 = data_LowOrder2(:,4);

x_Adaptive3 = data_Adaptive3(:,1);
rho_Adaptive3 = data_Adaptive3(:,2);
u_Adaptive3 = data_Adaptive3(:,3);
T_Adaptive3 = data_Adaptive3(:,4);

x_Reference3 = data_Reference3(:,1);
rho_Reference3 = data_Reference3(:,2);
u_Reference3 = data_Reference3(:,3);
T_Reference3 = data_Reference3(:,4);

x_LowOrder3 = data_LowOrder3(:,1);
rho_LowOrder3 = data_LowOrder3(:,2);
u_LowOrder3 = data_LowOrder3(:,3);
T_LowOrder3 = data_LowOrder3(:,4);

%% Error computing

rho_Diff_LowHigh1 = norm(rho_LowOrder1-rho_Reference1);
rho_Diff_LowAdaptive1 = norm(rho_LowOrder1-rho_Adaptive1);
rho_Diff_AdaptiveHigh1 = norm(rho_Adaptive1-rho_Reference1);

u_Diff_LowHigh1 = norm(u_LowOrder1-u_Reference1);
u_Diff_LowAdaptive1 = norm(u_LowOrder1-u_Adaptive1);
u_Diff_AdaptiveHigh1 = norm(u_Adaptive1-u_Reference1);

T_Diff_LowHigh1 = norm(T_LowOrder1-T_Reference1);
T_Diff_LowAdaptive1 = norm(T_LowOrder1-T_Adaptive1);
T_Diff_AdaptiveHigh1 = norm(T_Adaptive1-T_Reference1);

rho_Diff_LowHigh_ref1 = norm(rho_LowOrder1-rho_Reference1)/norm(rho_Reference1);
rho_Diff_LowAdaptive_ref1 = norm(rho_LowOrder1-rho_Adaptive1)/norm(rho_Reference1);
rho_Diff_AdaptiveHigh_ref1 = norm(rho_Adaptive1-rho_Reference1)/norm(rho_Reference1);

u_Diff_LowHigh_ref1 = norm(u_LowOrder1-u_Reference1)/norm(u_Reference1);
u_Diff_LowAdaptive_ref1 = norm(u_LowOrder1-u_Adaptive1)/norm(u_Reference1);
u_Diff_AdaptiveHigh_ref1 = norm(u_Adaptive1-u_Reference1)/norm(u_Reference1);

T_Diff_LowHigh_ref1 = norm(T_LowOrder1-T_Reference1)/norm(T_Reference1);
T_Diff_LowAdaptive_ref1 = norm(T_LowOrder1-T_Adaptive1)/norm(T_Reference1);
T_Diff_AdaptiveHigh_ref1 = norm(T_Adaptive1-T_Reference1)/norm(T_Reference1);

rho_Diff_LowHigh2 = norm(rho_LowOrder2-rho_Reference2);
rho_Diff_LowAdaptive2 = norm(rho_LowOrder2-rho_Adaptive2);
rho_Diff_AdaptiveHigh2 = norm(rho_Adaptive2-rho_Reference2);

u_Diff_LowHigh2 = norm(u_LowOrder2-u_Reference2);
u_Diff_LowAdaptive2 = norm(u_LowOrder2-u_Adaptive2);
u_Diff_AdaptiveHigh2 = norm(u_Adaptive2-u_Reference2);

T_Diff_LowHigh2 = norm(T_LowOrder2-T_Reference2);
T_Diff_LowAdaptive2 = norm(T_LowOrder2-T_Adaptive2);
T_Diff_AdaptiveHigh2 = norm(T_Adaptive2-T_Reference2);

rho_Diff_LowHigh_ref2 = norm(rho_LowOrder2-rho_Reference2)/norm(rho_Reference2);
rho_Diff_LowAdaptive_ref2 = norm(rho_LowOrder2-rho_Adaptive2)/norm(rho_Reference2);
rho_Diff_AdaptiveHigh_ref2 = norm(rho_Adaptive2-rho_Reference2)/norm(rho_Reference2);

u_Diff_LowHigh_ref2 = norm(u_LowOrder2-u_Reference2)/norm(u_Reference2);
u_Diff_LowAdaptive_ref2 = norm(u_LowOrder2-u_Adaptive2)/norm(u_Reference2);
u_Diff_AdaptiveHigh_ref2 = norm(u_Adaptive2-u_Reference2)/norm(u_Reference2);

T_Diff_LowHigh_ref2 = norm(T_LowOrder2-T_Reference2)/norm(T_Reference2);
T_Diff_LowAdaptive_ref2 = norm(T_LowOrder2-T_Adaptive2)/norm(T_Reference2);
T_Diff_AdaptiveHigh_ref2 = norm(T_Adaptive2-T_Reference2)/norm(T_Reference2);

rho_Diff_LowHigh3 = norm(rho_LowOrder3-rho_Reference3);
rho_Diff_LowAdaptive3 = norm(rho_LowOrder3-rho_Adaptive3);
rho_Diff_AdaptiveHigh3 = norm(rho_Adaptive3-rho_Reference3);

u_Diff_LowHigh3 = norm(u_LowOrder3-u_Reference3);
u_Diff_LowAdaptive3 = norm(u_LowOrder3-u_Adaptive3);
u_Diff_AdaptiveHigh3 = norm(u_Adaptive3-u_Reference3);

T_Diff_LowHigh3 = norm(T_LowOrder3-T_Reference3);
T_Diff_LowAdaptive3 = norm(T_LowOrder3-T_Adaptive3);
T_Diff_AdaptiveHigh3 = norm(T_Adaptive3-T_Reference3);

rho_Diff_LowHigh_ref3 = norm(rho_LowOrder3-rho_Reference3)/norm(rho_Reference3);
rho_Diff_LowAdaptive_ref3 = norm(rho_LowOrder3-rho_Adaptive3)/norm(rho_Reference3);
rho_Diff_AdaptiveHigh_ref3 = norm(rho_Adaptive3-rho_Reference3)/norm(rho_Reference3);

u_Diff_LowHigh_ref3 = norm(u_LowOrder3-u_Reference3)/norm(u_Reference3);
u_Diff_LowAdaptive_ref3 = norm(u_LowOrder3-u_Adaptive3)/norm(u_Reference3);
u_Diff_AdaptiveHigh_ref3 = norm(u_Adaptive3-u_Reference3)/norm(u_Reference3);

T_Diff_LowHigh_ref3 = norm(T_LowOrder3-T_Reference3)/norm(T_Reference3);
T_Diff_LowAdaptive_ref3 = norm(T_LowOrder3-T_Adaptive3)/norm(T_Reference3);
T_Diff_AdaptiveHigh_ref3 = norm(T_Adaptive3-T_Reference3)/norm(T_Reference3);
%% Plotting

rho_diffs1 = [rho_Diff_LowHigh_ref1,rho_Diff_AdaptiveHigh_ref1];
u_diffs1 = [u_Diff_LowHigh1, u_Diff_AdaptiveHigh1];
T_diffs1 = [T_Diff_LowHigh1, T_Diff_AdaptiveHigh1];

rho_diffs2 = [rho_Diff_LowHigh_ref2,rho_Diff_AdaptiveHigh_ref2];
u_diffs2 = [u_Diff_LowHigh2, u_Diff_AdaptiveHigh2];
T_diffs2 = [T_Diff_LowHigh2, T_Diff_AdaptiveHigh2];

rho_diffs3 = [rho_Diff_LowHigh_ref3,rho_Diff_AdaptiveHigh_ref3];
u_diffs3 = [u_Diff_LowHigh3, u_Diff_AdaptiveHigh3];
T_diffs3 = [T_Diff_LowHigh3, T_Diff_AdaptiveHigh3];

names = categorical({'HSM4-HSM6','SAMM46-HSM6'});
names = reordercats(names,{'HSM4-HSM6','SAMM46-HSM6'});
names2 = categorical({'\tau = 10.0','\tau = 1.0','\tau = 0.1'});
names2 = reordercats(names2,{'\tau = 10.0','\tau = 1.0','\tau = 0.1'});
x = 1:3; 

y=[u_diffs1;u_diffs2;u_diffs3];
%y=[rho_diffs1+u_diffs1+T_diffs1;rho_diffs2+u_diffs2+T_diffs2;rho_diffs3+u_diffs3+T_diffs3];

b=bar(x,y,0.72);
hAx=gca;
xticks(1:numel(names2))
xticklabels(names2)
hAx.XAxis.FontSize=14;

%ylabel('$\textsf{X-Label} \: x_{label,LaTeX} \textsf{(W)}$','Interpreter','latex')
ylabel('Relative distance of velocity $u$','Interpreter','latex','FontSize',14);

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels1 = string(b(1).YData);
labels2 = string(b(2).YData);
for i = 1:3
    labels1(i)=sprintf('%.1e',labels1(i));
    labels2(i)=sprintf('%.1e',labels2(i));
end
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',10)
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',10)

legend=legend(names,'Location','northeast');
set(legend,'FontSize',14);
ylim([0,1.2])

%% Export figure
 %addpath('C:\Users\rikve\Gitlab\PhD-RUG\spatially-adaptive-moment-models\RGD proceedings\Data\Export_fig\', '-end')
 %cd export_fig
 %export_fig('linearProfile+shock_mesh100+200_Time0.3_DifferenceBarGraph_u.pdf', '-pdf','-transparent');
