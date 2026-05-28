(* ::Package:: *)

(* ::Chapter:: *)
(*1D Depth-Projected Shallow Flow *)


(* ::Section:: *)
(*Function Definitions*)


(* ::Subsection::Closed:: *)
(*Tools*)


(* ::Input::Initialization:: *)
Clear["Global`*"]

interpolate[pts_,val_,dim_]:=Block[{len=Length[Flatten[val]]},
Interpolation[ArrayReshape[Flatten[Thread[{ArrayReshape[Flatten[pts],{len,dim}],Flatten[val]}]],{len,dim+1}],InterpolationOrder->1]
]
iSize=400;
color={Red,Green,Blue,Purple};


(* ::Subsection::Closed:: *)
(*Matrix Assemblation for 2D Data*)


(* ::Input::Initialization:: *)
Clear[H,dt,nx,ny,LS,A,id,R,\[Chi]]

Assemble[dt_,H_]:=Block[{nx,ny,id},
{nx,ny}=Dimensions[H];

A=SparseArray`SparseBlockMatrix[
Table[{ix,ix}->SparseArray[{
{1,1}->-(R/(H[[ix,1]]^2 dy^2)) (4 (3 H[[ix,1]]dy+2 \[Chi]))/(3 H[[ix,1]]dy+8 \[Chi]),{1,2}->-(R/(H[[ix,1]]^2 dy^2)) (4 (-H[[ix,1]]dy-2 \[Chi]))/(3 H[[ix,1]]dy+8 \[Chi]),
{ny,ny}->-(R/(H[[ix,ny]]^2 dy^2)),{ny,ny-1}->R/(H[[ix,ny]]^2 dy^2),
{i_,i_}:>-2 R/(H[[ix,i]]^2 dy^2),{i_,j_}/;i-j==1:>R/(H[[ix,i]]^2 dy^2),{i_,j_}/;j-i==1:>R/(H[[ix,i]]^2 dy^2)},{ny,ny}]
,{ix,1,nx}]
];

id=SparseArray[{{i_,i_}->1},{nx ny,nx ny}];

LS=LinearSolve[id-dt A];
];

impEulerY[rhs_]:=Block[{nx,ny},
{nx,ny}=Dimensions[rhs];
ArrayReshape[LS[Flatten[rhs]],{nx,ny}]
]



(* ::Subsection::Closed:: *)
(*Cell Interface Reconstruction*)


(* ::Input::Initialization:: *)
rc[dm_,dp_]=If[dp dm<=0,0,Sign[dp]Min[Abs[2dm],Abs[(2dp+dm)/3],Abs[3/2dp]]]; (*Limiter*)

PlusRecon[list_,opt_]:=Transpose[Map[Map[Function[{um,u,up},u+1/2 rc[u-um,up-u]]@@#&,Partition[ArrayPad[#,2,opt,InterpolationOrder->1],3,1]]&,list]]
MinusRecon[list_,opt_]:=Transpose[Map[Map[Function[{um,u,up},u-1/2 rc[up-u,u-um]]@@#&,Partition[ArrayPad[#,2,opt,InterpolationOrder->1],3,1]]&,list]]




(* ::Subsection::Closed:: *)
(*Full Grid Flux Calculation*)


(* ::Input::Initialization:: *)
kappa[list_]:=ArrayPad[Accumulate[list-Mean[list]],{1,0}] (* vertical averaging operator *)

Residuum[H_,HU_]:=Block[{nx,ny,Hleft,Hright,HUleft,HUright,U,HUbottom,HUtop,Hbottom,Htop,dxHU,W,CMax,delta=0.95},
{nx,ny}=Dimensions[HU];

{HUleft,Hleft,HUright,Hright,HUbottom,Hbottom,HUtop,Htop}=Parallelize[{
MinusRecon[Transpose[HU],"Fixed"], (*JK: boundary values now extrapolated constantly from interior*)
MinusRecon[Transpose[H],"Fixed"],(*JK: boundary values now extrapolated constantly from interior*)
PlusRecon[Transpose[HU],"Fixed"],(*JK: boundary values now extrapolated constantly from interior*)
PlusRecon[Transpose[H],"Fixed"],(*JK: boundary values now extrapolated constantly from interior*)
MinusRecon[HU,"Extrapolated"],
MinusRecon[H,"Extrapolated"],
PlusRecon[HU,"Extrapolated"],
PlusRecon[H,"Extrapolated"]}];

CMax=ListConvolve[{{1/2},{1/2}},H^-1 Abs[HU]+Sqrt[H],{1,-1},"Fixed"];(*JK: boundary values now extrapolated constantly from interior*)
dxHU=Differences[-(1/dx)Table[1/2 (HUright[[i]]+HUleft[[i+1]])+delta/2 CMax[[i]](Hright[[i]]-Hleft[[i+1]]),{i,1,nx+1}]] ;
W=Table[(Htop[[i]]+Hbottom[[i+1]])/2,{i,1,ny+1}]^-1 Transpose[Map[kappa,dy dxHU]];
cMax=Max[CMax]+Max[Abs[W]];


{Differences[-(1/dx)Table[1/2 (HUright[[i]]+HUleft[[i+1]])+delta/2 CMax[[i]](Hright[[i]]-Hleft[[i+1]]),{i,1,nx+1}]]+
Transpose[Differences[-(1/dy)Table[1/2 W[[i]](Htop[[i]]+Hbottom[[i+1]])+1/2 Abs[W[[i]]](Htop[[i]]-Hbottom[[i+1]]),{i,1,ny+1}]]],Differences[-(1/dx)Table[1/2 (HUright[[i]]^2/Hright[[i]]+Hright[[i]]^2/2+HUleft[[i+1]]^2/Hleft[[i+1]]+Hleft[[i+1]]^2/2)+delta/2 CMax[[i]](HUright[[i]]-HUleft[[i+1]]),{i,1,nx+1}]]+
Transpose[Differences[-(1/dy)Table[1/2 W[[i]](HUtop[[i]]+HUbottom[[i+1]])+1/2 Abs[W[[i]]](HUtop[[i]]-HUbottom[[i+1]]),{i,1,ny+1}]]]}
];



(* ::Subsection:: *)
(*Simulation Setup and Time Integration*)


(* ::Input::Initialization:: *)
\[Gamma]=1/2 (1+1/Sqrt[3]);
\[Phi]1[y_]=(1-2y);
\[Phi]2[y_]=(1-6 y+6 y^2);

FiniteVolumeRun[nx_,ny_,tend_]:=Block[{H1,HU1,H2,HU2},
dx=(x2-x1)/nx;
xs[i_]=x1+(i-0.5)dx;
dy=(y2-y1)/ny;
ys[j_]=y1+(j-0.5)dy;
pts=Table[{xs[i],ys[j]},{i,1,nx},{j,1,ny}];

Hraw=Table[h0[xs[i]],{i,1,nx},{j,1,ny}];
HUraw=Table[h0[xs[i]]u0[xs[i],ys[j]],{i,1,nx},{j,1,ny}];
alpha1Mean=Table[dy \[Phi]1[ys[i]],{i,1,ny}];
alpha2Mean=Table[dy (\[Phi]2[ys[i-1/2]]+4*\[Phi]2[ys[i]]+\[Phi]2[ys[i+1/2]])/6,{i,1,ny}];
alpha3Mean=Table[dy (\[Phi]3[ys[i-1/2]]+4*\[Phi]3[ys[i]]+\[Phi]3[ys[i+1/2]])/6,{i,1,ny}];
alpha4Mean=Table[dy (7*\[Phi]4[ys[i-1/2]]+32*\[Phi]4[ys[i-1/4]]+12*\[Phi]4[ys[i]]+32*\[Phi]4[ys[i+1/4]]+7*\[Phi]4[ys[i+1/2]])/90,{i,1,ny}];
alpha5Mean=Table[dy (7*\[Phi]5[ys[i-1/2]]+32*\[Phi]5[ys[i-1/4]]+12*\[Phi]5[ys[i]]+32*\[Phi]5[ys[i+1/4]]+7*\[Phi]5[ys[i+1/2]])/90,{i,1,ny}];

cMax=0.0;
dt=0.1dx;
time=0.0;
step=0;

While[time<tend,

{H1,HU1}={Hraw,HUraw}+dt Residuum[Hraw,HUraw];
{H2,HU2}=3/4{Hraw,HUraw}+1/4({H1,HU1}+dt Residuum[H1,HU1]);
{Hraw,HUraw}=1/3{Hraw,HUraw}+2/3({H2,HU2}+dt Residuum[H2,HU2]);

Assemble[\[Gamma] dt,Hraw];
Uraw=Hraw^-1 HUraw;
U1=impEulerY[Uraw];
U2=impEulerY[(3\[Gamma]-1)/\[Gamma] Uraw+(1-2\[Gamma])/\[Gamma] U1];
Uraw=(1-2\[Gamma])/\[Gamma] Uraw+(6\[Gamma]-3)/(2\[Gamma]) U1+1/(2\[Gamma]) U2;
HUraw=Hraw*Uraw;

time+=dt;
step++;

dt=0.6dx/cMax;
CFL=cMax dt/dx;
If[time+dt>=tend&&time<tend,dt=tend-time+10^-8];
];

{time,step}
];



(* ::Subsection::Closed:: *)
(*Visualization for Dynamic Output*)


(* ::Input::Initialization:: *)
Visual[H_,HU_,iSize_]:=Block[{nx,ny,U},
{nx,ny}=Dimensions[HU];
U=H^-1 HU;

Grid[{{
Show[Plot[0,{x,x1,x2},PlotLabel->"h(x,t) at time: "<>ToString[time]<>" (step: "<>ToString[step]<>", \[CapitalDelta]t = "<>ToString[dt]<>", CFL = "<>ToString[CFL]<>")",PlotRange->{{x1,x2},{0.8,5.5}},Frame->True,ImageSize->iSize],
ListPlot[{
Thread[{pts[[All,1,1]],H[[All,1]]}],
Thread[{pts[[All,1,1]],H[[All,Floor[ny/2]]]}],
Thread[{pts[[All,1,1]],H[[All,ny]]}]
},ImageSize->iSize]],
Show[Plot[{1000(x+0.5),1000x,1000(x-0.5)},{x,x1,x2},PlotStyle->{Red,Blue,Green},PlotLabel->"u_Ave(x,t)",PlotRange->{{x1,x2},{0,1}},Frame->True,ImageSize->iSize],
ListPlot[{Thread[{pts[[All,1,1]],Map[Mean,U]}]},ImageSize->iSize]],
ListPlot[{
Thread[{U[[Floor[nx/4],All]],pts[[1,All,2]]}],
Thread[{U[[Floor[nx/2],All]],pts[[1,All,2]]}],
Thread[{U[[3Floor[nx/4],All]],pts[[1,All,2]]}]
},Frame->True,PlotLabel->"velocity profiles at 3 positions",PlotStyle->{{Red,\[FilledCircle]},{Blue,\[FilledCircle]},{Green,\[FilledCircle]}},PlotRange->{{-0.4,2.7},{0,1}},ImageSize->iSize]}}]
]


(* ::Section:: *)
(*Simulations*)


(* ::Input::Initialization:: *)
x1=-3.0;
x2=2.0;
y1=0.0;
y2=1.0;

xOffset=0.5;
h0[x_]=If[x<0,5,1];
(*h0[x_]=If[x<-7,4,3+Exp[-1.5*(x-7)^2]];*)
h1[x_]=1+Exp[3Cos[\[Pi] (x+xOffset)]]/Exp[4];
h2[x_]=If[x>5,1.0,2.0];(*JK: adjusted IC*)
nx=800;  
ny=200;
tend=1.0;


(* ::Subsection:: *)
(*Time Evolution Output*)


(* ::Input::Initialization:: *)
show=False;
Button[" Show/Hide ",show=!show,BaseStyle->{"GenericButton",12}]
Dynamic[
If[show,
If[Length[Hraw]>0,
Visual[Hraw,HUraw,400],
"no data"
],""],SynchronousUpdating->True]


(* ::Subsection:: *)
(*Button["Open/Close", SelectionMove[EvaluationCell[], All, CellGroup]; FrontEndExecute[FrontEndToken["OpenCloseGroup"]], Appearance -> Automatic, Evaluator -> Automatic, Method -> "Preemptive"] Depth-Projected Shallow Flow (friction, linear / moderate quadratic /small quadratic profiles)*)


(* ::Input::Initialization:: *)
\[Chi]=0.1;
R=1.0;
nx=800;  
ny=200;
u0[x_,y_]:=0;
tend = 1.0;

FiniteVolumeRun[nx,ny,tend];
(*JK: boundary values now extrapolated constantly from interior*)
Uraw=Hraw^-1 HUraw;
ix=1;
hLL[ix]=interpolate[ArrayPad[pts[[All,1,1]],1,"Extrapolated"],ArrayPad[Map[Mean,Hraw],1,"Fixed"],1];
uMeanLL[ix]=interpolate[ArrayPad[pts[[All,1,1]],1,"Extrapolated"],ArrayPad[Map[Mean,Uraw],1,"Fixed"],1];
sLL[ix]=interpolate[ArrayPad[pts[[All,1,1]],1,"Extrapolated"],ArrayPad[Map[sMean . #&,Uraw],1,"Fixed"],1];
\[Kappa]LL[ix]=interpolate[ArrayPad[pts[[All,1,1]],1,"Extrapolated"],ArrayPad[Map[\[Kappa]Mean . #&,Uraw],1,"Fixed"],1];
uFullLL[ix]=interpolate[ArrayPad[pts,{1,1,0},"Extrapolated"],Transpose[ArrayPad[Transpose[ArrayPad[Uraw,{1},"Fixed"]],{1},"Extrapolated"]],2];
Referencevalues=Table[{xs[i],Hraw[[i,1]],Map[Mean,Uraw][[i]],Map[alpha1Mean . #&,Uraw][[i]],Map[alpha2Mean . #&,Uraw][[i]],Map[alpha3Mean . #&,Uraw][[i]],Map[alpha4Mean . #&,Uraw][[i]],Map[alpha5Mean . #&,Uraw][[i]],Map[alpha6Mean . #&,Uraw][[i]]},{i,1,nx}]//MatrixForm;
dataset=Flatten[Referencevalues];
Export["damBreakLarge_lambda0p1_nu1p0_t1p0_800x200.txt",dataset,"CSV"];
