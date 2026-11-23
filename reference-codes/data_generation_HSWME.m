(* ::Package:: *)

(* ::Chapter:: *)
(*2D Depth-Projected Shallow Flow (with vertical transport, periodic in x) *)


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



(* ::Subsection::Closed:: *)
(*Simulation Setup and Time Integration*)


(* ::Input::Initialization:: *)
\[Gamma]=1/2 (1+1/Sqrt[3]);
\[Phi]1[y_]:=-LegendreP[1,2*y-1];
\[Phi]2[y_]:=LegendreP[2,2*y-1];
\[Phi]3[y_]:=-LegendreP[3,2*y-1];
\[Phi]4[y_]:=LegendreP[4,2*y-1];
\[Phi]5[y_]:=-LegendreP[5,2*y-1];
\[Phi]6[y_]:=LegendreP[6,2*y-1];

FiniteVolumeRun[nx_,ny_,tend_]:=Block[{H1,HU1,H2,HU2},
dx=(x2-x1)/nx;
xs[i_]:=x1+(i-0.5)dx;
dy=(y2-y1)/ny;
ys[j_]:=y1+(j-0.5)dy;
pts=Table[{xs[i],ys[j]},{i,1,nx},{j,1,ny}];

Hraw=Table[h0[xs[i]],{i,1,nx},{j,1,ny}];
HUraw=Table[h0[xs[i]]u0[xs[i],ys[j]],{i,1,nx},{j,1,ny}];
alpha1Mean=Table[dy \[Phi]1[ys[i]],{i,1,ny}];
alpha2Mean=Table[dy (\[Phi]2[ys[i-1/2]]+4*\[Phi]2[ys[i]]+\[Phi]2[ys[i+1/2]])/6,{i,1,ny}];
alpha3Mean=Table[dy (\[Phi]3[ys[i-1/2]]+4*\[Phi]3[ys[i]]+\[Phi]3[ys[i+1/2]])/6,{i,1,ny}];
alpha4Mean=Table[dy (7*\[Phi]4[ys[i-1/2]]+32*\[Phi]4[ys[i-1/4]]+12*\[Phi]4[ys[i]]+32*\[Phi]4[ys[i+1/4]]+7*\[Phi]4[ys[i+1/2]])/90,{i,1,ny}];
alpha5Mean=Table[dy (7*\[Phi]5[ys[i-1/2]]+32*\[Phi]5[ys[i-1/4]]+12*\[Phi]5[ys[i]]+32*\[Phi]5[ys[i+1/4]]+7*\[Phi]5[ys[i+1/2]])/90,{i,1,ny}];
alpha6Mean=Table[dy (41*\[Phi]6[ys[i-1/2]]+216*\[Phi]6[ys[i-1/3]]+27*\[Phi]6[ys[i-1/6]]+272*\[Phi]6[ys[i]]+27*\[Phi]6[ys[i+1/6]]+216*\[Phi]6[ys[i+1/3]]+41*\[Phi]6[ys[i+1/2]])/840,{i,1,ny}];

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

dt=0.8dx/cMax;
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
Show[Plot[0,{x,x1,x2},PlotLabel->"h(x,t) at time: "<>ToString[time]<>" (step: "<>ToString[step]<>", \[CapitalDelta]t = "<>ToString[dt]<>", CFL = "<>ToString[CFL]<>")",PlotRange->{{x1,x2},{0.8,2.0}},Frame->True,ImageSize->iSize],
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
x1=-1.0;
x2=1.0;
y1=0.0;
y2=1.0;

(* Initial height profiles *)
xOffset1=-0.2;
xOffset2=0.2;
hDamLarge[x_]:=If[xOffset1<x<xOffset2,1,3];
hDamMid[x_]:=If[xOffset1<x<xOffset2,1,2];
hDamSmall[x_]:=If[xOffset1<x<xOffset2,1,1.5];
hSmoothSteep[x_]:=1+1.0*Exp[-10*(x+0.1)*(x+0.1)];
hSmoothMid[x_]:=1+0.5*Exp[-5*(x+0.1)*(x+0.1)];
hSmoothSmall[x_]:=1+0.2*Exp[-2*(x+0.1)*(x+0.1)];
initialHeightNames={"hDamLarge","hDamMid","hDamSmall","hSmoothSteep","hSmoothMid","hSmoothSmall"};
initialVelocityNames={"uConst","uVarXnotvarZ","uVarZnotvarX","uVarXvarZ"};
initialVelocityMagnitudeNames={"uFast","uMid","uSlow"};
(* Initial velocity profiles *)
start=-0.4;
end=0.4;
rangeFactor=end-start;
phi1[\[Zeta]_]:=-LegendreP[1,2*\[Zeta]-1];phi2[\[Zeta]_]:=LegendreP[2,2*\[Zeta]-1];phi3[\[Zeta]_]:=-LegendreP[3,2*\[Zeta]-1];phi4[\[Zeta]_]:=LegendreP[4,2*\[Zeta]-1];phi5[\[Zeta]_]:=-LegendreP[5,2*\[Zeta]-1];phi6[\[Zeta]_]:=LegendreP[6,2*\[Zeta]-1];
uconstFast[x_,y_,mag_]:=1.0*mag;
uvarXnotvarZ[x_,y_,mag_]:=If[start<x<end,mag*(x-start)/(end-start),0];
uvarZnotvarX[x_,y_,mag_]:=(1-4 y+78y^2-500 y^3+1225 y^4-1260 y^5+462 y^6)*mag;
uvarXvarZ[x_,y_,mag_]:=uvarXnotvarZ[x,y,mag]-If[start<x<end,mag*(x-start)/(end-start),0]*0.5*phi1[y]+If[start<x<end,mag*(x-start)/(end-start),0]*0.5*phi2[y]+If[start<x<end,mag*(x-start)/(end-start),0]*0.5*phi3[y]-If[start<x<end,mag*(x-start)/(end-start),0]*0.5*phi4[y]-If[start<x<end,mag*(x-start)/(end-start),0]*0.5*phi5[y]+If[start<x<end,mag*(x-start)/(end-start),0]*0.5*phi6[y];

(* viscosities and slip lengths *)
viscosities={0.05,0.5,5.0};
slipLengths={0.05,0.5,5.0};

nx=400;  
ny=70;
tend=0.2;


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


(* ::Input::Initialization:: *)
For[i=1,i<=Length[viscosities],i++,
For[j=1,j<=Length[slipLengths],j++,
For[m=1,m<=Length[initialHeightNames],m++,
For[n=1,n<=Length[initialVelocityNames],n++,
For[o=1,o<=Length[initialVelocityMagnitudeNames],o++,
Which[initialHeightNames[[m]]=="hDamLarge",h0[x_]:=hDamLarge[x],
initialHeightNames[[m]]=="hDamMid",h0[x_]:=hDamMid[x],
initialHeightNames[[m]]=="hDamSmall",h0[x_]:=hDamSmall[x],
initialHeightNames[[m]]=="hSmoothSteep",h0[x_]:=hSmoothSteep[x],
initialHeightNames[[m]]=="hSmoothMid",h0[x_]:=hSmoothMid[x],
initialHeightNames[[m]]=="hSmoothSmall",h0[x_]:=hSmoothSmall[x]];
mag=0;
Which[initialVelocityMagnitudeNames[[o]]=="uFast",mag=1.0,
initialVelocityMagnitudeNames[[o]]=="uMid",mag=0.5,initialVelocityMagnitudeNames[[o]]=="uFast",mag=0.2];
Which[initialVelocityNames[[n]]=="uConst",u0[x_,y_]:=uconstFast[x,y,mag],
initialHeightNames[[n]]=="uVarXnotvarZ",u0[x_,y_]:=uvarXnotvarZ[x,y,mag],
initialHeightNames[[n]]=="uVarZnotvarX",u0[x_,y_]:=uvarZnotvarX[x,y,mag],
initialHeightNames[[n]]=="uVarXvarZ",u0[x_,y_]:=uvarXvarZ[x,y,mag]];
\[Chi]=slipLengths[[j]];
R=viscosities[[i]];
FiniteVolumeRun[nx,ny,tend];
Uraw=Hraw^-1 HUraw;
Referencevalues=Table[{xs[i],Hraw[[i,1]],Map[Mean,Uraw][[i]],Map[alpha1Mean . #&,Uraw][[i]],Map[alpha2Mean . #&,Uraw][[i]],Map[alpha3Mean . #&,Uraw][[i]],Map[alpha4Mean . #&,Uraw][[i]],Map[alpha5Mean . #&,Uraw][[i]],Map[alpha6Mean . #&,Uraw][[i]]},{i,1,nx}]//MatrixForm;
dataset=Flatten[Referencevalues];
foldername="error_trainingData/"<>initialHeightNames[[m]]<>"_"<>initialVelocityNames[[n]]<>"_"<>initialVelocityMagnitudeNames[[o]]<>"/";
filename="lambda"<>ToString[slipLengths[[j]]]<>"_"<>"viscosity"<>ToString[viscosities[[i]]]<>"_ref.csv";
outputname=foldername<>filename;
Export[outputname,dataset,"CSV"];
];
];
];
];
];
