(* ::Package:: *)

(* ::Input::Initialization:: *)
maxRecursion=500;


(* ::Input::Initialization:: *)
{u0,u1}={218,218/1.75};


(* ::Input::Initialization:: *)
{r0,r1}={1.5,1.5};


(* ::Input::Initialization:: *)
a=3.3;


(* ::Input::Initialization:: *)
I0=Assuming[{r\[Lambda]>0,k>0},Simplify@Integrate[r^2*1/(2*\[Pi])^(3/2)*Exp[-(r^2/(2*r\[Lambda]^2))]*SphericalBesselJ[0,k*r],{r,0,Infinity}]];


(* ::Input::Initialization:: *)
I1=Assuming[{r\[Lambda]>0,k>0},Simplify@Integrate[r^2*1/(2*\[Pi])^(3/2)*Exp[-(r^2/(2*r\[Lambda]^2))]*SphericalBesselJ[1,k*r],{r,0,Infinity}]];


(* ::Input::Initialization:: *)
I2=Assuming[{r\[Lambda]>0,k>0},Simplify@Integrate[r^2*1/(2*\[Pi])^(3/2)*Exp[-(r^2/(2*r\[Lambda]^2))]*SphericalBesselJ[2,k*r],{r,0,Infinity}]];


(* ::Input::Initialization:: *)
u::impl="Currently implemented: \[Lambda]=0,1";


(* ::Input::Initialization:: *)
u[\[Lambda]_]:=Switch[\[Lambda],0,u0,1,u1,_,Message[u::impl];]


(* ::Input::Initialization:: *)
I\[Lambda]::impl="Currently implemented: \[Lambda]=0,1,2";


(* ::Input::Initialization:: *)
I\[Lambda][\[Lambda]_,kk_]:=Switch[\[Lambda],0,I0/.{r\[Lambda]->r0,k->kk},1,I1/.{r\[Lambda]->r1,k->kk},2,I2/.{r\[Lambda]->r2,k->kk},_,Message[I\[Lambda]::impl];]


(* ::Input::Initialization:: *)
\[Epsilon][k_]:=k^2/2


(* ::Input::Initialization:: *)
\[Omega]k[k_,{n_}]:=Sqrt[\[Epsilon][k]*(\[Epsilon][k]+8*\[Pi]*a*n)]


(* ::Input::Initialization:: *)
U[\[Lambda]_,k_,{n_}]:=u[\[Lambda]]*Sqrt[(8*n*k^2*\[Epsilon][k])/(\[Omega]k[k,{n}]*(2*\[Lambda]+1))]*I\[Lambda][\[Lambda],k]
