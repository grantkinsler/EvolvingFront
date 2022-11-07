(* ::Package:: *)

(* ::Section::Initialization::Closed:: *)
(**)
(*1. Read trajectories, read depths, number trajectories*)


(* ::Subsection::Initialization::Closed:: *)
(**)
(*i. Read trajectories*)


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/CarbonData"];

AllTrajectories2M3=Import["2M3_Paper_unique_counts_corrected.txt", "Table"];
AllTrajectories4M3=Import["4M3_Paper_unique_counts_corrected.txt", "Table"];

AllReadTrajectories2M3=Delete[AllTrajectories2M3,1];
AllReadTrajectories4M3=Delete[AllTrajectories4M3,1];



(* ::Text::Initialization:: *)
(*The time points that are used in this analysis are:*)


(* ::Input::Initialization:: *)
TimePoints2M3={0,16,32, 40, 48, 64, 72, 80, 88, 96, 104, 112};
TimePoints4M3={0,8,24, 40, 48,56, 64, 72,88, 96};


(* ::Text::Initialization:: *)
(*This function combined all the reads at a given time point. One has to be careful in doing this because each is prepped separatedly so the error associated with a DNA prep (some fraction of \[Kappa])  will add up for each prep.*)


(* ::Input::Initialization:: *)
CombineTimePoints2M3[Trajectory_]:=Module[{t0, t8, t16, t24, t32, t40, t48, t64, t72, t80, t88, t96, t104, t112},
t0={{0,If[Sum[Trajectory[[i]], {i, 1, 5}]==0,1,Sum[Trajectory[[i]], {i, 1, 5}]]}};
t16={{16,Sum[Trajectory[[i]], {i, 6,6+1}]}};
t32={{32,Sum[Trajectory[[i]], {i, 8,8+1}]}};
t40={{40,Sum[Trajectory[[i]], {i, 10,10+1}]}};
t48={{48,Sum[Trajectory[[i]], {i, 12,12+1}]}};
t64={{64,Sum[Trajectory[[i]], {i, 14,14+1}]}};
t72={{72,Sum[Trajectory[[i]], {i, 16,16+1}]}};
t80={{80,Sum[Trajectory[[i]], {i, 18,18+1}]}};
t88={{88,Sum[Trajectory[[i]], {i, 20,20+1}]}};
t96={{96,Sum[Trajectory[[i]], {i, 22,22+1}]}};
t104={{104,Sum[Trajectory[[i]], {i, 24,24+1}]}};
t112={{112,Sum[Trajectory[[i]], {i, 26,26+1}]}};

Join[t0,t16,t32, t40, t48, t64, t72, t80, t88, t96, t104, t112]
];

CombineTimePoints4M3[Trajectory_]:=Module[{t0, t8, t16, t24, t32, t40, t48,t56, t64, t72, t80, t88, t96},
t0={{0,If[Sum[Trajectory[[i]], {i, 1, 5}]==0,1,Sum[Trajectory[[i]], {i, 1, 5}]]}};
t8={{8,Sum[Trajectory[[i]], {i, 6,6+1}]}};
t24={{24,Sum[Trajectory[[i]], {i, 8,8}]}};
t40={{40,Sum[Trajectory[[i]], {i, 9,9}]}};
t48={{48,Sum[Trajectory[[i]], {i, 10,10}]}};
t56={{56,Sum[Trajectory[[i]], {i, 11,11}]}};
t64={{64,Sum[Trajectory[[i]], {i, 12,12}]}};
t72={{72,Sum[Trajectory[[i]], {i, 13,13}]}};
t88={{88,Sum[Trajectory[[i]], {i, 14,14}]}};
t96={{96,Sum[Trajectory[[i]], {i, 15,15}]}};

Join[t0,t8,t24, t40, t48,t56, t64, t72, t88, t96]
];


(* ::Text::Initialization:: *)
(*List of combined reads for each barcode. Export to a file. *)


(* ::Input::Initialization:: *)
ReadTrajectories2M3=Map[CombineTimePoints2M3, AllReadTrajectories2M3];
Export["ReadTrajectories2M3.csv",ReadTrajectories2M3,"CSV"];
TableForm[ReadTrajectories2M3[[1;;10]], TableHeadings -> {Table[StringJoin["Barcode",ToString[i]], {i, 1, 10}],Table[StringJoin["t=",ToString[TimePoints2M3[[i]]]], {i, 1, Length[TimePoints2M3]}]}]

ReadTrajectories4M3=Map[CombineTimePoints4M3, AllReadTrajectories4M3];
Export["ReadTrajectories4M3.csv",ReadTrajectories4M3,"CSV"];
TableForm[ReadTrajectories4M3[[1;;10]], TableHeadings -> {Table[StringJoin["Barcode",ToString[i]], {i, 1, 10}],Table[StringJoin["t=",ToString[TimePoints4M3[[i]]]], {i, 1, Length[TimePoints4M3]}]}]



(* ::Subsection::Initialization::Closed:: *)
(**)
(*ii. Cell number trajectories*)


(* ::Text::Initialization:: *)
(*List of total read depths. Export to file.*)


(* ::Input::Initialization:: *)
ReadDepths2M3=Table[
ReadsByTime=Transpose[ReadTrajectories2M3];
ReadsAddedUp=Total[ReadsByTime[[i]]];
{TimePoints2M3[[i]],ReadsAddedUp[[2]]}
, {i, 1, Length[TimePoints2M3]}];
Export["ReadDepths2M3.csv",ReadDepths2M3,"CSV"];
TableForm[ReadDepths2M3, TableHeadings -> {{},{"t", "ReadDepth"}}]

ReadDepths4M3=Table[
ReadsByTime=Transpose[ReadTrajectories4M3];
ReadsAddedUp=Total[ReadsByTime[[i]]];
{TimePoints4M3[[i]],ReadsAddedUp[[2]]}
, {i, 1, Length[TimePoints4M3]}];
Export["ReadDepths4M3.csv",ReadDepths4M3,"CSV"];
TableForm[ReadDepths4M3, TableHeadings -> {{},{"t", "ReadDepth"}}]


(* ::Input::Initialization:: *)

ListLogPlot[{ReadDepths2M3, ReadDepths4M3}, Joined -> True]


(* ::Text::Initialization:: *)
(*Convert read trajectories to approximate number trajectories for plotting purposes. Export to file.*)


(* ::Input::Initialization:: *)

PopSize=5*10^8;



(* ::Input::Initialization:: *)
ExtinctLineages={};
CellNumberTrajectories2M3=Table[
If[Mod[Barcode, 50000]==0, Print[Barcode];];
NumberTraj1=Table[{ReadTrajectories2M3[[Barcode]][[i]][[1]],If[ReadTrajectories2M3[[Barcode]][[i]][[2]]==0,0,N[(PopSize/ReadDepths2M3[[i]][[2]])*(ReadTrajectories2M3[[Barcode]][[i]][[2]])] ]}, {i, 1, Length[ReadDepths2M3]}];

TimePointsBarcodeUnRead=Sort[Map[First,Select[NumberTraj1, #[[2]]==0&]], #1>#2&];
(*Print["TimePointsBarcodeUnRead=", TimePointsBarcodeUnRead];*)
TimePointsBarcodeRead=Complement[TimePoints2M3,TimePointsBarcodeUnRead];
(*Print["TimePointsBarcodeRead=", TimePointsBarcodeRead];*)

If[Length[TimePointsBarcodeUnRead]>0&&Max[TimePointsBarcodeRead]<Min[TimePointsBarcodeUnRead], ExtinctLineage=1; AppendTo[ExtinctLineages, Barcode];, ExtinctLineage=0];

If[ExtinctLineage==1, NumberTraj2=Replace[NumberTraj1, 0-> 0.1, {2}], NumberTraj2=NumberTraj1];

NumberTraj2

, {Barcode, 1,Length[ReadTrajectories2M3]}];

TableForm[CellNumberTrajectories2M3[[1;;10]], TableHeadings -> {Table[StringJoin["Barcode",ToString[i]], {i, 1, 10}],Table[StringJoin["t=",ToString[TimePoints2M3[[i]]]], {i, 1, Length[TimePoints2M3]}]}]
Export["CellNumberTrajectories2M3.csv", CellNumberTrajectories2M3, "CSV"];
Export["ExtinctLineages2M3.csv", ExtinctLineages, "CSV"];


(* ::Input::Initialization:: *)

ExtinctLineages={};
CellNumberTrajectories4M3=Table[
If[Mod[Barcode, 50000]==0, Print[Barcode];];
NumberTraj1=Table[{ReadTrajectories4M3[[Barcode]][[i]][[1]],If[ReadTrajectories4M3[[Barcode]][[i]][[2]]==0,0,N[(PopSize/ReadDepths4M3[[i]][[2]])*(ReadTrajectories4M3[[Barcode]][[i]][[2]])] ]}, {i, 1, Length[ReadDepths4M3]}];

TimePointsBarcodeUnRead=Sort[Map[First,Select[NumberTraj1, #[[2]]==0&]], #1>#2&];
(*Print["TimePointsBarcodeUnRead=", TimePointsBarcodeUnRead];*)
TimePointsBarcodeRead=Complement[TimePoints4M3,TimePointsBarcodeUnRead];
(*Print["TimePointsBarcodeRead=", TimePointsBarcodeRead];*)

If[Length[TimePointsBarcodeUnRead]>0&&Max[TimePointsBarcodeRead]<Min[TimePointsBarcodeUnRead], ExtinctLineage=1; AppendTo[ExtinctLineages, Barcode];, ExtinctLineage=0];

If[ExtinctLineage==1, NumberTraj2=Replace[NumberTraj1, 0-> 0.1, {2}], NumberTraj2=NumberTraj1];

NumberTraj2

, {Barcode, 1,Length[ReadTrajectories4M3]}];

TableForm[CellNumberTrajectories4M3[[1;;10]], TableHeadings -> {Table[StringJoin["Barcode",ToString[i]], {i, 1, 10}],Table[StringJoin["t=",ToString[TimePoints4M3[[i]]]], {i, 1, Length[TimePoints4M3]}]}]
Export["CellNumberTrajectories4M3.csv", CellNumberTrajectories4M3, "CSV"];
Export["ExtinctLineages4M3.csv", ExtinctLineages, "CSV"];


(* ::Section::Initialization::Closed:: *)
(**)
(*2. Noise model*)


(* ::Subsection::Initialization::Closed:: *)
(**)
(*Notes on the model and importing data*)


(* ::Text::Initialization:: *)
(*This distribution comes up in many of the processes because it is the approximation of the inverse laplace transform of *)
(**)
(*M(\[Phi]) = exp(- a \[Phi] / (b \[Phi] +1))*)
(**)
(*which is a form that appears in birth death processes. The function is peaked around a but has an exponential tail rather than a Gaussian one. *)


(* ::Input::Initialization:: *)

LogProb[a_, b_, n_]:=-((Sqrt[n]-Sqrt[a])^2/b)+Log[Sqrt[Sqrt[a]/(4 \[Pi] b n^(3/2))]]


(* ::Text::Initialization:: *)
(*Use this form to estimate: \[Kappa], and xbar between neighboring time points. We condition on a number of reads r1 then fit the measured distribution of reads r2 conditioned on r1 i.e. P(r2|r1). We perform a two parameter fit for (xbar, \[Kappa]) using the above form for the distribution with*)
(**)
(*a=(R2/R1)*r1*Exp[-(t2-t1)*xbar  --- the mean decline*)
(*b=\[Kappa]*)
(**)
(*We do this conditioning on reads between 20 - 30 then take the mean of the fitted parameters. 20-30 is chosen because it is small enouht that we assume neutral, but not too small than integer number effects and fact we do not model things perfectly at small n become important.*)
(**)
(*The outputted mean fit values between the neighboring time points is written in to a file as is the series of plots showing how well the fit works across all r1 from LowerReadLimit (20 normally) to UpperReadLimit (40 usually).*)


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/CarbonData"];

ReadTrajectories2M3=ToExpression[Import["ReadTrajectories2M3.csv"]];
ReadDepths2M3=ToExpression[Import["ReadDepths2M3.csv"]];

ReadTrajectories4M3=ToExpression[Import["ReadTrajectories4M3.csv"]];
ReadDepths4M3=ToExpression[Import["ReadDepths4M3.csv"]];


(* ::Subsection::Initialization::Closed:: *)
(**)
(*Inferring the mean fitness and unadaptive fraction*)


(* ::Input::Initialization:: *)
TimePoints2M3={0,16,32, 40, 48, 64, 72, 80, 88, 96, 104, 112};
TimePoints4M3={0,8,24, 40, 48,56, 64, 72,88, 96};


(* ::Input::Initialization:: *)
LowerReadLimit=20;
UpperReadLimit=40;
GridOfFittedPlots=Table[Table[,{r, LowerReadLimit,UpperReadLimit}], {t, 1, Length[TimePoints2M3]-1}];
AllxBarKappaFits=Table[
t1=TimePoints2M3[[t]];
t2=TimePoints2M3[[t+1]];
R1=ReadDepths2M3[[t]][[2]];
R2=ReadDepths2M3[[t+1]][[2]];

ReadNumbersToUse=Table[i, {i, LowerReadLimit, UpperReadLimit}];

Allr1Fits=Table[
r1=r0;
pm=1;

BarcodesInRange={};

Do[
Measuredr1=ReadTrajectories2M3[[Barcode]][[t]][[2]];
If[Abs[Measuredr1-r1]<= pm, AppendTo[BarcodesInRange, Barcode]];, {Barcode, 1,Length[ReadTrajectories2M3]}];
Setofr2Reads=DeleteCases[Table[ReadTrajectories2M3[[Barcode]][[t+1]][[2]], {Barcode, BarcodesInRange}], 0];
NumberOfBarcodesInRange=Length[BarcodesInRange];
Minr2=Max[{Min[Setofr2Reads], 5}];Maxr2=Max[Setofr2Reads];
r2Distribution=BinCounts[Setofr2Reads, {Minr2, Maxr2, 1}];
(*Print[TableForm[r2Distribution]];*)
Predictedr2Reads[xbar_, \[Kappa]_]:=Table[NumberOfBarcodesInRange*Exp[LogProb[(R2/R1)*r1*Exp[-(t2-t1)*xbar],\[Kappa],r2]], {r2, Minr2, Maxr2-1, 1}];
GoodnessOfFit[xbar_, \[Kappa]_]:=Total[N[(r2Distribution-Predictedr2Reads[xbar, \[Kappa]])^2]];

FittedxBarKappa=NMinimize[{GoodnessOfFit[x, y], 0.00001<x<0.5, 0.5<y<6}, {x,y}];
xBar=FittedxBarKappa[[2]][[1]][[2]];
Kappa=FittedxBarKappa[[2]][[2]][[2]];
Print["xbar=",xBar,"Kappa=",Kappa];


Fittedr2Reads=Table[{r2+0.5,NumberOfBarcodesInRange*Exp[LogProb[(R2/R1)*r1*Exp[-(t2-t1)*xBar],Kappa,r2]]}, {r2, Minr2, 2*Maxr2, 1}];

Fittedr2ReadsNoMean=Table[{r2+0.5,NumberOfBarcodesInRange*Exp[LogProb[(R2/R1)*r1,Kappa,r2]]}, {r2, Minr2, 2*Maxr2, 1}];

LogData=Histogram[Setofr2Reads, {1}, "LogCount",PlotRange ->{{1, 80},All},Frame -> True,AxesOrigin -> {0, -1},FrameLabel -> {Text[Style["Number of reads", 12,FontFamily-> "Helvetica"]],Text[Style["Number of barcodes", 12,FontFamily-> "Helvetica"]]}, FrameTicksStyle->Directive[Black,12, FontFamily-> "Helvetica"], ChartStyle ->Directive[{c4, EdgeForm[None]}]];
LogTheory=ListLogPlot[{Fittedr2Reads, Fittedr2ReadsNoMean},Joined -> True, AxesOrigin -> {0, -1},PlotStyle-> {{Black, Thickness[0.01]},{Black, Dashing[0.03]} }, PlotRange ->{{1, 80},{-1,Log[Max[r2Distribution]]}}];
LOGPLOT=Show[{LogData, LogTheory}, PlotLabel-> Text[Style[StringJoin["\[Kappa]=",ToString[Kappa],"    \!\(\*OverscriptBox[\(x\), \(_\)]\)=",ToString[xBar],"   r=",ToString[Round[r0]]],12, FontFamily-> "Helvetica"]], ImageSize -> 250, PlotRange -> {{1, 80},{-1,Log[Max[r2Distribution]]}}, PlotRangeClipping -> True, AxesOrigin -> {0, -1}];
GridOfFittedPlots[[t]][[r0-LowerReadLimit+1]]=LOGPLOT;

{xBar, Kappa}

, {r0,ReadNumbersToUse}];
Mean[Allr1Fits]

, {t, 1,Length[TimePoints2M3]-1}];

Export["InferredxBarAndKappa2M3.csv", AllxBarKappaFits, "CSV"];
Export["PlotsOfInferredxBarAndKappa2M3.pdf",TableForm[GridOfFittedPlots, TableHeadings -> {Table[StringJoin[ToString[TimePoints2M3[[i]]],"-",ToString[TimePoints2M3[[i+1]]]], {i, 1, Length[TimePoints2M3]-1}],Table[r, {r, LowerReadLimit,UpperReadLimit}]}], "PDF"];
Speak["Finished evaluating"];



(* ::Input::Initialization:: *)

LowerReadLimit=20;
UpperReadLimit=40;
GridOfFittedPlots=Table[Table[,{r, LowerReadLimit,UpperReadLimit}], {t, 1, Length[TimePoints4M3]-1}];
AllxBarKappaFits=Table[
t1=TimePoints4M3[[t]];
t2=TimePoints4M3[[t+1]];
R1=ReadDepths4M3[[t]][[2]];
R2=ReadDepths4M3[[t+1]][[2]];

ReadNumbersToUse=Table[i, {i, LowerReadLimit, UpperReadLimit}];

Allr1Fits=Table[
r1=r0;
pm=1;

BarcodesInRange={};

Do[
Measuredr1=ReadTrajectories4M3[[Barcode]][[t]][[2]];
If[Abs[Measuredr1-r1]<= pm, AppendTo[BarcodesInRange, Barcode]];, {Barcode, 1,Length[ReadTrajectories4M3]}];
Setofr2Reads=DeleteCases[Table[ReadTrajectories4M3[[Barcode]][[t+1]][[2]], {Barcode, BarcodesInRange}], 0];
NumberOfBarcodesInRange=Length[BarcodesInRange];
Minr2=Max[{Min[Setofr2Reads], 5}];Maxr2=Max[Setofr2Reads];
r2Distribution=BinCounts[Setofr2Reads, {Minr2, Maxr2, 1}];
(*Print[TableForm[r2Distribution]];*)
Predictedr2Reads[xbar_, \[Kappa]_]:=Table[NumberOfBarcodesInRange*Exp[LogProb[(R2/R1)*r1*Exp[-(t2-t1)*xbar],\[Kappa],r2]], {r2, Minr2, Maxr2-1, 1}];
GoodnessOfFit[xbar_, \[Kappa]_]:=Total[N[(r2Distribution-Predictedr2Reads[xbar, \[Kappa]])^2]];

FittedxBarKappa=NMinimize[{GoodnessOfFit[x, y], 0.00001<x<0.5, 0.5<y<6}, {x,y}];
xBar=FittedxBarKappa[[2]][[1]][[2]];
Kappa=FittedxBarKappa[[2]][[2]][[2]];
Print["xbar=",xBar,"Kappa=",Kappa];


Fittedr2Reads=Table[{r2+0.5,NumberOfBarcodesInRange*Exp[LogProb[(R2/R1)*r1*Exp[-(t2-t1)*xBar],Kappa,r2]]}, {r2, Minr2, 2*Maxr2, 1}];

Fittedr2ReadsNoMean=Table[{r2+0.5,NumberOfBarcodesInRange*Exp[LogProb[(R2/R1)*r1,Kappa,r2]]}, {r2, Minr2, 2*Maxr2, 1}];

LogData=Histogram[Setofr2Reads, {1}, "LogCount",PlotRange ->{{1, 80},All},Frame -> True,AxesOrigin -> {0, -1},FrameLabel -> {Text[Style["Number of reads", 12,FontFamily-> "Helvetica"]],Text[Style["Number of barcodes", 12,FontFamily-> "Helvetica"]]}, FrameTicksStyle->Directive[Black,12, FontFamily-> "Helvetica"], ChartStyle ->Directive[{c4, EdgeForm[None]}]];
LogTheory=ListLogPlot[{Fittedr2Reads, Fittedr2ReadsNoMean},Joined -> True, AxesOrigin -> {0, -1},PlotStyle-> {{Black, Thickness[0.01]},{Black, Dashing[0.03]} }, PlotRange ->{{1, 80},{-1,Log[Max[r2Distribution]]}}];
LOGPLOT=Show[{LogData, LogTheory}, PlotLabel-> Text[Style[StringJoin["\[Kappa]=",ToString[Kappa],"    \!\(\*OverscriptBox[\(x\), \(_\)]\)=",ToString[xBar],"   r=",ToString[Round[r0]]],12, FontFamily-> "Helvetica"]], ImageSize -> 250, PlotRange -> {{1, 80},{-1,Log[Max[r2Distribution]]}}, PlotRangeClipping -> True, AxesOrigin -> {0, -1}];
GridOfFittedPlots[[t]][[r0-LowerReadLimit+1]]=LOGPLOT;

{xBar, Kappa}

, {r0,ReadNumbersToUse}];
Mean[Allr1Fits]

, {t, 1,Length[TimePoints4M3]-1}];

Export["InferredxBarAndKappa4M3.csv", AllxBarKappaFits, "CSV"];
Export["PlotsOfInferredxBarAndKappa4M3.pdf",TableForm[GridOfFittedPlots, TableHeadings -> {Table[StringJoin[ToString[TimePoints2M3[[i]]],"-",ToString[TimePoints2M3[[i+1]]]], {i, 1, Length[TimePoints2M3]-1}],Table[r, {r, LowerReadLimit,UpperReadLimit}]}], "PDF"];
Speak["Finished evaluating"];


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/EarlyTime"];
xBarKappa2M3=ToExpression[Import["InferredxBarAndKappa2M3.csv"]];
TimePoints2M3={0,16,32, 40, 48, 64, 72, 80, 88};
MeanFitnessList2M3=Table[{TimePoints2M3[[i+1]], xBarKappa2M3[[i]][[1]]},{i, 1, Length[xBarKappa2M3]}];
KappaList2M3=Table[{TimePoints2M3[[i+1]], xBarKappa2M3[[i]][[2]]},{i, 1, Length[xBarKappa2M3]}];

SetDirectory["~/Dropbox/BayesianAnalysis/EarlyTime"];
xBarKappa4M3=ToExpression[Import["InferredxBarAndKappa4M3.csv"]];
TimePoints4M3={0,8,24, 40, 48,56, 64, 72,88};
MeanFitnessList4M3=Table[{TimePoints4M3[[i+1]], xBarKappa4M3[[i]][[1]]},{i, 1, Length[xBarKappa4M3]}];
KappaList4M3=Table[{TimePoints4M3[[i+1]], xBarKappa4M3[[i]][[2]]},{i, 1, Length[xBarKappa4M3]}];


(* ::Input::Initialization:: *)

FitForm[U_, s_]:=Table[{t, U/s (Exp[s t]-1)}, {t, 0, 88,8}]


(* ::Input::Initialization:: *)

ListPlot[{MeanFitnessList2M3,MeanFitnessList4M3, FitForm[10^-5.3, 0.062]}, Joined -> True]


(* ::Input::Initialization:: *)

MeanFitnessList2M3=Delete[Table[{t,10^-5.3/ 0.062 (Exp[ 0.062 t]-1)},{t, TimePoints2M3}],1]
MeanFitnessList4M3=Delete[Table[{t,10^-5.3/ 0.062 (Exp[ 0.062 t]-1)},{t, TimePoints2M3}],1]


(* ::Input::Initialization:: *)




(* ::Section::Initialization:: *)
(**)
(*3. Calculating probabilities*)


(* ::Subsection::Initialization:: *)
(*Inference*)


(* ::Text::Initialization:: *)
(*Importing the data for both 2M3 and 4M3*)


(* ::Input::Initialization:: *)

$HistoryLength=0;


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/2M3"];
ReadTrajectories2M3=ToExpression[Import["ReadTrajectories2M3.csv"]];
ReadDepths2M3=ToExpression[Import["ReadDepths2M3.csv"]];
CellNumberTrajectories2M3=ToExpression[Import["CellNumberTrajectories2M3.csv"]];
ClearMemory


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/4M3"];
ReadTrajectories4M3=ToExpression[Import["ReadTrajectories4M3.csv"]];
ReadDepths4M3=ToExpression[Import["ReadDepths4M3.csv"]];
CellNumberTrajectories4M3=ToExpression[Import["CellNumberTrajectories4M3.csv"]];
ClearMemory


(* ::Input::Initialization:: *)




(* ::Text::Initialization:: *)
(*Parameters that are useful to define*)


(* ::Input::Initialization:: *)
T=8;
PopSize=5*10^8;
NumberOfBarcodes=Length[ReadTrajectories2M3];


(* ::Text::Initialization:: *)
(*Mean fitness list defined so that the mean fitness at 16 corresponds to the mean fitness accounting for reduction in freqs from 8-16.*)


(* ::Input::Initialization:: *)
SetDirectory["~/Dropbox/BayesianAnalysis/2M3"];
xBarKappa2M3=ToExpression[Import["InferredxBArAndKappa2M3.csv"]];
TimePoints2M3={0,16,32, 40, 48, 64, 72, 80, 88, 96, 104, 112};
MeanFitnessList2M3=Table[{TimePoints2M3[[i+1]], xBarKappa2M3[[i]][[1]]},{i, 1, Length[xBarKappa2M3]}];
KappaList2M3=Table[{TimePoints2M3[[i+1]], xBarKappa2M3[[i]][[2]]},{i, 1, Length[xBarKappa2M3]}];

SetDirectory["~/Dropbox/BayesianAnalysis/4M3"];
xBarKappa4M3=ToExpression[Import["InferredxBArAndKappa4M3.csv"]];
TimePoints4M3={0,8,24, 40, 48,56, 64, 72,88, 96};
MeanFitnessList4M3=Table[{TimePoints4M3[[i+1]], xBarKappa4M3[[i]][[1]]},{i, 1, Length[xBarKappa4M3]}];
KappaList4M3=Table[{TimePoints4M3[[i+1]], xBarKappa4M3[[i]][[2]]},{i, 1, Length[xBarKappa4M3]}];


(* ::Text::Initialization:: *)
(*Determine the interpolation formulae for mean fitness and the unadaptive fraction (integral of mean fitness over time)*)


(* ::Input::Initialization:: *)
MeanFitnessBeforeZero2M3=Table[ {8i, 10^-7}, {i, -25, TimePoints2M3[[1]]/T}];
MeanFitnessExtended2M3=Join[ MeanFitnessBeforeZero2M3,MeanFitnessList2M3];
UnAdaptedFraction2M3=Table[{MeanFitnessExtended2M3[[j]][[1]],Exp[-Sum[(MeanFitnessExtended2M3[[t]][[1]]-MeanFitnessExtended2M3[[t-1]][[1]])*MeanFitnessExtended2M3[[t]][[2]], {t, 2,j}]]}, {j, 2, Length[MeanFitnessExtended2M3]}];
InterpolatedMeanFitness2M3=Interpolation[MeanFitnessExtended2M3, InterpolationOrder->1];
InterpolatedUnAdaptiveFraction2M3=Interpolation[UnAdaptedFraction2M3, InterpolationOrder->1];

MeanFitnessBeforeZero4M3=Table[ {8i, 10^-7}, {i, -25, TimePoints4M3[[1]]/T}];
MeanFitnessExtended4M3=Join[ MeanFitnessBeforeZero4M3,MeanFitnessList4M3];
UnAdaptedFraction4M3=Table[{MeanFitnessExtended4M3[[j]][[1]],Exp[-Sum[(MeanFitnessExtended4M3[[t]][[1]]-MeanFitnessExtended4M3[[t-1]][[1]])*MeanFitnessExtended4M3[[t]][[2]], {t, 2,j}]]}, {j, 2, Length[MeanFitnessExtended4M3]}];
InterpolatedMeanFitness4M3=Interpolation[MeanFitnessExtended4M3, InterpolationOrder->1];
InterpolatedUnAdaptiveFraction4M3=Interpolation[UnAdaptedFraction4M3, InterpolationOrder->1];


(* ::Input::Initialization:: *)




(* ::Input::Initialization:: *)

TableForm[{{

Show[{Plot[InterpolatedMeanFitness2M3[x], {x, 0, 112}, PlotLabel -> "InterpolatedMeanFitness", PlotStyle -> Black],Plot[InterpolatedMeanFitness4M3[x], {x, 0, 98}, PlotLabel -> "InterpolatedMeanFitness", PlotStyle -> Directive[{Black, Dashing[0.02]}]]},ImageSize -> 400],

Show[{Plot[InterpolatedUnAdaptiveFraction2M3[x], {x, 0, 112}, PlotLabel -> "InterpolatedUnAdaptedFraction", PlotStyle -> Black],Plot[InterpolatedUnAdaptiveFraction4M3[x], {x, 0, 98}, PlotLabel ->"InterpolatedUnAdaptedFraction", PlotStyle -> Directive[{Black, Dashing[0.02]}]] },ImageSize -> 400],

Show[{ListPlot[KappaList2M3, PlotMarkers -> Markers[0.025, Black, Thickness[0.005]], PlotLabel -> "Kappa for each time"],ListPlot[KappaList4M3, PlotMarkers -> Markers[0.025,White, Thickness[0.005]], PlotLabel -> "Kappa for each time"]}, ImageSize -> 400]}}]


(* ::Text::Initialization:: *)
(*Establishment time as a function of s and \[Tau]. True establishment size is 1/(s-\!\(\*OverscriptBox[\(x\), \(_\)]\)(\[Tau]))*)


(* ::Input::Initialization:: *)

EstablishmentSize2M3[s_, \[Tau]_]:=1/(s-InterpolatedMeanFitness2M3[\[Tau]]) *UnitStep[s-(InterpolatedMeanFitness2M3[\[Tau]]+0.005)]+1/0.005*UnitStep[-(s-(InterpolatedMeanFitness2M3[\[Tau]]+0.005))];

EstablishmentSize4M3[s_, \[Tau]_]:=1/(s-InterpolatedMeanFitness4M3[\[Tau]]) *UnitStep[s-(InterpolatedMeanFitness4M3[\[Tau]]+0.005)]+1/0.005*UnitStep[-(s-(InterpolatedMeanFitness4M3[\[Tau]]+0.005))];


(* ::Text::Initialization:: *)
(*Given this the r2 is fixed by the deterministic growth of the beneficial cells emerging from \[Tau] with fitness s.*)


(* ::Input::Initialization:: *)
Meanr2Givenr12M3[r1_,R1_,R2_,t1_,t2_,s_,\[Tau]_, PopSize_]:=
Module[{TotalNumberOfCells,NumberOfBeneficialCells,NumberOfNeutralCells,TotalNumberOfCellsAtt2,NumberOfReadsr2},
TotalNumberOfCells=PopSize(r1/R1);
NumberOfBeneficialCells=(*Exp[-s \[Tau]](Exp[s t1]-1)*)Exp[s (t1-\[Tau])]*EstablishmentSize2M3[s, \[Tau]]*InterpolatedUnAdaptiveFraction2M3[t1]/InterpolatedUnAdaptiveFraction2M3[\[Tau]];
NumberOfNeutralCells=TotalNumberOfCells-NumberOfBeneficialCells;

If[NumberOfBeneficialCells>=TotalNumberOfCells, NumberOfNeutralCells=0;NumberOfBeneficialCells=TotalNumberOfCells;];

TotalNumberOfCellsAtt2=(NumberOfNeutralCells+Exp[s(t2-t1)]*NumberOfBeneficialCells)*InterpolatedUnAdaptiveFraction2M3[t2]/InterpolatedUnAdaptiveFraction2M3[t1];

NumberOfReadsr2=(TotalNumberOfCellsAtt2/PopSize)*R2;
NumberOfReadsr2];

Meanr2Givenr14M3[r1_,R1_,R2_,t1_,t2_,s_,\[Tau]_, PopSize_]:=
Module[{TotalNumberOfCells,NumberOfBeneficialCells,NumberOfNeutralCells,TotalNumberOfCellsAtt2,NumberOfReadsr2},
TotalNumberOfCells=PopSize(r1/R1);
NumberOfBeneficialCells=(Exp[s (t1-\[Tau])])*EstablishmentSize4M3[s, \[Tau]]*InterpolatedUnAdaptiveFraction4M3[t1]/InterpolatedUnAdaptiveFraction4M3[\[Tau]];
NumberOfNeutralCells=TotalNumberOfCells-NumberOfBeneficialCells;

If[NumberOfBeneficialCells>=TotalNumberOfCells, NumberOfNeutralCells=0;NumberOfBeneficialCells=TotalNumberOfCells;];

TotalNumberOfCellsAtt2=(NumberOfNeutralCells+Exp[s(t2-t1)]*NumberOfBeneficialCells)*InterpolatedUnAdaptiveFraction4M3[t2]/InterpolatedUnAdaptiveFraction4M3[t1];

NumberOfReadsr2=(TotalNumberOfCellsAtt2/PopSize)*R2;
NumberOfReadsr2];


(* ::Text::Initialization:: *)
(*The propagator is then simply the standard error function with \[Kappa] as the noise and the mean set by the above function*)
(**)
(*Dealing with logs is better than the probability itself or numbers get too large and computational overflow occurs.*)


(* ::Input::Initialization:: *)

LogPropagator2M3[r1_,r2_,R1_,R2_,t1_,t2_,PopSize_,\[Kappa]_,s_,\[Tau]_]:=LogProb[
Meanr2Givenr12M3[r1,R1,R2,t1,t2,s,\[Tau], PopSize], \[Kappa], r2]

LogPropagator4M3[r1_,r2_,R1_,R2_,t1_,t2_,PopSize_,\[Kappa]_,s_,\[Tau]_]:=LogProb[
Meanr2Givenr14M3[r1,R1,R2,t1,t2,s,\[Tau], PopSize], \[Kappa], r2]


(* ::Text::Initialization:: *)
(*The prior distriution on (s, \[Tau]) is set by the theory of a constant feeding process. If a constant population of size n0 is feeding mutants at rate \[Mu](s)\[Delta]s into fitness range (s, s+\[Delta]s) then the statistics of the establishment times follow from solving the backwards time generating function equations. *)
(**)
(*The distribution is \[Rho](\[Tau]) = (s d\[Tau] / \[CapitalGamma](n U))*exp(-n U s \[Tau] - e^(-s \[Tau]))*)
(**)
(*To do this we need to have a prior distribution over \[Mu](s)ds. This we have to guess initially and then can improve upon once *)


(* ::Input::Initialization:: *)

DeltaS=0.1;
\[Mu]Trial[s_]:=10^-5 1/DeltaS Exp[-s/DeltaS];
LogPriorDensityTau[\[Tau]_,s_,n_, U_]:=Log[n U s]


(* ::Input::Initialization:: *)
LogProbFirstTimePoint2M3[r_, s_, \[Tau]_, Neutral_]:=Module[{NumberOfBeneficialCellsAtZero,NumberOfCellsInferredAtZero, ExpectedReadsAtZero,TotalNumberOfCellsAtZero, ProbOft0},
If[Neutral==1,ProbOft0=LogProb[r, 2, r],
NumberOfBeneficialCellsAtZero=Exp[-s \[Tau]]*EstablishmentSize2M3[s, \[Tau]]*1/InterpolatedUnAdaptiveFraction2M3[\[Tau]];
NumberOfCellsInferredAtZero=(r*PopSize)/ReadDepths2M3[[1]][[2]];
 TotalNumberOfCellsAtZero=NumberOfCellsInferredAtZero;

If[NumberOfBeneficialCellsAtZero>NumberOfCellsInferredAtZero, TotalNumberOfCellsAtZero=NumberOfBeneficialCellsAtZero];
ExpectedReadsAtZero=TotalNumberOfCellsAtZero/PopSize*ReadDepths2M3[[1]][[2]];
ProbOft0=LogProb[ExpectedReadsAtZero, 2, r]];
ProbOft0]

LogProbFirstTimePoint4M3[r_, s_, \[Tau]_, Neutral_]:=Module[{NumberOfBeneficialCellsAtZero,NumberOfCellsInferredAtZero, ExpectedReadsAtZero,TotalNumberOfCellsAtZero, ProbOft0},
If[Neutral==1,ProbOft0=LogProb[r, 2, r],
NumberOfBeneficialCellsAtZero=Exp[-s \[Tau]]*EstablishmentSize4M3[s, \[Tau]]*1/InterpolatedUnAdaptiveFraction4M3[\[Tau]];
NumberOfCellsInferredAtZero=(r*PopSize)/ReadDepths4M3[[1]][[2]];
 TotalNumberOfCellsAtZero=NumberOfCellsInferredAtZero;

If[NumberOfBeneficialCellsAtZero>NumberOfCellsInferredAtZero, TotalNumberOfCellsAtZero=NumberOfBeneficialCellsAtZero];
ExpectedReadsAtZero=TotalNumberOfCellsAtZero/PopSize*ReadDepths4M3[[1]][[2]];
ProbOft0=LogProb[ExpectedReadsAtZero, 2, r]];
ProbOft0]




(* ::Input::Initialization:: *)
t1=SessionTime[];
TableOfPlots={};
TableOfData={};
Do[
Traj2M3=Replace[ReadTrajectories2M3[[Barcode]], 0-> 1, 2];(*Get rid of zeros since prob does not model them.*)
Traj4M3=Replace[ReadTrajectories4M3[[Barcode]], 0-> 1, 2];(*Get rid of zeros since prob does not model them.*)


(*DETERMINE INITIAL SIZE OF EACH LINEAGE*)

InitialLineageSize2M3=PopSize/ReadDepths2M3[[1]][[2]]*Traj2M3[[1]][[2]];
InitialLineageSize4M3=PopSize/ReadDepths4M3[[1]][[2]]*Traj4M3[[1]][[2]];
(*Print["n0_2M3=",N[InitialLineageSize2M3]];
Print["n0_4M3=",N[InitialLineageSize4M3]];*)

\[Delta]s=0.005;
\[Delta]\[Tau]=1;


(*PROBABILITY OF THE NULL HYPOTHESIS*)

WeightOfT0For2M3=-LogProbFirstTimePoint2M3[Traj2M3[[1]][[2]], 0, 0, 1];
WeightOfLaterTFor2M3=Sum[-LogPropagator2M3[Traj2M3[[t]][[2]],Traj2M3[[t+1]][[2]],ReadDepths2M3[[t]][[2]],ReadDepths2M3[[t+1]][[2]],TimePoints2M3[[t]],TimePoints2M3[[t+1]],10^5 PopSize,KappaList2M3[[t]][[2]],0.00001,0], {t, 1,  Length[TimePoints2M3]-1}];
WeightNotBeneficial2M3=WeightOfT0For2M3+WeightOfLaterTFor2M3;
(*Print["T0 2M3 =",N[WeightOfT0For2M3],"Later T 2M3 = ",WeightOfLaterTFor2M3,"Total Weight2M3 = ",WeightNotBeneficial2M3];*)

WeightOfT0For4M3=-LogProbFirstTimePoint4M3[Traj4M3[[1]][[2]], 0, 0, 1];
WeightOfLaterTFor4M3=Sum[-LogPropagator4M3[Traj4M3[[t]][[2]],Traj4M3[[t+1]][[2]],ReadDepths4M3[[t]][[2]],ReadDepths4M3[[t+1]][[2]],TimePoints4M3[[t]],TimePoints4M3[[t+1]],10^5 PopSize,KappaList4M3[[t]][[2]],0.00001,0], {t, 1,  Length[TimePoints4M3]-1}];
WeightNotBeneficial4M3=WeightOfT0For4M3+WeightOfLaterTFor4M3;
(*Print["T0 4M3 =",N[WeightOfT0For4M3],"Later T 4M3 = ",WeightOfLaterTFor4M3,"Total Weight 4M3 = ",WeightNotBeneficial4M3];*)

If[ Or[WeightNotBeneficial2M3>60,WeightNotBeneficial4M3>50],

t2=SessionTime[];
If[Mod[Barcode, 1]==0,Print[Barcode, "Memory=", MemoryInUse[], "  max",MaxMemoryUsed[]]; Print[t2-t1];];
t1=t2;

(*PROBABILITY OF THE S, T HYPOTHESIS*)

PriorWeight2M3[s_, \[Tau]_]:=-LogPriorDensityTau[\[Tau], s,InitialLineageSize2M3, \[Mu]Trial[s]*\[Delta]s]-Log[\[Delta]\[Tau]];
WeightT02M3[s_, \[Tau]_]:=-LogProbFirstTimePoint2M3[Traj2M3[[1]][[2]], s, \[Tau], 0];
WeightOfLaterT2M3[s_, \[Tau]_]:=Sum[-LogPropagator2M3[Traj2M3[[t]][[2]],Traj2M3[[t+1]][[2]],ReadDepths2M3[[t]][[2]],ReadDepths2M3[[t+1]][[2]],TimePoints2M3[[t]],TimePoints2M3[[t+1]],PopSize,KappaList2M3[[t]][[2]],s,\[Tau]], {t, 1, Length[TimePoints2M3]-1}];

ProbSAndTauGivenData2M3[s_, \[Tau]_]:=PriorWeight2M3[s, \[Tau]]+WeightT02M3[s, \[Tau]]+WeightOfLaterT2M3[s, \[Tau]];


PriorWeight4M3[s_, \[Tau]_]:=-LogPriorDensityTau[\[Tau], s,InitialLineageSize4M3, \[Mu]Trial[s]*\[Delta]s]-Log[\[Delta]\[Tau]];
WeightT04M3[s_, \[Tau]_]:=-LogProbFirstTimePoint4M3[Traj4M3[[1]][[2]], s, \[Tau], 0];
WeightOfLaterT4M3[s_, \[Tau]_]:=Sum[-LogPropagator4M3[Traj4M3[[t]][[2]],Traj4M3[[t+1]][[2]],ReadDepths4M3[[t]][[2]],ReadDepths4M3[[t+1]][[2]],TimePoints4M3[[t]],TimePoints4M3[[t+1]],PopSize,KappaList4M3[[t]][[2]],s,\[Tau]], {t, 1, Length[TimePoints4M3]-1}];

ProbSAndTauGivenData4M3[s_, \[Tau]_]:=PriorWeight4M3[s, \[Tau]]+WeightT04M3[s, \[Tau]]+WeightOfLaterT4M3[s, \[Tau]];


(*FIND MOST LIKELY S AND T FOR BOTH EXPERIMENTS*)

MostProbableSAndTau2M3=NMinimize[{ProbSAndTauGivenData2M3[s, \[Tau]], 0.005<s<0.4 && -150<\[Tau]<90}, {s, \[Tau]}, AccuracyGoal->6,MaxIterations->20,Method -> "NelderMead"(*{"SimulatedAnnealing", "SearchPoints"\[Rule] 100}*)];
MostProbableSAndTau4M3=NMinimize[{ProbSAndTauGivenData4M3[s, \[Tau]], 0.005<s<0.4 && -150<\[Tau]<90}, {s, \[Tau]}, AccuracyGoal->6,MaxIterations->20,Method -> "NelderMead"(*{"SimulatedAnnealing", "SearchPoints"\[Rule] 100}*)];
(*Print[Barcode,"   (s, \[Tau])_2M3 =  ", MostProbableSAndTau2M3,"  (s, \[Tau])_4M3 =  ", MostProbableSAndTau4M3];*)


(*IF PR[BENEFICIAL] > PR[NEUTRAL] THEN ASSIGN BEST FIT VALUES FOR S AND T, ELSE ASSIGN ZERO*)

InferredS2M3=MostProbableSAndTau2M3[[2]][[1]][[2]];
InferredTau2M3=MostProbableSAndTau2M3[[2]][[2]][[2]];
If[WeightNotBeneficial2M3>MostProbableSAndTau2M3[[1]],
InferredS2M3=MostProbableSAndTau2M3[[2]][[1]][[2]];
InferredTau2M3=MostProbableSAndTau2M3[[2]][[2]][[2]];
, 
InferredS2M3=0.00001; 
InferredTau2M3=0;];

InferredS4M3=MostProbableSAndTau4M3[[2]][[1]][[2]];
InferredTau4M3=MostProbableSAndTau4M3[[2]][[2]][[2]];
If[WeightNotBeneficial4M3>MostProbableSAndTau4M3[[1]],
InferredS4M3=MostProbableSAndTau4M3[[2]][[1]][[2]];
InferredTau4M3=MostProbableSAndTau4M3[[2]][[2]][[2]];
, 
InferredS4M3=0.00001; 
InferredTau4M3=0;];




(*LOOK AT CONTRIBUTIONS TO THE WEIGHTS*)
(*
Print["ML T0 2M3 =",N[WeightT02M3[InferredS2M3, InferredTau2M3]],"ML Later T 2M3 = ",WeightOfLaterT2M3[InferredS2M3, InferredTau2M3],"Total Weight ML 2M3 = ",ProbSAndTauGivenData2M3[InferredS2M3, InferredTau2M3]];
Print["ML T0 4M3 =",N[WeightT04M3[InferredS4M3, InferredTau4M3]],"ML Later T 4M3 = ",WeightOfLaterT4M3[InferredS4M3, InferredTau4M3],"Total Weight ML 4M3 = ",ProbSAndTauGivenData4M3[InferredS4M3, InferredTau4M3]];*)




(*ERRORS ON INFERRED VALUES FROM CURVATURE*)

If[InferredS2M3>0.001,
GradientInS2M3=-\!\(
\*SubscriptBox[\(\[PartialD]\), \(s\)]\ \(ProbSAndTauGivenData2M3[s, \ InferredTau2M3]\)\)/. s -> InferredS2M3;
GradientInTau2M3=-\!\(
\*SubscriptBox[\(\[PartialD]\), \(\[Tau]\)]\ \(ProbSAndTauGivenData2M3[InferredS2M3, \ \[Tau]]\)\)/. \[Tau] -> InferredTau2M3;

CurvatureSS2M3=\!\(
\*SubscriptBox[\(\[PartialD]\), \(s\)]\ \(
\*SubscriptBox[\(\[PartialD]\), \(s\)]\ ProbSAndTauGivenData2M3[s, \ InferredTau2M3]\)\)/. s -> InferredS2M3;
CurvatureST2M3=\!\(
\*SubscriptBox[\(\[PartialD]\), \(\[Tau]\)]\ \(
\*SubscriptBox[\(\[PartialD]\), \(s\)]\ ProbSAndTauGivenData2M3[s, \ \[Tau]]\)\)/. {s -> InferredS2M3, \[Tau] -> InferredTau2M3};
CurvatureTT2M3=\!\(
\*SubscriptBox[\(\[PartialD]\), \(\[Tau]\)]\ \(
\*SubscriptBox[\(\[PartialD]\), \(\[Tau]\)]\ ProbSAndTauGivenData2M3[InferredS2M3, \ \[Tau]]\)\)/. \[Tau] -> InferredTau2M3;

Kmatrix2M3={{CurvatureSS2M3,CurvatureST2M3},{CurvatureST2M3,CurvatureTT2M3}};
EigenValsAndVects2M3=Eigensystem[Kmatrix2M3];
PrincipleCurvature2M3=EigenValsAndVects2M3[[1]][[2]];
PrincipleCurvatureDirection2M3=EigenValsAndVects2M3[[2]][[2]];

ErrorInS2M3=Abs[PrincipleCurvatureDirection2M3[[1]]]/Sqrt[PrincipleCurvature2M3];ErrorInTau2M3=Abs[PrincipleCurvatureDirection2M3[[2]]]/Sqrt[PrincipleCurvature2M3];
,ErrorInS2M3=0.01;ErrorInTau2M3=5];


If[InferredS4M3>0.001,
GradientInS4M3=-\!\(
\*SubscriptBox[\(\[PartialD]\), \(s\)]\ \(ProbSAndTauGivenData4M3[s, \ InferredTau4M3]\)\)/. s -> InferredS4M3;
GradientInTau4M3=-\!\(
\*SubscriptBox[\(\[PartialD]\), \(\[Tau]\)]\ \(ProbSAndTauGivenData4M3[InferredS4M3, \ \[Tau]]\)\)/. \[Tau] -> InferredTau4M3;

CurvatureSS4M3=\!\(
\*SubscriptBox[\(\[PartialD]\), \(s\)]\ \(
\*SubscriptBox[\(\[PartialD]\), \(s\)]\ ProbSAndTauGivenData4M3[s, \ InferredTau4M3]\)\)/. s -> InferredS4M3;
CurvatureST4M3=\!\(
\*SubscriptBox[\(\[PartialD]\), \(\[Tau]\)]\ \(
\*SubscriptBox[\(\[PartialD]\), \(s\)]\ ProbSAndTauGivenData4M3[s, \ \[Tau]]\)\)/. {s -> InferredS4M3, \[Tau] -> InferredTau4M3};
CurvatureTT4M3=\!\(
\*SubscriptBox[\(\[PartialD]\), \(\[Tau]\)]\ \(
\*SubscriptBox[\(\[PartialD]\), \(\[Tau]\)]\ ProbSAndTauGivenData4M3[InferredS4M3, \ \[Tau]]\)\)/. \[Tau] -> InferredTau4M3;

Kmatrix4M3={{CurvatureSS4M3,CurvatureST4M3},{CurvatureST4M3,CurvatureTT4M3}};
EigenValsAndVects4M3=Eigensystem[Kmatrix4M3];
PrincipleCurvature4M3=EigenValsAndVects4M3[[1]][[2]];
PrincipleCurvatureDirection4M3=EigenValsAndVects4M3[[2]][[2]];

ErrorInS4M3=Abs[PrincipleCurvatureDirection4M3[[1]]]/Sqrt[PrincipleCurvature4M3];ErrorInTau4M3=Abs[PrincipleCurvatureDirection4M3[[2]]]/Sqrt[PrincipleCurvature4M3];
,ErrorInS4M3=0.01;ErrorInTau4M3=5];





(*PLOT THE RESULTS TO LOOK AT IF ADAPTIVE IN BOTH EXPERIMENTS AND OCCUR EARLY*)


If[Or[MostProbableSAndTau2M3[[1]]<WeightNotBeneficial2M3,MostProbableSAndTau4M3[[1]]<WeightNotBeneficial4M3] && RandomReal[]<1.0,

(*PosteriorPlot2M3=Show[{Plot3D[Exp[-ProbSAndTauGivenData2M3[s, \[Tau]]+ProbSAndTauGivenData2M3[InferredS2M3, InferredTau2M3]], {s,Max[InferredS2M3-5*ErrorInS2M3,0],InferredS2M3+5*ErrorInS2M3}, {\[Tau],InferredTau2M3-5*ErrorInTau2M3,InferredTau2M3+5*ErrorInTau2M3},PlotRange \[Rule] {{0.001,0.3}, {-100,80},All},Boxed \[Rule] False,BoxRatios \[Rule] {1,1,0.6}, MaxRecursion\[Rule]4, PlotPoints\[Rule] 40,ColorFunction\[Rule]Function[{x,y,z},Scheme3[1+Log[10,Clip[z, {10^-10, 1}]]]],ColorFunctionScaling \[Rule] False, Mesh \[Rule] None,AxesLabel \[Rule]  {Text[Style["s", FontFamily\[Rule] "Helvetica", 11]],Text[Style["\[Tau]", FontFamily\[Rule] "Helvetica", 11]]},AxesEdge\[Rule]{{-1,-1},{-1,-1},{-1,-1}}, ImageSize \[Rule] 200]}];

PosteriorPlot4M3=Show[{Plot3D[Exp[-ProbSAndTauGivenData4M3[s, \[Tau]]+ProbSAndTauGivenData4M3[InferredS4M3, InferredTau4M3]], {s,Max[InferredS4M3-5*ErrorInS4M3,0],InferredS4M3+5*ErrorInS4M3}, {\[Tau],InferredTau4M3-5*ErrorInTau4M3,InferredTau4M3+5*ErrorInTau4M3},PlotRange \[Rule] {{0.001,0.3}, {-100,80},All},Boxed \[Rule] False,BoxRatios \[Rule] {1,1,0.6}, MaxRecursion\[Rule]4, PlotPoints\[Rule] 40,ColorFunction\[Rule]Function[{x,y,z},Scheme3[1+Log[10,Clip[z, {10^-10, 1}]]]],ColorFunctionScaling \[Rule] False, Mesh \[Rule] None,AxesLabel \[Rule]  {Text[Style["s", FontFamily\[Rule] "Helvetica", 11]],Text[Style["\[Tau]", FontFamily\[Rule] "Helvetica", 11]]},AxesEdge\[Rule]{{-1,-1},{-1,-1},{-1,-1}}, ImageSize \[Rule] 200]}];*)

PosteriorPlot2M3=Show[{ArrayPlot[Table[-ProbSAndTauGivenData2M3[s, \[Tau]]+ProbSAndTauGivenData2M3[InferredS2M3, InferredTau2M3], {s,0.01,.15, .005}, {\[Tau],-150,50, 1}], FrameLabel->{"s","\[Tau]"}, PlotRange -> {-100, 50}, AspectRatio->0.25]}];

PosteriorPlot4M3=Show[{ArrayPlot[Table[-ProbSAndTauGivenData4M3[s, \[Tau]]+ProbSAndTauGivenData4M3[InferredS4M3, InferredTau4M3], {s,0.01,.15, .005}, {\[Tau],-150,50, 1}], FrameLabel->{"s","\[Tau]"}, PlotRange -> {-100, 50}, AspectRatio->0.25]}];

(*UnAdaptiveTrajectoryPlot2M3=ListLogPlot[Table[{t,CellNumberTrajectories2M3[[Barcode]][[1]][[2]]*InterpolatedUnAdaptiveFraction2M3[t]}, {t,0, 112}], Joined \[Rule] True, PlotStyle \[Rule] Directive[{c0, Thickness[0.007], Dashing[0.02]}], PlotRange \[Rule] {{0,115},{10, 10^6}},ImageSize \[Rule] 300];
UnAdaptiveTrajectoryPlot4M3=ListLogPlot[Table[{t,CellNumberTrajectories4M3[[Barcode]][[1]][[2]]*InterpolatedUnAdaptiveFraction4M3[t]}, {t,0, 96}], Joined \[Rule] True, PlotStyle \[Rule] Directive[{c3, Thickness[0.007], Dashing[0.02]}], PlotRange \[Rule] {{0,115},{10, 10^6}},ImageSize \[Rule] 300];*)

FittedTrajectoryPlot2M3=ListLogPlot[Table[{t,(CellNumberTrajectories2M3[[Barcode]][[1]][[2]]-Exp[InferredS2M3(0-InferredTau2M3)]*(InterpolatedUnAdaptiveFraction2M3[0]/InferredS2M3))*InterpolatedUnAdaptiveFraction2M3[t]+Exp[InferredS2M3(t-InferredTau2M3)]*(InterpolatedUnAdaptiveFraction2M3[t]/InferredS2M3)}, {t,0, 112}], Joined -> True, PlotStyle -> Directive[{c10, Thickness[0.007]}], PlotRange -> {{0,115},{10, 10^6}},ImageSize -> 300];

FittedTrajectoryPlot4M3=ListLogPlot[Table[{t,(CellNumberTrajectories4M3[[Barcode]][[1]][[2]]-Exp[InferredS4M3(0-InferredTau4M3)]*(InterpolatedUnAdaptiveFraction4M3[0]/InferredS4M3))*InterpolatedUnAdaptiveFraction4M3[t]+Exp[InferredS4M3(t-InferredTau4M3)]*(InterpolatedUnAdaptiveFraction4M3[t]/InferredS4M3)}, {t,0, 104}], Joined -> True, PlotStyle -> Directive[{RGBColor[0.77,0.04,0.], Thickness[0.007]}], PlotRange -> {{0,115},{10, 10^6}},ImageSize -> 300];

TrajectoryPlot2M3=ListLogPlot[CellNumberTrajectories2M3[[Barcode]], PlotMarkers ->Markers[0.02, c10, Thickness[0.003]],PlotRange -> {{0,115},{10, 10^6}}, ImageSize -> 300];
TrajectoryPlot4M3=ListLogPlot[CellNumberTrajectories4M3[[Barcode]], PlotMarkers ->Markers[0.02, RGBColor[0.77,0.04,0.], Thickness[0.003]],PlotRange -> {{0,115},{10, 10^6}}, ImageSize -> 300];

CompareTrajectoriesPlot=Show[{(*UnAdaptiveTrajectoryPlot2M3,UnAdaptiveTrajectoryPlot4M3,*)FittedTrajectoryPlot2M3,FittedTrajectoryPlot4M3, TrajectoryPlot2M3,TrajectoryPlot4M3}
, Frame -> True,FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Number of cells", FontFamily-> "Helvetica", 10]]}, FrameTicks->{{Table[{Log[10^k], Superscript[10,k]},{k, 0, 6}], None}, {Table[{8k, 8k}, {k, 0, 14,2}],None}}, ImageSize -> 200];

(*ADDING THE PLOT TO A TABLE*)

AppendTo[TableOfPlots, {Barcode,-MostProbableSAndTau2M3[[1]]+WeightNotBeneficial2M3,-MostProbableSAndTau4M3[[1]]+WeightNotBeneficial4M3,InferredS2M3,ErrorInS2M3, InferredS4M3,ErrorInS4M3, InferredTau2M3,ErrorInTau2M3, InferredTau4M3,ErrorInTau4M3,CompareTrajectoriesPlot, PosteriorPlot2M3,PosteriorPlot4M3}];
];



(*EXPORTING TABLE OF PLOTS IF MORE THAN 2 IN LENGTH*)


NumberInPlotList=Length[TableOfPlots];
Print["NumberInPlotList=", NumberInPlotList];

If[NumberInPlotList==2,
Export[StringJoin["BarcodesPLOTS",ToString[TableOfPlots[[1]][[1]]],"-",ToString[TableOfPlots[[2]][[1]]],".pdf"],TableForm[TableOfPlots, TableHeadings->{{},{"BC","Log(Pr[b]/Pr[n])_2M3","Log(Pr[b]/Pr[n])_4M3","s_2M3","\[Delta]s_2M3","s_4M3","\[Delta]s_4M3","\[Tau]_2M3","\[Delta]\[Tau]_2M3", "\[Tau]_4M3","\[Delta]\[Tau]_4M3", "Trajectory Plots", "Posterior Plot_2M3","Posterior Plot_4M3"}}],"AllowRasterization"-> True, ImageSize -> 300, ImageResolution -> 500];
TableOfPlots={};
];




(*ADDING THE DATA TO A TABLE, EXPORTING IF MORE THAN 500 IN LENGTH*)

AppendTo[TableOfData, {Barcode,-MostProbableSAndTau2M3[[1]]+WeightNotBeneficial2M3,-MostProbableSAndTau4M3[[1]]+WeightNotBeneficial4M3,InferredS2M3,ErrorInS2M3, InferredS4M3,ErrorInS4M3, InferredTau2M3,ErrorInTau2M3, InferredTau4M3,ErrorInTau4M3}];
NumberInDataList=Length[TableOfData];
(*Print["NumberInDataList=", NumberInDataList];*)

If[NumberInDataList==250,
Export[StringJoin["BarcodeDATA",ToString[TableOfData[[1]][[1]]],"-",ToString[TableOfData[[250]][[1]]],".CSV"],TableOfData,"CSV"];
TableOfData={};
];
ClearMemory

];


(*
{Barcode,CompareTrajectoriesPlot,-MostProbableSAndTau2M3[[1]]+WeightNotBeneficial2M3,-MostProbableSAndTau4M3[[1]]+WeightNotBeneficial4M3,InferredS2M3,ErrorInS2M3, InferredS4M3,ErrorInS4M3, InferredTau2M3,ErrorInTau2M3, InferredTau4M3,ErrorInTau4M3(*, PosteriorPlot2M3,PosteriorPlot4M3*)}*)

,{Barcode,{273092}(*ListOfClonesSequenced*)}];
(*Export[FileNames1[[j]],TableForm[PosteriorPlots,TableAlignments\[Rule]{Bottom, Bottom}],"AllowRasterization"\[Rule]True,ImageSize\[Rule]360,ImageResolution\[Rule]600];*)

Speak["Finished, evaluating"];



(* ::Input::Initialization:: *)

TableOfPlots


(* ::Input::Initialization:: *)

Export[StringJoin["BarcodeDATA",ToString[TableOfData[[1]][[1]]],"-",ToString[TableOfData[[Length[TableOfData]]][[1]]],".CSV"],TableOfData,"CSV"];


(* ::Input::Initialization:: *)

Export["ExamplePosteriorSurfacePlot.pdf",TableForm[TableOfPlots, TableHeadings->{{},{"BC","Log(Pr[b]/Pr[n])_2M3","Log(Pr[b]/Pr[n])_4M3","s_2M3","\[Delta]s_2M3","s_4M3","\[Delta]s_4M3","\[Tau]_2M3","\[Delta]\[Tau]_2M3", "\[Tau]_4M3","\[Delta]\[Tau]_4M3", "Trajectory Plots", "Posterior Plot_2M3","Posterior Plot_4M3"}}],"AllowRasterization"-> True, ImageSize -> 300, ImageResolution -> 500];



(* ::Subsection::Initialization::Closed:: *)
(**)
(*Uniting the data into one large table of adaptive lineages*)


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/AllProcessedBarcodeData"];

FilesToImport=FileNames[][[2;;434]]


(* ::Input::Initialization:: *)


Do[TableOfDat[i]=ToExpression[Import[FilesToImport[[i]]]], {i, 1, Length[FilesToImport]}];



(* ::Input::Initialization:: *)

UnitedData=Sort[Join[Flatten[Table[TableOfDat[i], {i, 1, Length[FilesToImport]}],1]], #1[[1]]<#2[[1]]&];


Export["AdaptiveLineagesE1&E2.csv",UnitedData,"CSV"]


(* ::Section::Initialization::Closed:: *)
(**)
(*5. MeanFitness comparison*)


(* ::Text::Initialization:: *)
(*Infer the mean fitness by adding up contribution of each lineage uring teh inferred (s, \[Tau]) while accounting for the eventual decline of some lineages using InterpolatedUnadaptiveFraction.*)


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/AllProcessedBarcodeData_Carbon"];
AllAdaptiveLineages=ToExpression[Import["AdaptiveLineagesE1&E2.csv"]];


(* ::Input::Initialization:: *)

TableForm[AllAdaptiveLineages[[1;;20]], TableHeadings->{{},{"BC","Log(Pr[b]/Pr[n])_1","Log(Pr[b]/Pr[n])_2","s_1","\[Delta]s_1","s_2","\[Delta]s_2","\[Tau]_1","\[Delta]\[Tau]_1", "\[Tau]_2","\[Delta]\[Tau]_2"}}]


(* ::Input::Initialization:: *)
SetDirectory["~/Dropbox/BayesianAnalysis/2M3"];
xBarKappa2M3=ToExpression[Import["InferredxBArAndKappa2M3.csv"]];
TimePoints2M3={0,16,32, 40, 48, 64, 72, 80, 88, 96, 104, 112};
MeanFitnessList2M3=Table[{TimePoints2M3[[i+1]], xBarKappa2M3[[i]][[1]]},{i, 1, Length[xBarKappa2M3]}];
KappaList2M3=Table[{TimePoints2M3[[i+1]], xBarKappa2M3[[i]][[2]]},{i, 1, Length[xBarKappa2M3]}];

SetDirectory["~/Dropbox/BayesianAnalysis/4M3"];
xBarKappa4M3=ToExpression[Import["InferredxBArAndKappa4M3.csv"]];
TimePoints4M3={0,8,24, 40, 48,56, 64, 72,88, 96};
MeanFitnessList4M3=Table[{TimePoints4M3[[i+1]], xBarKappa4M3[[i]][[1]]},{i, 1, Length[xBarKappa4M3]}];
KappaList4M3=Table[{TimePoints4M3[[i+1]], xBarKappa4M3[[i]][[2]]},{i, 1, Length[xBarKappa4M3]}];


(* ::Input::Initialization:: *)

Adaptive2M3=DeleteDuplicates[Select[AllAdaptiveLineages, #[[2]]>3 &]];
Adaptive4M3=Select[AllAdaptiveLineages, #[[3]]>0 &&#[[6]]<0.16 &];


(* ::Input::Initialization:: *)

Length[Adaptive2M3]
Length[Adaptive4M3]


(* ::Text::Initialization:: *)
(*Infer the mean fitness by adding up contribution of each lineage uring teh inferred (s, \[Tau]) while accounting for the eventual decline of some lineages using InterpolatedUnadaptiveFraction.*)


(* ::Input::Initialization:: *)

PopSize=5*10^8;
MeanF[0]=0;
UnadaptedFrac=1;
Do[
MeanF[j]=1/PopSize Sum[Exp[Adaptive2M3[[i]][[4]](j-Adaptive2M3[[i]][[8]])]*UnadaptedFrac,{i, 1, Length[Adaptive2M3]}];
UnadaptedFrac=UnadaptedFrac*Exp[-MeanF[j]];
(*Print[MeanF[j]];*)
j=j+1;, {j, 1, 130}];
MeanFit=Table[{i, MeanF[i]}, {i, 0, 130}];
Export["MeanFitnessFromBeneficials2M3.csv",MeanFit,"CSV"]


(* ::Input::Initialization:: *)
PopSize=5*10^8;
MeanF[0]=0;
UnadaptedFrac=1;
Do[
MeanF[j]=1/PopSize Sum[Exp[Adaptive4M3[[i]][[6]](j-Adaptive4M3[[i]][[10]])]*UnadaptedFrac,{i, 1, Length[Adaptive4M3]}];
UnadaptedFrac=UnadaptedFrac*Exp[-MeanF[j]];
(*Print[MeanF[j]];*)
j=j+1;, {j, 1, 132}];
MeanFit=Table[{i, MeanF[i]}, {i, 0, 132}];
Export["MeanFitnessFromBeneficials4M3.csv",MeanFit,"CSV"]


(* ::Input::Initialization:: *)

MeanFit2M3=ToExpression[Import["MeanFitnessFromBeneficials2M3.csv"]];
MeanFit4M3=ToExpression[Import["MeanFitnessFromBeneficials4M3.csv"]];


(* ::Input::Initialization:: *)

MeanFitnessCompare2M3=Show[{ListPlot[MeanFit2M3, PlotStyle -> c2, Joined -> True, PlotRange -> All],ListPlot[MeanFitnessList2M3, PlotMarkers -> Markers[0.02, c8, Thickness[0.005]]]}, ImageSize -> 170,Frame -> True,FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Mean Fitness", FontFamily-> "Helvetica", 10]]}, PlotRange -> {{0,120},{-0.001, 0.08}}]

MeanFitnessCompare4M3=Show[{ListPlot[MeanFit4M3, PlotStyle -> c2, Joined -> True, PlotRange -> All],ListPlot[MeanFitnessList4M3, PlotMarkers -> Markers[0.02, c8, Thickness[0.005]]]}, ImageSize -> 170,Frame -> True,FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Mean Fitness", FontFamily-> "Helvetica", 10]]}, PlotRange ->{{0,120}, {-0.001, 0.08}}]


(* ::Input::Initialization:: *)

UpperMeanFitness4M3=Delete[Table[{t,10^-5.2/ 0.062 (Exp[ 0.067t]-1)},{t, 1, 120}],1];
LowerMeanFitness4M3=Delete[Table[{t,10^-5.3/ 0.062 (Exp[ 0.058t]-1)},{t, 1,120}],1];
MeanFitnessBounds=ListPlot[{UpperMeanFitness4M3, LowerMeanFitness4M3},PlotStyle -> Gray,Filling -> 1-> {{2}, RGBColor[0.6,0.6,0.6,0.3]},ImageSize -> 250,Frame -> True,Joined -> True,FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Mean Fitness", FontFamily-> "Helvetica", 10]]}, PlotRange ->{{0,120}, {-0.001, 0.08}}];


(* ::Input::Initialization:: *)

MeanFitnessWithBounds=Show[{MeanFitnessBounds,MeanFitnessCompare4M3},PlotRange ->{{0,120}, {-0.001, 0.08}}, AspectRatio->1]
Export["MeanFitnessWithBounds.pdf",MeanFitnessWithBounds,"PDF"]


(* ::Input::Initialization:: *)



(* ::Input::Initialization:: *)
Export["MeanFitnessFromNeutralsE1.csv",MeanFitnessList2M3,"CSV"]
Export["MeanFitnessFromNeutralsE2.csv",MeanFitnessList4M3,"CSV"]


(* ::Text::Initialization:: *)
(*Comparison shows that we must be capturing all that matters to the mean fitness. *)


(* ::Input::Initialization:: *)

Export["MeanFitnessCompare2M3.eps",MeanFitnessCompare2M3,"EPS"]
Export["MeanFitnessCompare4M3.eps",MeanFitnessCompare4M3,"EPS"]


(* ::Section::Initialization::Closed:: *)
(**)
(*6. Plots of trajectories*)


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/2M3"];
AllAdaptiveLineages=ToExpression[Import["AdaptiveLineagesExperiment.csv"]];


(* ::Input::Initialization:: *)

TableForm[AllAdaptiveLineages[[1;;20]], TableHeadings->{{},{"BC","Log(Pr[b]/Pr[n])_1","Log(Pr[b]/Pr[n])_2","s_1","\[Delta]s_1","s_2","\[Delta]s_2","\[Tau]_1","\[Delta]\[Tau]_1", "\[Tau]_2","\[Delta]\[Tau]_2"}}]


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/2M3"];
CellNumberTrajectories2M3=ToExpression[Import["CellNumberTrajectories2M3.csv"]];

SetDirectory["~/Dropbox/BayesianAnalysis/4M3"];
CellNumberTrajectories4M3=ToExpression[Import["CellNumberTrajectories4M3.csv"]];


(* ::Input::Initialization:: *)
PlotOfTrajectory2M3[bc_]:=Module[{bcData, bcProbBene1, bcFitness1},
If[MemberQ[Map[First,AllAdaptiveLineages],bc],
bcData=Select[AllAdaptiveLineages, #[[1]]==bc&][[1]];
bcProbBene1=bcData[[2]];
bcFitness1=bcData[[4]];, bcProbBene1=-19.0;bcFitness1=0.0;];
ListLogPlot[CellNumberTrajectories2M3[[bc]], Joined -> True, PlotStyle -> {{Scheme3[(bcFitness1/0.13)], Thickness[0.0004]}}(*{{Scheme3[Clip[Log[bcProbBene1+20]/5, {0,1}]], Thickness[0.0004(1+bcFitness1)]}}*),ImageSize ->400, PlotRange -> {{0, 168},{0.03, 5*10^8}}, Frame -> True,FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Number of cells", FontFamily-> "Helvetica", 10]]}, FrameTicks->{{Prepend[Table[{10^k, Superscript[10,k]},{k, 0, 8}], {0.1, "Extinct"}], None}, {Table[{16k, 16k}, {k, 0, 25}],None}}]]


PlotOfTrajectory4M3[bc_]:=Module[{bcData, bcProbBene1, bcFitness1},
If[MemberQ[Map[First,AllAdaptiveLineages],bc],
bcData=Select[AllAdaptiveLineages, #[[1]]==bc&][[1]];
bcProbBene1=bcData[[3]];
bcFitness1=bcData[[6]];, bcProbBene1=-19.0;bcFitness1=0.0;];
ListLogPlot[CellNumberTrajectories4M3[[bc]], Joined -> True, PlotStyle -> {{Scheme3[(bcFitness1/0.13)], Thickness[0.0004]}}(*{{Scheme3[Clip[Log[bcProbBene1+20]/5, {0,1}]], Thickness[0.0004(1+bcFitness1)]}}*),ImageSize ->400, PlotRange -> {{0, 168},{0.03, 5*10^8}}, Frame -> True,FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Number of cells", FontFamily-> "Helvetica", 10]]}, FrameTicks->{{Prepend[Table[{10^k, Superscript[10,k]},{k, 0, 8}], {0.1, "Extinct"}], None}, {Table[{16k, 16k}, {k, 0, 25}],None}}]]



(* ::Input::Initialization:: *)
SetDirectory["~/Dropbox/BayesianAnalysis/2M3"];
AdaptiveBarcodes2M3=Map[First,Sort[Select[AllAdaptiveLineages, (#[[2]]>30&& #[[4]]>0.05)&], #1[[4]]>#2[[4]]&]][[1;;1000]];
BorderLineBarcodes2M3=RandomSample[Map[First,Select[AllAdaptiveLineages, 1<#[[2]]<30& ]], 500];
NeutralBarcodes2M3=Map[First,Sort[AllAdaptiveLineages, #1[[4]]<#2[[4]]& ][[1;;500]]];
ExtinctBarcodes2M3=Map[First,ToExpression[Import["ExtinctLineages2M3.csv"]]][[1;;500]];

SetDirectory["~/Dropbox/BayesianAnalysis/4M3"];
AdaptiveBarcodes4M3=Map[First,Sort[Select[AllAdaptiveLineages, (#[[3]]>30&& #[[6]]>0.05)&], #1[[6]]>#2[[6]]&]][[1;;1000]];
BorderLineBarcodes4M3=RandomSample[Map[First,Select[AllAdaptiveLineages, -3<#[[3]]<26& ]], 500];
NeutralBarcodes4M3=Map[First,Sort[AllAdaptiveLineages, #1[[6]]<#2[[6]]& ][[1;;500]]];
ExtinctBarcodes4M3=Map[First,ToExpression[Import["ExtinctLineages4M3.csv"]]][[1;;500]];


(* ::Input::Initialization:: *)
BarcodesToPlot2M3=RandomSample[Join[AdaptiveBarcodes2M3,BorderLineBarcodes2M3, NeutralBarcodes2M3, ExtinctBarcodes2M3], Length[AdaptiveBarcodes2M3]+Length[NeutralBarcodes2M3]+Length[BorderLineBarcodes2M3]+Length[ExtinctBarcodes2M3]];
BarcodesToPlot4M3=RandomSample[Join[AdaptiveBarcodes4M3,BorderLineBarcodes4M3, NeutralBarcodes4M3, ExtinctBarcodes4M3], Length[AdaptiveBarcodes4M3]+Length[NeutralBarcodes4M3]+Length[BorderLineBarcodes4M3]+Length[ExtinctBarcodes4M3]];


(* ::Input::Initialization:: *)

PlotOfTrajectory2M3[1]


(* ::Input::Initialization:: *)

Timing[LineagesSample2M3=Show[Table[PlotOfTrajectory2M3[bc], {bc, BarcodesToPlot2M3}],ImageSize -> 600]]
Timing[LineagesSample4M3=Show[Table[PlotOfTrajectory4M3[bc], {bc, BarcodesToPlot4M3}],ImageSize -> 600]]

Speak["Finished, evaluating"];



(* ::Input::Initialization:: *)


Export["AdaptiveLineagesPlot2M3_FitnessColored.eps",LineagesSample2M3,"EPS"]
Export["AdaptiveLineagesPlot4M3_FitnessColored.eps",LineagesSample4M3,"EPS"]


(* ::Input::Initialization:: *)

Length[Intersection[Map[First,Adaptive2M3], Map[First,Adaptive4M3]]]


(* ::Section::Initialization::Closed:: *)
(**)
(*7. Distribution of s, \[Tau] and \[Mu](s)*)


(* ::Subsection::Initialization::Closed:: *)
(**)
(*Importing data*)


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/AllProcessedBarcodeData_Carbon"];
AllAdaptiveLineages=ToExpression[Import["AdaptiveLineagesE1&E2.csv"]];


(* ::Input::Initialization:: *)

TableForm[AllAdaptiveLineages[[1;;20]], TableHeadings->{{},{"BC","Log(Pr[b]/Pr[n])_1","Log(Pr[b]/Pr[n])_2","s_1","\[Delta]s_1","s_2","\[Delta]s_2","\[Tau]_1","\[Delta]\[Tau]_1", "\[Tau]_2","\[Delta]\[Tau]_2"}}]


(* ::Input::Initialization:: *)
SetDirectory["~/Dropbox/BayesianAnalysis/2M3"];
xBarKappa2M3=ToExpression[Import["InferredxBArAndKappa2M3.csv"]];
TimePoints2M3={0,16,32, 40, 48, 64, 72, 80, 88, 96, 104, 112};
MeanFitnessList2M3=Table[{TimePoints2M3[[i+1]], xBarKappa2M3[[i]][[1]]},{i, 1, Length[xBarKappa2M3]}];
KappaList2M3=Table[{TimePoints2M3[[i+1]], xBarKappa2M3[[i]][[2]]},{i, 1, Length[xBarKappa2M3]}];

SetDirectory["~/Dropbox/BayesianAnalysis/4M3"];
xBarKappa4M3=ToExpression[Import["InferredxBArAndKappa4M3.csv"]];
TimePoints4M3={0,8,24, 40, 48,56, 64, 72,88, 96};
MeanFitnessList4M3=Table[{TimePoints4M3[[i+1]], xBarKappa4M3[[i]][[1]]},{i, 1, Length[xBarKappa4M3]}];
KappaList4M3=Table[{TimePoints4M3[[i+1]], xBarKappa4M3[[i]][[2]]},{i, 1, Length[xBarKappa4M3]}];


(* ::Input::Initialization:: *)

Adaptive2M3=Select[AllAdaptiveLineages, #[[2]]>8 &];
Adaptive4M3=Select[AllAdaptiveLineages, #[[3]]>0.1 &];
Length[Adaptive2M3]
Length[Adaptive4M3]

MeanFit2M3=ToExpression[Import["MeanFitnessFromBeneficials2M3.csv"]];
MeanFit4M3=ToExpression[Import["MeanFitnessFromBeneficials4M3.csv"]];


(* ::Input::Initialization:: *)

AdaptiveOnly2M3=Select[AllAdaptiveLineages, (#[[2]]>1 &&#[[3]]<0&& #[[8]]>-3/#[[4]])&];
AdaptiveOnly4M3=Select[AllAdaptiveLineages, (#[[3]]>1 &&#[[2]]<0 && #[[10]]>-3/#[[6]])&];
Print["ExclusivelyAdatpive2M3=",Length[AdaptiveOnly2M3]];
Print["ExclusivelyAdatpive4M3=",Length[AdaptiveOnly4M3]];



(* ::Subsection::Initialization::Closed:: *)
(**)
(*s-\[Tau] plots *)


(* ::Text::Initialization:: *)
(*First make all the -150 ones not quiite all on -150. Then add the log(c)/s correction.*)


(* ::Input::Initialization:: *)
InferredSvsTau2M3=Map[If[-160<#[[1]]<-140, {#[[1]]-RandomReal[{0,10}](*+(Log[3.5]/#[[2]])*), #[[2]]}, {#[[1]]+Log[3.5]/#[[2]], #[[2]]}]&,Table[{Adaptive2M3[[i]][[8]],Adaptive2M3[[i]][[4]]},{i, 1, Length[Adaptive2M3]}]];
TauList2M3=Map[#[[1]]&, InferredSvsTau2M3];
SList2M3=Map[#[[2]]&, InferredSvsTau2M3];

InferredSvsTau4M3=Map[If[#[[1]]==-150., {#[[1]]-RandomReal[{0,10}]+Log[3.5]/#[[2]], #[[2]]}, {#[[1]]+Log[3.5]/#[[2]], #[[2]]}]&,Table[{Adaptive4M3[[i]][[10]],Adaptive4M3[[i]][[6]]},{i, 1, Length[Adaptive4M3]}]];
TauList4M3=Map[#[[1]]&, InferredSvsTau4M3];
SList4M3=Map[#[[2]]&, InferredSvsTau4M3];


(* ::Input::Initialization:: *)
MarkerSizes2M3[i_,t_]:=0.3*10^-4 Sqrt[Exp[Adaptive2M3[[i]][[4]]*(t-Adaptive2M3[[i]][[8]])]/Adaptive2M3[[i]][[4]]]
MarkerSizes4M3[i_,t_]:=0.3*10^-4 Sqrt[Exp[Adaptive4M3[[i]][[6]]*(t-Adaptive4M3[[i]][[10]])]/Adaptive4M3[[i]][[6]]]


(* ::Text::Initialization:: *)
(*Function the colors lineages based on whether they are likely pre-existing or not. *)
(**)
(*1. If adaptive in other replicate AND establishment time earlier than -(c/s)*)


(* ::Input::Initialization:: *)
STauPlot2M3=Show[RandomSample[Table[ListPlot[{InferredSvsTau2M3[[i]]}, PlotMarkers -> Markers[(*0.003*)MarkerSizes2M3[i, 88], (*Black*)(*Scheme3[Adaptive2M3[[i]][[3]]/(10*Log[10])]*)If[Adaptive2M3[[i]][[3]]>3 , c9,c5],Thickness[0.0003]], PlotRange -> {{-100, 112},{0.0, 0.16}}, Frame -> True,FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Fitness effect, s", FontFamily-> "Helvetica", 10]]}, FrameTicks->{{Automatic, None}, {Table[{16k, 16k}, {k, -10, 7}],None}},AxesOrigin -> {-150,0}], {i, 1, Length[Adaptive2M3]}], Length[Adaptive2M3]]];
STauPlot4M3=Show[RandomSample[Table[ListPlot[{InferredSvsTau4M3[[i]]}, PlotMarkers -> Markers[(*0.003*)MarkerSizes4M3[i, 88], (*Black*)(*Scheme3[Adaptive4M3[[i]][[2]]/(10*Log[10])]*)If[Adaptive4M3[[i]][[2]]>3, c9,c5],Thickness[0.0003]], PlotRange -> {{-100, 112},{0.0, 0.16}}, Frame -> True,FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Fitness effect, s", FontFamily-> "Helvetica", 10]]}, FrameTicks->{{Automatic, None}, {Table[{16k, 16k}, {k, -10, 7}],None}},AxesOrigin -> {-150,0}], {i, 1, Length[Adaptive4M3]}], Length[Adaptive4M3]]];


(* ::Input::Initialization:: *)
n0=1000;

MeanFit2M3=ToExpression[Import["MeanFitnessFromBeneficials2M3.csv"]];
MeanFit4M3=ToExpression[Import["MeanFitnessFromBeneficials4M3.csv"]];

ObservableLimitTime2M3=Table[{MeanFit2M3[[i]][[1]]-(1/MeanFit2M3[[i]][[2]] Log[n0 MeanFit2M3[[i]][[2]]]), MeanFit2M3[[i]][[2]]}, {i, 71, Length[MeanFit2M3]}];
ObservableLimitTime4M3=Table[{MeanFit4M3[[i]][[1]]-(1.07/MeanFit4M3[[i]][[2]] Log[n0 MeanFit4M3[[i]][[2]]]), MeanFit4M3[[i]][[2]]}, {i, 71, Length[MeanFit4M3]}];

IdealLimitTime2M3=Table[{MeanFit2M3[[i]][[1]]-(2/MeanFit2M3[[i]][[2]]), MeanFit2M3[[i]][[2]]}, {i, 71, Length[MeanFit2M3]}];
IdealLimitTime4M3=Table[{MeanFit4M3[[i]][[1]]-(2/MeanFit4M3[[i]][[2]]), MeanFit4M3[[i]][[2]]}, {i, 71, Length[MeanFit4M3]}];


(* ::Input::Initialization:: *)

ObservableTimeLimitPlot2M3=ListPlot[ObservableLimitTime2M3, Frame -> True,PlotRange -> {{-100, 112},{0.0, 0.16}},Joined -> True,AxesOrigin -> {-150,0},FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Fitness effect, s", FontFamily-> "Helvetica", 10]]}, FrameTicks->{{Automatic, None}, {Table[{16k, 16k}, {k, -10, 7}],None}}, PlotStyle -> Directive[{Thickness[0.001],Black, Dashing[0.01]}]];
ObservableTimeLimitPlot4M3=ListPlot[ObservableLimitTime4M3, Frame -> True,PlotRange -> {{-100, 112},{0.0, 0.16}},Joined -> True,AxesOrigin -> {-150,0},FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Fitness effect, s", FontFamily-> "Helvetica", 10]]}, FrameTicks->{{Automatic, None}, {Table[{16k, 16k}, {k, -10, 7}],None}}, PlotStyle -> Directive[{Thickness[0.001],Black, Dashing[0.01]}]];

IdealTimeLimitPlot2M3=ListPlot[IdealLimitTime2M3, Frame -> True,PlotRange -> {{-100, 112},{0.0, 0.16}},Joined -> True,AxesOrigin -> {-150,0},FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Fitness effect, s", FontFamily-> "Helvetica", 10]]}, FrameTicks->{{Automatic, None}, {Table[{16k, 16k}, {k, -10, 7}],None}}, PlotStyle -> Directive[{Thickness[0.001],Black, Dashing[0.03]}]];
IdealTimeLimitPlot4M3=ListPlot[IdealLimitTime4M3, Frame -> True,PlotRange -> {{-100, 112},{0.0, 0.16}},Joined -> True,AxesOrigin -> {-150,0},FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Fitness effect, s", FontFamily-> "Helvetica", 10]]}, FrameTicks->{{Automatic, None}, {Table[{16k, 16k}, {k, -10, 7}],None}}, PlotStyle -> Directive[{Thickness[0.001],Black, Dashing[0.03]}]];


(* ::Input::Initialization:: *)

MeanFitnessPlot2M3=ListPlot[MeanFit2M3, Frame -> True,PlotRange -> {{-100, 112},{0.0, 0.16}},Joined -> True,AxesOrigin -> {-150,0},FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Fitness effect, s", FontFamily-> "Helvetica", 10]]}, FrameTicks->{{Automatic, None}, {Table[{16k, 16k}, {k, -10, 7}],None}}, PlotStyle -> Directive[{Thickness[0.001],Black}]];
MeanFitnessPlot4M3=ListPlot[MeanFit4M3, Frame -> True,PlotRange -> {{-100, 112},{0.0, 0.16}},Joined -> True,AxesOrigin -> {-150,0},FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Fitness effect, s", FontFamily-> "Helvetica", 10]]}, FrameTicks->{{Automatic, None}, {Table[{16k, 16k}, {k, -10, 7}],None}}, PlotStyle -> Directive[{Thickness[0.001],Black}]];

FitnessTauScatterPlot2M3=Show[{STauPlot2M3,MeanFitnessPlot2M3,ObservableTimeLimitPlot2M3,IdealTimeLimitPlot2M3}, ImageSize-> 400];
FitnessTauScatterPlot4M3=Show[{STauPlot4M3,MeanFitnessPlot4M3,ObservableTimeLimitPlot4M3,IdealTimeLimitPlot4M3}, ImageSize-> 400];

Export["FitnessTauScatterPlot2M3_Reverse.eps",FitnessTauScatterPlot2M3,"EPS"];
Export["FitnessTauScatterPlot4M3_Reverse.eps",FitnessTauScatterPlot4M3,"EPS"];


(* ::Input::Initialization:: *)
Speak["Finished, evaluating"]


(* ::Input::Initialization:: *)

DeltaS=0.001;

GatheredByFitness=GatherBy[SList, Round[#/DeltaS]&];

SDistribution=Table[
FitnessOfClass=DeltaS*Round[GatheredByFitness[[k]][[1]]/DeltaS];

NumberInBin=Length[GatheredByFitness[[k]]];

{N[FitnessOfClass],NumberInBin}, {k,1, Length[GatheredByFitness]}];


DeltaTau=1;

GatheredByTau=GatherBy[TauList, Round[#/DeltaTau]&];

TauDistribution=Table[
TauOfClass=DeltaTau*Round[GatheredByTau[[k]][[1]]/DeltaTau];

NumberInBin=Length[GatheredByTau[[k]]];

{N[TauOfClass],NumberInBin}, {k,1, Length[GatheredByTau]}];


(* ::Input::Initialization:: *)

FitnessDistribution=Show[{ListPlot[SDistribution,PlotMarkers -> Markers[0.0002, Black, Thin], PlotRange -> {{0.0, 0.14},All},Filling -> {1-> {Axis, Black}}, ImageSize->200, Ticks -> {None,{{0,0}, {3000,3000}}}, AspectRatio -> 3]/.List[x_,y_]->List[y,x]}, ImageSize -> 70]

EstablishmentTimeDistribution=ListPlot[TauDistribution,PlotMarkers -> Markers[0.0002, Black, Thin], PlotRange -> {{-80,112},All}, AxesOrigin-> {-80,0},Ticks -> {None, {{0,0},{600,600}}},Filling -> {1-> {Axis, Black}}, ImageSize->250, AspectRatio -> 0.1]


(* ::Input::Initialization:: *)

Export["S_dist.eps",FitnessDistribution,"EPS"]
Export["Tau_dist.eps",EstablishmentTimeDistribution,"EPS"]


(* ::Subsection::Initialization::Closed:: *)
(**)
(*Fitness correlation for the pre-existing class*)


(* ::Input::Initialization:: *)

PreExisting=Select[AllAdaptiveLineages, #[[2]]>-Log[0.01]&&#[[3]]>-Log[0.01]&];
Exclusive2M3=Select[AllAdaptiveLineages, #[[2]]>-Log[0.1]&&#[[3]]<-Log[0.1]&];


(* ::Input::Initialization:: *)

Length[PreExisting]
Length[Exclusive2M3]


(* ::Input::Initialization:: *)

PreExistingFitnessesIn2M3=Map[#[[4]]&,PreExisting];
PreExistingTauIn2M3=Map[#[[8]]&,PreExisting];

ExclusiveFitnessesIn2M3=Map[#[[4]]&,Exclusive2M3];
ExclusiveTauIn2M3=Map[#[[8]]&,Exclusive2M3];


(* ::Input::Initialization:: *)

PreExistingS=Histogram[{PreExistingFitnessesIn2M3}, {0.001},PlotRange -> {{0, 0.15},Automatic}, AxesOrigin->{0,0}, ChartStyle->Directive[{EdgeForm[None], c9}],Frame -> True,FrameLabel -> {Text[Style["s", 12,FontFamily-> "Helvetica"]],Text[Style["Number of lineages", 12,FontFamily-> "Helvetica"]]}, FrameTicksStyle->Directive[Black,12, FontFamily-> "Helvetica"]];
ExclusiveS=Histogram[{ExclusiveFitnessesIn2M3}, {0.001},PlotRange -> {{0, 0.15},Automatic}, AxesOrigin->{0,0}, ChartStyle->Directive[{EdgeForm[None], c5}],Frame -> True,FrameLabel -> {Text[Style["s", 12,FontFamily-> "Helvetica"]],Text[Style["Number of lineages", 12,FontFamily-> "Helvetica"]]}, FrameTicksStyle->Directive[Black,12, FontFamily-> "Helvetica"]];

PreExistingS2=Histogram[{PreExistingFitnessesIn2M3}, {0.001},PlotRange -> {{0.05, 0.15},{0,100}}, AxesOrigin->{0,0}, ChartStyle->Directive[{EdgeForm[None], c9}],Frame -> True,FrameLabel -> {Text[Style["s", 12,FontFamily-> "Helvetica"]],Text[Style["Number of lineages", 12,FontFamily-> "Helvetica"]]}, FrameTicksStyle->Directive[Black,12, FontFamily-> "Helvetica"], PlotRangeClipping->True];
ExclusiveS2=Histogram[{ExclusiveFitnessesIn2M3}, {0.001},PlotRange -> {{0.05, 0.15},{0,100}}, AxesOrigin->{0,0}, ChartStyle->Directive[{EdgeForm[None], c5}],Frame -> True,FrameLabel -> {Text[Style["s", 12,FontFamily-> "Helvetica"]],Text[Style["Number of lineages", 12,FontFamily-> "Helvetica"]]}, FrameTicksStyle->Directive[Black,12, FontFamily-> "Helvetica"],PlotRangeClipping->True];

PreExistingTau=Histogram[PreExistingTauIn2M3, {2}, AxesOrigin->{0,0}, ChartStyle->Directive[{EdgeForm[None], c9}],Frame -> True,FrameLabel -> {Text[Style["Establishment time, \[Tau]", 12,FontFamily-> "Helvetica"]],Text[Style["Number of lineages", 12,FontFamily-> "Helvetica"]]}, FrameTicksStyle->Directive[Black,12, FontFamily-> "Helvetica"]];
ExclusiveTau=Histogram[ExclusiveTauIn2M3, {2}, AxesOrigin->{0,0}, ChartStyle->Directive[{EdgeForm[None], c5}],Frame -> True,FrameLabel -> {Text[Style["Establishment time, \[Tau]", 12,FontFamily-> "Helvetica"]],Text[Style["Number of lineages", 12,FontFamily-> "Helvetica"]]}, FrameTicksStyle->Directive[Black,12, FontFamily-> "Helvetica"]];
PreS1=Show[{ExclusiveS,PreExistingS}]
PreS2=Show[{ExclusiveS2,PreExistingS2}, ImageSize -> 200]
Tau=Show[{ExclusiveTau,PreExistingTau}]


(* ::Input::Initialization:: *)

Export["PreExistingS.pdf",PreS1,"PDF"]
Export["PreExistingS2.pdf",PreS2,"PDF"]
Export["PreExistingTau.pdf",Tau,"PDF"]


(* ::Input::Initialization:: *)

sCorrelation=Map[{#[[4]], #[[6]]}&,PreExisting];


(* ::Input::Initialization:: *)

PointWithErrorBars[s1_,s2_,Error1_,Error2_, col_, Thic_, size_]:=Graphics[{Black, Thickness[0.00002], Line[{{s1-Error1,s2},{s1+Error1,s2}}], Line[{{s1,s2-Error2},{s1,s2+Error2}}],EdgeForm[Thickness[Thic]],col,Disk[{s1, s2},Scaled[size]]}];


(* ::Input::Initialization:: *)

sCorrelationPlot=Show[Table[PointWithErrorBars[PreExisting[[i]][[4]],PreExisting[[i]][[6]],0,0, Black(*Scheme3[PreExisting[[i]][[3]]/(10*Log[10])]*), 0.0000002, 0.000001], {i, 1, Length[PreExisting]}],Frame -> True,AxesOrigin -> {0,0},FrameLabel -> {Text[Style["Fitness in E1", FontFamily-> "Helvetica", 10]],Text[Style["Fitness in E2", FontFamily-> "Helvetica", 10]]}, AspectRatio -> 1, PlotRange -> {{0,0.15},{0,0.15}}, ImageSize -> 300];

sLinePlot=Plot[{x}, {x, 0,0.15},Frame -> True,AxesOrigin -> {0,0},FrameLabel -> {Text[Style["Fitness in E1", FontFamily-> "Helvetica", 10]],Text[Style["Fitness in E2", FontFamily-> "Helvetica", 10]]}, AspectRatio -> 1, PlotRange -> {{0,0.15},{0,0.15}},PlotStyle -> Directive[{Black, Thickness[0.003]}], ImageSize -> 300];


(* ::Input::Initialization:: *)

CorrelationPreExisting=Show[{(* sLinePlot,*)sCorrelationPlot}, ImageSize -> 500,Frame -> True,FrameLabel -> {Text[Style["Fitness effect in E1", 18,FontFamily-> "Helvetica"]],Text[Style["Fitness effect in E2", 18,FontFamily-> "Helvetica"]]}, FrameTicksStyle->Directive[Black,18, FontFamily-> "Helvetica"]]


(* ::Input::Initialization:: *)

Export["sCorrelationPlot.eps",CorrelationPreExisting,"EPS"]


(* ::Subsection::Initialization::Closed:: *)
(**)
(*Inferring \[Mu](s) with deterministic approx*)


(* ::Text::Initialization:: *)
(*An initial attempt to infer mu assuming all beneficial lineages are legitamately included*)


(* ::Input::Initialization:: *)
DeltaS=0.001;

ScaleUp=1/DeltaS;(*FITNESS BIN WIDTH =1/THIS NUMBER*)
TauOffSet=0; (*Where is the true zero of time? Really ~30 gens earlier so Tau should have 30 added to it*)

InferredAsBeneficial=Select[Select[PutativeBeneficials, #[[2]]<1&], #[[4]]<100&];

InferredByFitness=GatherBy[InferredAsBeneficial, Round[#[[3]]/DeltaS]&];

EffectiveNo=1000;
EffectiveTime=80;

PopSize=5*10^8;

InferredMutationRatesToClasses=Sort[Table[
FitnessOfClass=DeltaS*Round[InferredByFitness[[k]][[1]][[3]]/DeltaS];

ProbabilityOfClassClimbingAboven0=(Exp[FitnessOfClass*EffectiveTime -(EffectiveNo*FitnessOfClass)/(Exp[FitnessOfClass*EffectiveTime]-1)] 1/(Exp[FitnessOfClass*EffectiveTime]-1));
(*Print[ProbabilityOfClassClimbingAboven0];*)

MutationRateToBin=(*1/ProbabilityOfClassClimbingAboven0*)1/PopSize Sum[Exp[-InferredByFitness[[k]][[i]][[3]]*(InferredByFitness[[k]][[i]][[4]]+TauOffSet)], {i, 1, Length[InferredByFitness[[k]]]}];

{N[FitnessOfClass], MutationRateToBin}, {k,1, Length[InferredByFitness]}]];

CumulativeAboveS=Table[{InferredMutationRatesToClasses[[j]][[1]],Sum[InferredMutationRatesToClasses[[i]][[2]], {i, j, Length[InferredMutationRatesToClasses]}]}, {j, 1, Length[InferredMutationRatesToClasses]}];


(* ::Input::Initialization:: *)

LogMuOfsPlot1=Show[{ListLogPlot[InferredMutationRatesToClasses,PlotMarkers -> Markers[0.0002, Black, Thin], Filling -> {1-> {Axis, Black}}]}, PlotRange -> {{0.035, 0.14},{-32, -12}}, Frame -> True,(*FrameTicks\[Rule] {{Automatic, All},{Automatic, All}},*)FrameLabel -> {Text[Style["Fitness class, s", FontFamily-> "Helvetica", 10]],Text[Style["Mutation rate \[Mu](s)\[Delta]s", FontFamily-> "Helvetica", 10]]}, ImageSize->200]

Cumulative1=Show[{(*PredictedNoThreshold[100, DeltaS, Black],*)ListLogPlot[CumulativeAboveS,PlotMarkers -> Markers[0.0002, Black,Thin],PlotStyle -> c0,(*Joined\[Rule] True,*) Filling -> {1-> {Axis, Black}},FillingStyle -> Thick]}, PlotRange -> {{0.035, 0.14},{-27, -9}}, Frame -> True,FrameTicks-> {{Table[{Log[10^i], Superscript[10,i]}, {i, -11, -4}],None(*Table[{Log[10^i], Superscript[10,i+9]}, {i, -11, -4}]*)},{Automatic, None}},FrameLabel -> {{Text[Style["Mutation rate > s", FontFamily-> "Helvetica", 10]],None(*ext[Style["Target size, bp", FontFamily\[Rule] "Helvetica", 10]]*)}, {Text[Style["Fitness class, s", FontFamily-> "Helvetica", 10]],}}, ImageSize ->200 ]


(* ::Input::Initialization:: *)

Export["CumulativeMuOfSPlot.eps",Cumulative1,"EPS"]
Export["MuOfSPlot.eps",LogMuOfsPlot1,"EPS"]


(* ::Text::Initialization:: *)
(*Total rate comes out to about E-4, consistent with previous analysis.*)


(* ::Subsection::Initialization::Closed:: *)
(**)
(*\[Mu](s) by counting number of adaptive lineages at each s*)


(* ::Text::Initialization:: *)
(*Now infer mu(s) by counting number of lineages at each s*)


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/2M3"];


(* ::Input::Initialization:: *)
MeanFit2M3=ToExpression[Import["MeanFitnessFromBeneficials2M3.csv"]];
MeanFit4M3=ToExpression[Import["MeanFitnessFromBeneficials4M3.csv"]];


(* ::Input::Initialization:: *)

UnAdaptedFraction2M3=Interpolation[Table[{MeanFit2M3[[i]][[1]], Exp[-Sum[MeanFit2M3[[j]][[2]], {j, 1, i}]]}, {i, 1, Length[MeanFit2M3]}], InterpolationOrder -> 1];
UnAdaptedFraction4M3=Interpolation[Table[{MeanFit4M3[[i]][[1]], Exp[-Sum[MeanFit4M3[[j]][[2]], {j, 1, i}]]}, {i, 1, Length[MeanFit4M3]}], InterpolationOrder -> 1];


(* ::Text::Initialization:: *)
(*tstar is the time window a mutation has to occur in if it is to be observed *)
(**)
(*tstar ideal is equivalent window if lienages were all started with a single cell*)


(* ::Input::Initialization:: *)
n0=1000;
tstar2M3=Interpolation[Table[{MeanFit2M3[[i]][[2]],MeanFit2M3[[i]][[1]]- Log[n0*MeanFit2M3[[i]][[2]]]/MeanFit2M3[[i]][[2]]}, {i, 40, Length[MeanFit2M3]}], InterpolationOrder -> 1];
tstar4M3=Interpolation[Table[{MeanFit4M3[[i]][[2]],MeanFit4M3[[i]][[1]]- Log[n0*MeanFit4M3[[i]][[2]]]/MeanFit4M3[[i]][[2]]}, {i, 40, Length[MeanFit4M3]}], InterpolationOrder -> 1];


tstarIdeal2M3=Interpolation[Table[{MeanFit2M3[[i]][[2]],MeanFit2M3[[i]][[1]]-2/MeanFit2M3[[i]][[2]]}, {i, 50, Length[MeanFit2M3]}], InterpolationOrder -> 1];
tstarIdeal4M3=Interpolation[Table[{MeanFit4M3[[i]][[2]],MeanFit4M3[[i]][[1]]-2/MeanFit4M3[[i]][[2]]}, {i, 50, Length[MeanFit4M3]}], InterpolationOrder -> 1];

Plot[{tstar2M3[s], tstarIdeal2M3[s]}, {s, 0.02, 0.12}]

TimeToOccur2M3[s_]:=NIntegrate[Exp[-Exp[-s \[Tau]]]*UnAdaptedFraction2M3[\[Tau]+2/s], {\[Tau],-\[Infinity], tstar2M3[s]}]
TimeToOccur4M3[s_]:=NIntegrate[Exp[-Exp[-s \[Tau]]]*UnAdaptedFraction4M3[\[Tau]+2/s], {\[Tau],-\[Infinity], tstar4M3[s]}]

TimeToOccurIdeal2M3[s_]:=NIntegrate[Exp[-Exp[-s \[Tau]]]*UnAdaptedFraction2M3[\[Tau]+2/s], {\[Tau],-\[Infinity], tstarIdeal2M3[s]}]
TimeToOccurIdeal4M3[s_]:=NIntegrate[Exp[-Exp[-s \[Tau]]]*UnAdaptedFraction4M3[\[Tau]+2/s], {\[Tau],-\[Infinity], tstarIdeal4M3[s]}]


(* ::Text::Initialization:: *)
(*The above plot is tstar as a function of s and tstar ideal as a function of s. The time window is zero for mutations of effect size ~2.5%*)
(**)
(*Time to occur accounts for the fact that the feeding population is reduced as time goes on which will also affect the window in which a mutation can be observed. Integrate this effect up to the time tstar and it gives you the effective number of generations that the feeding population has to form a mutant.*)


(* ::Text::Initialization:: *)
(*1. Number of lineages in a given fitness bin*)
(**)
(*2. The number of mutations expected is (N  \[Mu](s)ds  s t) so can back out \[Mu](s)ds from this. *)
(**)
(*3. MutationRatePlot is simply the number of lineages that are adaptive in fitness bins*)
(**)
(*4. The DFE plots the estimated \[Mu](s) from the above equation*)


(* ::Input::Initialization:: *)

PopSize=5*10^8;
AdaptiveOnly2M3a=Select[AllAdaptiveLineages, (#[[2]]>8 &&#[[3]]<0&& #[[8]]>-1/#[[4]])&];
AdaptiveOnly2M3b=Select[AllAdaptiveLineages, (#[[2]]>8 &&#[[3]]<0&& #[[8]]>-2/#[[4]])&];
AdaptiveOnly2M3c=Select[AllAdaptiveLineages, (#[[2]]>8 &&#[[3]]<0&& #[[8]]>-4/#[[4]])&];


AdaptiveOnly4M3a=Select[AllAdaptiveLineages, (#[[3]]>5 &&#[[2]]<0 && #[[10]]>-1/#[[6]])&];
AdaptiveOnly4M3b=Select[AllAdaptiveLineages, (#[[3]]>5 &&#[[2]]<0 && #[[10]]>-2/#[[6]])&];
AdaptiveOnly4M3c=Select[AllAdaptiveLineages, (#[[3]]>5 &&#[[2]]<0 && #[[10]]>-4/#[[6]])&];
Print["ExclusivelyAdatpive2M3=",Length[AdaptiveOnly2M3a],",",Length[AdaptiveOnly2M3b],",",Length[AdaptiveOnly2M3c]];
Print["ExclusivelyAdatpive4M3=",Length[AdaptiveOnly4M3a],",",Length[AdaptiveOnly4M3b],",",Length[AdaptiveOnly4M3c]];


(* ::Input::Initialization:: *)

ListOfs2M3a=Map[#[[4]]&,AdaptiveOnly2M3a];
ListOfs2M3b=Map[#[[4]]&,AdaptiveOnly2M3b];
ListOfs2M3c=Map[#[[4]]&,AdaptiveOnly2M3c];

ListOfs4M3a=Map[#[[6]]&,AdaptiveOnly4M3a];
ListOfs4M3b=Map[#[[6]]&,AdaptiveOnly4M3b];
ListOfs4M3c=Map[#[[6]]&,AdaptiveOnly4M3c];


(* ::Input::Initialization:: *)
BinWidth=0.0025;

RawBinCounts1a=BinCounts[ListOfs2M3a, {0, 0.16, BinWidth}];
RawBinCounts1b=BinCounts[ListOfs2M3b, {0, 0.16, BinWidth}];
RawBinCounts1c=BinCounts[ListOfs2M3c, {0, 0.16, BinWidth}];

MutationRateByFitnessBin1a=Table[{i*BinWidth, RawBinCounts1a[[i]]/(PopSize*i*BinWidth*TimeToOccur2M3[i*BinWidth])}, {i, 1, Length[RawBinCounts1a]}];
MutationRateByFitnessBin1b=Table[{i*BinWidth, RawBinCounts1b[[i]]/(PopSize*i*BinWidth*TimeToOccur2M3[i*BinWidth])}, {i, 1, Length[RawBinCounts1b]}];
MutationRateByFitnessBin1c=Table[{i*BinWidth, RawBinCounts1c[[i]]/(PopSize*i*BinWidth*TimeToOccur2M3[i*BinWidth])}, {i, 1, Length[RawBinCounts1c]}];


DFE1=ListLogPlot[{MutationRateByFitnessBin1a, MutationRateByFitnessBin1b, MutationRateByFitnessBin1c}, Frame -> {True, True, True, True},Filling -> {1 -> {3}},PlotStyle -> Black, PlotMarkers -> {,Markers[0.015, Black, Thickness[0.005]], },FrameLabel -> {{Text[Style["Beneficial mutation rate", FontFamily-> "Helvetica", 10]],},{Text[Style["Fitness class, s", FontFamily-> "Helvetica", 10]],}},PlotRange -> {{0.04, 0.165},{10^-10, 10^-5}},FrameTicks->{{All,None(*Table[{10^i,Superscript[10,i-4]}, {i, 0, 3}]*)},{All,None}}, ImageSize -> 250]

TotalRate1a=Total[Map[#[[2]]&,MutationRateByFitnessBin1a]];
Print["2M3TotalRate=",TotalRate1a];
TotalRate1b=Total[Map[#[[2]]&,MutationRateByFitnessBin1b]];
Print["2M3TotalRate=",TotalRate1b];
TotalRate1c=Total[Map[#[[2]]&,MutationRateByFitnessBin1c]];
Print["2M3TotalRate=",TotalRate1c];


(* ::Input::Initialization:: *)

BinWidth=0.0025;

RawBinCounts2a=BinCounts[ListOfs4M3a, {0, 0.16, BinWidth}];
RawBinCounts2b=BinCounts[ListOfs4M3b, {0, 0.16, BinWidth}];
RawBinCounts2c=BinCounts[ListOfs4M3c, {0, 0.16, BinWidth}];

MutationRateByFitnessBin2a=Table[{i*BinWidth, RawBinCounts2a[[i]]/(PopSize*i*BinWidth*TimeToOccur4M3[i*BinWidth])}, {i, 1, Length[RawBinCounts2a]}];
MutationRateByFitnessBin2b=Table[{i*BinWidth, RawBinCounts2b[[i]]/(PopSize*i*BinWidth*TimeToOccur4M3[i*BinWidth])}, {i, 1, Length[RawBinCounts2b]}];
MutationRateByFitnessBin2c=Table[{i*BinWidth, RawBinCounts2c[[i]]/(PopSize*i*BinWidth*TimeToOccur4M3[i*BinWidth])}, {i, 1, Length[RawBinCounts2c]}];


DFE2=ListLogPlot[{MutationRateByFitnessBin2a, MutationRateByFitnessBin2b, MutationRateByFitnessBin2c}, Frame -> {True, True, True, True},Filling -> {1 -> {3}},PlotStyle -> Black, PlotMarkers -> {,Markers[0.015, Black, Thickness[0.005]], },FrameLabel -> {{Text[Style["Beneficial mutation rate", FontFamily-> "Helvetica", 10]],},{Text[Style["Fitness class, s", FontFamily-> "Helvetica", 10]],}},PlotRange -> {{0.04, 0.165},{10^-10, 10^-5}},FrameTicks->{{All,None(*Table[{10^i,Superscript[10,i-4]}, {i, 0, 3}]*)},{All,None}}, ImageSize -> 250]

TotalRate2a=Total[Map[#[[2]]&,MutationRateByFitnessBin2a]];
Print["4M3TotalRate=",TotalRate2a];
TotalRate2b=Total[Map[#[[2]]&,MutationRateByFitnessBin2b]];
Print["4M3TotalRate=",TotalRate2b];
TotalRate2c=Total[Map[#[[2]]&,MutationRateByFitnessBin2c]];
Print["4M3TotalRate=",TotalRate2c];


(* ::Input::Initialization:: *)

Export["dfe2M3.pdf",DFE1,"PDF"]
Export["dfe4M3.pdf",DFE2,"PDF"]


(* ::Input::Initialization:: *)

Export["DFE_2.eps",DFE,"EPS"]


(* ::Subsection::Initialization::Closed:: *)
(**)
(*Inferring \[Mu](s) with deterministic approx - EXCLUDING THE "pre-existing" ONES*)


(* ::Text::Initialization:: *)
(*An initial attempt to infer mu assuming all beneficial lineages are legitamately included*)


(* ::Input::Initialization:: *)
DeltaS=0.001;
InferredByFitness2M3=GatherBy[AdaptiveOnly2M3, Round[#[[4]]/DeltaS]&];
InferredByFitness4M3=GatherBy[AdaptiveOnly4M3, Round[#[[6]]/DeltaS]&];
TauOffset=0;

PopSize=5*10^8;

InferredMutationRatesToClasses2M3=Sort[Table[
FitnessOfClass=DeltaS*Round[InferredByFitness2M3[[k]][[1]][[4]]/DeltaS]+0.008;


MutationRateToBin=1/(PopSize DeltaS(1+FitnessOfClass *Log[10^-5*10^12])) Sum[Exp[-InferredByFitness2M3[[k]][[i]][[4]]*(InferredByFitness2M3[[k]][[i]][[8]]+TauOffset)], {i, 1, Length[InferredByFitness2M3[[k]]]}];

MutationBound= MutationRateToBin (*-2Sqrt[MutationRateToBin/(5*10^8)]*);
If[ MutationBound<0, MutationBound=10^-12, ];

{N[FitnessOfClass], MutationBound}, {k,1, Length[InferredByFitness2M3]}]];


InferredMutationRatesToClasses4M3=Sort[Table[
FitnessOfClass=DeltaS*Round[InferredByFitness4M3[[k]][[1]][[6]]/DeltaS];


MutationRateToBin=1/(PopSize DeltaS(1+FitnessOfClass *Log[10^-5*10^12])) Sum[Exp[-InferredByFitness4M3[[k]][[i]][[6]]*(InferredByFitness4M3[[k]][[i]][[10]]+TauOffset)], {i, 1, Length[InferredByFitness4M3[[k]]]}];
MutationBound= MutationRateToBin(* -2Sqrt[MutationRateToBin/(5*10^8)]*);
If[ MutationBound<0, MutationBound=10^-12, ];

{N[FitnessOfClass], MutationBound}, {k,1, Length[InferredByFitness4M3]}]];




CumulativeAboveS2M3=Table[{InferredMutationRatesToClasses2M3[[j]][[1]],DeltaS*Sum[InferredMutationRatesToClasses2M3[[i]][[2]], {i, j, Length[InferredMutationRatesToClasses2M3]}]}, {j, 1, Length[InferredMutationRatesToClasses2M3]}];

CumulativeAboveS4M3=Table[{InferredMutationRatesToClasses4M3[[j]][[1]],DeltaS*Sum[InferredMutationRatesToClasses4M3[[i]][[2]], {i, j, Length[InferredMutationRatesToClasses4M3]}]}, {j, 1, Length[InferredMutationRatesToClasses4M3]}];


(* ::Input::Initialization:: *)

LogMuOfsPlot1=Show[{ListLogPlot[InferredMutationRatesToClasses2M3,PlotMarkers -> Markers[0.0002, Black, Thin], Filling -> {1-> {Axis, Black}}]}, PlotRange -> {{0.03, 0.14},{-26, -5}}, Frame -> True,(*FrameTicks\[Rule] {{Automatic, All},{Automatic, All}},*)FrameLabel -> {Text[Style["Fitness class, s", FontFamily-> "Helvetica", 10]],Text[Style["Density \[Mu](s)", FontFamily-> "Helvetica", 10]]}, ImageSize->180]


LogMuOfsPlot2=Show[{ListLogPlot[InferredMutationRatesToClasses4M3,PlotMarkers -> Markers[0.0002, Black, Thin], Filling -> {1-> {Axis, Black}}]}, PlotRange -> {{0.03, 0.14},{-26, -5}}, Frame -> True,(*FrameTicks\[Rule] {{Automatic, All},{Automatic, All}},*)FrameLabel -> {Text[Style["Fitness class, s", FontFamily-> "Helvetica", 10]],Text[Style["Density \[Mu](s)", FontFamily-> "Helvetica", 10]]}, ImageSize->180]

Cumulative1=Show[{(*PredictedNoThreshold[100, DeltaS, Black],*)ListLogPlot[CumulativeAboveS2M3,PlotMarkers -> Markers[0.0002, Black,Thin],PlotStyle -> c0,(*Joined\[Rule] True,*) Filling -> {1-> {Axis, Black}},FillingStyle -> Thick]}, PlotRange -> {{0.015, 0.14},{-27, -9}}, Frame -> True,FrameTicks-> {{Table[{Log[10^i], Superscript[10,i]}, {i, -11, -4, 2}],None(*Table[{Log[10^i], Superscript[10,i+9]}, {i, -11, -4}]*)},{Automatic, None}},FrameLabel -> {{Text[Style["Mutation rate > s", FontFamily-> "Helvetica", 10]],None(*ext[Style["Target size, bp", FontFamily\[Rule] "Helvetica", 10]]*)}, {Text[Style["Fitness class, s", FontFamily-> "Helvetica", 10]],}}, ImageSize ->250 ]

Cumulative2=Show[{(*PredictedNoThreshold[100, DeltaS, Black],*)ListLogPlot[CumulativeAboveS4M3,PlotMarkers -> Markers[0.0002, Black,Thin],PlotStyle -> c0,(*Joined\[Rule] True,*) Filling -> {1-> {Axis, Black}},FillingStyle -> Thick]}, PlotRange -> {{0.015, 0.14},{-27, -9}}, Frame -> True,FrameTicks-> {{Table[{Log[10^i], Superscript[10,i]}, {i, -11, -4, 2}],None(*Table[{Log[10^i], Superscript[10,i+9]}, {i, -11, -4}]*)},{Automatic, None}},FrameLabel -> {{Text[Style["Mutation rate > s", FontFamily-> "Helvetica", 10]],None(*ext[Style["Target size, bp", FontFamily\[Rule] "Helvetica", 10]]*)}, {Text[Style["Fitness class, s", FontFamily-> "Helvetica", 10]],}}, ImageSize ->250 ]


(* ::Text::Initialization:: *)
(*Plots comparing the inferred distribution to an exponential and estimating the errors on the inferences.*)


(* ::Input::Initialization:: *)
LowerMuE1=InferredMutationRatesToClasses2M3;
LowerMuE2=InferredMutationRatesToClasses4M3;


(* ::Input::Initialization:: *)
LogMuOfsPlot1=Show[{ListLogPlot[{HigherMuE1,LowerMuE1},Joined -> True,PlotStyle -> Black,PlotMarkers -> Markers[0.0002, Black, Thickness[0.01]], Filling -> {1-> {{2}, Black}}]}, PlotRange -> {{0.02, 0.15},{-26, -4}}, Frame -> True,FrameTicks->{{Table[{Log[10^k], Superscript[10,k]},{k, -11, -2, 3}], None}, {Table[{k, k}, {k, 0, 0.15, 0.03}],None}},FrameLabel -> {Text[Style["Fitness class, s", FontFamily-> "Helvetica", 10]],Text[Style["Density \[Mu](s)", FontFamily-> "Helvetica", 10]]}, ImageSize->220,AspectRatio -> 1]
LogMuOfsPlot2=Show[{ListLogPlot[{HigherMuE2,LowerMuE2},Joined -> True,PlotStyle -> Black,PlotMarkers -> Markers[0.0002, Black, Thickness[0.01]], Filling -> {1-> {{2}, Black}}]}, PlotRange -> {{0.02, 0.15},{-26, -4}}, Frame -> True,FrameTicks->{{Table[{Log[10^k], Superscript[10,k]},{k, -11, -2, 3}], None}, {Table[{k, k}, {k, 0, 0.15, 0.03}],None}},FrameLabel -> {Text[Style["Fitness class, s", FontFamily-> "Helvetica", 10]],Text[Style["Density \[Mu](s)", FontFamily-> "Helvetica", 10]]}, ImageSize->220, AspectRatio -> 1]


(* ::Input::Initialization:: *)

DifferenceE1E2=Show[{ListLogPlot[{InferredMutationRatesToClasses2M3,InferredMutationRatesToClasses4M3},Joined -> True,PlotStyle -> Black,PlotMarkers -> Markers[0.0002, Black, Thickness[0.01]], Filling -> {1-> {{2}, Black}}]}, PlotRange -> {{0.03, 0.14},{-26, -4}}, Frame -> True,(*FrameTicks\[Rule] {{Automatic, All},{Automatic, All}},*)FrameLabel -> {Text[Style["Fitness class, s", FontFamily-> "Helvetica", 10]],Text[Style["Density \[Mu](s)", FontFamily-> "Helvetica", 10]]}, ImageSize->250, AspectRatio -> 1];
DifferenceE1E2WithE=Show[{DifferenceE1E2, ExponentialFit2M3}]


(* ::Input::Initialization:: *)

Export["MuOfsDifferenceE1E2.pdf",DifferenceE1E2WithE,"PDF"]


(* ::Input::Initialization:: *)

R=210;
Amp=40;
ExponentialFit2M3=LogPlot[{Amp Exp[-R x],Amp Exp[-R x]+(Amp Exp[-R x]/10^8)^(1/2),Amp Exp[-R x]-(Amp Exp[-R x]/10^8)^(1/2)+((Amp Exp[-R x]/10^8)^(1/2)-Amp Exp[-R x]+10^-12)UnitStep[x-0.105]}, {x, 0.0,0.14}, PlotRange -> {{0.02, 0.14},{Exp[-26], Exp[-3]}}, PlotStyle -> Directive[{RGBColor[0.7,0.7,0.7,0.5], Dashing[0.0]}], Filling -> {2 -> {{3}, RGBColor[0.7,0.7,0.7,0.5]}}, AspectRatio -> 1, Background->None, Frame -> True,FrameTicks->{{Table[{10^k, Superscript[10,k]},{k, -11, -2, 3}], None}, {Table[{k, k}, {k, 0, 0.15, 0.03}],None}},FrameLabel -> {Text[Style["Fitness class, s", FontFamily-> "Helvetica", 10]],Text[Style["Density \[Mu](s)", FontFamily-> "Helvetica", 10]]}, ImageSize->220]


(* ::Input::Initialization:: *)
R=170;
Amp=4;
ExponentialFit4M3=LogPlot[{Amp Exp[-R x],Amp Exp[-R x]+(Amp Exp[-R x]/10^8)^(1/2),Amp Exp[-R x]-(Amp Exp[-R x]/10^8)^(1/2)+((Amp Exp[-R x]/10^8)^(1/2)-Amp Exp[-R x]+10^-12)UnitStep[x-0.116]}, {x, 0.0,0.14}, PlotRange -> {{0.02, 0.14},{Exp[-26], Exp[-3]}}, PlotStyle -> Directive[{RGBColor[0.7,0.7,0.7,0.5], Dashing[0.0]}], Filling -> {2 -> {{3}, RGBColor[0.7,0.7,0.7,0.5]}}, AspectRatio -> 1, Background->None, Frame -> True,FrameTicks->{{Table[{10^k, Superscript[10,k]},{k, -11, -2, 3}], None}, {Table[{k, k}, {k, 0, 0.15, 0.03}],None}},FrameLabel -> {Text[Style["Fitness class, s", FontFamily-> "Helvetica", 10]],Text[Style["Density \[Mu](s)", FontFamily-> "Helvetica", 10]]}, ImageSize->220]


(* ::Input::Initialization:: *)

LogMuOfsPlot1=Show[{ExponentialFit2M3,LogMuOfsPlot1}]
LogMuOfsPlot2=Show[{ExponentialFit4M3,LogMuOfsPlot2}]


(* ::Input::Initialization:: *)

Export["ExponentialFit.pdf",ExponentialFit,"PDF"]


(* ::Input::Initialization:: *)

Export["MuOfSDifference.pdf",DifferencesE1E2,"PDF"]


(* ::Input::Initialization:: *)

Export["MutationSpectrumE1.csv",InferredMutationRatesToClasses2M3,"CSV"]
Export["MutationSpectrumE2.csv",InferredMutationRatesToClasses4M3,"CSV"]


(* ::Input::Initialization:: *)

Export["MuOfsDeterministic2M3.pdf",LogMuOfsPlot1,"PDF"]


(* ::Input::Initialization:: *)

Export["MuOfsDeterministic4M3.pdf",LogMuOfsPlot2,"PDF"]


(* ::Text::Initialization:: *)
(*Total rate comes out to about E-4, consistent with previous analysis.*)


(* ::Subsection::Initialization::Closed:: *)
(**)
(*Errors in s and tau binned by fitness.*)


(* ::Input::Initialization:: *)
TableForm[Adaptive2M3[[1;;20]], TableHeadings->{{},{"BC","Log(Pr[b]/Pr[n])_1","Log(Pr[b]/Pr[n])_2","s_1","\[Delta]s_1","s_2","\[Delta]s_2","\[Tau]_1","\[Delta]\[Tau]_1", "\[Tau]_2","\[Delta]\[Tau]_2"}}]



(* ::Input::Initialization:: *)




(* ::Input::Initialization:: *)
BinWidth=0.01;
BinWidthTau=20;
data=
Table[
ErrorsForS=Map[#[[5]]&,Select[Adaptive2M3, (s<#[[4]]<s+BinWidth &&\[Tau]<#[[8]]<\[Tau]+BinWidthTau && Im[#[[5]]]==0)&]];
ErrorsForTau=Map[#[[9]]&,Select[Adaptive2M3, (s<#[[4]]<s+BinWidth  &&\[Tau]<#[[8]]<\[Tau]+BinWidthTau&& Im[#[[9]]]==0)&]];

If[Length[ErrorsForS]>2,
MeanErrorS=Median[ErrorsForS];
ErrorErrorS=StandardDeviation[ErrorsForS];,MeanErrorS=0.0;
ErrorErrorS=0.0;];

If[Length[ErrorsForTau]>2,
MeanErrorTau=Median[ErrorsForTau];
ErrorErrorTau=StandardDeviation[ErrorsForTau];,MeanErrorTau=0.0;
ErrorErrorTau=0.0;];

{s,\[Tau], MeanErrorS,ErrorErrorS, MeanErrorTau,ErrorErrorTau}, {s, 0.03, 0.14,BinWidth},{\[Tau], -100, 60, BinWidthTau}];



(* ::Input::Initialization:: *)

sData=Map[{#[[2]], #[[1]], #[[3]]}&,Flatten[data, 1]];
TauData=Map[{#[[2]], #[[1]], #[[5]]}&,Flatten[data, 1]];


(* ::Input::Initialization:: *)

sErrorsPlot=Show[Table[ListPlot[{{sData[[i]][[1]],sData[[i]][[2]]}}, PlotMarkers -> Markers[(*0.003*)0.25sData[[i]][[3]]^0.5, Black,Thickness[0.0003]], PlotRange -> {{-105, 70},{0.02, 0.15}}, Frame -> True,FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Fitness effect, s", FontFamily-> "Helvetica", 10]]}, FrameTicks->{{Automatic, None}, {Table[{16k, 16k}, {k, -10, 7}],None}},AxesOrigin -> {-150,0}], {i, 1, Length[sData]}], ImageSize-> 400 ]
TauErrorsPlot=Show[Table[ListPlot[{{TauData[[i]][[1]],TauData[[i]][[2]]}}, PlotMarkers -> Markers[(*0.003*)0.009TauData[[i]][[3]]^0.5, Black,Thickness[0.0003]], PlotRange -> {{-105, 70},{0.02, 0.15}}, Frame -> True,FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Fitness effect, s", FontFamily-> "Helvetica", 10]]}, FrameTicks->{{Automatic, None}, {Table[{16k, 16k}, {k, -10, 7}],None}},AxesOrigin -> {-150,0}], {i, 1, Length[sData]}], ImageSize-> 400 ]


(* ::Input::Initialization:: *)

Export["sErrors.pdf",sErrorsPlot,"PDF"]
Export["TauErrors.pdf",TauErrorsPlot,"PDF"]


(* ::Input::Initialization:: *)




(* ::Input::Initialization:: *)




(* ::Input::Initialization:: *)

ErrorsForTau[[1;;5]]


(* ::Section::Initialization::Closed:: *)
(**)
(*8. Joy-division plot and movie*)


(* ::Subsection::Initialization::Closed:: *)
(**)
(*Importing data*)


(* ::Input::Initialization:: *)
SetDirectory["~/Dropbox/BayesianAnalysis/AllProcessedBarcodeData"];
AllAdaptiveLineages=ToExpression[Import["AdaptiveLineagesE1&E2.csv"]];


(* ::Input::Initialization:: *)

TableForm[AllAdaptiveLineages[[1;;20]], TableHeadings->{{},{"BC","Log(Pr[b]/Pr[n])_1","Log(Pr[b]/Pr[n])_2","s_1","\[Delta]s_1","s_2","\[Delta]s_2","\[Tau]_1","\[Delta]\[Tau]_1", "\[Tau]_2","\[Delta]\[Tau]_2"}}]


(* ::Input::Initialization:: *)

Adaptive2M3=Select[AllAdaptiveLineages, #[[2]]>8 &];
Print["NumberAdaptive2M3=",Length[Adaptive2M3]];
Adaptive4M3=Select[AllAdaptiveLineages, #[[3]]>2 &];
Print["NumberAdaptive4M3=",Length[Adaptive4M3]];

SetDirectory["~/Dropbox/BayesianAnalysis/4M3"];
MeanFit2M3=ToExpression[Import["MeanFitnessFromBeneficials2M3.csv"]];
MeanFit4M3=ToExpression[Import["MeanFitnessFromBeneficials4M3.csv"]];



(* ::Input::Initialization:: *)



(* ::Input::Initialization:: *)

UnAdaptedFraction4M3=Interpolation[Table[{MeanFit4M3[[i]][[1]],Exp[-Sum[MeanFit4M3[[j]][[2]], {j, 1, i}]]}, {i, 1, Length[MeanFit4M3]}]];

UnAdaptedFraction2M3=Interpolation[Table[{MeanFit2M3[[i]][[1]],Exp[-Sum[MeanFit2M3[[j]][[2]], {j, 1, i}]]}, {i, 1, Length[MeanFit2M3]}]];


(* ::Subsection::Initialization::Closed:: *)
(**)
(*Plotting with full diversity inside*)


(* ::Input::Initialization:: *)


AllRectangles={};
CrazyPlots=Table[

Print["Time point=",t];

(*break down adaptive set by barcode fitness into fitness classes*)
FitnessSizes=Sort[Table[
\[Tau]=Adaptive4M3[[BC]][[10]];
s=Adaptive4M3[[BC]][[6]];
Size=Exp[s (t-\[Tau])]/s*UnAdaptedFraction4M3[t];
{BC,s, Size}, {BC, 1, Length[Adaptive4M3]}], #1[[3]]<#2[[3]]&];
Sizes=Table[FitnessSizes[[i]][[3]], {i, 1, Length[FitnessSizes]}];
Fitnesses=Table[FitnessSizes[[i]][[2]], {i, 1, Length[FitnessSizes]}];

MaxFitness=Max[Fitnesses];
Print["Max Fitness",MaxFitness];
FitnessBinWidth=0.002;
BarcodesInEachFitnessBin={};


Do[

NumberOfCellsInFitnessBin=0.0;

BarcodesInRange={};
Do[If[FitnessBin<Fitnesses[[i]]<FitnessBin+FitnessBinWidth &&Sizes[[i]]>0, AppendTo[BarcodesInRange,i]];
, {i, 1, Length[FitnessSizes]}];

 AppendTo[BarcodesInEachFitnessBin, {FitnessBin,BarcodesInRange}];

,{FitnessBin, 0, MaxFitness, FitnessBinWidth}];

(*Put all rectangles in table*)
i=1;
PreviousHeight=0.0;
Do[
BCs=BarcodesInEachFitnessBin[[Bin]][[2]];
(*Print[BCs];*)
Binx=BarcodesInEachFitnessBin[[Bin]][[1]];
(*Print[Binx];*)
StartingN=0.0;
LowX=Binx-FitnessBinWidth;
HighX=Binx+FitnessBinWidth;

Do[

i=i+1;
If[EvenQ[i], Color=White,Color=ColorList[[((t-8)/8)]]];

RectangleCoords={EdgeForm[{Thickness[0.001],Black}],Specularity[0.2,1],Polygon[{{Binx-FitnessBinWidth/2,t,StartingN},{Binx-FitnessBinWidth/2,t,StartingN+Sizes[[Barcode]]},{Binx+FitnessBinWidth/2,t,StartingN+Sizes[[Barcode]]},{Binx+FitnessBinWidth/2,t,StartingN}}]};
StartingN=StartingN+Sizes[[Barcode]];
AppendTo[AllRectangles, Color];
AppendTo[AllRectangles,RectangleCoords];
, {Barcode, BCs}];
AppendTo[AllRectangles,{Black,Thickness[0.005],Line[{{Binx-FitnessBinWidth/2,t,PreviousHeight},{Binx-FitnessBinWidth/2,t,StartingN}, {Binx+FitnessBinWidth/2,t,StartingN}}]}];
(*Print[Bin];*)
PreviousHeight=StartingN;

, {Bin, 1,Length[BarcodesInEachFitnessBin]}];
, {t, 112, 112, 8}];
ClearMemory

FloorCoords=Table[{MeanFit4M3[[i]][[2]], MeanFit4M3[[i]][[1]],0}, {i, 1, Length[MeanFit4M3]}];
AppendTo[FloorCoords, {0,112,0}];
FloorPolygon={EdgeForm[{Thickness[0.001],Black}],Specularity[0.2,1],Polygon[FloorCoords]};
AppendTo[AllRectangles, Gray];
AppendTo[AllRectangles, FloorPolygon];


(* ::Input::Initialization:: *)
CrazyPlot3D=Graphics3D[AllRectangles,BoxRatios-> {3,3,3}, Background -> None, PlotRange -> {{0.0,0.165}, {23, 114}, {0, 20000000}},  ViewPoint-> {0., -25,15},AxesEdge->{{-1,-1},{-1,-1},{-1,1}},Lighting->"Neutral",Axes->True,AxesLabel->{Text[Style["Fitness, s", FontFamily-> "Helvetica", 10]],,Rotate[Text[Style["Cell Number", FontFamily-> "Helvetica", 10]], \[Pi]/2]},Ticks->{Automatic,Table[i, {i, 0, 112, 8}],{{20000000, "2\!\(\*SuperscriptBox[\(x10\), \(7\)]\)"}}}, TicksStyle ->  Directive[FontFamily-> "Helvetica", 10], ImageSize-> 500, Boxed-> False];


(* ::Input::Initialization:: *)

Export["JoyDivisionPlot12.eps",CrazyPlot3D,"EPS"]


(* ::Input::Initialization:: *)

Export["JoyDivisionPlot.pdf",CrazyPlot3D,"AllowRasterization"-> True, ImageResolution-> 300, ImageSize -> 400]


(* ::Subsection::Initialization::Closed:: *)
(**)
(*Plotting only the top *)


(* ::Input::Initialization:: *)


AllRectangles={};
CrazyPlots=Table[

Print["Time point=",t];

(*break down adaptive set by barcode fitness into fitness classes*)
FitnessSizes=Sort[Table[
\[Tau]=Adaptive4M3[[BC]][[10]];
s=Adaptive4M3[[BC]][[6]];
Size=Exp[s (t-\[Tau])]/s*UnAdaptedFraction4M3[t];
{BC,s, Size}, {BC, 1, Length[Adaptive4M3]}], #1[[3]]<#2[[3]]&];
Sizes=Table[FitnessSizes[[i]][[3]], {i, 1, Length[FitnessSizes]}];
Fitnesses=Table[FitnessSizes[[i]][[2]], {i, 1, Length[FitnessSizes]}];

MaxFitness=Max[Fitnesses];
Print["Max Fitness",MaxFitness];
FitnessBinWidth=0.002;
BarcodesInEachFitnessBin={};


Do[

NumberOfCellsInFitnessBin=0.0;

BarcodesInRange={};
Do[If[FitnessBin<Fitnesses[[i]]<FitnessBin+FitnessBinWidth &&Sizes[[i]]>0, AppendTo[BarcodesInRange,i]];
, {i, 1, Length[FitnessSizes]}];

 AppendTo[BarcodesInEachFitnessBin, {FitnessBin,BarcodesInRange}];

,{FitnessBin, 0, MaxFitness, FitnessBinWidth}];

(*Put all rectangles in table*)
i=1;
PreviousHeight=0.0;
Do[
BCs=BarcodesInEachFitnessBin[[Bin]][[2]];
(*Print[BCs];*)
Binx=BarcodesInEachFitnessBin[[Bin]][[1]];
(*Print[Binx];*)
StartingN=0.0;
LowX=Binx-FitnessBinWidth;
HighX=Binx+FitnessBinWidth;

BinHeight=Sum[Sizes[[Barcode]], {Barcode, BCs}];

RectangleCoords={EdgeForm[None],Specularity[0.5,1],Polygon[{{Binx-FitnessBinWidth/2,t,StartingN},{Binx-FitnessBinWidth/2,t,StartingN+BinHeight},{Binx+FitnessBinWidth/2,t,StartingN+BinHeight},{Binx+FitnessBinWidth/2,t,StartingN}}]};
StartingN=StartingN+BinHeight;
AppendTo[AllRectangles,Scheme4[(Bin-15)/30]];
AppendTo[AllRectangles,RectangleCoords];

AppendTo[AllRectangles,{Black,Thickness[0.003],Line[{{Binx-FitnessBinWidth/2,t,PreviousHeight},{Binx-FitnessBinWidth/2,t,StartingN}, {Binx+FitnessBinWidth/2,t,StartingN}}]}];
(*Print[Bin];*)
PreviousHeight=StartingN;

, {Bin, 1,Length[BarcodesInEachFitnessBin]}];
, {t, 0, 132, 4}];
ClearMemory



(* ::Input::Initialization:: *)

CrazyPlot3D=Graphics3D[AllRectangles,BoxRatios-> {3,6.5,4}, Background ->None, PlotRange -> {{0,0.165}, {0, 132}, {0, 50000000}},  ViewPoint-> {5, -45,15},AxesEdge->{{-1,-1},{1,-1},{1,1}},Lighting->"Neutral",Axes->True,AxesLabel->{Text[Style["Fitness, s", FontFamily-> "Helvetica", 10]],Text[Style["   Time", FontFamily-> "Helvetica", 10]],Rotate[Text[Style["Cell Number", FontFamily-> "Helvetica", 10]], \[Pi]/2]},Ticks->{Automatic,Table[i, {i, 0, 132, 24}],{{50000000, "5*\!\(\*SuperscriptBox[\(10\), \(7\)]\)"}}}, TicksStyle ->  Directive[FontFamily-> "Helvetica", 10], ImageSize-> 200, Boxed-> False]


(* ::Input::Initialization:: *)

MeanFitnessCurtain2M3=Join[Table[{MeanFit2M3[[i]][[2]],MeanFit2M3[[i]][[1]],0}, {i, 1, Length[MeanFit2M3]}], Table[{MeanFit2M3[[Length[MeanFit2M3]-i+1]][[2]],MeanFit2M3[[Length[MeanFit2M3]-i+1]][[1]],50000000}, {i, 1, Length[MeanFit2M3]}]];

MeanFitnessCurtain4M3=Join[Table[{MeanFit4M3[[i]][[2]],MeanFit4M3[[i]][[1]],0}, {i, 1, Length[MeanFit4M3]}], Table[{MeanFit4M3[[Length[MeanFit4M3]-i+1]][[2]],MeanFit4M3[[Length[MeanFit4M3]-i+1]][[1]],50000000}, {i, 1, Length[MeanFit4M3]}]];


(* ::Input::Initialization:: *)
Curtain2M3=Graphics3D[{Opacity[0.25],Gray,Polygon[MeanFitnessCurtain2M3]},BoxRatios-> {3,6.5,4}, Background ->None, PlotRange -> {{0.00,0.165}, {0, 132}, {0, 50000000}},  ViewPoint-> {5, -45,15},AxesEdge->{{-1,-1},{1,-1},{1,1}},Lighting->"Neutral",Axes->True,AxesLabel->{Text[Style["Fitness, s", FontFamily-> "Helvetica", 10]],,Rotate[Text[Style["Cell Number", FontFamily-> "Helvetica", 10]], \[Pi]/2]},Ticks->{Automatic,Table[i, {i, 0, 132, 16}],{{100000000, "\!\(\*SuperscriptBox[\(10\), \(8\)]\)"}}}, TicksStyle ->  Directive[FontFamily-> "Helvetica", 10], ImageSize-> 200];

Curtain4M3=Graphics3D[{Opacity[0.25],Gray,Polygon[MeanFitnessCurtain4M3]},BoxRatios-> {3,6.5,4}, Background ->None, PlotRange -> {{0.00,0.165}, {0, 132}, {0, 50000000}},  ViewPoint-> {5, -45,15},AxesEdge->{{-1,-1},{1,-1},{1,1}},Lighting->"Neutral",Axes->True,AxesLabel->{Text[Style["Fitness, s", FontFamily-> "Helvetica", 10]],Text[Style["  Time", FontFamily-> "Helvetica", 10]],Rotate[Text[Style["Cell Number", FontFamily-> "Helvetica", 10]], \[Pi]/2]},Ticks->{Automatic,Table[i, {i, 0, 132, 16}],{{100000000, "\!\(\*SuperscriptBox[\(10\), \(8\)]\)"}}}, TicksStyle ->  Directive[FontFamily-> "Helvetica", 10], ImageSize-> 200];


(* ::Input::Initialization:: *)

CellsByFitness=Show[{CrazyPlot3D, Curtain4M3}]


(* ::Input::Initialization:: *)


Export["JoyDivision4M3.pdf",CellsByFitness,"PDF"]


(* ::Subsection::Initialization::Closed:: *)
(**)
(*The plot in 3D*)


(* ::Text::Initialization:: *)
(*Determine the interpolation formulae for mean fitness and the unadaptive fraction (integral of mean fitness over time)*)


(* ::Input::Initialization:: *)
MeanFitnessBeforeZero=Table[ {8i, 10^-7}, {i, -15, TimePoints[[1]]/T}];
MeanFitnessExtended=Join[ MeanFitnessBeforeZero,MeanFitnessList];
UnAdaptedFraction=Table[{MeanFitnessExtended[[j]][[1]],Exp[-Sum[(MeanFitnessExtended[[t]][[1]]-MeanFitnessExtended[[t-1]][[1]])*MeanFitnessExtended[[t]][[2]], {t, 2,j}]]}, {j, 2, Length[MeanFitnessExtended]}];
InterpolatedMeanFitness=Interpolation[MeanFitnessExtended, InterpolationOrder->1];
InterpolatedUnAdaptiveFraction=Interpolation[UnAdaptedFraction, InterpolationOrder->1];


(* ::Input::Initialization:: *)
BeneficialLineages=Select[PutativeBeneficials, #[[2]]<0.1&];
Print["NumberBeneficial=", Length[BeneficialLineages]];


(* ::Input::Initialization:: *)

ColorOfEachBarcode=Table[{i,RandomChoice[{c1, c2, c3, c4, c5, c6, c7, c8, c9}]}, {i, 1, Length[BeneficialLineages]}];


(* ::Input::Initialization:: *)

AllRectangles={};
CrazyPlots=Table[
Print["Time point=",t];

(*break down adaptive set by barcode fitness into fitness classes*)
FitnessSizes=Sort[Table[
\[Tau]=BeneficialLineages[[BC]][[4]];
s=BeneficialLineages[[BC]][[3]];
Size=Exp[s (t-\[Tau])]/s*InterpolatedUnAdaptiveFraction[t];
{BC,s, Size}, {BC, 1, Length[BeneficialLineages]}], #1[[3]]<#2[[3]]&];
Sizes=Table[FitnessSizes[[i]][[3]], {i, 1, Length[FitnessSizes]}];
Fitnesses=Table[FitnessSizes[[i]][[2]], {i, 1, Length[FitnessSizes]}];

MaxFitness=Max[Fitnesses];
Print["Max Fitness",MaxFitness];
FitnessBinWidth=0.001;
BarcodesInEachFitnessBin={};


Do[

NumberOfCellsInFitnessBin=0.0;

BarcodesInRange={};
Do[If[FitnessBin<Fitnesses[[i]]<FitnessBin+FitnessBinWidth &&Sizes[[i]]>2000, AppendTo[BarcodesInRange,i]];
, {i, 1, Length[FitnessSizes]}];

 AppendTo[BarcodesInEachFitnessBin, {FitnessBin,BarcodesInRange}];

,{FitnessBin, 0, MaxFitness, FitnessBinWidth}];

(*Put all rectangles in table*)

PreviousHeight=0.0;
Do[
BCs=BarcodesInEachFitnessBin[[Bin]][[2]];
(*Print[BCs];*)
Binx=BarcodesInEachFitnessBin[[Bin]][[1]];
(*Print[Binx];*)
StartingN=0.0;
LowX=Binx-FitnessBinWidth;
HighX=Binx+FitnessBinWidth;

Do[
Color=ColorOfEachBarcode[[Position[ColorOfEachBarcode,Barcode][[1]][[1]]]][[2]];
RectangleCoords={EdgeForm[Color],Polygon[{{Binx-FitnessBinWidth/2,t,StartingN},{Binx-FitnessBinWidth/2,t,StartingN+Sizes[[Barcode]]},{Binx+FitnessBinWidth/2,t,StartingN+Sizes[[Barcode]]},{Binx+FitnessBinWidth/2,t,StartingN}}]};
StartingN=StartingN+Sizes[[Barcode]];
AppendTo[AllRectangles, Color];
AppendTo[AllRectangles,RectangleCoords];
, {Barcode, BCs}];
AppendTo[AllRectangles,{Black,Thickness[0.002],Line[{{Binx-FitnessBinWidth/2,t,PreviousHeight},{Binx-FitnessBinWidth/2,t,StartingN}, {Binx+FitnessBinWidth/2,t,StartingN}}]}];
(*Print[Bin];*)
PreviousHeight=StartingN;

, {Bin, 1,Length[BarcodesInEachFitnessBin]}];
, {t, 72, 88, 8}];


(* ::Input::Initialization:: *)
CrazyPlot3D=(*Rasterize[*)Graphics3D[AllRectangles,BoxRatios-> {3,5,3}, Background -> None, Axes -> True, PlotRange -> {{-0.01,0.15}, {31, 112}, {0, 15000000}},  ViewPoint-> {-0.3, -15,9},AxesEdge->{{-1,-1},{-1,-1},{-1,1}},Lighting->"Neutral",Axes->True,AxesLabel->{Text[Style["Fitness, s", FontFamily-> "Helvetica", 10]],,Rotate[Text[Style["Cell Number", FontFamily-> "Helvetica", 10]], \[Pi]/2]},Ticks->{Automatic,Table[i, {i, 0, 112, 16}],{{15000000, ""}}}, TicksStyle ->  Directive[FontFamily-> "Helvetica", 10], ImageSize-> 500, Boxed-> False](*, RasterSize\[Rule]1000]*);


(* ::Input::Initialization:: *)

CrazyPlot3D


(* ::Input::Initialization:: *)

Export[StringJoin["CrazyPlot2M3.pdf"],CrazyPlot3D,"PDF"];


(* ::Subsection::Initialization::Closed:: *)
(**)
(*The movie*)


(* ::Input::Initialization:: *)
SetDirectory["~/Dropbox/BayesianAnalysis/2M3"];
(*ReadTrajectories2M3=ToExpression[Import["ReadTrajectories2M3.csv"]];
ReadDepths2M3=ToExpression[Import["ReadDepths2M3.csv"]];*)
PutativeBeneficials=ToExpression[Import["BayesianAnalysis2M3.csv"]];
CellNumberTrajectories2M3=ToExpression[Import["CellNumberTrajectories2M3.csv"]];



(* ::Input::Initialization:: *)
xBarKappa=ToExpression[Import["InferredxBArAndKappa2M3.csv"]];
TimePoints={0,16,32, 40, 48, 64, 72, 80, 88, 96, 104, 112};
MeanFitnessList=Table[{TimePoints[[i+1]], xBarKappa[[i]][[1]]},{i, 1, Length[xBarKappa]}];
KappaList=Table[{TimePoints[[i]], xBarKappa[[i]][[2]]},{i, 1, Length[xBarKappa]}];


(* ::Text::Initialization:: *)
(*Import the true beneficial mutatation data and true mean fitness. SIMULATION only.*)


(* ::Input::Initialization:: *)
(*TrueFitnessesAndSizeTraj=ToExpression[Import["AllMutationTrajSimStepFunc.csv"]];
TrueMeanFitness=ToExpression[Import["MeanFitnessesSimStepFunc.csv"]];*)



(* ::Text::Initialization:: *)
(*Determine the interpolation formulae for mean fitness and the unadaptive fraction (integral of mean fitness over time)*)


(* ::Input::Initialization:: *)
MeanFitnessBeforeZero=Table[ {8i, 10^-7}, {i, -15, TimePoints[[1]]/T}];
MeanFitnessExtended=Join[ MeanFitnessBeforeZero,MeanFitnessList];
UnAdaptedFraction=Table[{MeanFitnessExtended[[j]][[1]],Exp[-Sum[(MeanFitnessExtended[[t]][[1]]-MeanFitnessExtended[[t-1]][[1]])*MeanFitnessExtended[[t]][[2]], {t, 2,j}]]}, {j, 2, Length[MeanFitnessExtended]}];
InterpolatedMeanFitness=Interpolation[MeanFitnessExtended, InterpolationOrder->1];
InterpolatedUnAdaptiveFraction=Interpolation[UnAdaptedFraction, InterpolationOrder->1];


(* ::Input::Initialization:: *)
BeneficialLineages=Select[PutativeBeneficials, #[[2]]<0.1&];
Print["NumberBeneficial=", Length[BeneficialLineages]];


(* ::Input::Initialization:: *)

ColorOfEachBarcode=Table[{i,RandomChoice[{c1, c2, c3, c4, c5, c6, c7, c8, c9}]}, {i, 1, Length[BeneficialLineages]}];


(* ::Input::Initialization:: *)

AllRectangles={};
CrazyPlots=Table[
Print["Time point=",t];

(*break down adaptive set by barcode fitness into fitness classes*)
FitnessSizes=Sort[Table[
\[Tau]=BeneficialLineages[[BC]][[4]];
s=BeneficialLineages[[BC]][[3]];
Size=Exp[s (t-\[Tau])]/s*InterpolatedUnAdaptiveFraction[t];
{BC,s, Size}, {BC, 1, Length[BeneficialLineages]}], #1[[3]]<#2[[3]]&];
Sizes=Table[FitnessSizes[[i]][[3]], {i, 1, Length[FitnessSizes]}];
Fitnesses=Table[FitnessSizes[[i]][[2]], {i, 1, Length[FitnessSizes]}];

MaxFitness=Max[Fitnesses];
Print["Max Fitness",MaxFitness];
FitnessBinWidth=0.001;
BarcodesInEachFitnessBin={};


Do[

NumberOfCellsInFitnessBin=0.0;

BarcodesInRange={};
Do[If[FitnessBin<Fitnesses[[i]]<FitnessBin+FitnessBinWidth &&Sizes[[i]]>3000, AppendTo[BarcodesInRange,i]];
, {i, 1, Length[FitnessSizes]}];

 AppendTo[BarcodesInEachFitnessBin, {FitnessBin,BarcodesInRange}];

,{FitnessBin, 0, MaxFitness, FitnessBinWidth}];

(*Put all rectangles in table*)

PreviousHeight=0.0;
Do[
BCs=BarcodesInEachFitnessBin[[Bin]][[2]];
(*Print[BCs];*)
Binx=BarcodesInEachFitnessBin[[Bin]][[1]];
(*Print[Binx];*)
StartingN=0.0;
LowX=Binx-FitnessBinWidth;
HighX=Binx+FitnessBinWidth;

Do[
Color=ColorOfEachBarcode[[Position[ColorOfEachBarcode,Barcode][[1]][[1]]]][[2]];
RectangleCoords={EdgeForm[Color],Polygon[{{Binx-FitnessBinWidth/2,StartingN},{Binx-FitnessBinWidth/2,StartingN+Sizes[[Barcode]]},{Binx+FitnessBinWidth/2,StartingN+Sizes[[Barcode]]},{Binx+FitnessBinWidth/2,StartingN}}]};
StartingN=StartingN+Sizes[[Barcode]];
AppendTo[AllRectangles, Color];
AppendTo[AllRectangles,RectangleCoords];
, {Barcode, BCs}];
AppendTo[AllRectangles,{Black,Thickness[0.002],Line[{{Binx-FitnessBinWidth/2,PreviousHeight},{Binx-FitnessBinWidth/2,StartingN}, {Binx+FitnessBinWidth/2,StartingN}}]}];
(*Print[Bin];*)
PreviousHeight=StartingN;

, {Bin, 1,Length[BarcodesInEachFitnessBin]}];
, {t, 90, 90, 1}];


(* ::Input::Initialization:: *)
CrazyPlot=(*Rasterize[*)Graphics[AllRectangles, Background -> None, Axes -> True, PlotRange -> {{-0.01,0.15},{0, 15000000}},AspectRatio-> 1,FrameLabel->{Text[Style["Fitness, s", FontFamily-> "Helvetica", 10]],Text[Style["Cell Number", FontFamily-> "Helvetica", 10]]}, Frame -> True,FrameTicksStyle ->  Directive[FontFamily-> "Helvetica", 10], ImageSize-> {500}](*, RasterSize\[Rule]1000]*)


(* ::Input::Initialization:: *)

Export[StringJoin["CrazyPlot2M3.pdf"],CrazyPlot3D,"PDF"];


(* ::Section::Initialization::Closed:: *)
(**)
(*9. Comparing to fluorescence based fitness measured by Sandeep*)


(* ::Input::Initialization:: *)
SetDirectory["~/Dropbox/BayesianAnalysis/4M3"];
YFPFitnesses=Delete[ToExpression[Import["YFPFitnessMeasurements.csv"]], {1}];


(* ::Input::Initialization:: *)
ListOfYFPBarcodes=Map[First, YFPFitnesses];


(* ::Input::Initialization:: *)

YFPComparison=Table[
AllAdatpiveBarcodesList=Map[First,AllAdaptiveLineages];
BarcodeID=ListOfYFPBarcodes[[i]];
(*Print[BarcodeID];*)
If[MemberQ[AllAdatpiveBarcodesList, BarcodeID],PositionInList=Position[AllAdatpiveBarcodesList,BarcodeID][[1]][[1]];
(*Print[PositionInList];*)
 SequenceFitness=AllAdaptiveLineages[[PositionInList]][[6]];SequenceFitnessError=AllAdaptiveLineages[[PositionInList]][[7]];
ProbabilityAdaptive=AllAdaptiveLineages[[PositionInList]][[3]];, SequenceFitness=0.0; SequenceFitnessError=0.01;ProbabilityAdaptive=0.0;];
{BarcodeID,YFPFitnesses[[i]][[2]],YFPFitnesses[[i]][[3]], SequenceFitness, SequenceFitnessError,ProbabilityAdaptive}, {i, 1, Length[ListOfYFPBarcodes]}];
TableForm[YFPComparison]


(* ::Input::Initialization:: *)

PointWithErrorBars[s1_,s2_,Error1_,Error2_, col_, Thic_, size_]:=Graphics[{Black, Thickness[0.002], Line[{{s1-Error1,s2},{s1+Error1,s2}}], Line[{{s1,s2-Error2},{s1,s2+Error2}}],EdgeForm[Thickness[Thic]],col,Disk[{s1, s2},Scaled[size]]}];

FitnessVerifyYFP=Show[{Table[
If[YFPComparison[[i]][[1]]==8825 || YFPComparison[[i]][[1]]==29375, Color=LightGray, Color=Scheme3[Log[YFPComparison[[i]][[6]]+20]/6]];
PointWithErrorBars[YFPComparison[[i]][[2]],YFPComparison[[i]][[4]],YFPComparison[[i]][[3]],YFPComparison[[i]][[5]], Color,0.001, 0.015], {i, 1, Length[YFPComparison]}], Plot[x, {x, -0.02, 0.15}, PlotStyle -> Directive[{Black, Thickness[0.003]}] ]}, Frame -> True, ImageSize -> 200,FrameLabel -> {Text[Style["Fitness, fluorescence-based assay ", FontFamily-> "Helvetica", 10]],Text[Style["Fitness, sequencing barcodes assay", FontFamily-> "Helvetica", 10]]}, PlotRange -> {{-0.02,0.15},{-0.02,0.15}}]



(* ::Input::Initialization:: *)

Export["FitnessVerify.eps",FitnessVerifyYFP,"EPS"]


(* ::Input::Initialization:: *)
Correlation[DeleteCases[Table[
If[YFPComparison[[i]][[1]]==8825 || YFPComparison[[i]][[1]]==29375|| YFPComparison[[i]][[1]]==39674 ,,YFPComparison[[i]][[2]]], {i, 1, Length[YFPComparison]}], Null], DeleteCases[Table[
If[YFPComparison[[i]][[1]]==8825 || YFPComparison[[i]][[1]]==29375|| YFPComparison[[i]][[1]]==39674 ,,YFPComparison[[i]][[4]]], {i, 1, Length[YFPComparison]}], Null]]^2


(* ::Input::Initialization:: *)

Length[YFPComparison]


(* ::Section::Initialization::Closed:: *)
(**)
(*10. Evolution of n(t)*)


(* ::Text::Initialization:: *)
(*Collect all barcode initially at an effective size of 400-480 cells i.e. *)


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/AllProcessedBarcodeData"];
AllAdaptiveLineages=ToExpression[Import["AdaptiveLineagesE1&E2.csv"]];
SetDirectory["~/Dropbox/BayesianAnalysis/2M3"];
ReadTrajectories2M3=ToExpression[Import["ReadTrajectories2M3.csv"]];
ReadDepths2M3=ToExpression[Import["ReadDepths2M3.csv"]];


(* ::Input::Initialization:: *)

AdaptiveBarcodes=Map[First, AllAdaptiveLineages];


(* ::Input::Initialization:: *)

SmallestDepth=Min[Map[#[[2]]&,ReadDepths2M3]];
MeansForEach=Table[SmallestDepth/ReadDepths2M3[[i]][[2]], {i, 1, Length[ReadDepths2M3]}];
EquallySampledReads2M3=Table[
If[Mod[bc, 1000]==0,Print[bc]];Table[{ReadTrajectories2M3[[bc]][[t]][[1]],If[ReadTrajectories2M3[[bc]][[t]][[2]]>0,RandomVariate[PoissonDistribution[MeansForEach[[t]]*ReadTrajectories2M3[[bc]][[t]][[2]]]],0]}, {t, 1, Length[ReadDepths2M3]}], {bc, 1, Length[ReadTrajectories2M3]}];
Export["EquallySampledReads2M3.csv",EquallySampledReads2M3,"CSV"]


(* ::Input::Initialization:: *)
AdaptiveSet={};
RangeOfInitialCellNumber=Interval[{98,102}];
Do[
If[IntervalMemberQ[ RangeOfInitialCellNumber,EquallySampledReads2M3[[Barcode]][[1]][[2]]], AppendTo[AdaptiveSet, Barcode]]
, {Barcode, 1, Length[EquallySampledReads2M3]}]
Length[AdaptiveSet]


(* ::Input::Initialization:: *)

Distributiont[t_]:=Table[EquallySampledReads2M3[[i]][[t]][[2]], {i, AdaptiveSet}]


(* ::Input::Initialization:: *)


AllRectangles={};
LowerReadLimit=0;
UpperReadLimit=500;
BinWidth=2;

CrazyPlots=Table[
Print["Time point=",t];
BarcodesInEachReadBin={};
AllRectangles={};

Do[
(*Print["ReadBin=", ReadBin];*)

NumberOfCellsInFitnessBin=0.0;

BarcodesInRange={};
Do[If[ReadBin<EquallySampledReads2M3[[i]][[t]][[2]]<ReadBin+BinWidth, AppendTo[BarcodesInRange,i]];
, {i, AdaptiveSet}];
(*Print["InRange=",TableForm[BarcodesInRange]];*)
(*Print["NumberInRange=",Length[BarcodesInRange]];*)

 AppendTo[BarcodesInEachReadBin, {ReadBin,BarcodesInRange}];

,{ReadBin, LowerReadLimit,UpperReadLimit, BinWidth}];

(*Put all rectangles in table*)

PreviousHeight=0.0;
Do[
BCs=BarcodesInEachReadBin[[Bin]][[2]];
(*Print[BCs];*)
Binx=BarcodesInEachReadBin[[Bin]][[1]];
(*Print[Binx];*)
StartingN=0.0;
LowX=Binx-BinWidth;
HighX=Binx+BinWidth;

Do[
If[MemberQ[AdaptiveBarcodes, Barcode], PositionInList=Position[AdaptiveBarcodes, Barcode][[1]][[1]];If[AllAdaptiveLineages[[PositionInList]][[1]]>10, Color=c2, Color=c6],Color=c6];
RectangleCoords={EdgeForm[Color],Polygon[{{Binx-BinWidth/2,StartingN},{Binx-BinWidth/2,StartingN+1},{Binx+BinWidth/2,StartingN+1},{Binx+BinWidth/2,StartingN}}]};
StartingN=StartingN+2;
AppendTo[AllRectangles, Color];
AppendTo[AllRectangles,RectangleCoords];
, {Barcode, BCs}];
(*AppendTo[AllRectangles,{Black,Thickness[0.002],Line[{{Binx-BinWidth/2,PreviousHeight},{Binx-BinWidth/2,StartingN}, {Binx+BinWidth/2,StartingN}}]}];*)
(*Print[Bin];*)
PreviousHeight=StartingN;

, {Bin, 1,Length[BarcodesInEachReadBin]}];




SetOfHeightsAt[t_]:=HistogramList[Table[EquallySampledReads2M3[[i]][[t]][[2]], {i,AdaptiveSet}], {BinWidth}][[2]];
SetOfHeightsForThisTime=SetOfHeightsAt[t];
DistanceOfFit[a_, b_, HeightScale_]:=Total[Table[N[(SetOfHeightsForThisTime[[i]]-HeightScale*Exp[LogProb[a, b, BinWidth*i]])^2], {i, 1, Length[SetOfHeightsForThisTime]}]];
BestFit=FindMinimum[{DistanceOfFit[a, b, cc], 50<a<110, 1<b<10, 1000<cc<30000}, {{a, 100}, {b, 4}, {cc, 10000}}];
Bestcc=BestFit[[2]][[3]][[2]];
Bestb=BestFit[[2]][[2]][[2]];
Besta=BestFit[[2]][[1]][[2]];

BestFitPlot=ListPlot[Table[{BinWidth*i,Bestcc*Exp[LogProb[Besta, Bestb, BinWidth*(i+1)]]}, {i, 1, Length[pp]}], Joined -> True, PlotStyle -> Black];

CrazyPlot=Show[{Graphics[AllRectangles, Background -> None, Axes -> True, PlotRange -> {{0, 500},{0, 300}},AspectRatio-> 1,FrameLabel->{Text[Style["Lineage size", FontFamily-> "Helvetica", 10]],Text[Style["Number of barcodes", FontFamily-> "Helvetica", 10]]}, Frame -> True,FrameTicksStyle ->  Directive[FontFamily-> "Helvetica", 10], ImageSize-> 200], BestFitPlot}];

Export[StringJoin["PorpagatedDistribution",ToString[t],".pdf"],CrazyPlot,"PDF"]


, {t, 8, 8, 1}];



(* ::Input::Initialization:: *)

CrazyPlot


(* ::Section::Initialization::Closed:: *)
(**)
(*11. Simulation of mean fitnes with inferred \[Mu](s)*)


(* ::Text::Initialization:: *)
(*Importing the data*)


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/AllProcessedBarcodeData_Carbon"];
AllAdaptiveLineages=ToExpression[Import["AdaptiveLineagesE1&E2.csv"]];


(* ::Input::Initialization:: *)

TableForm[AllAdaptiveLineages[[1;;20]], TableHeadings->{{},{"BC","Log(Pr[b]/Pr[n])_1","Log(Pr[b]/Pr[n])_2","s_1","\[Delta]s_1","s_2","\[Delta]s_2","\[Tau]_1","\[Delta]\[Tau]_1", "\[Tau]_2","\[Delta]\[Tau]_2"}}]


(* ::Input::Initialization:: *)

PreExistingBarcodes=Select[AllAdaptiveLineages, (#[[2]]>1 &&#[[3]]>1&& #[[8]]<-0/#[[4]]&&#[[10]]<-0/#[[6]])&];
Length[PreExistingBarcodes]
FitnessTimeOfPreX=Map[{#[[4]], #[[8]]}&,PreExistingBarcodes];


(* ::Input::Initialization:: *)

DeltaS=0.01;


(* ::Input::Initialization:: *)

InitialFrequencies[PreExistingFitnessTimes_, DeltaS_]:=Module[{FitnessNewTimes,InitialFreqs,FreqsByFitness,TotalFreqByFitness, FitnessOfClass, FreqOfBin},
FitnessNewTimes=Map[{#[[1]], #[[2]]+RandomVariate[NormalDistribution[0,1]]*1.5/#[[1]]}&,PreExistingFitnessTimes];
InitialFreqs=Map[{#[[1]], Exp[#[[1]]*(-#[[2]])]/(#[[1]]*5*10^8)}&,FitnessNewTimes];
FreqsByFitness=GatherBy[InitialFreqs, Round[#[[1]]/DeltaS]&];
TotalFreqByFitness=Sort[Table[
FitnessOfClass=DeltaS*Round[FreqsByFitness[[k]][[1]][[1]]/DeltaS];

FreqOfBin=Sum[FreqsByFitness[[k]][[i]][[2]], {i, 1, Length[FreqsByFitness[[k]]]}];

{N[FitnessOfClass],FreqOfBin}, {k,1, Length[FreqsByFitness]}]];
TotalFreqByFitness]


(* ::Text::Initialization:: *)
(*The inferred distribution*)


(* ::Input::Initialization:: *)
AdaptiveOnly2M3=Select[AllAdaptiveLineages, (#[[2]]>1 &&#[[3]]<0&& #[[8]]>-4/#[[4]])&];

AdaptiveOnly4M3=Select[AllAdaptiveLineages, (#[[3]]>1 &&#[[2]]<0&& #[[10]]>-4/#[[6]])&];

InferredByFitness2M3=GatherBy[AdaptiveOnly2M3, Round[#[[4]]/DeltaS]&];
InferredByFitness4M3=GatherBy[AdaptiveOnly4M3, Round[#[[6]]/DeltaS]&];
TauOffset=0;

PopSize=5*10^8;

InferredMutationRatesToClasses2M3=Sort[Table[
FitnessOfClass=DeltaS*Round[InferredByFitness2M3[[k]][[1]][[4]]/DeltaS];


MutationRateToBin=1/(PopSize (1+FitnessOfClass *Log[10^-5*10^12])) Sum[Exp[-InferredByFitness2M3[[k]][[i]][[4]]*(InferredByFitness2M3[[k]][[i]][[8]]+TauOffset)], {i, 1, Length[InferredByFitness2M3[[k]]]}];

MutationBound= MutationRateToBin ;

{N[FitnessOfClass], MutationBound}, {k,1, Length[InferredByFitness2M3]}]][[1;;Round[0.20/DeltaS]]];


InferredMutationRatesToClasses4M3=Sort[Table[
FitnessOfClass=DeltaS*Round[InferredByFitness4M3[[k]][[1]][[6]]/DeltaS];


MutationRateToBin=1/(PopSize(1+FitnessOfClass *Log[10^-5*10^12])) Sum[Exp[-InferredByFitness4M3[[k]][[i]][[6]]*(InferredByFitness4M3[[k]][[i]][[10]]+TauOffset)], {i, 1, Length[InferredByFitness4M3[[k]]]}];
MutationBound= MutationRateToBin;

{N[FitnessOfClass], MutationBound}, {k,1, Length[InferredByFitness4M3]}]][[1;;Round[0.20/DeltaS]]];


(* ::Text::Initialization:: *)
(*Frequencies of new mutants by fitness at time t*)


(* ::Input::Initialization:: *)

Tau[R_, s_, k_]:=Module[{MedianTau, LowerTau, UpperTau, NumberOfTauBins, dTau, RhoDeltaTau, CumulativeRho, NormalizedCumulative, TauList, NuList, NearestValueInCumulative, PositionInCumulative, TauSample},

If[R<10,
MedianTau=-(1/s)Log[R];
LowerTau=MedianTau-(4/s);
UpperTau=MedianTau+(10/s);
NumberOfTauBins=1000;
dTau=(UpperTau-LowerTau)/NumberOfTauBins;

RhoDeltaTau=Table[(s*dTau)/Gamma[R]*Exp[-R s tau-Exp[-s tau]], {tau, LowerTau,UpperTau, dTau}];
CumulativeRho=Table[Sum[RhoDeltaTau[[i]], {i, 1, j}], {j,1,Length[RhoDeltaTau]}];

NormalizedCumulative=CumulativeRho/Max[CumulativeRho];

TauList=Table[
beta=RandomReal[];
(*Print[beta];*)
NearestValueInCumulative=Nearest[NormalizedCumulative, beta][[1]];
(*Print[NearestValueInCumulative];*)
PositionInCumulative=Position[NormalizedCumulative,NearestValueInCumulative][[1]][[1]];
TauSample=LowerTau+PositionInCumulative*dTau;
TauSample, {j, 1, k}];,

NuList=Table[RandomVariate[GammaDistribution[R,1]], {j, 1, k}];
TauList=Table[-(1/s)Log[NuList[[i]]], {i, 1, Length[NuList]}];
];
TauList
]


(* ::Input::Initialization:: *)

MeanFitnessTrajectory[FitnessRates_,InitialFitnessFreqList_]:=Module[{R, s, \[Tau], InitialSize, InitialFreq, FitnessFreqListPreX,FitnessFreqList, MeanAttPreX, MeanAttNotPreX, MeanAtt, InitialFitnessFreqsPreX, OldFreq, NewFreq},

InitialFitnessFreqsPreX=Table[
R=6*10^8*FitnessRates[[i]][[2]];
s= FitnessRates[[i]][[1]];
\[Tau]=Tau[R, s,1][[1]];
InitialSize=Exp[s(-\[Tau])]/s;
InitialFreq=InitialSize/(5*10^8);
{s,InitialFreq}, {i, 1, Length[FitnessRates]}];
(*Print[FitnessFreqs];*)

MeanAtt=0;

FitnessFreqListPreX=InitialFitnessFreqsPreX;
FitnessFreqList=InitialFitnessFreqList;

Table[

FitnessFreqListPreX=Table[
s=FitnessFreqListPreX[[i]][[1]];
OldFreq=FitnessFreqListPreX[[i]][[2]];
NewFreq=OldFreq*Exp[(s-MeanAtt)];
{s,NewFreq}
, {i, 1, Length[FitnessFreqListPreX]}];
MeanAttPreX=Sum[FitnessFreqListPreX[[i]][[1]]*FitnessFreqListPreX[[i]][[2]],{i, 1, Length[FitnessFreqListPreX]}];

FitnessFreqList=Table[
s=FitnessFreqList[[i]][[1]];
OldFreq=FitnessFreqList[[i]][[2]];
NewFreq=OldFreq Exp[(s-MeanAtt)];
{s,NewFreq}, {i, 1, Length[FitnessFreqList]}];
MeanAttNotPreX=Sum[FitnessFreqList[[i]][[1]]*FitnessFreqList[[i]][[2]],{i, 1, Length[FitnessFreqList]}];

MeanAtt=MeanAttPreX+MeanAttNotPreX;

MeanAtt
, {t, 1, 120}]
]



(* ::Input::Initialization:: *)

SimulatedMeans2M3=Table[MeanFitnessTrajectory[InferredMutationRatesToClasses2M3,InitialFrequencies[FitnessTimeOfPreX, DeltaS]], {k, 1, 500}];
SimulatedMeans4M3=Table[MeanFitnessTrajectory[InferredMutationRatesToClasses4M3,InitialFrequencies[FitnessTimeOfPreX, DeltaS]], {k, 1, 500}];


(* ::Input::Initialization:: *)

SimulatedTrajectoryPlot2M3[traj_]:=ListPlot[traj, Joined -> True,ImageSize -> 250,Frame -> True,PlotStyle -> RGBColor[1.0,0.75,0.75,0.1],FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Mean Fitness", FontFamily-> "Helvetica", 10]]}, PlotRange -> {{0, 120},{0, 0.1}}];

SimulatedTrajectoryPlot4M3[traj_]:=ListPlot[traj, Joined -> True,ImageSize -> 250,Frame -> True,PlotStyle ->RGBColor[64/255, 204/255, 255/255, 0.1],FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Mean Fitness", FontFamily-> "Helvetica", 10]]}, PlotRange -> {{0, 120},{0, 0.1}}];


(* ::Input::Initialization:: *)

BothSimulatedMeans=RandomSample[Join[Table[SimulatedTrajectoryPlot2M3[SimulatedMeans2M3[[i]]], {i, 1, Length[SimulatedMeans2M3]}],Table[SimulatedTrajectoryPlot4M3[SimulatedMeans4M3[[i]]], {i, 1, Length[SimulatedMeans4M3]}]], Length[SimulatedMeans4M3]+Length[SimulatedMeans2M3]];


(* ::Input::Initialization:: *)
SimulatedTrajectories=Show[BothSimulatedMeans];
MeasuredTrajectories=Show[{ListPlot[{MeanFit2M3, MeanFit4M3}, PlotStyle -> {Darker[c2, 0.1],Darker[c7, 0.2]}, Joined -> True, PlotRange -> {{0,120},{-0.001, 0.1}}, ImageSize -> 170,Frame -> True,FrameLabel -> {Text[Style["Time (generations)", FontFamily-> "Helvetica", 10]],Text[Style["Mean Fitness", FontFamily-> "Helvetica", 10]]}]}];


(* ::Input::Initialization:: *)

MeansSimulated=Show[{SimulatedTrajectories, MeasuredTrajectories}, ImageSize -> 250]


(* ::Input::Initialization:: *)

Export["MeansSimulated.pdf",MeansSimulated,"PDF"]


(* ::Input::Initialization:: *)




(* ::Section::Initialization::Closed:: *)
(**)
(*End: Color palette*)


(* ::Input::Initialization:: *)

c0=RGBColor[142/255, 142/255, 147/255, 1];
c1a=RGBColor[0.6*255/255, 0.6*45/255, 0.6*85/255, 1];
c1=RGBColor[255/255, 45/255, 85/255, 1];
c2=RGBColor[255/255, 59/255, 48/255, 1];
c3=RGBColor[255/255, 149/255, 0/255, 1];
c4=RGBColor[255/255, 204/255, 0/255, 1];
c5=RGBColor[76/255, 217/255, 100/255, 1];
c6=RGBColor[90/255, 200/255, 250/255, 1];
c7=RGBColor[52/255, 170/255, 220/255, 1];
c8=RGBColor[0/255, 122/255, 255/255, 1];
c9=RGBColor[88/255, 86/255, 214/255, 1];
c10=RGBColor[0/255, 80/255, 145/255, 1];
c11=RGBColor[30/255, 75/255, 175/255, 1];
c12=RGBColor[0.5*88/255,0.5* 86/255,0.5* 214/255, 1];
c13=RGBColor[1,1,1,1]

ColorList=Append[Reverse[Table[ToExpression[StringJoin["c",ToString[i]]], {i, 1, 12}]], c1a];

Markers[s_, col_, Edge_]:=Graphics[{EdgeForm[Edge],col,Disk[{0,0},Scaled[s]]}];
Markers3[s_, col_, Edge_]:=Graphics[{EdgeForm[Edge],col,Rectangle[Scaled[{s,s}],Scaled[{0,0}]]}];
Markers2[s_,Error1_,Error2_, col_, Thic_]:=Graphics[{Thickness[Thic],col,Rectangle[{{-Error1,-Error2},{Error1,Error2}},Scaled[s]]}];

ClearMemory:=Module[{},Unprotect[In,Out];
Clear[In,Out];
Protect[In,Out];
ClearSystemCache[];];

Scheme1[x_]:=Blend[{c0,White}, x]
Scheme2[x_]:=Blend[{c1,c5, c0}, x]
Scheme3[x_]:=Blend[{c10,c8,c13,c2, RGBColor[0.77,0.04,0.]}, x]
Scheme4[x_]:=Blend[{c1,c3,c5, c7, c9}, x]

TableForm[{Table[Graphics[{ToExpression[StringJoin["c",ToString[i]]],Disk[{0.5,0.5},0.05]}], {i, 0, 12}]}, TableHeadings -> {{},Table[StringJoin["c", ToString[i]], {i, 0, 12}]}, TableSpacing -> {0, -3}]


(* ::Input::Initialization:: *)

Pos={{0, -20},{0.2, -16}, {0.4, -8}, {0.5, 0}, {0.6, 15}, {0.8, 100}, {1.0, 300}};

Pos2=Table[{i/0.14, i}, {i, 0, 0.14, 0.02}];


(* ::Input::Initialization:: *)
ScaleBar=DensityPlot[y,{x,0,1},{y,0,1},ColorFunction->Scheme3,ColorFunctionScaling->False,AspectRatio->15, FrameTicks->{{Table[{Pos[[i]][[1]], Superscript[10, Pos[[i]][[2]]]}, {i, 1, 7}],None},{None,None}} (*,FrameLabel \[Rule] {,Text[Style["P(beneficial)/P(neutral)", FontFamily\[Rule] "Helvetica", 13]]}*), PlotRangePadding ->None, ImageSize -> 36]

ScaleBar2=DensityPlot[y,{x,0,1},{y,0,1},ColorFunction->Scheme3,ColorFunctionScaling->False,AspectRatio->15, FrameTicks->{{Table[{Pos2[[i]][[1]], Pos2[[i]][[2]]}, {i, 1, Length[Pos2]}],None},{None,None}} (*,FrameLabel \[Rule] {,Text[Style["P(beneficial)/P(neutral)", FontFamily\[Rule] "Helvetica", 13]]}*), PlotRangePadding ->None, ImageSize -> 36]




(* ::Input::Initialization:: *)

Export["ScaleBar2.pdf",ScaleBar2,"AllowRasterization"->True,ImageResolution->400]


(* ::Input::Initialization:: *)

SetDirectory["~/Dropbox/BayesianAnalysis/2M3"];
pb1=ToExpression[Import["BayesianAnalysis2M3_1-100000.csv"]];
pb2=ToExpression[Import["BayesianAnalysis2M3_100001-300000.csv"]];
pb3=ToExpression[Import["BayesianAnalysis2M3_300001-End.csv"]];
Combined=Flatten[Join[pb1, pb2, pb3], {1}];


(* ::Input::Initialization:: *)

Export["BayesianAnalysis2M3.csv",Combined,"CSV"]
