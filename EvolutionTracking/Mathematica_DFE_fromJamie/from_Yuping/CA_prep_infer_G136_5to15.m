(* ::Package:: *)

(* ::Section:: *)
(*1. Read trajectories, read depths, number trajectories (shift + enter to run current chunk)*)

SetDirectory["/labs/dpetrov/yupingli/DA_Mathematica/CA"];
AllTrajectories=Import["CA_CountsTimecourse_GoodTimepoints.txt", "Table"];
AllTrajectories =AllTrajectories[[ All, 2;;]];
PrepIDs=AllTrajectories[[1]];
AllReadTrajectories=Delete[AllTrajectories,1];
NumberOfBarcodes=Length[AllReadTrajectories]

FactorInteger[NumberOfBarcodes]


(*The time points that are used in this analysis are:*)
TimePoints={0,16,24,32,40,48,56,64,72,80,88,96,104,112,120,128,136};


(*This function combined all the reads at a given time point. One has to be careful in doing this because each is prepped separatedly so the error associated with a DNA prep (some fraction of \[Kappa])  will add up for each prep.*)
(*(* before join, join[ {{0,xx}}, {{16,xx}}, ....]; after join, {0,xx},{16,xx}, ...; thus each row is a [length(time_points) * 2] matrix*)*)

CombineTimePoints[Trajectory_]:=Module[{t0,t16,t24,t32,t40,t48,t56,t64,t72,t80,t88,t96,t104,t112,t120,t128,t136},
	t0={{    0,If[   Trajectory[[1]]==0,1, Trajectory[[1]]  ]     }};
	t16={{  16,Trajectory[[2]]    }};
	t24={{  24,Trajectory[[3]]     }};
	t32={{  32,Trajectory[[4]]     }};
	t40={{  40,Trajectory[[5]]     }};
	t48={{  48,Trajectory[[6]]     }};
	t56={{  56,Trajectory[[7]]     }};
	t64={{  64,Trajectory[[8]]     }};
	t72={{  72,Trajectory[[9]]     }};
	t80={{  80,Trajectory[[10]]   }};
	t88={{  88,Trajectory[[11]]    }};
	t96={{  96,Trajectory[[12]]    }};
	t104={{  104,Trajectory[[13]]    }};
	t112={{  112,Trajectory[[14]]    }};
	t120={{  120,Trajectory[[15]]    }};
	t128={{  128,Trajectory[[16]]  }};
	t136={{  136,Trajectory[[17]]    }};

	Join[t0,t16,t24,t32,t40,t48,t56,t64,t72,t80,t88,t96,t104,t112,t120,t128,t136]
]; (* before join, join[ {{0,xx}}, {{16,xx}}, ....]; after join, {0,xx},{16,xx}, ...; thus each row is a [length(time_points) * 2] matrix *)


(*List of combined reads for each barcode. Export to a file. *)
ReadTrajectories2M3=Map[CombineTimePoints, AllReadTrajectories]; (* map apply to function to every row *)
Export["ReadTrajectories_CA_CountsTimecourse_GoodTimepoints.csv",ReadTrajectories2M3,"CSV"];


(*List of total read depths. Export to file.*)
ReadDepths2M3=Table[
	ReadsByTime=Transpose[ReadTrajectories2M3];
	ReadsAddedUp=Total[ReadsByTime[[i]]]; (* total for each column *)
	{TimePoints[[i]],ReadsAddedUp[[2]]}
	, {i, 1, Length[TimePoints]}];
Export["ReadDepths_CA_CountsTimecourse_GoodTimepoints.csv",ReadDepths2M3,"CSV"];


(*Convert read trajectories to approximate number trajectories for plotting purposes. Export to file.*)
ReadTrajectories2M3=ToExpression[Import["ReadTrajectories_CA_CountsTimecourse_GoodTimepoints.csv"]];
ReadDepths2M3=ToExpression[Import["ReadDepths_CA_CountsTimecourse_GoodTimepoints.csv"]];

PopSize=5*10^8;
ExtinctLineages={};


(*add random noise from RandomReal; go through every TIMEPOINT {i,1,length} and every BARCODE *)
CellNumberTrajectories2M3=Table[
	If[Mod[Barcode, 50000]==0, Print[Barcode];];
	NumberTraj1=Table[{ReadTrajectories2M3[[Barcode]][[i]][[1]],If[ReadTrajectories2M3[[Barcode]][[i]][[2]]==0,0,N[(PopSize/ReadDepths2M3[[i]][[2]])*(ReadTrajectories2M3[[Barcode]][[i]][[2]])]+RandomReal[{-0.5,0.5}] ]}, {i, 1, Length[ReadDepths2M3]}]; 

	TimePointsBarcodeUnRead=Sort[Map[First,Select[NumberTraj1, #[[2]]==0&]], #1>#2&];
	TimePointsBarcodeRead=Complement[TimePoints,TimePointsBarcodeUnRead];

	If[Length[TimePointsBarcodeUnRead]>0&&Max[TimePointsBarcodeRead]<Min[TimePointsBarcodeUnRead], ExtinctLineage=1; AppendTo[ExtinctLineages, Barcode];, ExtinctLineage=0];

	If[ExtinctLineage==1, NumberTraj2=Replace[NumberTraj1, 0-> 0.1, {2}], NumberTraj2=NumberTraj1];

	NumberTraj2
	, {Barcode, 1,Length[ReadTrajectories2M3]}];

Export["CellNumberTrajectories_CA_CountsTimecourse_GoodTimepoints.csv", CellNumberTrajectories2M3, "CSV"];
Export["ExtinctLineages_CA_CountsTimecourse_GoodTimepoints.csv", ExtinctLineages, "CSV"];

(*Does not seem that we need Extinct Lineages later, let us move forward for now *)


(* ::Section:: *)
(*2. Noise model (Supp P20)*)

(*This distribution comes up in many of the processes because it is the approximation of the inverse laplace transform of *)
(**)
(* M(\[Phi]) = exp(- a \[Phi] / (b \[Phi] +1)) *)
(**)
(*which is a form that appears in birth death processes. The function is peaked around a but has an exponential tail rather than a Gaussian one. *)


LogProb[a_, b_, n_]:=-((Sqrt[n]-Sqrt[a])^2/b)+Log[Sqrt[Sqrt[a]/(4 \[Pi] b n^(3/2))]] (* Equal 33 *)


(*Use this form to estimate: \[Kappa], and xbar between neighboring time points. We condition on a number of reads r1 then fit the measured distribution of reads r2 conditioned on r1 i.e. P(r2|r1). We perform a two parameter fit for (xbar, \[Kappa]) using the above form for the distribution with*)
(*b=\[Kappa]*)
(*We do this conditioning on reads between 20 - 30 then take the mean of the fitted parameters. 20-30 is chosen because it is small enough that we assume neutral, but not too small than integer number effects and fact we do not model things perfectly at small n become important. *)
(***YL:  lineages with 20-30 reads are considered to be neutral. It might be too high for DA data, but still use 20-30 (Check fitness trajectories to visualize whether it is reasonable)*)
(*The outputted mean fit values between the neighboring time points is written in to a file as is the series of plots showing how well the fit works across all r1 from LowerReadLimit (20 normally) to UpperReadLimit (40 usually).*)



ReadTrajectories2M3=ToExpression[Import["ReadTrajectories_CA_CountsTimecourse_GoodTimepoints.csv"]];
ReadDepths2M3=ToExpression[Import["ReadDepths_CA_CountsTimecourse_GoodTimepoints.csv"]];


(*The noise in the frequency of a barcode read r2 times conditioned on r1 is captures by \[Kappa]. *)
(*With a \[Kappa]=1 the variance in reads would be simply given by  ~ r2 i.e. only sequencing noise. \[Kappa] >1 indicatse how much additional noise there is from other sources e.g. ampliciation / extraction and drift. From time point 0-8 in 2M3, for example, \[Kappa]~4. This means that the variance in r2 reads is about 4r2. Four times larger than expected from sequencing alone. Why? The reads for 0 are actually pooled together from 5 separate preps and sequencing runs while the reads in 8 are from two separate preps / runs. For each separate prep / run variances associated with preps ADD. If the true noise from a SINGLE sequencing + prep was \[Kappa]r2=(1+\[Epsilon])r2 [\[Epsilon] here the noise from the prep in units of the sequcing noise, which in principle is probably frequency dependent, but not v strongly], then pooling n such preps and counting total reads, we should expect \[Epsilon] -> n\[Epsilon].  From previous analysis, we know that additional noise due to prep is \[Epsilon]~0.3 - 0.5 (a reasonable fraction of the sequecning noise) so given n~7 for 0-8 this is roughly consistent with the \[Kappa]~4 that we see.*)
(*I want to make sure that \[Kappa] really is not very frequency dependent. From previous analysis it should have a dependence like (1+R/\[Beta] N) where \[Beta] is the fraction of the N cells at saturation that yield DNA and contribute to PCR. This would predict that it should be a linearly increasing function of read depth.  *)


ReadsToUse={};
ListOfr1Reads=Table[ReadTrajectories2M3[[i]][[1]][[2]], {i, 1, Length[ReadTrajectories2M3]}]; (* Only reads at T1 *)
BinSize=4;
NumbersInBins=BinCounts[ListOfr1Reads, {0, 2000,BinSize}]; (* YL: JB bin range-{0,2000}, width-BinSize and NumbersInBinds>1000 can be modified *) 
Do[If[NumbersInBins[[Bin]]>1000&&(Bin-0.5)*BinSize>20, AppendTo[ReadsToUse, Round[(Bin-0.5)*BinSize]]],{Bin, 1, Length[NumbersInBins]}]; (* YL: JB Use the bin with more than 1000 barcodes and with a read number > 24 *)

Table[
	t1=TimePoints[[t]];
	t2=TimePoints[[t+1]];
	R1=ReadDepths2M3[[t]][[2]];
	R2=ReadDepths2M3[[t+1]][[2]];

	KappaWithr0=Table[
	r1=r0;
	pm=2;
	
	BarcodesInRange={};

	Do[
		Measuredr1=ReadTrajectories2M3[[Barcode]][[t]][[2]]; (* Reads at T1 *)
		If[Abs[Measuredr1-r1]<= pm, AppendTo[BarcodesInRange, Barcode]];, {Barcode, 1,Length[ReadTrajectories2M3]}]; (* YL: JB e.g. A barcode with read 20 will be added when r0=20 from ReadsToUse *)
		Setofr2Reads=DeleteCases[   Table[ ReadTrajectories2M3[[Barcode]][[t+1]][[2]], {Barcode, BarcodesInRange} ], 0    ];
		NumberOfBarcodesInRange=Length[BarcodesInRange];
		Minr2=Max[{Min[Setofr2Reads], 10}];Maxr2=Max[Setofr2Reads];(*Do not fit to dist below 10 reads since discrete effects come in*)
		r2Distribution=BinCounts[Setofr2Reads, {Minr2, Maxr2, 1}];(* YL: CHANGE step from 1 to longer? *)
		Predictedr2Reads[xbar_, \[Kappa]_]:=Table[NumberOfBarcodesInRange*Exp[LogProb[(R2/R1)*r1*Exp[-(t2-t1)*xbar],\[Kappa],r2]], {r2, Minr2, Maxr2-1, 1}];
		GoodnessOfFit[xbar_, \[Kappa]_]:=Total[N[(r2Distribution-Predictedr2Reads[xbar, \[Kappa]])^2]];

		FittedxBarKappa=NMinimize[{GoodnessOfFit[0.00001, y], 0.5<y<6}, {y},AccuracyGoal->6,MaxIterations->50,Method -> {"RandomSearch", "SearchPoints"-> 200}];
		Kappa=FittedxBarKappa[[2]][[1]][[2]];
		Print["Kappa=",Kappa];
		{r0,Kappa}
		, {r0,ReadsToUse}
	];

, {t, 1,1} (* YL: JB kappa calculated based on t1 ---commonly defined as t0 -- & t2 *)
];


Kappa_plot = ListPlot[KappaWithr0, PlotRange -> {0, 6}]


(* ::Subsection:: *)
(**)
(*Inferring the mean fitness and unadaptive fraction (Section 6.1: Use Neutrals to Estimate Mean Fitness)*)



LowerReadLimit=5; (* YL: Consider change lower/upper ReadLimit *)
UpperReadLimit=15;
GridOfFittedPlots=Table[Table[,{r, LowerReadLimit,UpperReadLimit}], {t, 1, Length[TimePoints]-1}];
AllxBarKappaFits=Table[
	t1=TimePoints[[t]];
	t2=TimePoints[[t+1]];
	R1=ReadDepths2M3[[t]][[2]];
	R2=ReadDepths2M3[[t+1]][[2]];

	ReadNumbersToUse=Table[i, {i, LowerReadLimit, UpperReadLimit}];

	Allr1Fits=Table[
		r1=r0;
		pm=1;
		BarcodesInRange={};

		Do[
			Measuredr1=ReadTrajectories2M3[[Barcode]][[t]][[2]];
			If[Abs[Measuredr1-r1]<= pm, AppendTo[BarcodesInRange, Barcode]];, {Barcode, 1,Length[ReadTrajectories2M3]}
		];
		
		Setofr2Reads=DeleteCases[Table[ReadTrajectories2M3[[Barcode]][[t+1]][[2]], {Barcode, BarcodesInRange}], 0];
		NumberOfBarcodesInRange=Length[BarcodesInRange];
		Minr2=Max[{Min[Setofr2Reads], 10}];Maxr2=Max[Setofr2Reads];
		r2Distribution=BinCounts[Setofr2Reads, {Minr2, Maxr2, 1}];
		Predictedr2Reads[xbar_, \[Kappa]_]:=Table[NumberOfBarcodesInRange*Exp[LogProb[(R2/R1)*r1*Exp[-(t2-t1)*xbar],\[Kappa],r2]], {r2, Minr2, Maxr2-1, 1}];
		GoodnessOfFit[xbar_, \[Kappa]_]:=Total[N[(r2Distribution-Predictedr2Reads[xbar, \[Kappa]])^2]];

		FittedxBarKappa=NMinimize[{GoodnessOfFit[x, y], 0.00001<x<0.5, 0.5<y<6}, {x,y}];
		xBar=FittedxBarKappa[[2]][[1]][[2]];
		Kappa=FittedxBarKappa[[2]][[2]][[2]];

		Fittedr2Reads=Table[{r2+0.5,NumberOfBarcodesInRange*Exp[LogProb[(R2/R1)*r1*Exp[-(t2-t1)*xBar],Kappa,r2]]}, {r2, Minr2, 2*Maxr2, 1}];

		Fittedr2ReadsNoMean=Table[{r2+0.5,NumberOfBarcodesInRange*Exp[LogProb[(R2/R1)*r1,Kappa,r2]]}, {r2, Minr2, 2*Maxr2, 1}];

		LogData=Histogram[Setofr2Reads, {1}, "LogCount",PlotRange ->{{1, 80},All},Frame -> True,AxesOrigin -> {0, -1},FrameLabel -> {Text[Style["Number of reads", 12,FontFamily-> "Helvetica"]],Text[Style["Number of barcodes", 12,FontFamily-> "Helvetica"]]}, FrameTicksStyle->Directive[Black,12, FontFamily-> "Helvetica"], ChartStyle ->Directive[{c4, EdgeForm[None]}]];
		LogTheory=ListLogPlot[{Fittedr2Reads, Fittedr2ReadsNoMean},Joined -> True, AxesOrigin -> {0, -1},PlotStyle-> {{Black, Thickness[0.01]},{Black, Dashing[0.03]} }, PlotRange ->{{1, 80},{-1,Log[Max[r2Distribution]]}}];
		LOGPLOT=Show[{LogData, LogTheory}, PlotLabel-> Text[Style[StringJoin["Kappa=",ToString[Kappa],"   xBar=",ToString[xBar],"   r=",ToString[Round[r0]]],12, FontFamily-> "Helvetica"]], ImageSize -> 250, PlotRange -> {{1, 80},{-1,Log[Max[r2Distribution]]}}, PlotRangeClipping -> True, AxesOrigin -> {0, -1}];
		GridOfFittedPlots[[t]][[r0-LowerReadLimit+1]]=LOGPLOT;

		{xBar, Kappa}

		, {r0,ReadNumbersToUse}];
		
	Mean[Allr1Fits]
		
, {t, 1,Length[TimePoints]-1}];

Export["InferredxBarAndKappa_CA_CountsTimecourse_GoodTimepoints.csv", AllxBarKappaFits, "CSV"];
Export["PlotsOfInferredxBarAndKappa_CA_CountsTimecourse_GoodTimepoints.pdf", TableForm[GridOfFittedPlots, TableHeadings -> {Table[TimePoints[[i]], {i, 1, Length[TimePoints]}],Table[r, {r, LowerReadLimit,UpperReadLimit}]}], "PDF"];
(* Figure 17 & 19 *)





(* ::Section:: *)
(**)
(*3. Calculating probabilities (Bayesian P36-P39: Infer s and tau Using Mean Fitness Estimated from Neutrals)*)


(* ::Subsection:: *)
(*Inference*)


(* ::Text:: *)
(*Importing the data for both CA_CountsTimecourse_GoodTimepoints*)



$HistoryLength=0;


ReadDepths2M3=ToExpression[Import["ReadDepths_CA_CountsTimecourse_GoodTimepoints.csv"]];
ReadTrajectories2M3=ToExpression[Import["ReadTrajectories_CA_CountsTimecourse_GoodTimepoints.csv"]]; (*all lines*)
CellNumberTrajectories2M3=ToExpression[Import["CellNumberTrajectories_CA_CountsTimecourse_GoodTimepoints.csv"]]; (*all lines*)
ClearMemory


(* ::Text:: *)
(*Parameters that are useful to define*)

SetDirectory["/labs/dpetrov/yupingli/DA_Mathematica/CA"];
T=8;
PopSize=5*10^8;
NumberOfBarcodes=Length[ReadTrajectories2M3];
LogProb[a_, b_, n_]:=-((Sqrt[n]-Sqrt[a])^2/b)+Log[Sqrt[Sqrt[a]/(4 \[Pi] b n^(3/2))]]; (* Equal 33; Need THIS from earlier *)


(* ::Text:: *)
(*Mean fitness list defined so that the mean fitness at 16 corresponds to the mean fitness accounting for reduction in freqs from 8-16.*)


xBarKappa2M3=ToExpression[Import["InferredxBarAndKappa_CA_CountsTimecourse_GoodTimepoints.csv"]];
TimePoints2M3={0,16,24,32,40,48,56,64,72,80,88,96,104,112,120,128,136};
MeanFitnessList2M3=Table[{TimePoints2M3[[i+1]], xBarKappa2M3[[i]][[1]]},{i, 1, Length[xBarKappa2M3]}];
KappaList2M3=Table[{TimePoints2M3[[i+1]], xBarKappa2M3[[i]][[2]]},{i, 1, Length[xBarKappa2M3]}];


(* ::Text:: *)
(*Determine the interpolation formulae for mean fitness and the unadaptive fraction (integral of mean fitness over time)*)


MeanFitnessBeforeZero2M3=Table[ {8i, 10^-7}, {i, -25, TimePoints2M3[[1]]/T}];
MeanFitnessExtended2M3=Join[ MeanFitnessBeforeZero2M3,MeanFitnessList2M3];
UnAdaptedFraction2M3=Table[{MeanFitnessExtended2M3[[j]][[1]],Exp[-Sum[(MeanFitnessExtended2M3[[t]][[1]]-MeanFitnessExtended2M3[[t-1]][[1]])*MeanFitnessExtended2M3[[t]][[2]], {t, 2,j}]]}, {j, 2, Length[MeanFitnessExtended2M3]}];
InterpolatedMeanFitness2M3=Interpolation[MeanFitnessExtended2M3, InterpolationOrder->1];
InterpolatedUnAdaptiveFraction2M3=Interpolation[UnAdaptedFraction2M3, InterpolationOrder->1];


(* ::Text:: *)
(*Establishment time as a function of s and \[Tau]. True establishment size is 1/(s-\!\(\*OverscriptBox[\(x\), \(_\)]\)(\[Tau]))*)



EstablishmentSize2M3[s_, \[Tau]_]:=1/(s-InterpolatedMeanFitness2M3[\[Tau]]) *UnitStep[s-(InterpolatedMeanFitness2M3[\[Tau]]+0.005)]+1/0.005*UnitStep[-(s-(InterpolatedMeanFitness2M3[\[Tau]]+0.005))];



(* ::Text:: *)
(*Given this the r2 is fixed by the deterministic growth of the beneficial cells emerging from \[Tau] with fitness s.*)


Meanr2Givenr12M3[r1_,R1_,R2_,t1_,t2_,s_,\[Tau]_, PopSize_]:=
Module[{TotalNumberOfCells,NumberOfBeneficialCells,NumberOfNeutralCells,TotalNumberOfCellsAtt2,NumberOfReadsr2},
TotalNumberOfCells=PopSize(r1/R1);
NumberOfBeneficialCells=(*Exp[-s \[Tau]](Exp[s t1]-1)*)Exp[s (t1-\[Tau])]*EstablishmentSize2M3[s, \[Tau]]*InterpolatedUnAdaptiveFraction2M3[t1]/InterpolatedUnAdaptiveFraction2M3[\[Tau]];
NumberOfNeutralCells=TotalNumberOfCells-NumberOfBeneficialCells;

If[NumberOfBeneficialCells>=TotalNumberOfCells, NumberOfNeutralCells=0;NumberOfBeneficialCells=TotalNumberOfCells;];

TotalNumberOfCellsAtt2=(NumberOfNeutralCells+Exp[s(t2-t1)]*NumberOfBeneficialCells)*InterpolatedUnAdaptiveFraction2M3[t2]/InterpolatedUnAdaptiveFraction2M3[t1];

NumberOfReadsr2=(TotalNumberOfCellsAtt2/PopSize)*R2;
NumberOfReadsr2];



(* ::Text:: *)
(*The propagator is then simply the standard error function with \[Kappa] as the noise and the mean set by the above function*)
(**)
(*Dealing with logs is better than the probability itself or numbers get too large and computational overflow occurs.*)



LogPropagator2M3[r1_,r2_,R1_,R2_,t1_,t2_,PopSize_,\[Kappa]_,s_,\[Tau]_]:=LogProb[
Meanr2Givenr12M3[r1,R1,R2,t1,t2,s,\[Tau], PopSize], \[Kappa], r2]


(* ::Text:: *)
(*The prior distriution on (s, \[Tau]) is set by the theory of a constant feeding process. If a constant population of size n0 is feeding mutants at rate \[Mu](s)\[Delta]s into fitness range (s, s+\[Delta]s) then the statistics of the establishment times follow from solving the backwards time generating function equations. *)
(**)
(*The distribution is \[Rho](\[Tau]) = (s d\[Tau] / \[CapitalGamma](n U))*exp(-n U s \[Tau] - e^(-s \[Tau]))*)
(**)
(*To do this we need to have a prior distribution over \[Mu](s)ds. This we have to guess initially and then can improve upon once *)



DeltaS=0.1;
\[Mu]Trial[s_]:=10^-5 1/DeltaS Exp[-s/DeltaS];
LogPriorDensityTau[\[Tau]_,s_,n_, U_]:=Log[n U s]


LogProbFirstTimePoint2M3[r_, s_, \[Tau]_, Neutral_]:=Module[{NumberOfBeneficialCellsAtZero,NumberOfCellsInferredAtZero, ExpectedReadsAtZero,TotalNumberOfCellsAtZero, ProbOft0},
If[Neutral==1,ProbOft0=LogProb[r, 2, r],
NumberOfBeneficialCellsAtZero=Exp[-s \[Tau]]*EstablishmentSize2M3[s, \[Tau]]*1/InterpolatedUnAdaptiveFraction2M3[\[Tau]];
NumberOfCellsInferredAtZero=(r*PopSize)/ReadDepths2M3[[1]][[2]];
 TotalNumberOfCellsAtZero=NumberOfCellsInferredAtZero;

If[NumberOfBeneficialCellsAtZero>NumberOfCellsInferredAtZero, TotalNumberOfCellsAtZero=NumberOfBeneficialCellsAtZero];
ExpectedReadsAtZero=TotalNumberOfCellsAtZero/PopSize*ReadDepths2M3[[1]][[2]];
ProbOft0=LogProb[ExpectedReadsAtZero, 2, r]];
ProbOft0]


t1=SessionTime[];
TableOfData={};
Do[
	If[ Mod[Barcode, 20000]==0,Print[Barcode] ];
	Traj2M3=Replace[ReadTrajectories2M3[[Barcode]], 0-> 1, 2];
	(*Get rid of zeros since prob does not model them. YL NOTE: TIME 0 will be changed to TIME 1*)


	(*DETERMINE INITIAL SIZE OF EACH LINEAGE*)

	InitialLineageSize2M3=PopSize/ReadDepths2M3[[1]][[2]]*Traj2M3[[1]][[2]];

	\[Delta]s=0.005; (* YL: Increase this to accelerate computation *)
	\[Delta]\[Tau]=1; (* YL: Increase to accelerate computation *)


	(*PROBABILITY OF THE NULL HYPOTHESIS*)

	WeightOfT0For2M3=-LogProbFirstTimePoint2M3[Traj2M3[[1]][[2]], 0, 0, 1];
	WeightOfLaterTFor2M3=Sum[-LogPropagator2M3[Traj2M3[[t]][[2]],Traj2M3[[t+1]][[2]],ReadDepths2M3[[t]][[2]],ReadDepths2M3[[t+1]][[2]],TimePoints2M3[[t]],TimePoints2M3[[t+1]],10^5 PopSize,KappaList2M3[[t]][[2]],0.00001,0], {t, 1,  Length[TimePoints2M3]-1}];
	WeightNotBeneficial2M3=WeightOfT0For2M3+WeightOfLaterTFor2M3;

	If[ Mod[Barcode, 2000]==0,Print[WeightNotBeneficial2M3] ];

	If[ WeightNotBeneficial2M3>30, (* YL: how to decide the cutoff here? 60 seems more likely *)

		t2=SessionTime[];
		If[Mod[Barcode, 10000]==0,Print[Barcode, "Memory=", MemoryInUse[], "  max",MaxMemoryUsed[]]; Print[t2-t1];];
		t1=t2;

		(*PROBABILITY OF THE S, T HYPOTHESIS*)

		PriorWeight2M3[s_, \[Tau]_]:=-LogPriorDensityTau[\[Tau], s,InitialLineageSize2M3, \[Mu]Trial[s]*\[Delta]s]-Log[\[Delta]\[Tau]];
		WeightT02M3[s_, \[Tau]_]:=-LogProbFirstTimePoint2M3[Traj2M3[[1]][[2]], s, \[Tau], 0];
		WeightOfLaterT2M3[s_, \[Tau]_]:=Sum[-LogPropagator2M3[Traj2M3[[t]][[2]],Traj2M3[[t+1]][[2]],ReadDepths2M3[[t]][[2]],ReadDepths2M3[[t+1]][[2]],TimePoints2M3[[t]],TimePoints2M3[[t+1]],PopSize,KappaList2M3[[t]][[2]],s,\[Tau]], {t, 1, Length[TimePoints2M3]-1}];

		ProbSAndTauGivenData2M3[s_, \[Tau]_]:=PriorWeight2M3[s, \[Tau]]+WeightT02M3[s, \[Tau]]+WeightOfLaterT2M3[s, \[Tau]];

		(*FIND MOST LIKELY S AND T FOR BOTH EXPERIMENTS*)

		MostProbableSAndTau2M3=NMinimize[{ProbSAndTauGivenData2M3[s, \[Tau]], 0.005<s<0.4 && -150<\[Tau]<90}, {s, \[Tau]}, AccuracyGoal->6,MaxIterations->20,Method -> "NelderMead"];
		Print[Barcode,"   (s, \[Tau])_2M3 =  ", MostProbableSAndTau2M3];

		(*IF PR[BENEFICIAL] > PR[NEUTRAL] THEN ASSIGN BEST FIT VALUES FOR S AND T, ELSE ASSIGN ZERO*)

		InferredS2M3=MostProbableSAndTau2M3[[2]][[1]][[2]];
		InferredTau2M3=MostProbableSAndTau2M3[[2]][[2]][[2]];

		If[WeightNotBeneficial2M3>MostProbableSAndTau2M3[[1]],
			InferredS2M3=MostProbableSAndTau2M3[[2]][[1]][[2]];
			InferredTau2M3=MostProbableSAndTau2M3[[2]][[2]][[2]];
			, 
			InferredS2M3=0.00001; 
			InferredTau2M3=0;];


		(*LOOK AT CONTRIBUTIONS TO THE WEIGHTS*)
		Print["ML T0 2M3 =",N[WeightT02M3[InferredS2M3, InferredTau2M3]],"ML Later T 2M3 = ",WeightOfLaterT2M3[InferredS2M3, InferredTau2M3],"Total Weight ML 2M3 = ",ProbSAndTauGivenData2M3[InferredS2M3, InferredTau2M3]]; 



		(*ERRORS ON INFERRED VALUES FROM CURVATURE*)

		If[InferredS2M3>0.001,
			GradientInS2M3=-\!\(\*SubscriptBox[\(\[PartialD]\), \(s\)]\ \(ProbSAndTauGivenData2M3[s, \ InferredTau2M3]\)\)/. s -> InferredS2M3;
			GradientInTau2M3=-\!\(\*SubscriptBox[\(\[PartialD]\), \(\[Tau]\)]\ \(ProbSAndTauGivenData2M3[InferredS2M3, \ \[Tau]]\)\)/. \[Tau] -> InferredTau2M3;

			CurvatureSS2M3=\!\(\*SubscriptBox[\(\[PartialD]\), \(s\)]\ \(\*SubscriptBox[\(\[PartialD]\), \(s\)]\ ProbSAndTauGivenData2M3[s, \ InferredTau2M3]\)\)/. s -> InferredS2M3;
			CurvatureST2M3=\!\(\*SubscriptBox[\(\[PartialD]\), \(\[Tau]\)]\ \(\*SubscriptBox[\(\[PartialD]\), \(s\)]\ ProbSAndTauGivenData2M3[s, \ \[Tau]]\)\)/. {s -> InferredS2M3, \[Tau] -> InferredTau2M3};
			CurvatureTT2M3=\!\(\*SubscriptBox[\(\[PartialD]\), \(\[Tau]\)]\ \(\*SubscriptBox[\(\[PartialD]\), \(\[Tau]\)]\ ProbSAndTauGivenData2M3[InferredS2M3, \ \[Tau]]\)\)/. \[Tau] -> InferredTau2M3;

			Kmatrix2M3={{CurvatureSS2M3,CurvatureST2M3},{CurvatureST2M3,CurvatureTT2M3}};
			EigenValsAndVects2M3=Eigensystem[Kmatrix2M3];
			PrincipleCurvature2M3=EigenValsAndVects2M3[[1]][[2]];
			PrincipleCurvatureDirection2M3=EigenValsAndVects2M3[[2]][[2]];

			ErrorInS2M3=Abs[PrincipleCurvatureDirection2M3[[1]]]/Sqrt[PrincipleCurvature2M3];
			ErrorInTau2M3=Abs[PrincipleCurvatureDirection2M3[[2]]]/Sqrt[PrincipleCurvature2M3];
			,
			ErrorInS2M3=0.01;
			ErrorInTau2M3=5];


		AppendTo[TableOfData, {Barcode,-MostProbableSAndTau2M3[[1]]+WeightNotBeneficial2M3,InferredS2M3,ErrorInS2M3, InferredTau2M3,ErrorInTau2M3}];
		NumberInDataList=Length[TableOfData];
		(*Print["NumberInDataList=", NumberInDataList];*)

		If[NumberInDataList==50000,(*YL: Need to change this for output *)
		Export[StringJoin["CA_BarcodeDATA_G136_5to15_",ToString[TableOfData[[1]][[1]]],"-",ToString[TableOfData[[50000]][[1]]],".CSV"],TableOfData,"CSV"]; (* all lines *)
		TableOfData={};
		];

		ClearMemory

	];

	,{Barcode, Length[ReadTrajectories2M3]}
];

Print["Finished, evaluating"];

(*all lines *)
Export[StringJoin["CA_BarcodeDATA_G136_5to15_",ToString[TableOfData[[1]][[1]]],"-",ToString[TableOfData[[Length[TableOfData]]][[1]]],".CSV"],TableOfData,"CSV"];
Exit[]