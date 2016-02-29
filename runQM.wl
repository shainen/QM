(* ::Package:: *)

(*SetDirectory[Directory[]<>"/QM"];*)


SetDirectory[NotebookDirectory[]]


tscale=.2;


tmax=1;
steps=500;
times=Range[0,tmax,tmax/(steps-1)];


length1=3;
length2=3;
startnumferms=5;


(*finMu=10;*)


(*chemPot[t_]:=-finMu(1-2E^(-t^2/tscale^2))*)


(*chemPot[t_]:=Piecewise[{{finMu Cos[2 \[Pi] t/tmax],t<2 tscale},{-finMu,t\[GreaterEqual]2 tscale}}]*)


coup[t_]:=(1-E^(-t^2/tscale^2))


(*coup[t_]:=Piecewise[{{(1-Cos[\[Pi] t/tscale/2])/2,t<2 tscale},{1,t\[GreaterEqual]2 tscale}}]*)


<<makeStates.wl


thestates=makeStates[length1,length2,startnumferms];


dim=Length[thestates];


bosonchempot=SparseArray[{ii_,ii_}:>N[Total[thestates[[ii,1]],2]],{dim,dim}];


fermenergy=Table[N[-2(Cos[2\[Pi] ii/length1]+Cos[2\[Pi] jj/length2])],{ii,0,length1-1},{jj,0,length2-1}];


fermikinham=SparseArray[{ii_,ii_}:>Total[fermenergy (thestates[[ii,2]]+thestates[[ii,3]]),2],{dim,dim}];


Get["151008_1_mkham3b3s5.dat"];


hamndintfull=hamInt3b3s5;


(*hamTot[t_]:=chemPot[t]bosonchempot+fermikinham+coup[t]hamndintfull/Sqrt[length1 length2];*)


constHam=-10 bosonchempot+fermikinham;


changeHam=20 hamndintfull/Sqrt[length1 length2];


(*hamTotF[t_]:=-10 bosonchempot+fermikinham+20 coup[t] hamndintfull/Sqrt[length1 length2];*)


hamTot=constHam+changeHam;


init=SparseArray[Position[thestates,{({
 {0, 0, 0},
 {0, 0, 0},
 {0, 0, 0}
}),({
 {1, 1, 1},
 {1, 0, 0},
 {1, 0, 0}
}),({
 {1, 1, 1},
 {1, 0, 0},
 {1, 0, 0}
})}][[1,1]]->1,{dim}];


Timing[ket=NDSolveValue[{psi[0]==Normal[init],psi'[t]==-I (constHam.psi[t]+coup[t](changeHam.psi[t]))},psi,{t,0,tmax}]/@times;]


momNumsQM=Transpose[(Abs[ket]^2).thestates,{4,1,2,3}];


intObs=(#\[Conjugate].hamndintfull.#&/@ket)/Sqrt[length1 length2];


mmu=MaxMemoryUsed[]/10.^6;


SetDirectory[ParentDirectory[]];


(*Save["fullKet.dat",{mmu,ket}];*)


(*Save["momNumsQM.dat",{mmu,momNumsQM,intObs}];*)
