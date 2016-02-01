(* ::Package:: *)

SetDirectory[Directory[]<>"/QM"];


<<makeStates.wl


tmax=50;
steps=500;
times=Range[0,tmax,tmax/(steps-1)];


length1=3;
length2=3;
startnumferms=5;


thestates=makeStates[length1,length2,startnumferms];


dim=Length[thestates];


bosonchempot=SparseArray[{ii_,ii_}:>N[Total[thestates[[ii,1]],2]],{dim,dim}];


fermenergy=Table[N[-2(Cos[2\[Pi] ii/length1]+Cos[2\[Pi] jj/length2])],{ii,0,length1-1},{jj,0,length2-1}];


fermikinham=SparseArray[{ii_,ii_}:>Total[fermenergy (thestates[[ii,2]]+thestates[[ii,3]]),2],{dim,dim}];


Get["151008_1_mkham3b3s5.dat"];


hamndintfull=hamInt3b3s5


hamTot=-bosonchempot+fermikinham+hamndintfull/Sqrt[length1 length2];


diag=Eigensystem[N[hamTot]];


mmu=MaxMemoryUsed[]/10.^6


SetDirectory[ParentDirectory[]];


Save["diag3x3s5.dat",{mmu,diag}];
