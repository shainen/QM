(* ::Package:: *)

makeStates=Function[{length1,length2,startnumferms},
(*Block[{},*)
momentum[fs_?MatrixQ]:={Mod[Total[fs\[Transpose]].Range[0,length1-1],length1],Mod[Total[fs].Range[0,length2-1],length2]};
possiblenumparticles=
{#,startnumferms-#,startnumferms-#}&/@Range[0,startnumferms];
possiblemomenta={Select[Tuples[Range[0,length1-1],3],Mod[Total[#],length1]==0&],Select[Tuples[Range[0,length2-1],3],Mod[Total[#],length2]==0&]};
possibleposany[list_,num_,mom_]:=Select[ArrayReshape[#,{length1,length2}]&/@Permutations[list~Join~Table[0,{length1 length2-num}]],momentum[#]==mom&];
possibleposf[num_,mom_]:=possibleposany[Table[1,{num}],num,mom];
possibleposb[num_,mom_]:=Flatten[possibleposany[#,num,mom]&/@DeleteDuplicates[Sort/@Select[Tuples[Range[0,num],num],Total[#]==num&]],1];
makestate[part_?VectorQ,moms1_?VectorQ,moms2_?VectorQ]:=Tuples[{possibleposb[part[[1]],{moms1[[1]],moms2[[1]]}],possibleposf[part[[2]],{moms1[[2]],moms2[[2]]}],possibleposf[part[[3]],{moms1[[3]],moms2[[3]]}]}];
thestatesnotflat={};
Do[
Do[
Do[
AppendTo[thestatesnotflat,makestate[possiblenumparticles[[pp]],possiblemomenta[[1,m1]],possiblemomenta[[2,m2]]]];
,{m2,Length[possiblemomenta[[2]]]}]
,{m1,Length[possiblemomenta[[1]]]}]
,{pp,Length[possiblenumparticles]}];
thestates=Flatten[thestatesnotflat,1];
thestates
];
