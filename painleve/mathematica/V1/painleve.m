(* : Title : Painleve Integrabilty Test *)
(* : Author : Douglas Baldwin *)
(* : Summary : 
      This package test a system of PDEs for the Painleve property. *)
(* : Package Version : 1.0 *)
(* : Data : August 8, 2000 *)
(* : Mathematica Version : 4.0 *)
(* : Keywords : Painleve, integrability *)
(* : Warnings : None *)
(* : Limitations : Only deals with cases with positive resonances, as 
   all negative resonances would require the Taylor series rather than 
   the Laurent series. *)
(* : Example : PainleveTest[{D[u[x, t], t] + 6*u[x, t]*D[u[x, t], x] + 
   D[u[x, t], {x, 3}] == 0}, {u[x, t]}, {x, t}, Kruskal->x] *)

BeginPackage["Painleve`"]

PainleveTest::usage = 
"PainleveTest[{eqn1, eqn2, ...}, {funcs}, {vars}, opt___]
test for the Painleve property of a system of nonlinear ODEs and/or PDEs. 
\[NewLine]\[NewLine]
PainleveTest[] brings up a menu of test cases.
\[NewLine]\[NewLine]
The Painleve Test is made up of three steps: 1) the
dominant behavior is determined, 2) the resonances are computed, and 3) the
constants of integration are determined.  For the system of ODEs/PDEs to be
of the Painleve Property all three steps must passed unambiguously. 
If freedom or compatibility conditions exist, a message will be sent
to the user.  
\[NewLine]\[NewLine]
The optional argument are, i) Kruskal->var,
which will use the Kruskal Simplification.  For example Kruskal->x, for
varaibles {x,y,z,t} would assumes g[x,y,z,t]->x-h[y,z,t] and the partial
derivative with respect to x of u[i,j] are zero, ii) WorkLog->\"FileName\",
will write to a file in the current directory the equations, the alphas, the
resonances, the constants, and the net time.
\[NewLine]\[NewLine]
Example:  PainleveTest[{D[u[x, t], t] +
6*u[x, t]*D[u[x, t], x] + D[u[x, t], {x, 3}] == 0}, {u[x, t]}, {x, t},
Kruskal->x, WorkLog->\"worklog.m\"]"

alphaSolver::freedom = "Freedom exists in this system, as `1` are both
dominate powers in expression `2`."

alphaSolver::fail = "This system failed the Painleve test while
determinining the negative integer values of alpha[n]."

alphaSolver::noSoln = "Alpha Solver did not find any negative
integer values for alpha, and will test the values from -3 to -1
to see if the package missed a solution due to inequalities."

initialSolver::fail = "This system failed the Painleve test while
determining the expression(s) for u[i,0]."

resonanceSolver::fail = "This system failed the Painleve test while
determining the resonances of the system."

resonanceSolver::fail2 = "This system failed the Painleve test while
determining the resonances of the system.  As, `1`, are not
non-negative integers."

PainleveTest::symbol = "One of the equations had a power that was of 
type Symbol."

PainleveTest::exp = "All the terms in the equations must be polynomial."

constantsSolver::fail = "This system failed the Painleve test while
determining the constants of the system."

constantsSolver::true = "This system is free at resonance `1`."

Unprotect[alphaSolver, alphaSolve, resonanceSolver, resonanceSolve,
  constantsSolve, PainleveTest];

(* : Title: pdeForm *)
(* : Summary : Displays functions and expressions in standard
  mathematical notation. *)
(* : Authors : Initially done by Tracy Otto and Tony Miller (1995),
   later updated first by U. Goktas, then by Mike Colagrosso, and 
   finally rewritten by Douglas Baldwin *)

pdeForm[expr_] := 
  expr /. 
    {Derivative[n__][u_[a_,b_]][x__] :>
       SequenceForm[
         Subscript[u,
           SequenceForm["(",a,",",b,")"]
         ], 
         Subscript[
           SequenceForm @@ 
             Flatten[
               Table[#[[1]], 
                 {#[[2]]}
               ]& /@
                 Transpose[{{x}, {n}}]
             ]
         ]
       ],
    Derivative[n__][F[a_]][x__] :>
      SequenceForm[
        Subscript[u,a],
        Subscript[
          SequenceForm @@ 
            Flatten[
              Table[#[[1]], {#[[2]]}]& /@
                Transpose[{{x}, {n}}]
            ]
        ]
      ],
    Derivative[n__][F_[a_]][x__] :>
      SequenceForm[
        Subscript[F,a],
        Subscript[
          SequenceForm @@ 
            Flatten[
              Table[#[[1]], {#[[2]]}]& /@
                Transpose[{{x}, {n}}]
           ]
        ]
      ],
    Derivative[n__][g_][x__] :>
      SequenceForm[g,
        Subscript[
          SequenceForm @@ 
            Flatten[
              Table[#[[1]], {#[[2]]}]& /@
                Transpose[{{x}, {n}}]
            ]
          ]
        ],
    F[a_][__] :> Subscript[u,a],
    g[__] :> g, 
    u[a_,b_][__] :>
      Subscript[u,
        SequenceForm["(",a,",",b,")"]
      ],
    u[a_,b_] :> 
      Subscript[u,
        SequenceForm["(",a,",",b,")"]
      ]
    };

(* : Title: toPure *)
(* : Author : Douglas Baldwin *)
(* : Summary : Converts u[i,j][x,y,...]->expr to 
   u[i,j]:>Function[{x,y,z,t,...}, expr]. *)

toPure = (# /. (u_[a_,b_][var__]->temp__):>
  (u[a,b]:>Function[{var}, temp]))&;

(* : Title: alphaSolver: *)
(* : Summary : Determines the dominant behavior of the system. *)
(* : Input : system of equations, the dependent variables, and the
   independent variables. *)
(* : Output : A list of sub-list of alpha solutions, e.g. 
  {{alpha[1]->-2, alpha[2]->-2}, {alpha[1]->-1, alpha[2]->-1}} *)

debug = debug1 = debug2 = debug3 = debug4 = False;

alphaSolver[eqList_List, funcs_List, vars_List] := 
  Block[{newFuncs = Array[F, Length[funcs]], funcsRules, eqnList0, 
    numberOfEquations, replaceRules, eqnList, alphaList0, alphaList, 
    alphaListStructure, perm, comboList, alphaRules, comboList2, 
    alphaSoln2, freeVariable, freeVariableMin, alphaSolnFree, 
    alphaSoln, alphaSoln0, checkList2}, 
    Off[Solve::svars];
    Off[Solve::ifun]; 
    funcsRules = 
      Table[
        Head[funcs[[i]]] -> 
          newFuncs[[i]], 
        {i, Length[funcs]}
      ];
    eqnList0 = 
      (eqList /. 
        (a__ == b__) :> 
          a-b) /. 
      funcsRules;
    numberOfEquations = 
      Length[eqList]; 
    replaceRules = 
      F[n_] :> 
        Function[vars, 
          Evaluate[u0[n]*g[Sequence @@ vars]^alpha[n]]
        ];
    eqnList = 
      If[Head[#]===Plus, 
        List @@ #, 
        #
      ]& /@ 
        ExpandAll[eqnList0 /. replaceRules];
    alphaList0 = 
      Exponent[#, 
        g[Sequence @@ vars]
      ]& /@ 
        #& /@ 
          eqnList;
    If[debug, 
      Print["The possible lead terms are:"];  
      Print[alphaList0]];
    alphaList = 
      Map[
        If[Head[#]===Plus, 
          List @@ #, 
          If[Head[#]===alpha || Head[#]===Times, 
            {#}, 
            #
          ]
        ]&, 
      alphaList0, 
      2];
    If[debug, 
      Print["The lead terms where List has been applied:"];
      Print[alphaList]];
    alphaListStructure = 
      Table[
        Union[# /. 
          {a_Integer, b__} :> {b}& /@ 
            alphaList[[i]] ], 
        {i, numberOfEquations}
      ];
    If[debug, 
      Print["The underlying alphaListStructure is: "];
      Print[alphaListStructure]];
    alphaList = 
      Table[
        Flatten[
          Table[
            Cases[#, 
              Flatten[
                {_, alphaListStructure[[i,j]]}
              ] |
              alphaListStructure[[i,j]] 
            ]& /@ 
              {alphaList[[i]]}, 
            {j, Length[alphaListStructure[[i]] ]}
          ], 
          1
        ], 
        {i, numberOfEquations}
      ];
If[debug, Print["The alphaList sorted by alphaListStructure: "];
  Print[alphaList]];
alphaList = {Min[# /. {a_, ___}:>If[IntegerQ[a], a, 0]& /@ #]}& /@ 
  #& /@ alphaList;
alphaList = Plus @@ #& /@ #& /@ Table[Table[Join[alphaList[[i,j]], 
  alphaListStructure[[i,j]]], {j, Length[alphaList[[i]] ]}], 
  {i, numberOfEquations}];
If[debug, Print["The alphaList is: "];Print[alphaList]];
alphaList = 
  Table[If[Length[ alphaList[[i]] ]===1,
    PadRight[alphaList[[i]], 2],
    alphaList[[i]] ],
  {i, numberOfEquations}];
perm := Permutations[Range[Length[alphaList[[j]]]]];
comboList = 
  Flatten[Table[Table[Partition[Flatten[Union[
    Table[Prepend[Sort[Drop[perm[[i]], 
      -Length[ perm[[1]] ] + k]], j],
    {i, 1, Length[perm]}]]], k+1], 
   {k, 2, Length[alphaList[[j]] ]}], 
  {j, numberOfEquations}], 2];
If[debug1, Print["comboList1: "];Print[comboList]];
alphaRules = 
  Table[Union[Flatten[Table[
    Solve[Equal @@
      Flatten[
        Table[
          alphaList[[comboList[[i, 1]], comboList[[i, k]] ]], 
        {k, 2, Length[comboList[[i]] ]}]], 
      alpha[j]], 
    {i, 1, Length[comboList]}]]], 
   {j, 1, numberOfEquations}];
If[debug1, Print["alphaRules = "];Print[alphaRules]];
comboList2 = 
  Partition[Flatten[
    Fold[Table, 
      Array[beta, numberOfEquations],
      Table[{beta[i], Length[alphaRules[[i]] ]}, 
        {i,numberOfEquations}] ] ], 
    numberOfEquations];
If[debug1, Print["comboList2: "];Print[comboList2]];
alphaSoln = 
  Table[Solve[
    Table[Equal @@ alphaRules[[i, comboList2[[j, i]] ]], 
      {i, numberOfEquations}],  Array[alpha, numberOfEquations]], 
    {j, Length[comboList2]}];
alphaSoln = 
  Union[Table[Flatten[alphaSoln[[i]]], {i, Length[alphaSoln]}]];
If[debug2, Print["The resulting solutions are:"];Print[alphaSoln]];
alphaSoln2 = 
  Complement[alphaSoln, 
    Flatten[Union[
      Table[
       If[NumericQ[ Evaluate[alpha[j] /. alphaSoln[[i]] ]],
          alphaSoln[[i]], Null ],
       {i, Length[alphaSoln]}, {j, numberOfEquations}]], 1]];    
alphaSoln = 
  Complement[alphaSoln, 
    Flatten[Union[
      Table[
        If[IntegerQ[ Evaluate[alpha[j] /. alphaSoln[[i]] ]],
          Null, alphaSoln[[i]] ],
          {i, Length[alphaSoln]}, {j, numberOfEquations}]], 1]];
   alphaSoln = 
     Complement[alphaSoln,
       Flatten[Union[
         Table[
           If[Negative[ Evaluate[alpha[j] /. alphaSoln[[i]] ]],
             Null,alphaSoln[[i]] ],
         {i, Length[alphaSoln]}, {j, numberOfEquations}]], 1]];  
 alphaSoln0 = alphaSoln;
 If[debug2, Print["Of these, we are left with:"]; Print[alphaSoln]];
 checkList2=
   Table[DeleteCases[Union[Flatten[Table[
     Table[If[
       Evaluate[alphaList0[[i,k]] /. alphaSoln[[j]] ]==
         Min[alphaList0[[i]]/.alphaSoln[[j]]],
       alphaList0[[i,k]]],
     {k, Length[alphaList0[[i]]]}],
   {j, Length[alphaSoln]}]]], Null], 
   {i,numberOfEquations}];
 If[debug4, Print["The alpha terms that are of highest order are: "];
   Print[checkList2]];
freeVariable = 
  DeleteCases[Flatten[Table[
    If[Length[Cases[checkList2[[i]], 
      ___Integer+__*alpha[j]|___Integer+alpha[j] ] ]>1 && 
      numberOfEquations>1,
      alpha[j], Null],
  {j, numberOfEquations},
  {i, numberOfEquations}]], Null];
If[Length[freeVariable]!=0,
  freeVariableMin =   
    Min[Table[freeVariable[[j]] /. alphaSoln[[i]], {i, Length[alphaSoln]},
    {j, Length[freeVariable]}] ]];
If[debug4, Print["The free variable(s) is/are: "];Print[freeVariable];
  Print["The minimum value of either being: "];Print[freeVariableMin]];
If[Length[freeVariable]!=0,
  alphaSolnFree =    
      Partition[Flatten[
        Fold[Table, Table[alpha[i]->beta[i], {i, numberOfEquations}], 
            Table[{beta[k],freeVariableMin,-1},{k,numberOfEquations}]]],
          numberOfEquations]];
If[Length[freeVariable]!=0,
  alphaSolnFree =   
    Complement[DeleteCases[Table[If[Plus @@ Flatten[Table[
      If[Length[Cases[alphaList0[[i]] /. alphaSolnFree[[j]], 
      Min[alphaList0[[i]] /. alphaSolnFree[[j]]] ]]>=2,
      1, 0], {i, numberOfEquations}]]>=numberOfEquations,
      alphaSolnFree[[j]], Null], 
      {j, Length[alphaSolnFree]}], Null], alphaSoln]];
Do[If[Length[Cases[checkList2[[i]], 
  ___Integer+__*alpha[j]|___Integer+alpha[j] ] ]>1 && 
    numberOfEquations>1,
  Message[alphaSolver::freedom, Cases[checkList2[[i]],
  ___Integer+__*alpha[j]|___Integer+alpha[j] ], i] ],
  {j, numberOfEquations},
  {i, numberOfEquations}];
alphaSoln = 
    DeleteCases[
      Table[
        If[Plus @@ Flatten[
          Table[If[Length[Cases[alphaList0[[i]] /. alphaSoln[[j]], 
                     Min[alphaList0[[i]] /. alphaSoln[[j]] ] ]]>=2, 1, 0], 
            {i, numberOfEquations}]]>=numberOfEquations, 
          alphaSoln[[j]], Null ], 
      {j, Length[alphaSoln]}], Null];
 If[debug1, Print["alphaSoln at point ZZZZ is: "];Print[alphaSoln]];
 If[Length[alphaSolnFree]!=0,
     alphaSoln = Join[alphaSoln, alphaSolnFree]];
 On[Solve::svars];
 On[Solve::ifun]; 
 If[debug4, Print["AlphaSoln: "];Print[alphaSoln]];
 If[ Length[alphaSoln] === 0, 
   Message[alphaSolver::noSoln];
   alphaSoln = 
     Partition[Flatten[
       Fold[Table, 
         Table[alpha[i]->beta[i], {i, numberOfEquations}],
         Table[{beta[i], -3, -1}, 
           {i,numberOfEquations}] ] ], 
       numberOfEquations];
alphaSoln = 
    DeleteCases[
      Table[
        If[Plus @@ Flatten[
          Table[If[Length[Cases[alphaList0[[i]] /. alphaSoln[[j]], 
                     Min[alphaList0[[i]] /. alphaSoln[[j]] ] ]]>=2, 1, 0], 
            {i, numberOfEquations}]]>=numberOfEquations, 
          alphaSoln[[j]], Null ], 
      {j, Length[alphaSoln]}], Null],
   alphaSoln = alphaSoln];
 If[Length[alphaSoln] == 0, 
   Message[alphaSolver::fail]; Abort[ ], alphaSoln]
];

(* Initial Equation Solver: *)

debug2 = debug3 = False;

initialSolver[eqList_List, funcs_List, vars_List]:=
  Block[{alphaSoln, u0Soln},
    alphaSoln = alphaSolver[eqList, funcs, vars];
     u0Soln = 
       initialSolve[eqList, funcs, vars, #]& /@ alphaSoln;
     u0Soln = DeleteCases[u0Soln /. 
       {{__},{___,u[_,0][__]->0,___}}:>Null, Null, Infinity];

     If[Length[u0Soln] == 0, 
       Message[initialSolver::fail];Abort[ ], 
       Return[u0Soln]]
   ];

initialSolve[eqList_List, funcs_List, vars_List, 
  alphaSoln_List, prnt_:True]:=
Block[{newFuncs = Array[F, Length[funcs]], funcsRules, eqnList0, 
  numberOfEquations, replaceRules, eqnList, uSoln0eqnList, uSoln0},

Off[Solve::svars];

funcsRules = Table[Head[funcs[[i]]]->newFuncs[[i]], {i, Length[funcs]}];
eqnList0 = (eqList /. (a__ == b__):>a-b) /. funcsRules;

numberOfEquations = Length[eqList];

If[prnt, Print["Leading order analysis determines: "];
  Print[alphaSoln /. (alpha[a_]->b_):>(Subscript[\[Alpha], a]->b)]];

replaceRules = F[n_] :> 
  Function[vars, 
    Evaluate[u[n,0][Sequence @@ vars]*g[Sequence @@ vars]^alpha[n]]];
eqnList = (eqnList0 /. replaceRules) /. alphaSoln;
eqnList = 
  Expand[Divide[#, 
     (PolynomialGCD @@ If[Head[#] === Plus, List @@ #, {#}])]] & /@ eqnList;
eqnList0 = eqnList;

uSoln0eqnList = Coefficient[#, g[Sequence @@ vars], 0]& /@ eqnList0;
uSoln0eqnList = 
  Expand[Divide[#, 
    (PolynomialGCD @@ If[Head[#] === Plus, List @@ #, {#}])]] & /@ 
    uSoln0eqnList;
uSoln0eqnList = (# == 0)& /@ uSoln0eqnList;
uSoln0 = Solve[uSoln0eqnList,
  Sort[Table[u[i,0][Sequence @@ vars], {i, numberOfEquations}], Greater]];

On[Solve::svars];

Return[{alphaSoln, #}& /@ uSoln0]];


(* Resonance Finder *)

debug5 = debug6 = False;

resonanceSolver[eqList_List, funcs_List, vars_List, prnt_:True]:=
  Block[{alphaSoln, resSoln},
    alphaSoln = alphaSolver[eqList, funcs, vars];
    resSoln = 
      {#1, 
       #2, 
       resonanceSolve[eqList, funcs, vars, #1, #2, prnt]
      }& @@ #& 
        /@ initialSolve[eqList, funcs, vars, #, prnt]& /@ alphaSoln;
    resSoln = Flatten[resSoln, 1];
    resSoln = DeleteCases[resSoln //. {{__}, {__}, Null}->Null, Null];
    resSoln = DeleteCases[resSoln //. {___, List[], ___}->Null, Null];
    If[Length[resSoln] == 0,
      Message[resonanceSolver::fail];Abort[ ],
      resSoln]
   ];

resonanceSolve[eqList_List, funcs_List, vars_List, alphaSoln_, 
  u0Soln_, prnt_:True]:=
Block[{newFuncs = Array[F, Length[funcs]], funcsRules, eqnList,
  numberOfEquations, replaceRules, (* eqnReplaceRules,termList, *) rList0, 
  rList, rListMin, matrixQ, detQ, resSoln},

funcsRules = Table[Head[funcs[[i]]]->newFuncs[[i]], {i, Length[funcs]}];
eqnList0 = (eqList /. (a__ == b__):>a-b) /. funcsRules;

numberOfEquations = Length[eqList];

replaceRules = 
  Flatten[Join[alphaSoln,
    toPure[{u0Soln}],
    {F[n_] :> Function[vars, 
      Evaluate[u[n,0][Sequence @@ vars]*g[Sequence @@ vars]^alpha[n]+
      u[n,r][Sequence @@ vars]*g[Sequence @@ vars]^(alpha[n]+r)]]}]];
If[debug5, Print["The replacement rules for this function are:"];
  Print[replaceRules]];

eqnList = eqnList0 //. replaceRules;
eqnList = 
  Expand[Divide[#, 
     (PolynomialGCD @@ If[Head[#] === Plus, List @@ #, {#}])]]& 
  /@ eqnList;

matrixQ = Table[Table[
  Coefficient[
    Coefficient[eqnList[[i]], u[j,r][Sequence @@ vars]*g[Sequence @@ vars]^r],
    g[Sequence @@ vars], 0], 
  {i, numberOfEquations}],
  {j, numberOfEquations}];

If[prnt, StylePrint["The matrix for the resonances is:", "Text"];
  Print[MatrixForm[pdeForm /@ matrixQ]]];
detQ = Factor[ExpandAll[Det[matrixQ]]];
If[prnt, 
  StylePrint["The characteristic equation for the resonances is:", "Text"];
  Print[pdeForm[detQ]==0]];
resSoln = Solve[detQ==0, r];
If[FreeQ[resSoln, r->-1], Return[Null]];
resSoln = 
  DeleteCases[
    If[FreeQ[resSoln, (u[_,0][Sequence @@ vars]|Derivative[__][g][__])],
      resSoln, Null], Null];
resSoln = Flatten[{resSoln}];
If[prnt, 
  StylePrint["The resonances are: ", "Text"];
  Print[resSoln];
  If[Length[DeleteCases[Select[resSoln /. (r->a_):>a, Negative], -1]]>0,
   CellPrint[Cell["This case has negative resonances, other than r\[Rule]-1, 
which will be ignored.", "Message"]]]];
(* Removed on 8-8-00 to show which cases were ignored. *)
(* If[Length[resSoln]>=1, resSoln, Return[Null]];
resSoln = If[Length[resSoln] === 
  Length[Cases[#[[2]]& /@ (List @@ #& /@ resSoln), _Integer]],
  resSoln, If[prnt, Message[resonanceSolver::fail2, 
   DeleteCases[resSoln, r->_Integer]]];resSoln];*)
Return[resSoln]
];


(* : Title : Constants of Integration Solver *)
(* : Summary : Completes the last step of the ARS algorithm. *)

debug7 = debug8 = debug9 = debug10 = debug11 = debug12 = False;

PainleveTest[eqList_List, funcs_List, vars_List, opt___] := 
Block[{eqnList, newFuncs = Array[F, Length[funcs]],
    funcsRules, resSoln, uSoln, time},
  If[!FreeQ[eqList, Power[__, _Symbol]], 
    Message[PainleveTest::symbol];Abort[]];
  If[!FreeQ[eqList, Power[E, __]], 
    Message[PainleveTest::exp];Abort[]];
  funcsRules = Table[Head[funcs[[i]]]->newFuncs[[i]], {i, Length[funcs]}];
  eqnList  = (eqList /. (a__==b__):>(a-b)) /. funcsRules;
  time = TimeUsed[ ];
  StylePrint["Painleve analysis of the equation(s):", "Text"];
  Print[TableForm[pdeForm[(# == 0)&/@eqnList]]];
  CellPrint[Cell[
    TextData[{"Using the Laurent series of the solution ",
      Cell[BoxData[FormBox[RowBox[{SubscriptBox["u", "i"],
        ToBoxes[vars] /. {"{"->"(", "}"->")"}}],
        TraditionalForm]]], ":"}], 
    "Text"]];
  Print[SequenceForm[Subscript[u,i],
    StringReplace[ToString[vars], {"{" -> "(", "}" -> ")", ", " -> "," }],
    " = ",Power[g,Subscript[\[Alpha],i]],StringReplace[ToString[vars], 
    {"{"->"(","}"->")",", "->","}], HoldForm[Sum[" ", {k,0,\[Infinity]}]],
    HoldForm[Subscript[u,"(i,k)"]],StringReplace[ToString[vars], 
    {"{"->"(","}"->")",", "->","}],Power[g,k],
    StringReplace[ToString[vars], {"{"->"(","}"->")",", "->","}]]];
  resSoln = resonanceSolver[eqnList, funcs, vars, False];
  StylePrint["Leading order analysis determines "<>ToString[Length[resSoln]]<>
    " possible case(s).  ", 
    "Text"];
  Print[TableForm[
   {#[[1]] /. (alpha[a_]->b_):>(Subscript[\[Alpha],a]->b), 
     pdeForm[#[[2]]], #[[3]]}& /@ resSoln]];
  resSoln = 
    Select[resSoln, 
      (Length[#[[3]]]===Length[Select[#[[3]] /. (r->a_):>a, IntegerQ]] &&
        NonNegative[Max[Select[#[[3]] /. (r->a_):>a, IntegerQ]]])&];
  If[Length[resSoln]===0, Message[resonanceSolver::fail];Abort[]];
  StylePrint["The "<>ToString[Length[resSoln]]<>" case(s) whose resonances
are all integers and whose resonances are all non-negative.  These cases
will now be evaluated one at a time.", "Text"];
  Print[TableForm[
   {#[[1]] /. (alpha[a_]->b_):>(Subscript[\[Alpha],a]->b), 
     pdeForm[#[[2]]], #[[3]]}& /@ resSoln]];
  uSoln = 
    {#[[1]], #[[3]], 
      constantsSolve[eqnList, funcs, vars, #[[1]], #[[2]], #[[3]], opt]}& 
      /@ resSoln;
  time = Round[(TimeUsed[ ] - time)*100]/100.0;
  If[!FreeQ[{opt}||{}, WorkLog->_], Save[WorkLog /. {opt}, time]];
  If[Length[uSoln] == 0, 
    Message[constantsSolver::fail]; Abort[]];
  StylePrint["Painleve analysis of this system determined: ", "Text"];
  (StylePrint["The dominate behaviour is: ","Text"];
     Print[#[[1]] /. alpha[i_]:>Subscript[\[Alpha], i]];
   StylePrint["The resonances are: ", "Text"];Print[#[[2]]];
   StylePrint["The constants of integration are: ", "Text"];
   Print[pdeForm[#[[3]]]])& /@ uSoln; 
  Return[uSoln]];


PainleveTest[eqList_List, funcs_, vars_, opt___] := 
  PainleveTest[eqList, If[Head[funcs]===List, funcs, {funcs}],
    If[Head[vars]===List, vars, {vars}], opt];
  
constantsSolve[eqList_List, funcs_List, vars_List, alphaSoln_List, 
  u0Soln_List, resSoln_List, opt___] := 
Block[{numberOfEquations, maxRes, resList, fail, eqnList, 
        replaceRules, minGPower, uSoln, uSolnMatrix1, 
        uSolnMatrix2, uSolnEqnList, uSoln0eqnList, temp2,
        eqnReplaceRules, termList, minPhiPower, phiPower, 
        eqnSolnList, eqnSolnList0, uSolni, solution, ultUSoln,
        eqnList0, tempUSoln, kruskal, kruskalG, conditions}, 
Off[\[Infinity]::indet, Power::infy];
Unprotect[$MessageList];
Do[$MessageList = Append[$MessageList, HoldForm[Solve::svars]],{3}];

numberOfEquations = Length[eqList];

StylePrint["_______________________________________________________________",
  "Text"];
StylePrint["In this case, leading order analysis determined: ", "Text"];
Print[alphaSoln /. (alpha[a_]->b_):>(Subscript[\[Alpha],a]->b)];

resList = (r /. #)& /@ resSoln;
resList = Cases[resList, _Integer];
maxRes = Max[r /. #& /@ resSoln];
If[!NonNegative[maxRes], 
  StylePrint["This case does not have a positive integer resonance.", "Text"];
  Return[Null]];
If[debug7, Print["For which the maximum integer value is: ", maxRes]];
fail = {};

If[FreeQ[{opt}||{}, Kruskal],
  kruskalG = {},
  kruskalG = 
    {g :> Function[vars, Evaluate[(Kruskal /. {opt})-
       h[Sequence @@ (Complement[vars, {Kruskal /. {opt}}])]]],
     Derivative[var__][u[_,_]][Sequence @@ vars]:>0 /;
       {var}[[Position[vars, (Kruskal /. {opt})][[1,1]]]] != 0}];

uSoln = toPure[u0Soln /. kruskalG];

CellPrint[Cell[
    TextData[{"The values(s) of ",
      Cell[BoxData[FormBox[SubscriptBox["u", "(i,0)"], 
        TraditionalForm]]], " are: "}], 
    "Text"]];
Print[pdeForm[Table[u[i,0] -> (u[i,0] /. uSoln) /. 
  Function[vars, a_]:>a, {i,numberOfEquations}]]];

resonanceSolve[eqList, funcs, vars, alphaSoln, u0Soln];

replaceRules = F[n_] :> 
  Function[vars, 
    Sum[Evaluate[u[n,m][Sequence @@ vars]*g[Sequence @@ vars]^(m+alpha[n])],
      {m, 0, maxRes}]];
If[debug7, Print["The replacement rules for this function are:"];
  Print[replaceRules]];

eqnList = eqList /. replaceRules;
If[debug8, Print["The equation list, eqnList, after replacement is:"];
  Print[eqnList]];
eqnList = eqnList /. alphaSoln;
If[debug8, Print["The equation list after the alpha[i] are fixed is:"];
  Print[eqnList]];

eqnList = 
  Expand[Divide[#, 
     (PolynomialGCD @@ If[Head[#] === Plus, List @@ #, {#}])]] & /@ eqnList;
If[debug8, Print["The equation list after expanding and factoring is:"];
  Print[eqnList]];

eqnList0 = eqnList;

conditions = {};

Do[
uSolnEqnList = Coefficient[#, g[Sequence @@ vars], ii]& /@ eqnList0;
uSolnEqnList = uSolnEqnList //. uSoln;
If[debug12, Print["uSolnEqnList is: "];
  (Print[pdeForm[#]];Print[pdeForm[(#  /. kruskalG)]])&[uSolnEqnList]];
uSolnMatrix1 = 
  Table[Table[
     Expand[Coefficient[uSolnEqnList[[i]], u[j, ii][Sequence @@ vars]]],
       {j, numberOfEquations}],
       {i, numberOfEquations}];
uSolnMatrix2 = 
  Table[-(Plus @@ DeleteCases[
       If[Head[#]===Plus, List @@ #, {#}]&[Expand[uSolnEqnList[[i]] ]], 
       (__*u[_,ii][Sequence @@ vars]|u[_,ii][Sequence @@ vars])]),
       {i, numberOfEquations}];
uSolnMatrix1 = uSolnMatrix1 /. kruskalG;
uSolnMatrix2 = uSolnMatrix2 /. kruskalG;

If[MemberQ[resList, ii],
  StylePrint["The system at resonance level k="<>ToString[ii]<>" is:", "Text"],
  StylePrint["The system at non-resonance level k="<>ToString[ii]<>" is:", 
    "Text"]];
Print[MatrixForm[pdeForm[uSolnMatrix1]],
  MatrixForm[pdeForm[Table[{u[i,ii]}, {i, numberOfEquations}]]],
  " = "];Print[If[FreeQ[(Length[#]>7)&/@#, True], MatrixForm[#], #]&
  [pdeForm[uSolnMatrix2]]];

tempUSoln = 
  If[FreeQ[uSolnMatrix1, Table[0, {numberOfEquations}]],
    Solve[uSolnMatrix1.Table[u[i,ii][Sequence @@ vars], 
      {i, numberOfEquations}] 
    == uSolnMatrix2, Table[u[numberOfEquations-i+1,ii][Sequence @@ vars], 
    {i, numberOfEquations}]],
    If[FreeQ[Table[Append[#1[[i]], -#2[[i]]], 
      {i, numberOfEquations}]&[uSolnMatrix1,uSolnMatrix2], 
      Table[0, {numberOfEquations+1}]], 
      If[numberOfEquations==1,
        fail = Append[fail, ii];
        (uSoln = Flatten[Join[uSoln, 
          Flatten[CellPrint[Cell["Compatibility conditions:", "Message"]];
          StylePrint[pdeForm[#[[2]]], "Message"]; conditions = 
          Flatten[Join[conditions, {#[[2]]}]];#[[1]]]]];
          #[[1]])&[mySolve[
          uSolnMatrix1.Table[u[i,ii][Sequence @@ vars], {i, numberOfEquations}] 
          == uSolnMatrix2, ii]], {}],
      If[numberOfEquations==1,
        Print["The compatibility condition is satisfied."];{},
        Solve[uSolnMatrix1.Table[u[i,ii][Sequence @@ vars], 
        {i, numberOfEquations}]  == uSolnMatrix2, 
        Table[u[numberOfEquations-i+1,ii][Sequence @@ vars], 
        {i, numberOfEquations}]]]]];
tempUSoln = Map[Factor[Expand[#]]&, tempUSoln, Infinity];
If[Length[Flatten[tempUSoln]]>0,
  If[ii===1,
   StylePrint["The solutions at level "<>ToString[ii]<>" is: ", "Text"],
   StylePrint["The solutions at level "<>ToString[ii]<>" are: ", "Text"]];
  Print[pdeForm[tempUSoln]]];
If[debug11, Print["tempUSoln using Solve is:", toPure[tempUSoln]]];

If[FreeQ[resList, ii], 
  If[numberOfEquations === Length[Flatten[tempUSoln]],
    uSoln = Flatten[Join[uSoln, toPure[tempUSoln]]];
    If[debug11, Print["uSoln is: ",uSoln]],
    fail = Append[fail, ii];
    uSoln = Flatten[Join[uSoln, 
      Flatten[(CellPrint[Cell["Compatibility conditions:", "Message"]];
        StylePrint[pdeForm[#[[2]]], "Message"]; conditions = 
        Flatten[Join[conditions, {#[[2]]}]];#[[1]])&[mySolve[
        uSolnMatrix1.Table[u[i,ii][Sequence @@ vars], {i, numberOfEquations}] 
        == uSolnMatrix2, ii]]]]]],
  If[Det[uSolnMatrix1] === 0,
    If[numberOfEquations - Length[Cases[resList, ii]] === 
      Length[Flatten[tempUSoln]],
      uSoln = Flatten[Join[uSoln, toPure[tempUSoln]]];
      temp2 = Length[Cases[resList, ii]];
      If[temp2>1, 
        StylePrint["There are "<>ToString[temp2]<>" free functions.", "Text"],
        If[numberOfEquations != 1, 
        StylePrint["There is one free function.", "Text"]]],
      fail = Append[fail, ii];
      uSoln = Flatten[Join[uSoln, 
        Flatten[(CellPrint[Cell["Compatibility conditions:", "Message"]];
          StylePrint[pdeForm[#[[2]]], "Message"]; conditions = 
          Flatten[Join[conditions, {#[[2]]}]];#[[1]])&[mySolve[
          uSolnMatrix1.Table[u[i,ii][Sequence @@ vars], 
          {i, numberOfEquations}]  == uSolnMatrix2, ii]]]]]],
      fail = Append[fail, ii];
      uSoln = Flatten[Join[uSoln, 
        Flatten[(CellPrint[Cell["Compatibility conditions:", "Message"]];
          StylePrint[pdeForm[#[[2]]], "Message"]; conditions = 
          Flatten[Join[conditions, {#[[2]]}]];#[[1]])&[mySolve[
          uSolnMatrix1.Table[u[i,ii][Sequence @@ vars], 
          {i, numberOfEquations}] == uSolnMatrix2, ii]]]]]]],
  {ii, maxRes}];

ultUSoln = Join[
  Table[u[j,i] -> Evaluate[u[j,i][Sequence @@ vars] /. uSoln], 
    {i,0, maxRes}, {j, numberOfEquations}], {conditions}];
ultUSoln = DeleteCases[ultUSoln, List[]|{}, Infinity];
ii = 0;
ultUSoln = If[Length[ultUSoln]==1, If[Length[ultUSoln[[1]] ]==1, 
  Flatten[ultUSoln], Flatten[ultUSoln, 1]], ultUSoln] /. 
    (u[a_,b_]->u[_,_][Sequence @@ vars]):>If[MemberQ[fail, b],(u[a,b]->$Undetermined),
      (ii++;u[a,b]->C[ii][Sequence @@ vars])];

  Protect[$MessageList];
  On[\[Infinity]::indet, Power::infy];

  equations = eqList;
  alphas = alphaSoln;
  resonances = resSoln;
  constants = ultUSoln;
  If[!FreeQ[{opt}||{}, WorkLog->_], 
    Save[WorkLog /. {opt}, equations, alphas, resonances, 
      constants, ultUSoln]];

  Return[ultUSoln]];

(* : Title : mySolve *)
(* : Summary :
   This ``Solver'' attempts to determine the conditions when the Solve
   function fails to determine the constants of integration. *)

mySolveDebug = False;

mySolve[system_, level_]:=
  Block[{newSystem, numberOfEquations,allVar, param, iVar, 
    iVarSoln, iVarSystem, iVarSystem2},
    newSystem = DeleteCases[system, True];
    newSystem = (#[[1]]-#[[2]])& /@ Transpose[(system /. Equal->List)];
    newSystem = Factor[Expand[newSystem]];
      If[mySolveDebug, Print["The newSystem is: "];Print[newSystem]];
    iVarSoln = {};
    iVarSystem = newSystem;

    iVarSystem = FixedPoint[(
      allVar = Union[DeleteCases[
        DeleteCases[Flatten[# //. {Equal -> List, Times->List, 
        Plus->List, Power->List, Rational->List}], _Integer],_Rational]];
      iVar = Cases[allVar, u[_,level][__]];
      If[Length[iVar]>0 && Length[#]>1,
        iVarSoln = Join[iVarSoln,
          toPure[{Flatten[{(Solve[# == 0, iVar[[1]] ]& /@ 
          DeleteCases[If[FreeQ[#,iVar[[1]]], Null, #]& /@ 
            #, Null])[[1]]}]}]];
        iVarSystem = 
          DeleteCases[(# == 0)& /@ 
          Factor[Expand[Flatten[# /. iVarSoln]]], True];
        (iVarSystem /. Equal[a__,0]:>a),
        iVarSystem])&, newSystem];
     Return[{iVarSoln, 
     DeleteCases[(# == 0)& /@ (iVarSystem /. Equal[a__,b__]:>a-b), 
       True, 
       Infinity]}]
];

PainleveTest[]:=
Block[{choice},
  StylePrint["Welcome to the Painleve Test Package for Mathematica 4.0:\[NewLine]"<>
    "By Douglas Baldwin and Willy Hereman\[NewLine]"<>
    "\[Copyright] 2000, Colorado School of Mines",
    "Subtitle"];
  StylePrint["The built-in test cases are:", "Section"];
  StylePrint["Single Equations:", "Subsection"];
  StylePrint["  1) Fisher Equation;\[NewLine]"<>
    "  2) Cylindrical KdV Equation;\[NewLine]"<>
    "  3) Burgers' Equation;\[NewLine]"<>
    "  4) Korteweg de Vries Equation;\[NewLine]"<>
    "  5) Modified Korteweg de Vries Equation;\[NewLine]"<>
    "  6) Ito Equation;\[NewLine]"<>
    "  7) Boussinesq Equation;\[NewLine]"<>
    "  8) Non-Linear Sine-Gorden Equation\[NewLine]"<>
    "  9) Liouville Non-Linear Klein-Gorden Equation", "Text"];
  StylePrint["Systems of Equations:", "Subsection"];
  StylePrint[" 10) Hirota Satsuma w/ parameter;\[NewLine]"<>
    " 11) Hirota Satsuma w/o parameter;\[NewLine]"<>
    " 12) One-Dimensional Non-Linear Schrodinger Equations;\[NewLine]"<>
    " 13) Reduced Maxwell Bloch Equations;\[NewLine]"<>
    " 14) Lorenz Model of a Dissipative System;\[NewLine]"<>
    " 15) Hu System;\[NewLine]\[NewLine]"<>
    "  x) Exit", "Text"]; 
choice = Input["Enter your choice: "];
Switch[choice,
  x, Abort[],
  1, StylePrint["The Fisher Equation:", "Section"];
     Print[HoldForm[PainleveTest[{D[u[x,t], {x, 2}]+c*D[u[x,t], x]-u[x,t]^2+u[x,t]==0},
       {u[x,t]}, {x,t}, Kruskal->x]]];
     PainleveTest[{D[u[x,t], {x, 2}]+c*D[u[x,t], x]-u[x,t]^2+u[x,t]==0},
       {u[x,t]}, {x,t}, Kruskal->x],
  2, StylePrint["The Cylindrical KdV Equation:", "Section"];
     Print[HoldForm[PainleveTest[{a[t]*D[u[x,t], x]+D[u[x,t], {x, 4}]+
       6*u[x,t]*D[u[x,t], {x, 2}]+6*D[u[x,t], x]^2+
       D[u[x,t], {x, 1}, {t, 1}]}, {u[x,t]}, {x,t}]]];
     PainleveTest[{a[t]*D[u[x,t], x]+D[u[x,t], {x, 4}]+
       6*u[x,t]*D[u[x,t], {x, 2}]+6*D[u[x,t], x]^2+
       D[u[x,t], {x, 1}, {t, 1}]}, {u[x,t]}, {x,t}],
  3, StylePrint["The Burgers' Equation:", "Section"];
     Print[HoldForm[PainleveTest[{D[u[x,t], t]+u[x,t]*D[u[x,t], x]-
       sigma*D[u[x,t], {x,2}]}, {u[x,t]}, {x,t}]]];
     PainleveTest[{D[u[x,t], t]+u[x,t]*D[u[x,t], x]-
       sigma*D[u[x,t], {x,2}]}, {u[x,t]}, {x,t}],
  4, StylePrint["The Korteweg de Vries Equation:", "Section"];
     Print[HoldForm[PainleveTest[{D[u[x,t], t]+ 
       6*u[x,t]*D[u[x,t], x]+D[u[x,t], {x,3}]}, {u[x,t]}, {x,t}]]];
     PainleveTest[{D[u[x,t], t]+ 6*u[x,t]*D[u[x,t], x]+D[u[x,t], {x,3}]},
     {u[x,t]}, {x,t}],
  5, StylePrint["The Modified Korteweg de Vries Equation:", "Section"];
     Print[HoldForm[PainleveTest[{D[u[x,t], t]-3* u[x,t]^2*D[u[x,t], x]+
        2*sigma^2*D[u[x,t], {x,3}]},{u[x,t]}, {x,t}]]];
     PainleveTest[{D[u[x,t], t]-3* u[x,t]^2*D[u[x,t], x]+2*sigma^2*
        D[u[x,t], {x,3}]},{u[x,t]}, {x,t}],
  6, StylePrint["The Ito Equation:", "Section"];
     Print[HoldForm[PainleveTest[{D[u[x,t], t]+a*u[x,t]^2*D[u[x,t], x]+
        b*D[u[x,t], x]*D[u[x,t], {x,2}] + c*u[x,t]*D[u[x,t], {x, 3}]+
        D[u[x,t], {x,5}]} /. {a->2, b->6, c->3}, {u[x,t]}, {x,t},
        Kruskal->t]]];
     PainleveTest[{D[u[x,t], t]+a*u[x,t]^2*D[u[x,t], x]+
        b*D[u[x,t], x]*D[u[x,t], {x,2}] + c*u[x,t]*D[u[x,t], {x, 3}]+
        D[u[x,t], {x,5}]} /. {a->2, b->6, c->3}, {u[x,t]}, {x,t},
        Kruskal->t], 
  7, StylePrint["The Boussinesq Equation:", "Section"];
     Print[HoldForm[PainleveTest[{D[u[x,t], {t,2}]+2*u[x,t]*D[u[x,t], {x,2}]+
        2*D[u[x,t], x]^2+(1/3)*D[u[x,t], {x,4}]}, {u[x,t]}, {x,t}]]];
     PainleveTest[{D[u[x,t], {t,2}]+2*u[x,t]*D[u[x,t], {x,2}]+
        2*D[u[x,t], x]^2+(1/3)*D[u[x,t], {x,4}]}, {u[x,t]}, {x,t}],
  8, StylePrint["The Non-Linear Sine-Gorden Equation:", "Section"];
     Print[HoldForm[PainleveTest[{u[x,t]*(D[u[x,t], {t, 2}]-D[u[x,t], {x,2}])-
       D[u[x,t], t]^2+D[u[x,t], x]^2 - (1/2)*u[x,t]^3+(1/2)*u[x,t]}, 
       {u[x,t]}, {x,t}]]];
     PainleveTest[{u[x,t]*(D[u[x,t], {t, 2}]-D[u[x,t], {x,2}])-
       D[u[x,t], t]^2+D[u[x,t], x]^2 - (1/2)*u[x,t]^3+(1/2)*u[x,t]}, 
       {u[x,t]}, {x,t}],
  9, StylePrint["Liouville Non-Linear Klein-Gorden Equation:", "Section"];
     Print[HoldForm[PainleveTest[{u[x,t]*(D[u[x,t], {t, 2}]-D[u[x,t], {x,2}])-
       D[u[x,t], t]^2+D[u[x,t], x]^2 - u[x,t]^3}, {u[x,t]}, {x,t}]]];
     PainleveTest[{u[x,t]*(D[u[x,t], {t, 2}]-D[u[x,t], {x,2}])-
       D[u[x,t], t]^2+D[u[x,t], x]^2 - u[x,t]^3}, {u[x,t]}, {x,t}],
  10, StylePrint["The Hirota Satsuma w/ parameter:", "Section"];
      Print[HoldForm[PainleveTest[{-aa*D[u[x,t], {x, 3}]-
        6*aa*u[x,t]*D[u[x,t], x]+6*v[x,t]*D[v[x,t], x]+D[u[x,t], t],
        D[v[x,t], {x, 3}]+3*u[x,t]*D[v[x,t], x]+D[v[x,t], t]}, {u[x,t],v[x,t]},
        {x,t}, Kruskal->x]]];
      PainleveTest[{-aa*D[u[x,t], {x, 3}]-6*aa*u[x,t]*D[u[x,t], x]+6*
        v[x,t]*D[v[x,t], x]+D[u[x,t], t], D[v[x,t], {x, 3}]+
        3*u[x,t]*D[v[x,t], x]+D[v[x,t], t]}, {u[x,t],v[x,t]},
        {x,t}, Kruskal->x],
  11, StylePrint["The Hirota Satsuma w/o parameter:", "Section"];
      Print[HoldForm[PainleveTest[{-aa*D[u[x,t], {x, 3}]-
        6*aa*u[x,t]*D[u[x,t], x]+6*v[x,t]*D[v[x,t], x]+D[u[x,t], t],
        D[v[x,t], {x, 3}]+3*u[x,t]*D[v[x,t], x]+D[v[x,t], t]} /. aa->1/2, 
        {u[x,t],v[x,t]}, {x,t}, Kruskal->x]]];
      PainleveTest[{-aa*D[u[x,t], {x, 3}]-6*aa*u[x,t]*D[u[x,t], x]+
        6*v[x,t]*D[v[x,t], x]+D[u[x,t], t], D[v[x,t], {x, 3}]+
        3*u[x,t]*D[v[x,t], x]+D[v[x,t], t]} /. aa->1/2, {u[x,t],v[x,t]},
        {x,t}, Kruskal->x],
  12, StylePrint["The One-Dimensional Non-Linear Schrodinger Equations:", 
        "Section"];
      Print[HoldForm[PainleveTest[{D[u[x,t], t]+D[v[x,t], {x,2}]+
        a*v[x,t]*(u[x,t]^2+v[x,t]^2), D[v[x,t], t]-D[u[x,t], {x,2}]-
        a*u[x,t]*(u[x,t]^2+v[x,t]^2)}, {u[x,t], v[x,t]}, {x,t}, 
        Kruskal->x]]];
      PainleveTest[{D[u[x,t], t]+D[v[x,t], {x,2}]+
        a*v[x,t]*(u[x,t]^2+v[x,t]^2), D[v[x,t], t]-D[u[x,t], {x,2}]-
        a*u[x,t]*(u[x,t]^2+v[x,t]^2)}, {u[x,t], v[x,t]}, {x,t}, Kruskal->x],
  13, StylePrint["The Reduced Maxwell Bloch Equations:", "Section"];
      Print[HoldForm[PainleveTest[{D[u[x,t],x]+mu*v[x,t], D[v[x,t], x]-
        e[x,t]*w[x,t]-mu*u[x,t], D[w[x,t], x]+e[x,t]*v[x,t], D[e[x,t], t]+
        v[x,t]}, {u[x,t],v[x,t],w[x,t],e[x,t]}, {x,t}, Kruskal->x]]];
      PainleveTest[{D[u[x,t],x]+mu*v[x,t], D[v[x,t], x]-e[x,t]*w[x,t]-
        mu*u[x,t], D[w[x,t], x]+e[x,t]*v[x,t], D[e[x,t], t]+v[x,t]}, 
        {u[x,t],v[x,t],w[x,t],e[x,t]}, {x,t}, Kruskal->x],
  14, StylePrint["The Lorenz Model of a Dissipative System:", "Section"];
      Print[HoldForm[PainleveTest[{v[x,t]-D[u[x,t], t], u[x,t]-u[x,t]*w[x,t]-
        D[v[x,t], t], u[x,t]*v[x,t]-D[w[x,t], t]}, {u[x,t], v[x,t], w[x,t]}, 
        {x,t},  Kruskal->x]]];
      PainleveTest[{v[x,t]-D[u[x,t], t], u[x,t]-u[x,t]*w[x,t]-D[v[x,t], t], 
        u[x,t]*v[x,t]-D[w[x,t], t]}, {u[x,t], v[x,t], w[x,t]}, {x,t}, 
        Kruskal->x],
  15, StylePrint["The Hu System:", "SubSection"];
      Print[HoldForm[PainleveTest[{D[u[x,t], t]- D[u[x,t], {x, 2}]- 
        2*u[x,t]*D[u[x,t], x]-D[v[x,t], x], D[v[x,t], t]- D[v[x,t], {x,3}]-
        6*D[u[x,t], x]*D[v[x,t], x]}, {u[x,t], v[x,t]}, {x,t}, 
        Kruskal->x]]];
      PainleveTest[{D[u[x,t], t]- D[u[x,t], {x, 2}]- 2*u[x,t]*D[u[x,t], x]-
        D[v[x,t], x], D[v[x,t], t]- D[v[x,t], {x,3}]-
        6*D[u[x,t], x]*D[v[x,t], x]}, {u[x,t], v[x,t]}, {x,t}, Kruskal->x],
  ___, Abort[]]];

Protect[alphaSolver, alphaSolve, resonanceSolver, resonanceSolve,
  constantsSolve, PainleveTest];

EndPackage[]