(* ************************************************************************* *)

(* Last updated and adjusted for Mathematica version 13.1 *)
(* by Unal Goktas on July 10, 2023 in Houston, TX, USA *)

(* ::Package:: *)

(* :Title: PDESolitonSolutions.m *)

(* :Context: Integrability`PDESolitonSolutions` *)
(* :Author: Unal Goktas and Willy Hereman *)

(* :Summary:
    This package is used for the computation of 1, 2 and 3-soliton
    solutions of partial differential equations (PDEs) of up to
    (3+1) dimensions by transforming the given PDE by a suitable logarithmic
    derivative transformation, into a PDE that is homogeneous of degree
    in a new dependent variable.
*)

(* :Copyright: \[Copyright] 2023 by Unal Goktas and Willy Hereman *)
(* :Package Version: 1.0 *)
(* :Mathematica Version: 13.1 *)
(* :History: Software last updated by Unal Goktas on July 10, 2023 in Houston, TX, USA *)
(*           Previous update by Willy Hereman on June 23, 2023 in Boulder, CO, USA *)
(* :Keywords: soliton solutions *)

(* :Sources:

    W. Hereman and A. Nuseir,
    Symbolic Methods to Construct Exact Solutions of Nonlinear Partial
    Differential Equations,
    Mathematics and Computers in Simulation, 43 (1997) 13--27.

    Willy Hereman and Unal Goktas,
    Symbolic computation of solitary wave solutions and solitons through homogenization of degree,
    Proceedings Conference on Nonlinear and Modern Mathematical Physics (NMMP-2022),
    Springer Proceedings in Mathematics and Statistics, 57pp, submitted (2023).
*)

(* :Warnings: None *)
(* :Limitations: see on line documentation *)
(* :Discussion: see on line documentation *)
(* :Requirements: None *)
(* :Examples: see on line documentation *)

Print["Loading Package PDESolitonSolutions.m of July 10, 2023."]

BeginPackage["Integrability`PDESolitonSolutions`"]

PDESolitonSolutions::usage =
"PDESolitonSolutions[pde, u, vars: {x, ___}, t, R, opts] computes the
N-soliton solutions for N = 1, .., R, with integer R < = 3, of the given
partial differential equation (pde).
PDESolitonSolutions[pde, u, vars: {x, ___}, t, {N}, opts] computes only
the N-soliton solution where N < = 3, of the given partial differential
equation (pde).
PDESolitonSolutions[pde, u, vars: {x, ___}, t, {min, max}, opts] computes
the N-soliton solutions where 1 <= min <= N <= max <= 3, of the given
partial differential equation (pde).
x is understood as the space variable and t as the time variable.
u is understood as the dependent variable of the pde.
Equations with including and up to 3-dimensional space variables are allowed."

BoundForOrderOfLogarithmicDerivative::usage =
"BoundForOrderOfLogarithmicDerivative is an option to be used when finding
the logarithmic transformation u -> K*D[Log[f[x, ___, t], {x, n}] that
homogenizes the partial differential equation with minimal order n where
n = 1, ...., BoundForOrderOfLogarithmicDerivative.
The default is BoundForOrderOfLogarithmicDerivative -> 4."

BoundForDepthOfExpansion::usage =
"BoundForDepthOfExpansion is an option to be used when a solution
of type f -> 1 + Sum[epsilon^n*f[1][x, ___, t], {n, BoundForDepthOfExpansion}]
is seeked for the homogenized equation, after the logarithmic transformation
u -> K*D[Log[f[x, ___, t], {x, n}] that homogenizes the given partial
differential equation, with u as the dependent variable. In the expansion
f -> 1 + Sum[epsilon^n*f[1][x, ___, t], {n, BoundForDepthOfExpansion}],
epsilon serves as a book-keeping parameter (not a small quantity).
The default is BoundForDepthOfExpansion -> 8."

PrintInformation::usage =
"PrintInformation is an option to enable or disable printing of extra
information on the screen while the methods for computing soliton solutions
are applied. It has to be set to a boolean variable.
The default is PrintInformation -> False."

UseTruncatedPainleveExpansion::usage =
"UseTruncatedPainleveExpansion is an option to enable or disable of
using a truncated Painleve expansion (Laurent series solution)
for finding the transformation to homogenize the given PDE.
It has to be set to a boolean variable.
The default is UseTruncatedPainleveExpansion -> True."

k::usage = "k is the default wave number component in the x-axis direction."

l::usage = "l is the default wave number component in the y-axis direction."

m::usage = "m is the default wave number component in the z-axis direction."

omega::usage = "omega is the angular frequency
(gives the dispersion relation when expressed in terms of the wave number)."

delta::usage = "delta is the phase constant."

PDESolitonSolutionsParameters::usage =
"PDESolitonSolutionsConstants is an option to denote the wave number, dispersion and phase constant parameters."

Begin["`Private`"]

Options[PDESolitonSolutions] = {
                                BoundForOrderOfLogarithmicDerivative -> 4,
                                BoundForDepthOfExpansion -> 8,
                                PrintInformation -> False,
                                UseTruncatedPainleveExpansion -> True,
                                PDESolitonSolutionsParameters -> {k, l, m, omega, delta}
                               }

(* for printing a lesser amount of information (useful for debugging)        *)
(* set to true for testing *)
debugprintWH = False;
If[debugprintWH,
   debugprintWH = printinfo,
   debugprintWH = Hold
];

(* for printing a more detailed information of intermediate computations     *)
(* for debugging purposes (while printing includes contexts of the variables)*)
debugprint = False;
If[debugprint,
   debugprint = Print,
   debugprint = Hold
];

(* Main function which gets exported from the package is PDESolitonSolutions *)

PDESolitonSolutions[lhs_ == rhs_, u_/; !ListQ[u], vars:{x_, ___}, t_,
   n_/; (IntegerQ[n] && (1<= n <= 3)), opts___?OptionQ]:=
   With[{result = subPDESolitonSolutions[lhs == rhs, u, vars, t, {1, n}, opts]},
      result /; result =!= $Failed
   ]

PDESolitonSolutions[lhs_ == rhs_, u_/; !ListQ[u], vars:{x_, ___}, t_,
   {n_}/; (IntegerQ[n] && (1 <= n <= 3)), opts___?OptionQ]:=
   With[{result = subPDESolitonSolutions[lhs == rhs, u, vars, t, {n, n}, opts]},
      result /; result =!= $Failed
   ]

PDESolitonSolutions[lhs_ == rhs_, u_/; !ListQ[u], vars:{x_, ___}, t_,
   {min_, max_}/; (IntegerQ[min] && IntegerQ[max] && 1 <= min <= max <= 3), opts___?OptionQ]:=
   With[{result = subPDESolitonSolutions[lhs == rhs, u, vars, t, {min, max}, opts]},
      result /; result =!= $Failed
   ]

subPDESolitonSolutions[lhs_ == rhs_, u_, vars:{x_, ___}, t_, {min_, max_}, opts___] /;
   (
    !FreeQ[lhs-rhs, Derivative[__, q_/; q >= 1][u][Sequence @@ vars, t]] &&
    Length[vars] <= 3 &&
    FreeQ[lhs == rhs, k | l | m | omega | delta] &&
    IntegerQ[BoundForOrderOfLogarithmicDerivative /. {opts} /. Options[PDESolitonSolutions]] &&
    Positive[BoundForOrderOfLogarithmicDerivative /. {opts} /. Options[PDESolitonSolutions]] &&
    IntegerQ[BoundForDepthOfExpansion /. {opts} /. Options[PDESolitonSolutions]] &&
    Positive[BoundForDepthOfExpansion /. {opts} /. Options[PDESolitonSolutions]] &&
    BooleanQ[PrintInformation /. {opts} /. Options[PDESolitonSolutions]] &&
    BooleanQ[UseTruncatedPainleveExpansion /. {opts} /. Options[PDESolitonSolutions]] &&
    ListQ[PDESolitonSolutionsParameters /. {opts} /. Options[PDESolitonSolutions]]
   ):=
   Block[{result},
      If[PrintInformation /. {opts} /. Options[PDESolitonSolutions],
         infoprint = printinfo,
         infoprint = Hold
      ];
      infoprint["Constructing solitary wave or soliton solutions for the PDE: "];
      infoprint[lhs == rhs];
      result = pdeSolitonSolutions[lhs-rhs, u, vars, t, Range[min, max], opts];
      result /; result =!= $Failed
   ] (* end of Block subPDESolitonSolutions *)

subPDESolitonSolutions[___]:= $Failed

printinfo[expr__]:= Print[removeContextAnd$xxx[expr]]

removeContextAnd$xxx[expr__] := Sequence @@ (removeContextAnd$xxx /@ {expr})

removeContextAnd$xxx["expr_"] := expr

removeContextAnd$xxx[expr_] :=
   DeleteCases[
      (
       HoldForm[expr] /.
          With[{vars = Select[
                          Cases[Union[Level[expr, {-1}, Heads -> True]], _Symbol],
                          Not[MatchQ[Context[#], "System`" | "Global`"]]&
                       ]},
             (# -> ToExpression[
                      StringReplace[
                         StringReplace[
                            ToString[#],
                            StringExpression[Longest[__], "`", v__] :> v
                         ],
                         StringExpression[v_, "$", __] -> v
                      ],
                      StandardForm,
                      HoldForm
                   ]&) /@ vars
          ]
      ) /. Verbatim[Pattern][HoldForm[a_], Blank[]] ->
              Verbatim[Pattern][a, Blank[]],
      Verbatim | HoldForm,
      {0, Infinity},
      Heads -> True
   ]

pdeSolitonSolutions[pde_, u_, vars:{x_, ___}, t_, numsolitonsList_, opts___]:=
   Module[{f, K, logderivorder, localhomogeneouspde, savedseed = $RandomGeneratorState},
      logderivorder = BoundForOrderOfLogarithmicDerivative /. {opts} /. Options[PDESolitonSolutions];
      debugprint["logderivorder defined by BoundForOrderOfLogarithmicDerivative is: ", logderivorder];
      localhomogeneouspde = homogenizePDE[pde, u, vars, t, f, K, logderivorder, opts];
      debugprint["localhomogeneouspde after homogenizePDE is:"];
      debugprint[localhomogeneouspde];
      (
       If[Length[localhomogeneouspde[[1]]] == 3, homogenizedviatruncatedPainleve = True, homogenizedviatruncatedPainleve = False];
       localhomogeneouspde = subpdeSolitonSolutions[u, vars, t, numsolitonsList, f, K, #, homogenizedviatruncatedPainleve, opts]& /@ localhomogeneouspde;
       localhomogeneouspde = DeleteCases[localhomogeneouspde, $Failed];
       SeedRandom[savedseed];
       ({#}& /@ Flatten[localhomogeneouspde]) /; (Length[localhomogeneouspde] > 0)
      ) /; localhomogeneouspde =!= $Failed
   ]

pdeSolitonSolutions[___]:= $Failed

subpdeSolitonSolutions[u_, vars:{x_, ___}, t_, numsolitonsList_, f_, K_, localhomogeneouspde_, homogenizedviatruncatedPainleve_, opts___] :=
   Module[{homorder, homogeneouspde, transformation, logderivorder, kval, expansiondepth, k, r, s, w, d, a, phi, phirules, uinphibase, uinphi,
      homogeneouspdebase, homogeneouspdebasenew, homogeneouspdebasenewer, homogeneouspdebasenewest, homogeneouspdebaseultimate, testresult},
      If[homogenizedviatruncatedPainleve,
         {homorder, homogeneouspde, transformation} = localhomogeneouspde;
         logderivorder = 0;
         transformation = Factor[transformation];
         infoprint["The PDE can be homogenized with the transformation"];
         infoprint[u, " = ", transformation];
         infoprint["based on a Laurent series solution."];
         infoprint["The resulting homogeneous equation has degree ", homorder,":"];
         infoprint[homogeneouspde == 0],
         (* else *)
         homogeneouspde = localhomogeneouspde /. (K -> K) -> (K -> 1);
         debugprint["After setting the free constant K to 1: ", homogeneouspde];
         {logderivorder, kval, homorder, homogeneouspde} = homogeneouspde;
         transformation = Factor[kval[[2]]*D[Log[f @@ Flatten[{vars, t}]], {x, logderivorder}]];
         infoprint["The PDE can be homogenized with the transformation"];
         infoprint[u, " = (", kval[[2]],")*D[Log[", f @@ Flatten[{vars, t}], "], {x, ", logderivorder, "}]."];
         infoprint["based on a single-term logarithmic derivative with respect to x."];
         infoprint["The resulting homogeneous equation has degree ", homorder, ":"];
         infoprint[homogeneouspde == 0]
      ];
      expansiondepth = BoundForDepthOfExpansion /. {opts} /. Options[PDESolitonSolutions];
      debugprintWH["expansiondepth defined by BoundForDepthOfExpansion is: ", expansiondepth];
      {k, r, s, w, d} = PDESolitonSolutionsParameters /. {opts} /. Options[PDESolitonSolutions];
      homogeneouspde = computeSolitonSolutions[numsolitonsList, homogeneouspde,
                          f, a, phi, vars, t, {k, r, s, w}, expansiondepth, u, transformation];
      debugprintWH["At pt AA-1, inside pdeSolitonSolutions, after computeSolitonSolutions one gets, homogeneouspde ="];
      debugprintWH[homogeneouspde];
      (* Put below lines lower??? *)
      If[(numsolitonsList === {1} && homogeneouspde =!= $Failed),
         debugprintWH["At pt AA-2, inside if statement:"];
         infoprint["A candidate for a solitary (or one-soliton) solution is "];
         infoprint["f = ", Part[homogeneouspde,1][[1]]]
      ];
      If[((Min[numsolitonsList] === 1) && (Max[numsolitonsList] > 1)),
         debugprintWH["At pt AA-3, inside if statement:"];
         infoprint["A candidate for a solitary (or one-soliton) solution is "];
         infoprint["f = ", Part[homogeneouspde, 1][[1]]]
      ];
      (
       infoprint["While computing these solutions, the replacement rules ",
          DeleteCases[#, (Rule | RuleDelayed)[a[__], 0]]& /@ (Last /@ homogeneouspde)];
       infoprint["and the transformation: "];
       infoprint[u, " = ", transformation];
       infoprint["were used."];
       infoprint["With ", phi[i] @@ Flatten[{vars, t}], " = Exp[",
          Take[{k[i], r[i], s[i]}, Length[vars]].vars-w[i]*t+d[i], "],"];
       infoprint["one gets the following solution(s) for f:"];
       infoprint[(f == #)& /@ (First /@ homogeneouspde)];
       infoprint["corresponding to: "];
       debugprintWH[(u == (transformation /. f -> Function @@ {Flatten[{vars, t}], #}))& /@ (First /@ homogeneouspde)];
       debugprintWH["At pt A, apply phirules"];
       phirules = {
                   Derivative[p__, q_][phi[i_]][Sequence @@ vars, t] :>
                      Times @@ Thread[Take[{k[i], r[i], s[i]}, Length[{p}]]^{p}]*(-w[i])^q*phi[i][Sequence @@ vars, t]
                  };
       debugprintWH["At pt B, phirules = "];
       debugprintWH[phirules];
       debugprintWH["After replacing all the derivatives of the phi[i]: "];
       uinphibase = ((transformation /. f -> Function @@ {Flatten[{vars, t}], #}) /. phirules)& /@ (First /@ homogeneouspde);
       debugprintWH["At pt C, uinphibase = "];
       debugprintWH[uinphibase];
       uinphi = (u == #)& /@ uinphibase;
       debugprintWH["At pt D: uinphi = "];
       infoprint[uinphi];
       If[numsolitonsList === {1},
          (* || numsolitonsList === {2} *)
          infoprint["After factoring the (simplest) solution:"];
          debugprintWH["At pt D-1, in factor routine1, applying MapAll Factor to the simplest solution: "];
          uinphibase = MapAll[Factor, uinphibase];
          debugprintWH["At pt D-2, in factor routine1, uinphibase = ", uinphibase];
          debugprintWH["At pt D-3, in factor routine1, uinphibase = "];
          infoprint[(u == # &) /@ uinphibase]
       ];
       If[((Min[numsolitonsList] === 1) && (Max[numsolitonsList] > 1)),
          debugprintWH["After factoring the simplest solution: "];
          debugprintWH["At pt E-1, in factor routine2, applying MapAll Factor only to the simplest solution: "];
          uinphibase = MapAt[MapAll[Factor, #]&, uinphibase, {{1}}];
          debugprintWH["At pt E-2, in factor routine2, uinphibase = ", uinphibase];
          debugprintWH["At pt E-3, in factor routine2, uinphibase = "];
          debugprintWH[(u == # &) /@ uinphibase]
       ];
       debugprintWH["At pt H, old code, before replacing phi[i] in terms of exponentials, homogeneouspde = "];
       homogeneouspdebase = homogeneouspde;
       debugprintWH[homogeneouspdebase];
       homogeneouspde = ((#1 /. phi[i_][Sequence @@ vars, t] :>
                                    (Exp[Take[{k[i], r[i], s[i]}, Length[vars]].vars-w[i]*t+d[i]] //.
                                       #2)&) @@ # &) /@ homogeneouspde;
       debugprintWH["At pt I, old code, after replacing phi[i] in terms of exponentials, homogeneouspde = "];
       debugprintWH[homogeneouspde];
       debugprint["After returning (from phi) to Exp[k[i]x+___-w[i]*t], " <>
          "and applying the rules for ", w[i], ", a[i, j], a[i, j, p],..., and (if applicable) ", k[i], ", one gets:"];
       debugprint[homogeneouspde];
       infoprint["Substituting phi[i] = Exp[", Take[{k[i], r[i], s[i]}, Length[vars]].vars-w[i]*t+d[i],"]"];
       infoprint["as well as ", w[i], " and (if applicable), ",
          k[i], "..., ", d[i], ", one gets:"];
       debugprintWH["At pt J-1, old code, before the transformation from f to u, homogeneouspde = "];
       debugprintWH[homogeneouspde];
       debugprintWH["At pt J-2: old code, results for f in terms of exponentials:"];
       debugprintWH[(f == # &) /@ homogeneouspde];
       homogeneouspde = (transformation /. f -> Function @@ {Flatten[{vars, t}], #})& /@ homogeneouspde;
       homogeneouspde = DeleteCases[homogeneouspde, 0];
       homogeneouspde = DeleteCases[homogeneouspde, q_/; FreeQ[q, t]];
       debugprintWH["At pt K-1, old code, after the transformation from f to u, homogeneouspde = "];
       debugprintWH[homogeneouspde];
       debugprintWH["At pt K-2, new code, after swapping the f's with u's already (but in terms of phi[i]), homogeneouspdenew = "];
       homogeneouspdebasenew = Table[{Part[uinphibase, k], Part[homogeneouspdebase, k][[2]]}, {k, 1, Length[homogeneouspdebase]}];
       debugprintWH[homogeneouspdebasenew];
       homogeneouspdebasenewer = ((#1 /. phi[i_][Sequence @@ vars, t] :>
                                    (Exp[Take[{k[i], r[i], s[i]}, Length[vars]].vars-w[i]*t+d[i]] //.
                                       #2)&) @@ # &) /@ homogeneouspdebasenew;
       debugprintWH["At pt L-1, new code, after replacing phi[i] in terms of exponentials, homogeneouspdebasenewer = "];
       debugprintWH[homogeneouspdebasenewer];
       homogeneouspdebasenewest =
          Table[
             {
              Part[homogeneouspdebasenewer, k] //. Part[homogeneouspdebase, k][[2]] /.
                 ConditionalExpression[expr_, cond_] :> ConditionalExpression[expr, Simplify[cond]],
              Part[homogeneouspdebase, k][[2]]
             }, {k, 1, Length[homogeneouspdebase]}];
       debugprintWH["At pt L-2, new code, after applying rules for k[i], etc., homogeneouspdebasenewest = "];
       debugprintWH[homogeneouspdebasenewest];
       debugprintWH["Undoing the logarithmic transformation yields:"];
       debugprintWH[homogeneouspde];
       debugprintWH["At pt M-1, old code, Returning (from f) to ", u, ", yields:"];
       debugprintWH[(u == # &) /@ homogeneouspde];
       debugprintWH["At pt M-2, new code, Returning (from f) to ", u, ", yields, homogeneouspdebaseultimate:"];
       homogeneouspdebaseultimate =
          Table[Part[homogeneouspdebasenewest, k][[1]], {k, 1, Length[homogeneouspdebasenewest]}];
       homogeneouspdebaseultimate = DeleteCases[homogeneouspdebaseultimate, 0];
       homogeneouspdebaseultimate = DeleteCases[homogeneouspdebaseultimate, q_/; FreeQ[q, t]];
       infoprint[(u == # &) /@ homogeneouspdebaseultimate];
       If[numsolitonsList === {1},
          debugprintWH["At pt N-1, comparing old and new ways, testresult = "];
          testresult = MapAll[Factor,homogeneouspde-homogeneouspdebaseultimate];
          debugprintWH[testresult]
       ];
       debugprintWH["At pt N-2, renaming homogeneousbaseultimate to homogeneouspde = "];
       homogeneouspde = homogeneouspdebaseultimate;
       debugprintWH[homogeneouspde];
       debugprintWH["At pt 0, numsolitonsList = ", numsolitonsList];
       debugprintWH["At pt P, homogeneouspde = "];
       debugprintWH[homogeneouspde];
       If[(numsolitonsList === {1} || numsolitonsList === {2}),
          infoprint["After factoring the (simplest) solution(s): "];
          debugprintWH["At pt Q, in factor routine3, applying MapAll Factor to the solution(s): "];
          homogeneouspde = MapAll[Factor, homogeneouspde];
          debugprintWH["At pt Q, in factor routine3, homogeneouspde = ", homogeneouspde];
          infoprint[(u == # &) /@ homogeneouspde];
       ];
       If[((Min[numsolitonsList] === 1) && (Max[numsolitonsList] > 1)),
          infoprint["After factoring the (simplest) solution(s): "];
          debugprintWH["At pt R, in factor routine3, applying MapAll Factor only to the simplest solution(s): "];
          homogeneouspde = MapAt[MapAll[Factor, #]&, homogeneouspde, {{1}, {2}}];
          debugprintWH["At pt R, in factor routine3, homogeneouspde = ", homogeneouspde];
          infoprint[(u == # &) /@ homogeneouspde];
       ];
       (* QUESTION: homorder is used in two meanings *)
       (* homorder stands for the order of the homogeneous pde which makes *)
       (* sense but below it also refers to the replacement of homogeneouspde *)
       (* in terms of trigonometric functions which makes no sense. *)
       (* Replace homorder by homogeneouspdetrig or something else? *)
       (* ALSO: homogeneouspde seems to indicate that the result is a pde *)
       (* Replace homogeneouspde by solshomogeneouspde or something else? *)
       debugprint["Calling ExpToTrig on homogeneouspde"];
       homorder = TimeConstrained[
                      ExpToTrig[homogeneouspde],
                      30,
                      $Failed
                  ];
       If[homorder =!= $Failed,
          debugprint["Calling Together on ExpToTrig[homogeneouspde]"];
          homorder = TimeConstrained[
                        Together[
                           homorder,
                           Extension -> Automatic,
                           Trig -> True
                        ],
                        90,
                        $Failed
                     ];
          If[homorder =!= $Failed,
             debugprintWH["Conversion into trigonometric functions yields:"];
             debugprintWH["At pt S, homorder = "];
             debugprintWH[homorder];
             infoprint["In terms of hyperbolic (or other special) functions:"];
             infoprint[(u == # &) /@ homorder];
             If[(numsolitonsList === {1} || numsolitonsList === {2}),
                infoprint["After factoring the (simplest) solution(s): "];
                debugprintWH["At pt T-1, in factor routine4, applying MapAll Factor to the solution: "];
                homorder = MapAll[Factor, homorder];
                debugprintWH["At pt T-2, in factor routine4, homorder = ", homorder];
                infoprint[(u == # &) /@ homorder];
             ];
             If[((Min[numsolitonsList] === 1) && (Max[numsolitonsList] > 1)),
                infoprint["After factoring the (simplest) solution(s): "];
                debugprintWH["At pt U-1, in factor routine5, applying MapAll Factor only to the first two solutions: "];
                homorder = MapAt[MapAll[Factor, #]&, homorder, {{1}, {2}}];
                debugprintWH["At pt U-2, in factor routine5, homorder = ", homorder];
                infoprint[(u == # &) /@ homorder];
             ]
          ],
          debugprint["Trying without ExpToTrig"];
          homorder = TimeConstrained[
                        Together[
                           homogeneouspde
                        ],
                        90,
                        $Failed
                     ]
       ];
       If[homorder =!= $Failed, homogeneouspde = homorder];
       ({u[Sequence @@ vars, t] -> #}& /@ homogeneouspde)
      ) /; homogeneouspde =!= $Failed
   ]

subpdeSolitonSolutions[___]:= $Failed

(* The function homogenizePDE tries to homogenize the given PDE by trying  *)
(* a transformation either by using truncated Painleve expansion or        *)
(* or a transformation of the form u -> K*D[Log[f[x, ___]], {x, n}]        *)

homogenizePDE[pde_, u_, vars:{x_, ___}, t_, f_, K_, logderivorder_, opts___]:=
   Block[{check, expr},
      debugprint["Entered the homogenizePDE function with:" ];
      debugprint[{pde, u, vars, t, f, K, logderivorder}];
      check = checkIfHomogeneousAndItsDegree[pde, {}, u, vars, t];
      debugprint["Checked if the given PDE is already homogenized:"];
      debugprint[check];
      If[Length[check] == 0,
         expr = subhomogenizePDE[pde, u, vars, t, f, K, logderivorder, opts],
         expr = {0, K -> K, check[[2]], pde /. u -> f}
      ];
      expr /; expr =!= $Failed
   ] (* end of homogenizePDE *)

homogenizePDE[___]:= $Failed

subhomogenizePDE[pde_, u_, vars:{x_, ___}, t_, f_, K_, logderivorder_, opts___] /;
   TrueQ[UseTruncatedPainleveExpansion /. {opts} /. Options[PDESolitonSolutions]]:=
   Block[{variables = Flatten[{vars, t}], truncated, integratedpde, expr, min, a, phi, k, r, s, w},
      debugprint["Entering the subhomogenizePDE function via truncated Painleve expansion with:"];
      debugprint[{pde, u, vars, t, f, K, logderivorder}];
      truncated = Quiet[
                     Check[
                        TruncatedPainleveExpansion[{pde}, {u @@ variables}, variables, f],
                        $Failed
                     ]
                  ];
      (
       integratedpde = integratePDERepeatedly[pde, u, vars, t];
       debugprint["integratedpde in subhomogenizePDE:"];
       debugprint[integratedpde];
       (
        debugprint["Attempting to homogenize the PDE with the truncated Painleve expansion"];
        expr = ((integratedpde /. u -> (Function @@ {variables, #}))&) /@ truncated;
        expr = Function[{par}, Times @@ DeleteCases[First /@ FactorList[Numerator[Together[par]]], _?(FreeQ[#, f]&)]] /@ expr;
        expr = checkIfHomogeneousAndItsDegree[#, {}, f, vars, t]& /@ expr;
        integratedpde = DeleteCases[expr, {}];
        If[Length[integratedpde] > 0,
           min = Min[#[[2]]& /@ integratedpde];
           integratedpde = Select[integratedpde, (#[[2]] == min)&];
           integratedpde = Block[{infoprint},
                              Select[integratedpde,
                                 (
                                  computeSolitonSolutions[
                                    {1}, #[[3]], f, Unique[a, Temporary], Unique[phi, Temporary], vars, t, Unique[{k, r, s, w}, Temporary], 2, u
                                  ] =!= $Failed
                                 )&
                              ]
                           ]
        ];
        (
         (*min = Min[LeafCount /@ integratedpde];
         integratedpde = SelectFirst[integratedpde, ((LeafCount[#] == min)&)];*)
         truncated = Extract[truncated, Flatten[Position[expr, #]]& /@ integratedpde];
         MapThread[Flatten[{Take[#1, -2], #2}]&, {integratedpde, truncated}]
        ) /; Length[integratedpde] > 0
       ) /; integratedpde =!= $Failed
      ) /; truncated =!= $Failed
   ] (* end of subhomogenizePDE *)

subhomogenizePDE[pde_, u_, vars:{x_, ___}, t_, f_, K_, logderivorder_, opts___]:=
   Module[{integratedpde, i, expr, monomials, termswithK, solK, min, a, phi, k, r, s, w},
      debugprint["Entering the subhomogenizePDE function with:"];
      debugprint[{pde, u, vars, t, f, K, logderivorder}];
      integratedpde = integratePDERepeatedly[pde, u, vars, t];
      debugprint["integratedpde in subhomogenizePDE:"];
      debugprint[integratedpde];
      (
       Do[
          debugprint["Attempting to homogenize the PDE with the transformation"];
          debugprint["i = ", i, " and ", u[Sequence @@ vars, t], " = ",
             K*D[Log[f[Sequence @@ vars, t]], {x, i}]];
          expr = integratedpde /.
                    u -> Function[Evaluate[{Sequence @@ vars, t}], K*D[Log[f[Sequence @@ vars, t]], {x, i}]];
          expr = Times @@ DeleteCases[First /@ FactorList[Numerator[Together[expr]]], _?(FreeQ[#, f]&)];
          monomials = monomialList[expr, rootMonomialList[expr, {f}]];
          termswithK = DeleteCases[Coefficient[expr, #]& /@ monomials, _?(FreeQ[#, K]&)];
          solK = Union[Solve[Or @@ (#1 == 0 && K != 0 &) /@ termswithK, K]];
          expr = checkIfHomogeneousAndItsDegree[expr, #, f, vars, t]& /@ solK;
          expr = DeleteCases[expr, {}];
          If[Length[expr] > 0,
             min = Min[#[[2]]& /@ expr];
             expr = Select[expr, (#[[2]] == min)&];
             expr = Block[{infoprint},
                       Select[expr,
                          (
                           computeSolitonSolutions[
                             {1}, #[[3]], f, Unique[a, Temporary], Unique[phi, Temporary], vars, t, Unique[{k, r, s, w}, Temporary], 2, u
                           ] =!= $Failed
                          )&
                       ]
                    ];
             If[Length[expr] > 0,
                (*min = Min[LeafCount /@ expr];
                expr = SelectFirst[expr, ((LeafCount[#] == min)&)];*)
                expr = Flatten[{i, #}]& /@ expr;
                Break[]
             ]
          ],
          (* tries upto BoundForOrderOfLogarithmicDerivative the derivative of the logarithmic transformation *)
          {i, logderivorder}
       ];
       expr /; Length[expr] > 0
      ) /; integratedpde =!= $Failed
   ] (* end of subhomogenizePDE *)

subhomogenizePDE[___]:= $Failed

(* To make the homogenization easier, we attempt to integrate the PDE with   *)
(* respect to x as many times as possible. We continue doing this until we   *)
(* no longer can integrate the terms which DO NOT contain a t-derivative.    *)

integratePDERepeatedly[pde_, u_, vars:{x_, ___}, t_]:=
   Block[{xderivterms, tderivterms, integratedpde},
      debugprint["Entering the integratePDERepeatedly function with:"];
      debugprint[{pde, u, vars, x, t}];
      {xderivterms, tderivterms} = xtderivTerms[pde, u, vars, t];
      debugprint["Terms with x- and t-derivatives are separated as:"];
      debugprint[{xderivterms, tderivterms}];
      integratedpde = subintegratePDERepeatedly[xderivterms, tderivterms, x];
      debugprint["Integrated PDE (after integration wrt x) reads:"];
      debugprint[integratedpde];
      integratedpde /; integratedpde =!= $Failed
   ] (* end of integratePDERepeatedly *)

integratePDERepeatedly[___]:= $Failed

subintegratePDERepeatedly[xderivterms_, tderivterms_, x_]/; (!FreeQ[xderivterms, x]):=
   Block[{intxderivterms, integratedpde, numofints},
      intxderivterms = NestWhileList[Integrate[#, x]&, xderivterms, FreeQ[#, Integrate]&];
      numofints = Length[intxderivterms]-2;
      intxderivterms = intxderivterms[[-2]];
      integratedpde = intxderivterms+Nest[Integrate[#, x]&, tderivterms, numofints];
      infoprint["The given PDE can be integrated ", numofints, " time(s) with " <>
         "respect to ", x,"."];
      infoprint["After integration (if possible) the PDE becomes:"];
      infoprint[integratedpde == 0];
      integratedpde
   ] (* end of subintegratePDERepeatedly -- case 1 *)

subintegratePDERepeatedly[xderivterms_, tderivterms_, x_]/; (FreeQ[xderivterms, x] && !FreeQ[tderivterms, x]):=
   Block[{inttderivterms, integratedpde, numofints},
      inttderivterms = NestWhileList[Integrate[#, x]&, tderivterms, FreeQ[#, Integrate]&];
      numofints = Length[inttderivterms]-2;
      inttderivterms = inttderivterms[[-2]];
      integratedpde = inttderivterms + Nest[Integrate[#, x] &, xderivterms, numofints];
      infoprint["The given PDE can be integrated ", numofints, " time(s) with " <>
         "respect to ", x,"."];
      infoprint["After integration (if possible) the PDE becomes:"];
      infoprint[integratedpde];
      integratedpde
   ] (* end of subintegratePDERepeatedly -- case 2 *)

subintegratePDERepeatedly[___]:= $Failed

xtderivTerms[pde_, u_, vars:{x_, ___}, t_]:=
   Block[{rootmonomials, monomials, tmonomials, tderivterms, xderivterms},
      debugprint["Entering the xtderivTerms function with"];
      debugprint[{pde, u, vars, t}];
      rootmonomials = rootMonomialList[pde, {u}];
      monomials = monomialList[pde, rootmonomials];
      tmonomials = Cases[monomials, p_/; (!FreeQ[p, Derivative[__, q_/; q >= 1][u][Sequence @@ vars, t]])];
      tderivterms = Coefficient[pde, #] & /@ tmonomials;
      tderivterms = Plus @@ (tderivterms*tmonomials);
      xderivterms = pde-tderivterms;
      {xderivterms, tderivterms}
   ] (* end of xtderivTerms *)

rootMonomialList[expr_, funcs_List]:=
   Cases[{expr}, q_/; (Not[And @@ (FreeQ[q, #]& /@ funcs)] &&
      Head[q] =!= Power && Head[q] =!= Times && Head[q] =!= Plus), -1]

If[$VersionNumber > 4.,
   monomialList[expr_, rootmonomials_, opts___]:=
      Module[{s},
         If[$VersionNumber > 6.,
            s = GroebnerBasis`DistributedTermsList[expr, rootmonomials, opts],
            s = Internal`DistributedTermsList[expr, rootmonomials, opts]
         ];
         (Times @@ MapThread[Power, {s[[2]], #}] &) /@ (#[[1]] & /@ s[[1]])],
   monomialList = MonomialList
]

checkIfHomogeneousAndItsDegree[expr_, sol_, f_, vars: {x_, ___}, t_]:=
   Block[{solexpr, lambda = Unique[lambda, Temporary], lambdapower, fpower},
      debugprint["Entering the checkIfHomogeneousAndItsDegree function with:"];
      debugprint[{expr, sol, f, vars, t}];
      solexpr = FactorList[(expr /. sol) /.
                   {
                    f[Sequence @@ vars, t] -> lambda*f[Sequence @@ vars, t],
                    Derivative[p__][f][Sequence @@ vars, t] -> lambda*Derivative[p][f][Sequence @@ vars, t]
                   }
                ];
      lambdapower = Cases[solexpr, {lambda, p_} -> p];
      (
       solexpr = DeleteCases[solexpr, {lambda, _}];
       (
        lambdapower = lambdapower[[1]];
        fpower = Cases[solexpr, {f[Sequence @@ vars, t], p_} -> p];
        If[Length[fpower] == 1, fpower = fpower[[1]], fpower = 0];
        lambdapower = lambdapower-fpower;
        solexpr = First /@ solexpr;
        solexpr = DeleteCases[solexpr, q_/; FreeQ[q, f]];
        solexpr = DeleteCases[solexpr, f[Sequence @@ vars, t]];
        solexpr = Times @@ solexpr;
        debugprint["Finishing checkIfHomogeneousAndItsDegree with:"];
        debugprint[{sol, lambdapower, solexpr}];
        {sol, lambdapower, solexpr}
       ) /; FreeQ[solexpr, lambda]
      ) /; Length[lambdapower] == 1
   ] (* end of checkIfHomogeneousAndItsDegree *)

checkIfHomogeneousAndItsDegree[___]:= {}

(* computeSolitonSolutions computes the N-soliton solutions, where N is *)
(* in the list numsolitonsList.                                                 *)

computeSolitonSolutions[numsolitonsList_, homogeneouspde_, f_, a_, phi_, vars:{x_, ___}, t_,
   {k_, r_, s_, w_}, expansiondepth_, u_, transformation_:Undefined]:=
   Block[{epsilon = Unique[epsilon, Temporary], hompde, phirules, epsexp,
      minepsexp = 1, solw, sol, failpos, func},
      (* epsilon serves as a book-keeping parameter *)
      debugprint["Entering the function computeSolitonSolutions with"];
      debugprint[{numsolitonsList, homogeneouspde, f, a, phi, vars, t,
                  {k, r, s, w}, expansiondepth}];
      hompde = homogeneouspde /.
                  f -> Function[
                          Evaluate[{Sequence @@ vars, t}],
                          Evaluate[1+Sum[epsilon^i*f[i][Sequence @@ vars, t], {i, expansiondepth}]]
                       ];
      phirules = {
                  Derivative[p__, q_][phi[i_]][Sequence @@ vars, t] :>
                     Times @@
                        Thread[Take[{k[i], r[i], s[i]}, Length[{p}]]^{p}]*(-w[i])^q*phi[i][Sequence @@ vars, t]
                 };
      debugprint["After the expansion f -> 1+Sum[epsilon^i*f[i][x, ___, t],
         {i, expansiondepth}], where expansiondepth = ", expansiondepth];
      debugprint["The homogeneous equation becomes:"];
      debugprint[hompde];
      infoprint["After substitution of the expansion"];
      infoprint["f = 1+Sum[epsilon^i*", f[i] @@ Flatten[{vars, t}],
                ", {i, ", expansiondepth, "}]"];
      infoprint["the homogenized PDE becomes:"];
      infoprint[hompde == 0];
      epsexp = Exponent[hompde, epsilon];
      debugprint["epsexp in computeSolitonSolutions is: ", epsexp];
      While[
         PossibleZeroQ[Coefficient[hompde, epsilon, minepsexp]],
         minepsexp++
      ];
      debugprint["minepsexp in computeSolitonSolutions is: ", minepsexp];
      debugprintWH["At pt AAA, computation of solw using computeOneSolitonSolution, solw = "];
      solw = computeOneSolitonSolution[hompde, epsilon, epsexp, minepsexp,
                f, a, phi, phirules, vars, t, {k, r, s, w}, expansiondepth, u, transformation];
      debugprint["solw after computeOneSolitonSolution:"];
      debugprint[solw];
      debugprint[solw];
      (
       sol[1] = Last /@ solw;
       debugprint["sol[1] after computeOneSolitonSolution:"];
       debugprint[sol[1]];
       solw = First /@ solw;
       debugprint["solw after computeOneSolitonSolution:"];
       debugprint[solw];
       If[MemberQ[numsolitonsList, 2],
          solw = computeTwoAndThreeSolitonSolution[2, hompde, epsilon, epsexp,
                     minepsexp, f, a, phi, phirules, vars, t, {k, r, s, w}, #,
                     expansiondepth, u, transformation]& /@ solw;
          solw = Cases[{solw}, {{___, (Rule | RuleDelayed)[(k | r | s | w)[_],  _], ___}, _}, -1];
          debugprint["solw after computeTwoAndThreeSolitonSolution:"];
          debugprint[solw];
          failpos = Position[solw, $Failed];
          debugprint["failpos after computeTwoAndThreeSolitonSolution: ", failpos];
          sol[1] = Delete[sol[1], failpos];
          debugprint["sol[1] after computeTwoAndThreeSolitonSolution:"];
          debugprint[sol[1]];
          solw = Delete[solw, failpos];
          debugprint["solw after computeTwoAndThreeSolitonSolution:"];
          debugprint[solw];
          If[Length[solw] == 0, Return[$Failed]];
          sol[2] = Last /@ solw;
          solw = First /@ solw
       ];
       debugprint["sol[2] after computeTwoAndThreeSolitonSolution:"];
       debugprint[sol[2]];
       debugprint["solw after computeTwoAndThreeSolitonSolution:"];
       debugprint[solw];
       If[MemberQ[numsolitonsList, 3],
          solw = computeTwoAndThreeSolitonSolution[3, hompde, epsilon, epsexp,
                     minepsexp, f, a, phi, phirules, vars, t, {k, r, s, w}, #,
                     expansiondepth, u, transformation]& /@ solw;
          solw = Cases[{solw}, {{___, (Rule | RuleDelayed)[(k | r | s | w)[_],  _], ___}, _}, -1];
          debugprint["solw after computeTwoAndThreeSolitonSolution:"];
          debugprint[solw];
          failpos = Position[solw, $Failed];
          debugprint["failpos after computeTwoAndThreeSolitonSolution: ", failpos];
          sol[1] = Delete[sol[1], failpos];
          debugprint["sol[1] after computeTwoAndThreeSolitonSolution:"];
          debugprint[sol[1]];
          sol[2] = Delete[sol[2], failpos];
          debugprint["sol[2] after computeTwoAndThreeSolitonSolution:"];
          debugprint[sol[2]];
          solw = Delete[solw, failpos];
          debugprint["solw after computeTwoAndThreeSolitonSolution:"];
          debugprint[solw];
          If[Length[solw] == 0, Return[$Failed]];
          sol[3] = Last /@ solw;
          solw = First /@ solw
       ];
       debugprint["sol[3] after computeTwoAndThreeSolitonSolution:"];
       debugprint[sol[3]];
       debugprint["solw after computeTwoAndThreeSolitonSolution:"];
       debugprint[solw];
       debugprint["Returning the ", numsolitonsList, "-soliton solutions"];
       Flatten[MapThread[func, {#, solw}]& /@ (sol /@ numsolitonsList)] /. func -> List
      ) /; solw =!= $Failed
   ] (* end of computeSolitonSolutions *)

computeSolitonSolutions[___]:= $Failed

sortSolutionsPrioritizingPsuedoUpValues[sol_] :=
   Sort[sol, ((FreeQ[#1, Pattern] && (!FreeQ[#2, Pattern]))&)]

fromSolutionToGeneralRelation[sol_, unknowns_List, {k_, r_, s_}, params_List]:=
   sortSolutionsPrioritizingPsuedoUpValues[
      (
       sol /. {
               ((pp_/; (MemberQ[unknowns, pp] || MemberQ[{k, r, s}, pp]))[ii__] -> (rhs_?((!FreeQ[#, k | r | s])&))) :>
                  Block[{ruleslhs, rulesrhs, i, j, p, q, o, v, rr, ss, vv},
                     ruleslhs = Union[{ii}];
                     rulesrhs = Thread[ruleslhs -> Take[{i, j, p, q, o, v, rr, ss, vv}, Length[ruleslhs]]];
                     ruleslhs = Thread[ruleslhs -> Take[{i_, j_, p_, q_, o_, v_, rr_, ss_, vv_}, Length[ruleslhs]]];
                     ruleslhs = ((pp[ii] /. ruleslhs) -> (rhs /.
                                                            {
                                                             k[jj_] :> k[jj /. rulesrhs],
                                                             r[jj_] :> r[jj /. rulesrhs],
                                                             s[jj_] :> s[jj /. rulesrhs]
                                                            }
                                                          ));
                     ruleslhs
                  ],
               ((pp_/; MemberQ[params, pp])[ii__] -> rhs_) :>
                  Block[{ruleslhs, rulesrhs, i, j, p, q, o, v, rr, ss, vv},
                     ruleslhs = Union[{ii}];
                     rulesrhs = Thread[ruleslhs -> Take[{i, j, p, q, o, v, rr, ss, vv}, Length[ruleslhs]]];
                     ruleslhs = Thread[ruleslhs -> Take[{i_, j_, p_, q_, o_, v_, rr_, ss_, vv_}, Length[ruleslhs]]];
                     ruleslhs = ((pp[ii] /. ruleslhs) -> (rhs /.
                                                            {
                                                             k[jj_] :> k[jj /. rulesrhs],
                                                             r[jj_] :> r[jj /. rulesrhs],
                                                             s[jj_] :> s[jj /. rulesrhs]
                                                            }
                                                          ));
                     ruleslhs
                  ]
              }
      ) /. Rule -> RuleDelayed
   ]

computeOneSolitonSolution[hompde_, epsilon_, epsexp_, minepsexp_, f_, a_,
   phi_, phirules_, vars:{x_, ___}, t_, {k_, r_, s_, w_}, expansiondepth_, u_, transformation_]:=
   Block[{epscoeff, equation, parconds, solw, sequence},
      debugprintWH["At pt. BB-1, entering computeOneSolitonSolution!"];
      debugprint["Inside computeOneSolitonSolution"];
      epscoeff = Coefficient[hompde, epsilon, minepsexp];
      epscoeff = Expand[
                    epscoeff /.
                       {
                        f[1] -> Function[
                                   Evaluate[{Sequence @@ vars, t}],
                                   Evaluate[phi[1][Sequence @@ vars, t]]
                                ]
                       }
                 ];
      epscoeff = Expand[epscoeff /. phirules];
      equation = Coefficient[
                    epscoeff, phi[1][Sequence @@ vars, t],
                    minepsexp
                 ] == 0;
      parconds = And @@ (
                         ((# != 0) && Element[#, Reals])& /@
                             Union[Cases[{hompde}, q_ /; (FreeQ[q, Alternatives @@ Flatten[{f, vars, t, epsilon}]] && !NumericQ[q]), {-1}]]
                        );
      solw = Quiet[{ToRules[Reduce[equation && parconds && k[1] != 0 && r[1] != 0 && s[1] != 0, w[1]]]}];
      If[!FreeQ[solw, ToRules],
         solw = Quiet[{ToRules[Reduce[equation && k[1] != 0 && r[1] != 0 && s[1] != 0 , w[1]]]}]
      ];
      If[!FreeQ[solw, ToRules],
         solw = Solve[equation, w[1]]
      ];
      sequence = Unique[sequence, Temporary];
      solw = Quiet[
                solw /. (patt:{___, _?(FreeQ[#, w[1] | k[1] | r[1] | s[1]]&) -> _, ___}) :>
                           (
                            sequence @@ Solve[((Equal @@ #)& /@ patt), {w[1], k[1], r[1], s[1]}]
                           )
             ] /. sequence -> Sequence;
      (
       solw = Factor[Union[solw]];
       solw = DeleteCases[solw, {___, _ -> 0, ___}];
       (
        solw = fromSolutionToGeneralRelation[solw, {w}, {k, r, s}, {a}];
        debugprint["solw in computeOneSolitonSolution: ", solw];
        solw = subcomputeOneSolitonSolution[hompde, epsilon, epsexp, minepsexp,
                  f, a, phi, phirules, vars, t, {k, r, s, w}, #,
                  expansiondepth, u, transformation, parconds]& /@ solw;
        solw = DeleteCases[solw, $Failed];
        (
         debugprintWH["At pt. BB-2, leaving computeOneSolitonSolution!"];
         Union[Flatten[solw, 1, List]]
        ) /; Length[solw] > 0
       ) /; (Length[Flatten[solw]] > 0)
      ) /; FreeQ[solw, ToRules]
   ] (* end of computeOneSolitonSolution *)

computeOneSolitonSolution[___]:= $Failed

subcomputeOneSolitonSolution[hompde_, epsilon_, epsexp_, minepsexp_, f_, a_,
   phi_, phirules_, vars:{x_, ___}, t_, {k_, r_, s_, w_}, solw_, expansiondepth_, u_, transformation_, parconds_]:=
   Catch[
      Module[{frules, fterms, epscoeff, zeroQ, sol = solw, unknown, sola, check = True, ruleDelayed, i, j},
         debugprintWH["Inside subcomputeOneSolitonSolution"];
         frules = {
                   f[1] -> Function[
                              Evaluate[{Sequence @@ vars, t}],
                              Evaluate[phi[1][Sequence @@ vars, t]]
                           ]
                  };
         Do[
            debugprintWH["Inside subcomputeOneSolitonSolution: "];
            debugprintWH["At pt 0a, expansiondepth = ", expansiondepth];
            debugprintWH["At pt 0b, j = ", j];
            fterms = Coefficient[hompde, epsilon, j+minepsexp-1] /. frules;
            epscoeff = fterms /. f[_] -> Function[Evaluate[{Sequence @@ vars, t}], 0];
            fterms = Expand[fterms-epscoeff];
            epscoeff = Expand[epscoeff /. phirules];
            epscoeff = Expand[
                          (Coefficient[epscoeff, phi[1][Sequence @@ vars, t], j+minepsexp-1] //. sol) /. ConditionalExpression[expr_, cond_] :> expr
                       ];
            zeroQ = PossibleZeroQ[epscoeff];
            If[!zeroQ,
               epscoeff = Together[epscoeff];
               zeroQ = PossibleZeroQ[epscoeff]
            ];
            If[zeroQ,
               debugprintWH["Inside subcomputeOneSolitonSolution: "];
               debugprintWH["At pt 0c, rhs is zero!"];
               sol = Flatten[{sol, {a @@ ConstantArray[i_, j] -> 0}}];
               frules = Flatten[
                           {
                            frules,
                            f[j] -> Function[Evaluate[{Sequence @@ vars, t}], 0]
                           }
                        ];
               Continue[]
            ];
            unknown = a @@ ConstantArray[1, j];
            frules = Flatten[
                        {
                         frules,
                         f[j] -> Function[
                                    Evaluate[{Sequence @@ vars, t}],
                                    Evaluate[unknown*phi[1][Sequence @@ vars, t]^j]
                                 ]
                        }
                     ];
            fterms = Expand[fterms /. frules];
            fterms = Expand[fterms /. phirules];
            fterms = Together[(Coefficient[fterms, phi[1][Sequence @@ vars, t], j+minepsexp-1] //. sol) /. ConditionalExpression[expr_, cond_] :> expr];
            sola = Factor[Solve[fterms+epscoeff == 0, unknown]];
            debugprintWH["Inside subcomputeOneSolitonSolution: "];
            debugprintWH["At pt 0d, sola = ", sola];
            sola = fromSolutionToGeneralRelation[sola, {w}, {k, r, s}, {a}];
            sol = Flatten[{sol, sola}],
            {j, 2, expansiondepth}
         ];
         debugprintWH["At pt 1, check is ", check];
         infoprint["Verifying that the perturbation scheme terminates. Be patient!"];
         SeedRandom[10000];
         Do[
            epscoeff = (((Coefficient[hompde, epsilon, j+minepsexp-1] /. frules) /. phirules) //. sol) /. ConditionalExpression[expr_, cond_] :> expr;
            check = check &&
                       PossibleZeroQ[
                          Chop[Quiet[
                             epscoeff /. ((# -> RandomReal[{-2, 2}, WorkingPrecision -> 300]&) /@
                                                   {k[1], r[1], s[1], phi[1][Sequence @@ vars, t]})
                          ]]
                       ];
            debugprintWH["Inside subcomputeOneSolitonSolution: "];
            infoprint["Verifying that the perturbation scheme terminates. Be patient!"];
            debugprintWH["At j = ", j, ", at pt 2 (after numerical check), check is ", check];
            If[!check,
               Break[]
            ],
            {j, expansiondepth+1, epsexp-minepsexp+1}
         ];
         If[!check,
            frules = {
                      f[1] -> Function[
                                 Evaluate[{Sequence @@ vars, t}],
                                 Evaluate[phi[1][Sequence @@ vars, t]]
                              ]
                     };
            sol = solw;
            Do[
               debugprintWH["Inside subcomputeOneSolitonSolution: "];
               debugprintWH["At pt 3a, expansiondepth = ", expansiondepth];
               debugprintWH["At pt 3b, j = ", j];
               epscoeff = Coefficient[hompde, epsilon, j+minepsexp-1] /. frules;
               epscoeff = epscoeff /. f[_] -> Function[Evaluate[{Sequence @@ vars, t}], 0];
               epscoeff = Expand[epscoeff /. phirules];
               epscoeff = Expand[
                             (Coefficient[epscoeff, phi[1][Sequence @@ vars, t], j+minepsexp-1] //. sol) /. ConditionalExpression[expr_, cond_] :> expr
                          ];
               zeroQ = PossibleZeroQ[epscoeff];
               If[Not[And @@ zeroQ],
                  epscoeff = Together[Extract[epscoeff, Position[zeroQ, False]]];
                  zeroQ = PossibleZeroQ[epscoeff]
               ];
               If[And @@ zeroQ,
                  debugprintWH["Inside subcomputeOneSolitonSolution: "];
                  debugprintWH["At pt 4 rhs is zero!"];
                  sol = Flatten[{#, a @@ ConstantArray[i_, j] -> 0}]& /@ sol;
                  frules = Flatten[
                              {
                               frules,
                               f[j] -> Function[Evaluate[{Sequence @@ vars, t}], 0]
                              }
                           ];
                  Continue[]
               ];
               sola = Assuming[
                         parconds && k[1] != 0 && r[1] != 0 && s[1] != 0 && Element[k[1], Reals] && Element[r[1], Reals] && Element[s[1], Reals],
                         Simplify[
                            Quiet[
                               Solve[(And @@ Thread[epscoeff == 0]) && parconds && k[1] != 0 && r[1] != 0 && s[1] != 0 &&
                                  Element[k[1], Reals] && Element[r[1], Reals] && Element[s[1], Reals], {k[1], r[1], s[1], w[1]}],
                                  {Solve::svars, Solve::fulldim}
                            ]
                         ]
                      ];
               sola = sola /. (lhs_ -> ConditionalExpression[expr_, cond_]) :> (lhs -> ConditionalExpression[Assuming[cond, Simplify[expr]], cond]);
               sola = Union[sola];
               (*sola = {ToRules[Reduce[sola, {k[1], r[1], s[1]}]]};*)
               If[Length[Flatten[sola]] == 0, Throw[$Failed, "subcomputeOneSolitonSolution"]];
               sola = fromSolutionToGeneralRelation[sola, {w}, {k, r, s}, {a}];
               debugprintWH["Inside subcomputeOneSolitonSolution: "];
               debugprintWH["At pt 4, sola = ", sola];
               sol = Thread[
                        Flatten[
                           {
                            (((#1 //. #2) /. ConditionalExpression[expr_, cond_] :> expr) /.
                               ((Rule | RuleDelayed)[p_, q_] :> (ruleDelayed[p, Together[q]]))) /. ruleDelayed -> RuleDelayed,
                            #2,
                            a[Repeated[i_, {2, expansiondepth}]] -> 0
                           }
                        ]&[sol, #]
                     ]& /@ sola;
               sol = DeleteCases[sol, {___, (Rule | RuleDelayed)[w[_], 0], ___}];
               sol = Extract[sol, Position[((w[1] //. #)&) /@ sol, _?(((!PossibleZeroQ[#]))&), {1}, Heads -> False]];
               If[Length[Flatten[sol]] == 0, Throw[$Failed, "subcomputeOneSolitonSolution"]];
               frules = Flatten[
                           {
                            frules,
                            f[j] -> Function[Evaluate[{Sequence @@ vars, t}], 0]
                           }
                        ],
               {j, 2, expansiondepth}
            ];
            sola = {};
            frules = Flatten[
                        {
                         frules,
                         f[_] -> Function[Evaluate[{Sequence @@ vars, t}], 0]
                        }
                     ];
            SeedRandom[10000];
            Do[
               check = True;
               Do[
                  epscoeff = (((Coefficient[hompde, epsilon, j+minepsexp-1] /. frules) /. phirules) //. sol[[i]]) /. ConditionalExpression[expr_, cond_] :> expr;
                  check = check &&
                             PossibleZeroQ[
                                Chop[Quiet[
                                   epscoeff /. ((# -> RandomReal[{-2, 2}, WorkingPrecision -> 300]&) /@
                                                      {k[1], r[1], s[1], phi[1][Sequence @@ vars, t]})
                             ]]
                          ];
                  debugprintWH["Inside subcomputeOneSolitonSolution: "];
                  debugprintWH["At j = ", j, ", at pt 5, check is ", check];
                  If[!check,
                     Break[]
                  ],
                  {j, expansiondepth+1, epsexp-minepsexp+1}
               ];
               If[check,
                  sola = Append[sola, {sol[[i]], (1+(Sum[f[i][Sequence @@ vars, t], {i, expansiondepth}] /. frules) //. sol[[i]]) /.
                                                    ConditionalExpression[expr_, cond_] :> expr}]
               ],
               {i, Length[sol]}
            ],
            sola = {{sol, (1+(Sum[f[i][Sequence @@ vars, t], {i, expansiondepth}] /. frules) //. sol) /.
                             ConditionalExpression[expr_, cond_] :> expr}}
         ];
         Throw[Union[sola], "subcomputeOneSolitonSolution"]  /; Length[sola] > 0;
         debugprintWH["At pt 1sol pt 1, phirules = "];
         debugprintWH[phirules];
         infoprint["A candidate for a ", numsolitons, "-soliton solution is "];
         infoprint["f = ", monomials];
         infoprint["corresponding to: "];
         debugprintWH[(u == (transformation /. f -> Function @@ {Flatten[{vars, t}], monomials}))];
         debugprintWH["After replacing all derivatives of the phi[i]:"];
         (* infoprint[(u == Factor[(transformation /. f -> Function @@ {Flatten[{vars, t}], monomials}) /. phirules])] *)
         infoprint[(u == (transformation /. f -> Function @@ {Flatten[{vars, t}], monomials}) /. phirules)];
         debugprintWH["At pt CC-2, leaving subcomputeOneSolitonSolution"]
        ], (* end of Block subcomputeOneSolitonSolution *)
      "subcomputeOneSolitonSolution"
   ] (* end of Catch for subcomputeOneSolitonSolution *)

subcomputeOneSolitonSolution[___]:= $Failed

computeTwoAndThreeSolitonSolution[numsolitons_, hompde_, epsilon_, epsexp_,
   minepsexp_, f_, a_, phi_, phirules_, vars:{x_, ___}, t_, {k_, r_, s_, w_},
   solw_, expansiondepth_, u_, transformation_]:=
   Block[{epscoeff, sol},
      debugprintWH["At pt DD-1, entering computeTwoAndThreeSolitonSolution"];
      debugprint["Inside computeTwoAndThreeSolitonSolution"];
      debugprint["solw in computeTwoAndThreeSolitonSolution", solw];
      epscoeff = Coefficient[hompde, epsilon, minepsexp];
      debugprintWH["Inside computeTwoAndThreeSolitonSolution"];
      debugprintWH["At pt 1, epscoeff = ", epscoeff];
      epscoeff = Expand[
                    epscoeff /. {
                                 f[1] -> Function[
                                            Evaluate[{Sequence @@ vars, t}],
                                            Evaluate[Sum[phi[i][Sequence @@ vars, t], {i, numsolitons}]]
                                         ]
                                }
                 ];
      debugprintWH["At pt 2, after replacing f[1], epscoeff = ", epscoeff];
      epscoeff = Expand[epscoeff /. phirules];
      debugprintWH["At pt 3, after applying phirules, epscoeff = ", epscoeff];
      epscoeff = Together[epscoeff //. solw, Extension -> Automatic, Trig -> True];
      debugprintWH["At pt 4, after applying dispersion law, epscoeff = ", epscoeff];
      (
       debugprintWH["solw in computeTwoAndThreeSolitonSolution: ", solw];
       sol = subcomputeTwoAndThreeSolitonSolution[numsolitons, hompde, epsilon,
                epsexp, minepsexp, f, a, phi, phirules, vars, t, {k, r, s, w},
                solw, expansiondepth, u, transformation];
       If[!FreeQ[sol, $Failed],
          sol = subcomputeBiSolitonSolution[numsolitons, hompde, epsilon, f, phi, phirules, vars, t,
                   {k, r, s, w}, solw]
       ];
       (
        sol
       ) /; sol =!= $Failed
      ) /; PossibleZeroQ[epscoeff]
   ] (* end of computeTwoAndThreeSolitonSolution *)

computeTwoAndThreeSolitonSolution[___] := $Failed

subcomputeTwoAndThreeSolitonSolution[numsolitons_, hompde_, epsilon_, epsexp_,
   minepsexp_, f_, a_, phi_, phirules_, vars:{x_, ___}, t_, {k_, r_, s_, w_},
   solw_, expansiondepth_, u_, transformation_]:=
   Module[{frules, sol, ftermssaved, epscoeffsaved, ftermsmonomials, counterSuccessiveZeroQ = 0,
      zeroQ, retryZeroQ, numericalZeroQTest = True, monomials, coeffsrules, monomialscoeffs,
      frulessaved, unknown, monomial, fterms, epscoeff, sola, check = True, i, j},
      debugprintWH["At pt DD-1, entering subcomputeTwoAndThreeSolitonSolution"];
      debugprint["Inside subcomputeTwoAndThreeSolitonSolution:"];
      frules = {
                f[1] -> Function[
                           Evaluate[{Sequence @@ vars, t}],
                           Evaluate[Sum[phi[i][Sequence @@ vars, t], {i, numsolitons}]]
                        ]
               };
      infoprint["Verifying that the perturbation scheme terminates. Be patient!"];
      (*sol = sortSolutionsPrioritizingPsuedoUpValues[DeleteCases[solw, (Rule | RuleDelayed)[a[__], _]]];*)
      sol = sortSolutionsPrioritizingPsuedoUpValues[solw];
      Do[
         debugprintWH["Inside subcomputeTwoAndThreeSolitonSolution:"];
         debugprintWH["Verification of termination of the perturbation scheme can take many iterations. Be patient!"];
         debugprintWH["At pt 1a, j = ", j];
         ftermssaved = Coefficient[hompde, epsilon, j+minepsexp-1] /. frules;
         epscoeffsaved = ftermssaved /. f[_] -> Function[Evaluate[{Sequence @@ vars, t}], 0];
         ftermssaved = Expand[ftermssaved-epscoeffsaved];
         ftermssaved = Expand[ftermssaved /. phirules];
         ftermssaved = Expand[(ftermssaved //. sol) /. ConditionalExpression[expr_, cond_] :> expr];
         debugprintWH["At pt 1b, with j = ", j,", ftermssaved = ",ftermssaved];
         ftermsmonomials = monomialList[
                              ftermssaved,
                              Table[
                                 phi[i][Sequence @@ vars, t],
                                 {i, numsolitons}
                              ]
                           ];
         debugprintWH["Inside subcomputeTwoAndThreeSolitonSolution:"];
         debugprintWH["At pt 1c, with j = ", j,", ftermsmonomials = ", ftermsmonomials];
         If[Length[ftermsmonomials] > 1,
            ftermsmonomials = $Failed;
            Break[],
            ftermsmonomials = ftermsmonomials[[1]]
         ];
         debugprintWH["Inside subcomputeTwoAndThreeSolitonSolution:"];
         debugprintWH["At pt 2a, ftermsmonomials = ", ftermsmonomials];
         debugprintWH["At pt 2b, numericalZeroQTest = ", numericalZeroQTest];
         debugprintWH["At pt 2c, counterSuccessiveZeroQ = ", counterSuccessiveZeroQ];
         If[counterSuccessiveZeroQ < 2,
            epscoeffsaved = epscoeffsaved /. phirules;
            epscoeffsaved = TimeConstrained[Expand[epscoeffsaved], 90, epscoeffsaved];
            epscoeffsaved = (epscoeffsaved //. sol) /. ConditionalExpression[expr_, cond_] :> expr;
            epscoeffsaved = TimeConstrained[Expand[epscoeffsaved], 90, epscoeffsaved];
            zeroQ = PossibleZeroQ[epscoeffsaved];
            retryZeroQ = True;
            If[Not[zeroQ],
               epscoeffsaved = TimeConstrained[Together[epscoeffsaved], 90, retryZeroQ = False; epscoeffsaved];
               If[retryZeroQ, zeroQ = PossibleZeroQ[epscoeffsaved]]
            ];
            If[zeroQ, counterSuccessiveZeroQ++, counterSuccessiveZeroQ = 0],
            SeedRandom[10000];
            If[numericalZeroQTest && (counterSuccessiveZeroQ >= 2),
               debugprintWH["At pt 2d, using numerical zero testing"];
               zeroQ = PossibleZeroQ[
                          Chop[Quiet[
                             (((epscoeffsaved /. phirules) //. sol) /. ConditionalExpression[expr_, cond_] :> expr) /.
                                                   ((# -> RandomReal[{-2, 2}, WorkingPrecision -> 300]&) /@
                                                   {
                                                    k[1], k[2], k[3], r[1], r[2], r[3],
                                                    s[1], s[2], s[3], phi[1][Sequence @@ vars, t],
                                                    phi[2][Sequence @@ vars, t], phi[3][Sequence @@ vars, t]
                                                   })
                          ]]
                       ];
               If[Not[zeroQ], numericalZeroQTest = False]
            ]
         ];
         (* epscoeffsaved = Together[epscoeffsaved, Extension -> Automatic, Trig -> True]; *)
         If[zeroQ,
            debugprint["Inside subcomputeTwoAndThreeSolitonSolution:"];
            debugprint["epscoeffsaved is zero!"];
            frules = Flatten[
                       {
                        frules,
                        f[j] -> Function[Evaluate[{Sequence @@ vars, t}], 0]
                       }
                     ];
            Continue[]
         ];
         (*epscoeffsaved = epscoeffsaved /. ConditionalExpression[expr_, cond_] :> expr;*)
         monomials = Block[{lengthepscoeffsaved, rootvariables = Table[phi[i][Sequence @@ vars, t], {i, numsolitons}], partitioned,
                        partitionedmonomialslist = {}, tempmonomials},
                        If[(j >= 6) && (Head[epscoeffsaved] === Plus) && ((lengthepscoeffsaved = Length[epscoeffsaved]) > 600),
                           partitioned = Partition[Range[lengthepscoeffsaved], UpTo[75]];
                           Do[
                              tempmonomials = monomialList[Part[epscoeffsaved, partitioned[[jjj]]], rootvariables];
                              partitionedmonomialslist = Flatten[{partitionedmonomialslist, tempmonomials}],
                              {jjj, Length[partitioned]}
                           ];
                           Union[partitionedmonomialslist],
                           monomialList[epscoeffsaved, rootvariables]
                        ]
                     ];
         debugprintWH["Inside subcomputeTwoAndThreeSolitonSolution:"];
         monomials = Together[monomials/ftermsmonomials];
         debugprintWH["At pt 2e, monomials = ", monomials];
         ftermsmonomials = And @@ (PolynomialQ[#,
                                     Table[phi[i][Sequence @@ vars, t],
                                        {i, numsolitons}]
                                   ]& /@ monomials);
         debugprintWH["At pt 3, ftermsmonomials = ", ftermsmonomials];
         If[!ftermsmonomials,
            ftermsmonomials = $Failed;
            Break[]
         ];
         coeffsrules = {
                        Times[
                           Power[phi[i_][Sequence @@ vars, t], p_.],
                           Power[phi[j_][Sequence @@ vars, t], q_.],
                           Power[phi[m_][Sequence @@ vars, t], n_.]
                        ] :> Flatten[{Table[i, p], Table[j, q], Table[m, n]}],
                        Times[
                           Power[phi[i_][Sequence @@ vars, t], p_.],
                           Power[phi[j_][Sequence @@ vars, t], q_.]
                        ] :> Flatten[{Table[i, p], Table[j, q]}],
                        Power[phi[i_][Sequence @@ vars, t], q_] :> Table[i, q]
                       };
         (*monomialscoeffs = Sort[monomials /. coeffsrules];*)
         monomialscoeffs = monomials /. coeffsrules;
         frulessaved = frules;
         monomials = Plus @@ (((a @@ # &) /@ monomialscoeffs)*monomials);
         monomialscoeffs = DeleteCases[monomialscoeffs, q_/; FreeQ[(a @@ q) /. sol, a]];
         While[
            debugprint["monomialscoeffs in While:", monomialscoeffs];
            Length[monomialscoeffs] > 0,
            unknown = a @@ monomialscoeffs[[1]];
            monomial = Times @@ Table[phi[monomialscoeffs[[1, i]]][Sequence @@ vars, t],
                                   {i, Length[monomialscoeffs[[1]]]}];
            debugprint["monomial in While:", monomial];
            frules = Flatten[
                        {
                         frulessaved,
                         f[j] -> Function[
                                    Evaluate[{Sequence @@ vars, t}],
                                    Evaluate[unknown*monomial]
                                 ]
                        }
                     ];
            fterms = Expand[ftermssaved /. frules];
            fterms = Expand[fterms /. phirules];
            epscoeff = (Coefficient[epscoeffsaved, monomial] //. sol) /. ConditionalExpression[expr_, cond_] :> expr;
            epscoeff = TimeConstrained[Together[epscoeff], 90, epscoeff];
            fterms = Together[(Coefficient[fterms, monomial] //. sol) /. ConditionalExpression[expr_, cond_] :> expr];
            sola = Solve[fterms+epscoeff == 0, unknown];
            sola = TimeConstrained[Factor[sola], 90, TimeConstrained[Together[sola], 90, sola]];
            sola = fromSolutionToGeneralRelation[sola, {w}, {k, r, s}, {a}];
            sol = Flatten[{sol, sola}];
            monomialscoeffs = DeleteCases[monomialscoeffs, q_/; FreeQ[(a @@ q) /. sola, a]]
         ];
         infoprint["Verifying that the perturbation scheme terminates. Be patient!"];
         frules = Flatten[
                     {
                      frulessaved,
                      f[j] -> Function[
                                 Evaluate[{Sequence @@ vars, t}],
                                 Evaluate[monomials]
                              ]
                     }
                  ],
         {j, 2, expansiondepth}
      ];
      debugprintWH["Inside subcomputeTwoAndThreeSolitonSolution:"];
      debugprintWH["At pt 1-0, check is ", check];
      debugprintWH["frules: ", frules];
      debugprintWH["sol: ", sol];
      (
       infoprint["Verifying that the perturbation scheme terminates. Be patient!"];
       SeedRandom[10000];
       Do[
          debugprintWH["Inside subcomputeTwoAndThreeSolitonSolution:"];
          debugprintWH["Verification of termination of the perturbation scheme can take many iterations. Be patient!"];
          debugprintWH["At pt 2-0, j = ", j];
          epscoeff = (((Coefficient[hompde, epsilon, j+minepsexp-1] /. frules) /. phirules) //. sol) /. ConditionalExpression[expr_, cond_] :> expr;
          check = check &&
                     TimeConstrained[
                        PossibleZeroQ[
                           Chop[Quiet[
                              epscoeff /. ((# -> RandomReal[{-2, 2}, WorkingPrecision -> 300]&) /@
                                                    {
                                                     k[1], k[2], k[3], r[1], r[2], r[3],
                                                     s[1], s[2], s[3], phi[1][Sequence @@ vars, t],
                                                     phi[2][Sequence @@ vars, t], phi[3][Sequence @@ vars, t]
                                                    })
                           ]]
                        ],
                        60,
                        False
                     ];
          debugprintWH["Inside subcomputeTwoAndThreeSolitonSolution:"];
          debugprintWH["For j = ", j, ", at pt 2-1 (after numerical check), check is ", check];
          infoprint["Verifying that the perturbation scheme terminates. Be patient!"];
          If[!check,
             Break[]
          ],
          {j, expansiondepth+1, epsexp-minepsexp+1}
       ];
       debugprintWH["Inside subcomputeTwoAndThreeSolitonSolution:"];
       debugprintWH["At pt 3, after numerical checks, check is ", check];
       (
        monomials = 1+(Sum[f[i][Sequence @@ vars, t], {i, expansiondepth}] /. frules);
        debugprintWH["At pt 4, phirules = "];
        debugprintWH[phirules];
        debugprintWH["At pt CCC-0, numsolitons = "];
        debugprintWH[numsolitons];
        If[MemberQ[numsolitons, 1],
           debugprintWH["At pt CCC-1, inside if statement:"];
           debugprintWH["A candidate (at pt 5) for a solitary (or one-soliton) solution is "];
           debugprintWH["f = ", Part[homogeneouspde,1][[1]]];
           debugprintWH["At pt CCC-2, corresponding to: "];
           debugprintWH[(u == (transformation /. f -> Function @@ {Flatten[{vars, t}], monomials}))];
           debugprintWH["After replacing all derivatives of the phi[i]:"];
           debugprintWH["At pt CCC-3:"];
           debugprintWH[(u == (transformation /. f -> Function @@ {Flatten[{vars, t}], monomials}) /. phirules)]
        ];
        infoprint["A candidate for a ", numsolitons, "-soliton solution is "];
        infoprint["f = ", monomials];
        infoprint["corresponding to: "];
        debugprintWH[(u == (transformation /. f -> Function @@ {Flatten[{vars, t}], monomials}))];
        debugprintWH["After replacing all derivatives of the phi[i]:"];
        infoprint[(u == (transformation /. f -> Function @@ {Flatten[{vars, t}], monomials}) /. phirules)];
        {sol, (monomials //. sol) /. ConditionalExpression[expr_, cond_] :> expr}
       ) /; check
      ) /; ftermsmonomials =!= $Failed
      debugprintWH["At pt DD-2, leaving subcomputeTwoAndThreeSolitonSolution"]
   ] (* end of subcomputeTwoAndThreeSolitonSolution *)

subcomputeTwoAndThreeSolitonSolution[___]:= $Failed

subcomputeBiSolitonSolution[numsolitons_, hompde_, epsilon_, f_, phi_, phirules_, vars:{x_, ___}, t_,
   {k_, r_, s_, w_}, solw_]:=
   Catch[
      Module[{hompdelocal, frules, sol, monomials, unknowns, knownunknowns, conds, parconds, sola, i, j},
         debugprintWH["At pt 1, entering subcomputeBiSolitonSolution"];
         debugprint["Inside subcomputeBiSolitonSolution:"];
         hompdelocal = ((hompde /. epsilon -> 1) /. {f[_?((# >= 2)&)] -> Function[Evaluate[{Sequence @@ vars, t}], 0]});
         debugprintWH["At pt DD-1, entering subcomputeBiSolitonSolution"];
         debugprint["Inside subcomputeBiSolitonSolution:"];
         frules = {
                   f[1] -> Function[
                              Evaluate[{Sequence @@ vars, t}],
                              Evaluate[Sum[phi[i][Sequence @@ vars, t], {i, numsolitons}]]
                           ]
                  };
         infoprint["No ", numsolitons, "-soliton solutions are found. Seeking ", If[numsolitons === 2, "bi-", "tri-"], "soliton solutions!"];
         hompdelocal = hompdelocal /. frules;
         hompdelocal = Expand[hompdelocal /. phirules];
         sol = Select[solw, (!FreeQ[#, k | r | s | w])&];
         hompdelocal = Expand[hompdelocal //. sol];
         hompdelocal = hompdelocal /. ConditionalExpression[expr_, cond_] :> expr;
         monomials = monomialList[
                        hompdelocal,
                        Table[phi[i][Sequence @@ vars, t], {i, numsolitons}]
                     ];
         monomials = Union[Flatten[(Coefficient[hompdelocal, #]& /@ monomials) /. phi[_][__] -> 0]];
         unknowns = Table[{k[i], r[i], s[i]}[[j]], {j, Length[vars]}, {i, numsolitons}];
         knownunknowns = (And @@ (((Equal @@ #)&) /@ Union[Cases[{sol}, patt:((Rule | RuleDelayed)[(k | r | s | w)[_Integer], _]), -1]])) /.
                            ConditionalExpression[expr_, cond_] :> (expr && cond);
         conds = And @@ (((#[[1]] != #[[2]]) &) /@  Flatten[Subsets[#, {2}]& /@ unknowns, 1, List]);
         unknowns = Flatten[unknowns];
         conds = conds && (And @@ ((((# != 0) && Element[#, Reals])&) /@ unknowns));
         parconds = And @@ (
                         ((# != 0) && Element[#, Reals])& /@
                             Union[Cases[{hompde}, q_/; (FreeQ[q, Alternatives @@ Flatten[{f, vars, t, epsilon}]] && !NumericQ[q]), {-1}]]
                        );
         sola = Assuming[
                   parconds && conds,
                   Simplify[
                      Quiet[
                         Solve[(And @@ Thread[monomials == 0]) && parconds && conds && knownunknowns, Flatten[{unknowns, Table[w[i], {i, numsolitons}]}]] /.
                            ConditionalExpression[expr_, cond_] :> ConditionalExpression[expr, Simplify[cond]],
                         {Solve::svars, Solve::fulldim}
                      ]
                   ]
                ];
         If[Length[Flatten[sola]] == 0,
            infoprint["No ", If[numsolitons === 2, "bi-", "tri-"], "soliton solutions are found!"];
            Throw[$Failed, "subcomputeBiSolitonSolution"]
         ];
         sola = sola /. (lhs_ -> ConditionalExpression[expr_, cond_]) :> (lhs -> ConditionalExpression[Assuming[cond, Simplify[expr]], cond]);
         sola = Union[sola];
         sola = fromSolutionToGeneralRelation[sola, {w}, {k, r, s}, {}];
         sol = Select[#, ((!FreeQ[#, k | r | s | w])&)]& /@ (Union[Flatten[({#1 //. #2, #2} /.
                  ((Rule | RuleDelayed)[p_, q_] :> ruleDelayed[p, Together[q]])) /. ruleDelayed -> RuleDelayed]]&[sol, #]& /@ sola);
         sol = DeleteCases[sol, {___, (Rule | RuleDelayed)[w[_], 0], ___}];
         If[Length[Flatten[sol]] == 0,
            infoprint["No ", If[numsolitons === 2, "bi-", "tri-"], "soliton solutions are found!"];
            Throw[$Failed, "subcomputeBiSolitonSolution"]
         ];
         sola = {};
         Do[
            sola = Append[sola, {sol[[i]], 1+(f[1][Sequence @@ vars, t] /. frules) //. sol[[i]]}],
            {i, Length[sol]}
         ];
         (
          infoprint["Successfully found ", If[numsolitons === 2, "bi-", "bi/tri-"], "soliton solutions!"];
          Throw[Union[sola], "subcomputeBiSolitonSolution"]
         ) /; (Length[sola] > 0);
      ], (* end of Block subcomputeOneSolitonSolution *)
      "subcomputeBiSolitonSolution"
   ] (* end of subcomputeBiSolitonSolution *)

subcomputeBiSolitonSolution[___]:= $Failed

(* ************************************************************************* *)
(* Code for truncated Painleve expansion based on PainleveTestV5-2000.m      *)
(* by Douglas Baldwin and Willy Hereman                                      *)
(* modified for the purpose of this package                                  *)
(* ************************************************************************* *)

TruncatedPainleveExpansionSimplify[theSystem_List]:=
   Expand[Numerator[Together[#]]] & /@ theSystem

TruncatedPainleveExpansion[eqns_List, functions_List, variables_List, g_] /;
   (FreeQ[eqns, Power[__, _Symbol]] && FreeQ[eqns, Power[E, __]]) :=
   Block[{equations = TruncatedPainleveExpansionSimplify[eqns],
      myTrackingVariableMax, thedominantBehavior, theinitialConstantsOfIntegration,
      theResonances, theConstantsOfIntegration, alpha, gLowestDegree, r},
      alpha = Unique[alpha, Temporary];
      gLowestDegree = Unique[gLowestDegree, Temporary];
      thedominantBehavior = Internal`DeactivateMessages[
                               dominantBehavior[equations, functions, variables, g],
                               Solve::svars
                            ];
      debugprint["thedominantBehavior: ", thedominantBehavior];
      (
       theinitialConstantsOfIntegration = Internal`DeactivateMessages[
                                             Flatten[
                                                initialConstantsOfIntegration[
                                                   equations, functions, variables, #, g]& /@
                                                   thedominantBehavior, 1],
                                             Solve::svars
                                          ];
       debugprint["theinitialConstantsOfIntegration: ", theinitialConstantsOfIntegration];
       (
        r = Unique[r, Temporary];
        theResonances = Internal`DeactivateMessages[
                           (resonances[equations, functions, variables, #[[1]], #[[2]], g]& /@
                           theinitialConstantsOfIntegration) /. {} :> Sequence[],
                           Solve::svars
                        ];
        debugprint["theResonances: ", theResonances];
        (
         theResonances = (theResonances /. patt:{{alpha[1] -> alphaval_, gLowestDegree[1] -> _}, {u[1, 0][x, t] -> _},
                            {r_ -> -1, res__}} :> If[Abs[alphaval] > Min[r /. {res}], {}, patt]) /. {} :> Sequence[];
         debugprint["theResonances-2: ", theResonances];
         (
          theConstantsOfIntegration = Internal`DeactivateMessages[
                                         constantsOfIntegration[equations, functions, variables, #[[1]], #[[2]], #[[3]], g]& /@
                                            theResonances,
                                         Solve::svars
                                      ];
          With[{vec = Last /@ Last[#], alphaval = Abs[First[#][[1, 2]]]},
             vec.Table[(g @@ variables)^k, {k, 0, Length[vec]-1}]/(g @@ variables)^alphaval]& /@
                (theConstantsOfIntegration /. (gLowestDegree[_] -> _) :> Sequence[])
         ) /; (Length[theResonances] > 0)
        ) /; (Length[theResonances] > 0)
       ) /; (Length[theinitialConstantsOfIntegration] > 0)
      ) /; (thedominantBehavior =!= $Failed)
   ]

TruncatedPainleveExpansion[___] := $Failed

dominantBehavior[equations_List, functions_List, variables_List, g_] :=
   Block[{theSystem, myTrackingVariable, myTrackingVariableMax, alphaList0, alphaList, myAlphaList,
   alphaRules, alphaSoln0, alphaSoln},
      myTrackingVariable = Unique[myTrackingVariable, Temporary];
      {theSystem, myTrackingVariableMax} = dominantBehaviorattachTrackingVariables[equations];
      theSystem = dominantBehaviorGenerateSystem[theSystem, functions, variables, g];
      alphaList0 = dominantBehaviorListFormation[theSystem, functions, variables, myTrackingVariableMax, g];
      alphaList = Flatten[dominantBehaviorSimplification[#]& /@ #]& /@ alphaList0;
      myAlphaList = {};
      alphaList /. alpha[i_Integer] :> (myAlphaList = Append[myAlphaList, alpha[i]]; alpha[i]);
      myAlphaList = Union[myAlphaList];
      alphaSoln0 = dominantBehaviorSolveForAlpha[myAlphaList, alphaList];
      alphaRules = dominantBehaviorRulesSolver[#, myAlphaList]& /@ alphaList;
      alphaSoln = dominantBehaviorSystemCleanUp[alphaList, alphaSoln0, myAlphaList];
      (
       Join[#,
          Table[gLowestDegree[i] -> Min[alphaList[[i]] /. #], {i, Length[equations]}]]& /@ alphaSoln
      ) /; (Length[alphaSoln] != 0)
   ]

dominantBehavior[___] := $Failed

dominantBehaviorattachTrackingVariables[eqns_List] :=
   Block[{equations, i = 0},
      equations = If[Head[#] === Plus,
                     Plus @@ ((myTrackingVariable[++i]*#) & /@ List @@ #),
                     myTrackingVariable[++i]*#
                  ]& /@ Expand[eqns];
      {equations, i}
   ]

dominantBehaviorGenerateSystem[equations_List, functions_List, variables_List, g_] :=
   Block[{ansatzRules },
      ansatzRules = RuleDelayed @@ #& /@
                       Table[{Head[functions[[n]]],
                          Function @@ {variables, u[n, 0]*g[Sequence @@ variables]^alpha[n]}},
                          {n, Length[functions]}];
      MapAll[Expand, equations /. ansatzRules]
   ]

dominantBehaviorListFormation[equations_List, functions_List, variables_List, myTrackingVariableMax_Integer, g_] :=
   Block[{theSystem, alphaList},
      theSystem = (Table[Coefficient[#, myTrackingVariable[i]],
                     {i, 1, myTrackingVariableMax}]& /@ equations) /. 0 :> Sequence[];
      theSystem =
      If[Head[#] === Plus, List @@ #, {#}] & /@ # & /@ MapAll[Expand, theSystem];
      alphaList = (Union[Exponent[#, g[Sequence @@ variables]]] & /@ #) & /@ theSystem;
      alphaList /. {0} :> Sequence[]
   ]

dominantBehaviorSimplification[alphaList0_List] :=
   Block[{alphaList, alphaListStructure },
      alphaList = If[Head[#] === Plus,
                     List @@ #,
                     If[Head[#] === alpha || Head[#] === Times, {#}, {#, 0}]
                  ]& /@ alphaList0;
      alphaListStructure = Union[alphaList /. {a_Integer, b__} :> {b}];
      alphaList = Cases[alphaList, {_, Sequence @@ #} | #]& /@ alphaListStructure;
      alphaList = {Min[# /. {a_, ___} :> If[IntegerQ[a], a, 0]& /@ #]}& /@ alphaList;
      alphaList = (Plus @@ Flatten[#])& /@ Transpose[{alphaList, alphaListStructure}];
      alphaList
   ]

dominantBehaviorSolveForAlpha[myAlphaList_List, alphaList_List] :=
   Block[{alphaRules, alphaSoln0},
      alphaRules = dominantBehaviorRulesSolver[#, myAlphaList]& /@ alphaList;
      alphaSoln0 = dominantBehaviorPowerSolver[alphaRules, myAlphaList];
      alphaRules = dominantBehaviorFixFreeAlpha[alphaSoln0, myAlphaList];
      Union[Join[alphaSoln0, alphaRules]]
   ]

dominantBehaviorRulesSolver[alphaList0_, myAlphaList_] :=
   Block[{alphaList},
      alphaList = Union[dominantBehaviorSimplification[alphaList0]];
      alphaList = Flatten[Map[Thread[Table[#[[1]], {Length[#]-1}] == Drop[#, 1]]&,
                     Table[Drop[alphaList, i], {i, 0, Length[alphaList]-2}]], 1];
      Union[Flatten[Solve[#, myAlphaList]& /@ alphaList]]
   ]

dominantBehaviorPowerSolver[alphaRules0_, myAlphaList_] :=
   Block[{alphaRules = Sort[alphaRules0, (Length[#1] < Length[#2])&]},
      alphaRules = Outer[List, Sequence @@ (alphaRules /. Rule -> Equal)];
      alphaRules = Partition[Flatten[alphaRules], Length[myAlphaList]];
      Union[Sequence @@ Solve[#, myAlphaList] & /@ alphaRules]
   ]

dominantBehaviorFixFreeAlpha[alphaSoln_List, myAlphaList_List] :=
   Block[{alphaSoln0, alphaFree, alphaFreeValues, alphaFixedValues, alphaValues},
      alphaSoln0 = ({Rest /@ (# /. Rule -> List), #}& /@ alphaSoln) /.
                      {a_List, {(_Rule)..}} :> Sequence[] /; Not[And @@ (FreeQ[a, #]& /@ myAlphaList)];
      alphaSoln0 = #[[2]]& /@ alphaSoln0;
      alphaFree = Select[alphaSoln0, Length[#] < Length[myAlphaList]&];
      If[Length[alphaFree] >= 1,
         alphaFreeValues = (myAlphaList /. #)& /@ alphaFree;
         alphaFixedValues = (myAlphaList /. #)& /@ Complement[alphaSoln0, alphaFree];
         If[Length[alphaFixedValues] == 0,
            alphaFixedValues = alphaFreeValues /. alpha[_] :> -3
         ];
         alphaValues = Transpose[{#, Sequence @@
                          Cases[alphaFixedValues, # /. (alpha[i_] :> _)]}]& /@ alphaFreeValues;
         alphaValues = Sequence @@ Thread[If[Head[First[#]] === Integer, First[#],
                          Range[Min[Rest[#], {-3}], -1]]& /@ #]& /@ alphaValues;
         alphaValues = (Rule @@ # & /@ Transpose[{myAlphaList, #}])& /@ alphaValues;
         alphaValues
      ];
      alphaSoln0
   ]

dominantBehaviorSystemCleanUp[alphaList_List, alphaSoln0_List, myAlphaList_List] :=
   Block[{alphaSoln = Union[Sort /@ alphaSoln0]},
      alphaSoln = Transpose[{(alphaList //. #) & /@ alphaSoln, alphaSoln}];
      alphaSoln = alphaSoln /. {a_List, {(_Rule)..}} :>
                                  Sequence[] /; Not[And @@ (FreeQ[a, #]& /@ myAlphaList)];
      alphaSoln = alphaSoln /. {{a_List, b:{(_Rule)..}} :>
                                  b /; And @@ ((Length[Cases[#, Min[#]]] >= 2)& /@ a),
                                  {a_List, b:{(_Rule)..}} :> Sequence[]};
      alphaSoln = alphaSoln //. {a:{(_Rule)..} :> Sequence[] /; (Or @@ (!IntegerQ[#[[2]]]& /@ a)),
                                 a:{(_Rule)..} :> Sequence[] /; (Or @@ (#[[2]] > -1 & /@ a)),
                                 a:{(_Rule)..} :> Sequence[] /; (And @@ (Positive[#[[2]]]& /@ a))};
      alphaSoln
   ]

initialConstantsOfIntegration[equations_List, functions_List, variables_List, alphaSoln_List, g_] :=
   Block[{theSystem, uList, i, uSoln0},
      theSystem = initialConstantsOfIntegrationGenerateSystem[
                     equations, functions, variables, alphaSoln, g];
      theSystem = Table[g[Sequence @@ variables]^(-gLowestDegree[i])*theSystem[[i]] /.
                     alphaSoln, {i, Length[equations]}];
      theSystem = Coefficient[#, g[Sequence @@ variables], 0] & /@ theSystem;
      theSystem = (# == 0) & /@ TruncatedPainleveExpansionSimplify[theSystem];
      uList = Sort[Table[u[i, 0][Sequence @@ variables], {i, Length[equations]}], Greater];
      uSoln0 = Solve[Join[theSystem, (# != 0)& /@ uList], uList];
      uSoln0 = MapAll[Factor, uSoln0];
      If[Length[uSoln0] == 0, Return[{}]];
      {alphaSoln, #}& /@ uSoln0
   ]

initialConstantsOfIntegrationGenerateSystem[equations_List, functions_List,
   variables_List, alphaSoln_List, g_] :=
   Block[{ansatzRules},
      ansatzRules = RuleDelayed @@ # & /@
                       Table[{Head[functions[[n]]],
                          Function @@ {variables, u[n, 0][Sequence @@ variables]*
                          g[Sequence @@ variables]^alpha[n]}}, {n, Length[functions]}];
      MapAll[Expand, (equations /. ansatzRules) /. alphaSoln]
   ]

resonances[equations_List, functions_List, variables_List, thedominantBehavior_List,
   initialConstantsOfIntegration_List, g_] :=
   Block[{theSystem, epsilon, myCoefficient, matrixQ, detQ, resSoln},
      epsilon = Unique[epsilon, Temporary];
      theSystem = resonancesGenerateSystem[equations, functions, variables,
                     thedominantBehavior, initialConstantsOfIntegration, g];
      theSystem = Table[Expand[g[Sequence @@ variables]^(-gLowestDegree[i])*theSystem[[i]] /.
                     thedominantBehavior], {i, Length[equations]}];
      theSystem = TruncatedPainleveExpansionSimplify[theSystem];
      myCoefficient[p_Plus, q_, 0]:= Plus @@ Select[p, FreeQ[#, q]&];
      myCoefficient[p_, q_, 0]:= Plus @@ Select[{p}, FreeQ[#, q]&];
      myCoefficient[p_Plus, q_, r_:1]:= Plus @@ Cases[p, z_. q^r -> z];
      myCoefficient[p_, q_, r_:1]:= Plus @@ Cases[{p}, z_. q^r -> z];
      matrixQ = Table[Table[myCoefficient[
                   myCoefficient[myCoefficient[theSystem[[i]], epsilon[j]],
                      u[j, r][Sequence @@ variables]], g[Sequence @@ variables], r],
                      {i, Length[theSystem]}], {j, Length[theSystem]}];
      detQ = Factor[MapAll[Expand, Det[matrixQ]]];
      resSoln = Flatten[Solve[detQ == 0, r]];
      resSoln = MapAll[Factor, resSoln];
      If[Or[!FreeQ[resSoln, (u[_, 0][Sequence @@ variables] | Derivative[__][g][__])],
         Or @@ ((!IntegerQ[#[[2]]])& /@ resSoln),
         Max[Cases[(r /. #)& /@ resSoln, _Integer]] < 0],
         {}
      ];
      {thedominantBehavior, initialConstantsOfIntegration, resSoln}
   ]

resonancesGenerateSystem[equations_List, functions_List, variables_List, alphaSoln_List,
   initialConstantsOfIntegration_List, g_] :=
   Block[{ansatzRules , n},
      ansatzRules = RuleDelayed @@ # & /@ Table[{Head[functions[[n]]],
                       Function @@ {variables, Evaluate[u[n, 0][Sequence @@ variables]*
                                    g[Sequence @@ variables]^alpha[n]+
                                    epsilon[n]*u[n, r][Sequence @@ variables]*
                                    g[Sequence @@ variables]^(alpha[n]+r)]}},
                                    {n, Length[functions]}];
      MapAll[Expand, ((equations /. ansatzRules) /. alphaSoln) //. initialConstantsOfIntegration]
   ]

constantsOfIntegration[equations_List, functions_List, variables_List,
   thedominantBehavior_List, theInitialConstantsOfIntegration_List, theResonances_List, g_] :=
   Block[{alphaval, theSystem, theConstantsOfIntegration, i, j, pureRules},
      alphaval = alpha[1] /. thedominantBehavior;
      alphaval = Abs[alphaval]-1;
      theSystem = constantsOfIntegrationGenerateSystem[equations, functions,
                     variables, thedominantBehavior, alphaval, g];
      theConstantsOfIntegration = DSolve`DSolveToPureFunction[theInitialConstantsOfIntegration];
      theConstantsOfIntegration = Fold[constantsOfIntegrationNthResonance[#2, theSystem,
                                     variables, #1, g]&, theConstantsOfIntegration, Range[alphaval]];
      theConstantsOfIntegration = Flatten[Table[u[j, i] -> Evaluate[
                                                              u[j, i][Sequence @@ variables] /.
                                                              theConstantsOfIntegration],
                                            {i, 0, alphaval}, {j, Length[equations]}]];
      pureRules = theConstantsOfIntegration /. (c_ -> d_) :> (c -> (Function @@ {variables, d}));
      theConstantsOfIntegration = theConstantsOfIntegration /. Rule[a_, b_] :> Rule[a, b /. pureRules];
      {thedominantBehavior, theResonances, theConstantsOfIntegration}
   ]

constantsOfIntegrationGenerateSystem[equations_List, functions_List,
   variables_List, thedominantBehavior_List, theMaximumResonance_Integer, g_] :=
   Block[{m, n, seriesRules, theSystem},
      seriesRules = Table[functions[[n]] -> Sum[u[n, m][Sequence @@ variables]*
                                               g[Sequence @@ variables]^(m + alpha[n]),
                                               {m, 0, theMaximumResonance}] /. thedominantBehavior,
                                               {n, Length[functions]}];
      seriesRules = seriesRules /. (u_[var__] -> temp__) :> (u :> Function @@ {{var}, temp});
      theSystem = MapAll[Expand, equations /. seriesRules];
      theSystem = Table[Expand[g[Sequence @@ variables]^(-gLowestDegree[m])*theSystem[[m]] /.
                     thedominantBehavior], {m, Length[equations]}];
      TruncatedPainleveExpansionSimplify[theSystem]
   ]

constantsOfIntegrationNthResonance[n_Integer, equations_List,
   variables_List, theConstantsOfIntegration0_List, g_] :=
   Block[{theSystem, rightHandSideOfMatrix, leftHandSideOfMatrix, theConstantsOfIntegration},
      theSystem = Coefficient[#, g[Sequence @@ variables], n]& /@ equations;
      theSystem = Factor[MapAll[Expand, theSystem /. theConstantsOfIntegration0]];
      {rightHandSideOfMatrix, leftHandSideOfMatrix} = constantsOfIntegrationNthResonanceBuildMatrices[
                                                         n, theSystem, variables];
      theConstantsOfIntegration = constantsOfIntegrationNthResonanceSolve[n, rightHandSideOfMatrix,
                                     leftHandSideOfMatrix, variables];
      Flatten[Join[theConstantsOfIntegration0, theConstantsOfIntegration]]
   ]

constantsOfIntegrationNthResonanceBuildMatrices[n_Integer, theSystem_List, variables_List] :=
   Block[{i, j, rightHandSide, leftHandSide},
      rightHandSide = Table[Table[Expand[Coefficient[theSystem[[i]],
                         u[j, n][Sequence @@ variables]]], {j, Length[theSystem]}], {i,
                         Length[theSystem]}];
      leftHandSide = -(MapAll[Expand, theSystem] /. u[_, n] :> Function[Evaluate[variables], 0]);
      {rightHandSide, leftHandSide}
   ]

constantsOfIntegrationNthResonanceSolve[n_Integer, rightHandSideOfMatrix_List,
   leftHandSideOfMatrix_List, variables_List] :=
   Block[{theSystem, uList, theSolution},
      {theSystem, uList} = constantsOfIntegrationNthResonanceSolveBuildSystem[
                              n, rightHandSideOfMatrix, leftHandSideOfMatrix, variables];
      theSolution = Solve[theSystem, uList];
      theSolution = MapAll[Factor, theSolution];
      {DSolve`DSolveToPureFunction[MapAll[Factor, theSolution]], {}}
   ]

constantsOfIntegrationNthResonanceSolveBuildSystem[n_Integer, rightHandSideOfMatrix_List,
   leftHandSideOfMatrix_List, variables_List] :=
   Block[{i, theSystem, uList},
      theSystem = rightHandSideOfMatrix.Table[u[i, n][Sequence @@ variables],
                     {i, Length[rightHandSideOfMatrix]}] == leftHandSideOfMatrix;
      uList = Reverse[Table[u[i, n][Sequence @@ variables], {i, Length[rightHandSideOfMatrix]}]];
      {theSystem, uList}
   ]

End[]

SetAttributes[PDESolitonSolutions, ReadProtected]

Protect[
   {
    PDESolitonSolutions,
    BoundForOrderOfLogarithmicDerivative,
    BoundForDepthOfExpansionWeightedParameters,
    PrintInformation,
    PDESolitonSolutionsParameters
   }
];

EndPackage[]

Print["Package PDESolitonSolutions.m of July 10, 2023 is successfully loaded."]

(* ************************************************************************* *)