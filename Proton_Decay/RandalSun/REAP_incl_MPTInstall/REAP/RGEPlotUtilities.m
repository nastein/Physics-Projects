(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)



BeginPackage["REAP`RGEPlotUtilities`"];

$TextStyle = {FontFamily -> "Times", FontSize -> 10, SingleLetterItalics -> True};
SetOptions[Plot, Axes -> False, Frame -> True, PlotLabel -> ""];

ClearAll[RGELogTicks, RGELogTicksLabeled, RGELogTicksLabeledNegExp, RGEShadowEFT];

RGELogTicks::usage="LogTicks[min,max] returns a list of tickmarks between min
and max.";

RGELogTicksLabeled::usage="LogTicksLabeled[min,max] returns a list of labeled
tickmarks between min and max.";

RGELogTicksLabeledNegExp::usage="LogTicksLabeledNegExp[min,max] returns a list
of labeled tickmarks with negative exponent between min and max.";

RGEShadowEFT::usage="ShadowEFT[list of scales] produces from the given list of scales a list of rectangles with different gray levels.";

Begin["`Private`"];

RGELogTicks[xmin_?NumericQ, xmax_?NumericQ] := Block[{lj,lx},
    Return[Flatten[Table[
        Table[If[
            lj == 1, {lx, ""}, {Log[10, lj 10^lx], "", {0.00375, 
                0}}], {lj, 1, 9}], {lx, Floor[xmin], xmax}], 1]];
];

RGELogTicksLabeled[xmin_?NumericQ, xmax_?NumericQ] := Block[{lj,lx},
    Return[Flatten[Table[
        Table[If[
            lj == 1, {lx, lx}, {Log[10, lj 10^lx], "", {0.00375, 0}}], {lj, 
            1, 9}], {lx, Floor[xmin], xmax}], 1]];
];

RGELogTicksLabeledNegExp[xmin_?NumericQ, xmax_?NumericQ] := Block[{lj,lx},
    Return[Flatten[Table[
        Table[If[lj == 1, {lx, StringJoin["\!\(10\^\(-",ToString[-lx],"\)\)"]}, {Log[10, 
                lj 10^lx], "", {0.00375, 0}}], {lj, 1, 9}], {lx, Floor[xmin], xmax}], 1]];
];

RGEShadowEFT[tList___:{}] := Block[{plist,lList},
      plist = {};
      lList=Sort[Log[10,{Sequence[tList]}]];
      Do[plist = Append[plist, {GrayLevel[1 - 0.1 i], Rectangle[{lList[[i]], -10^4}, {lList[[i + 1]], 10^4}]}], {i, 1, Length[lList] - 1}];
      Return[plist];
];

End[];

EndPackage[];
