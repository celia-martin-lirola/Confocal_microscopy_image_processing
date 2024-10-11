#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Evaluaci\[OAcute]n entrada ATP en n\[UAcute]cleo*)


(* ::Section:: *)
(*Foto confocal (1 frame)*)


(* ::Text:: *)
(*Se importan im\[AAcute]genes que contengan parte de n\[UAcute]cleo y parte de citoplasma. Primero se deben alinear los canales de confocal con el canal de bright field, para poder seleccionar sobre el canal de ATP la regi\[OAcute]n citoplasm\[AAcute]tica. *)
(*Se seleccionan dos regiones:*)
(* 	- Regi\[OAcute]n nuclear: m\[AAcute]scara con el canal azul.*)
(* 	- Regi\[OAcute]n citoplasm\[AAcute]tica: seleccionar a mano.*)
(* Calcular media de ATP en ambos y comparar. *)


(* ::Subsubsection::Closed:: *)
(*Importar fotos*)


(* ::Text:: *)
(*Importamos el listado de todos los archivos de la carpeta y seleccionamos las fotos para hacer la evaluaci\[OAcute]n de ATP.*)


directory = "C:\\Users\\celia\\OneDrive\\Escritorio\\BioCompLab\\2024 - AECC\\nucleos_HeLa\\Nuevos resultados\\ATP_HO\\Confocal\\240729_Hela_ATP_Ho fijadas"; (*especificamos el directorio*)
fileList = FileNames["*.nd2",directory]; (*crea una lista con todos los archivos en formato .n2d*)


BFList = Select[fileList, (StringCases[FileNameTake[#],"BF"]//Length)>0 &];
BFList = SortBy[BFList, ToExpression@StringCases[FileNameTake[#], RegularExpression["\\d+"]]&];
nBF = BFList[[15;;]];

atpHoList = Select[fileList, (StringCases[FileNameTake[#],"ATP_Ho.nd2"]//Length)>0 &];
atpHoList = SortBy[atpHoList, ToExpression@StringCases[FileNameTake[#], RegularExpression["\\d+"]]&];
natpHo = atpHoList[[13;;]];


(* ::Text:: *)
(*Trabajamos  con  una  sola  imagen (x)  del  listado*)


Needs["BioFormatsLink`"];
x = 1;

(*Importar las fotos*)
bf = Import[nBF[[x]],"BioImageFormat"];
bf = bf[[1,1,1]];
file = FileNameTake[nBF[[x]]];

fluo = Import[natpHo[[x]],"BioImageFormat"];
{atp,ho}=ImageAdjust/@ColorSeparate[fluo[[1,1,1]]];

(*Reescalado ppara igualar dimensiones*)
hos = ImageResize[ho,Scaled[2]];
atps = ImageResize[atp,Scaled[2]];


(* ::Subsubsection::Closed:: *)
(*Alineamiento y m\[AAcute]scara*)


(* ::Text:: *)
(*M\[AAcute]scara nuclear: *)
(*	Existe una diferencia de calibraci\[OAcute]n entre el objetivo de campo claro y el de confocal, se realiza un alineamiento previo para ajustarlo:*)
(*	- Manipulate: se modifican las coordenadas para alineamiento y se suman a cada pixel, devuelve unas coordenadas.*)
(*	- newMask: se pasan de coordenadas a imagen*)
(*	- coord: listado de las coordenadas de cada imagen*)
(*M\[AAcute]scara extracelular: *)
(*	Se seleccionan a mano alzada las regiones extracelulares de la imagen BF y se exportan como imagen.*)


(*Alineamiento de la imagen*)
bin = Binarize[Blur[hos,3]//ImageAdjust]//FillingTransform//DeleteSmallComponents;
cont = ComponentMeasurements[bin,"PerimeterPositions"][[1,2,1]];

Manipulate[
	contf = ({x,y}+#)&/@cont;
	Show[bf^2//ImageAdjust,ListPlot[contf,Joined->True,PlotStyle->Directive[Thickness[0.002],Red]]],

	{{x,-9},-15,5,1},
	{{y,-8},-10,5,1}
]


(*M\[AAcute]scara nuclear*)
newMask = Rasterize[Show[0  bf,Graphics[{White,FilledCurve[Line[contf]]}]],RasterSize->512];
HighlightImage[bf^2//ImageAdjust,newMask]


(*Listado coordenadas*)
(*listado = {15, 16, 17, 18, 19, 20}*)
coord = {{-11, -9},{-11, -8}, {-12,-7}, {-12,-6}, {-11,-7}, {-9,-8}};


(*M\[AAcute]scara extracelular*)
bf^2//ImageAdjust


(*Exportar m\[AAcute]scara estracelular*)
folder = "C:\\Users\\celia\\OneDrive\\Escritorio\\BioCompLab\\2024 - AECC\\nucleos_HeLa\\Nuevos resultados\\ATP_HO\\Confocal\\240729_Hela_ATP_Ho fijadas\\mask\\";
file = FileNameTake[nBF[[x]]];
Export[folder<>StringReplace[file,".nd2"->"_ext_mask.jpeg"],extMask];


(* ::Subsection:: *)
(*Operaciones*)


(* ::Subsubsection::Closed:: *)
(*Importar imagen*)


processImport[fullpathATP_] := Module[
	{atp,ho,fluo,hos,atps},
	
	(*Importar las fotos*)
	Needs["BioFormatsLink`"];
	fluo = Import[fullpathATP,"BioImageFormat"];
	{atp,ho}=ImageAdjust/@ColorSeparate[fluo[[1,1,1]]];
	
	file = FileNameTake[fullpathATP];
	
	(*Reescalado ppara igualar dimensiones*)
	hos = ImageResize[ho,Scaled[2]];
	atps = ImageResize[atp,Scaled[2]];
	
	Return[{hos, atps}]
		
	];


(* ::Subsubsection::Closed:: *)
(*C\[AAcute]lculo intensidad media*)


averageInt[atps_, hos_, extMask_, coord_] := Module[ 
	{extMaskTrans,nucMask,b,c,citMask,meanNuc,meanCit,meanExt},
	
	(*Nuclear mask*)
	nucMask = Binarize[Blur[hos,3]//ImageAdjust]//FillingTransform//DeleteSmallComponents;
	
	(*Image alignment: extracellular mask*)
	extMaskTrans = ImageTransformation[extMask, TranslationTransform[coord], DataRange->Full];
	
	(*Citoplasm mask*)
	citMask = 1-nucMask-extMask;
	
	(*Average intensity*)
	meanNuc = ImageMeasurements[atps*nucMask, "Mean"];
	meanCit = ImageMeasurements[atps*citMask, "Mean"];
	meanExt = ImageMeasurements[atps*extMaskTrans, "Mean"];
	
	Return[{meanNuc, meanCit, meanExt}]
	
	];


(* ::Subsubsection::Closed:: *)
(*Iteraci\[OAcute]n sobre el directorio*)


processFile[nBF_, natpHo_, extMaskList_, coord_, listas1_] := Module[
	{extMask, hos, atps, bf, meanNuc, meanCit, meanExt, listas},
	
	(*Extracellular mask import*)
	extMask = Import[extMaskList,"BioImageFormat"];
	extMask = extMask[[1,1,1]];
	
	{hos, atps} = processImport[natpHo];
	{meanNuc, meanCit, meanExt} = averageInt[atps, hos, extMask, coord];
	
	listas = listas1;
	AppendTo[listas[[1]], meanNuc];
	AppendTo[listas[[2]], meanCit];
	AppendTo[listas[[3]], meanExt];
	
	Return[{listas}]
	
];


(*Lista m\[AAcute]scaras extracelulares*)
direct = "C:\\Users\\celia\\OneDrive\\Escritorio\\BioCompLab\\2024 - AECC\\nucleos_HeLa\\Nuevos resultados\\ATP_HO\\Confocal\\240729_Hela_ATP_Ho fijadas\\mask";
extMaskList = FileNames["*.jpeg",direct];

(*Lista coordenadas*)
coord = {{-11, -9},{-11, -8}, {-12,-7}, {-12,-6}, {-11,-7}, {-9,-8}};

(*Iteraci\[OAcute]n sobre todas las im\[AAcute]genes*)
num = {1, 2, 3, 4, 5, 6};
listas = {{},{},{}};

Do[
	listas = processFile[nBF[[x]], natpHo[[x]], extMaskList[[x]], coord[[x]], listas][[1]],
	{x, num}	
	]


(*Exportar lista medias*)
folder = "C:\\Users\\celia\\OneDrive\\Escritorio\\BioCompLab\\2024 - AECC\\nucleos_HeLa\\Nuevos resultados\\ATP_HO\\Confocal\\240729_Hela_ATP_Ho fijadas\\Resultados\\";
Export[folder<>"medias_fotos.csv",listas];


(* ::Subsubsection::Closed:: *)
(*Representaci\[OAcute]n: boxplot*)


listas = Import["C:\\Users\\celia\\OneDrive\\Escritorio\\BioCompLab\\2024 - AECC\\nucleos_HeLa\\Nuevos resultados\\ATP_HO\\Confocal\\240729_Hela_ATP_Ho fijadas\\Resultados\\medias_fotos.csv"]


(*Representaci\[OAcute]n*)
plot = BoxWhiskerChart[listas,
	"Outliers",
	ChartStyle->{Red, Green, Blue},
	Frame->True,
	FrameStyle->Directive[12,Black],
	FrameLabel->{"","\[LeftAngleBracket][ATP]\[RightAngleBracket]"},
	ImageSize->500,
	PlotLabel->"[ATP] fixed Hela",
	LabelStyle->Directive[15,Black],
	ChartLabels->{"Nuc","Cit","Ext"},
	PlotRange->{0,0.1}
	]


(*Exportar gr\[AAcute]fica*)
Export[folder<>"\\boxplot_fotos.jpeg",plot];


(* ::Subsubsection::Closed:: *)
(*Test diferencia estad\[IAcute]stica*)


nuc = listas[[1]];
cit = listas[[2]];
ext = listas[[3]];

nucCitMWT = MannWhitneyTest[{nuc,cit}];
nucExtMWT = MannWhitneyTest[{nuc,ext}];
citExtMWT = MannWhitneyTest[{cit,ext}];


MWT = Dataset[{{"Nuc/Cit", "Nuc/Ext", "Cit/Ext"},{nucCitMWT, nucExtMWT, citExtMWT}}]
Export[folder<>"\\Mann_Whitney_Test.xlsx", MWT, "XLSX"];


(* ::Subsubsection::Closed:: *)
(*Test bootstraping*)


bootstrapping[values_, n_] := Module[
	{muestreo, medianas, med, desv},	
	
	muestreo = Table[RandomChoice[values, Length[values]],n];
	
	medianas = Map[Median, muestreo];
	med = Mean[medianas];
	desv = StandardDeviation[medianas];
	
	Return[{med,desv}]
	
	];


n = 1000;

nucMed = bootstrapping[nuc,n];
citMed = bootstrapping[cit,n];
extMed = bootstrapping[ext,n];


Boots = Dataset[{{"Nuclei", "Citoplasm", "Extracellular"},
	{ToString[nucMed[[1]]]<>" +/- "<>ToString[nucMed[[2]]],
	ToString[citMed[[1]]]<>" +/- "<>ToString[citMed[[2]]],
	ToString[extMed[[1]]]<>" +/- "<>ToString[extMed[[2]]]}}
	]


Export[folder<>"\\Bootstrapping.xlsx", Transpose[Boots], "XLSX"];


(* ::Section:: *)
(*V\[IAcute]deo .nd2*)


(* ::Subsubsection::Closed:: *)
(*Lista  de archivos .nd2*)


(* ::Text:: *)
(*Seleccionamos las rutas de acceso de los v\[IAcute]deos utilizados para la evaluaci\[OAcute]n de ATP (del 15 al 20), junto con sus respectivas fotos de bright field.*)


directory = "C:\\Users\\celia\\OneDrive\\Escritorio\\BioCompLab\\2024 - AECC\\nucleos_HeLa\\Nuevos resultados\\ATP_HO\\Confocal\\240729_Hela_ATP_Ho fijadas"; (*especificamos el directorio*)
fileList = FileNames["*.nd2",directory]; (*crea una lista con todos los archivos en formato .n2d*)


(*Bright field images*)
BFList = Select[fileList, (StringCases[FileNameTake[#],"BF"]//Length)>0 &];
BFList = SortBy[BFList, ToExpression@StringCases[FileNameTake[#], RegularExpression["\\d+"]]&];
nBF = BFList[[15;;]];

(*Confocal videos*)
atpHoList = Select[fileList, (StringCases[FileNameTake[#],"ATP_Ho0"]//Length)>0 &];
atpHoList = SortBy[atpHoList, ToExpression@StringCases[FileNameTake[#], RegularExpression["\\d+"]]&];
natpHo = atpHoList[[15;;]];


(* ::Subsubsection::Closed:: *)
(*Importamos v\[IAcute]deo*)


Needs["BioFormatsLink`"];
x = 6;

(*Importar foto BF*)
bf = Import[nBF[[x]],"BioImageFormat"];
bf = bf[[1,1,1]];
file = FileNameTake[nBF[[x]]];

(*Importar v\[IAcute]deo atp_Ho*)
fluo = Import[natpHo[[x]],"BioImageFormat"];
cs = Map[ColorSeparate, fluo[[1,All,1]]];
{atp,ho} = Transpose[cs];

(*Reescalado para igualar dimensiones*)
hos = Map[ImageResize[#,Scaled[2]]&, ho];
atps = Map[ImageResize[#,Scaled[2]]&, atp];


(* ::Subsubsection::Closed:: *)
(*Alineamiento imagen (para m\[AAcute]scara extracel)*)


(*Alineamiento de la imagen*)
bin = Binarize[Blur[hos[[1]],3]//ImageAdjust]//FillingTransform//DeleteSmallComponents;
cont = ComponentMeasurements[bin,"PerimeterPositions"][[1,2,1]];

Manipulate[
	contf = ({x,y}+#)&/@cont;
	Show[bf^2//ImageAdjust,ListPlot[contf,Joined->True,PlotStyle->Directive[Thickness[0.002],Red]]],

	{{x,-9},-15,5,1},
	{{y,-8},-10,5,1}
]


coord = {{-12,-10}, {-14,-8}, {-11,-8}, {-10,-8}, {-12,-7}, {-11,-9}};


(* ::Subsection:: *)
(*Operaciones*)


(* ::Subsubsection::Closed:: *)
(*Video import*)


processImport[fullpathATP_] := Module[
	{atp,ho,fluo,hos,atps},
	
	(*Importar v\[IAcute]deo atp_Ho*)
	Needs["BioFormatsLink`"];
	fluo = Import[fullpathATP,"BioImageFormat"];
	cs = Map[ColorSeparate, fluo[[1,All,1]]];
	{atp,ho} = Transpose[cs];
	
	(*Reescalado para igualar dimensiones*)
	hos = Map[ImageResize[#,Scaled[2]]&, ho];
	atps = Map[ImageResize[#,Scaled[2]]&, atp];
	
	Return[{hos, atps}]
		
	];


(* ::Subsubsection::Closed:: *)
(*Average calculation*)


averageInt[atps_, hos_, extMask_, coord_] := Module[ 
	{atpsAdj, extMaskTrans, nucMask, citMask,meanNuc,meanCit,meanExt},
	
	(*Ajuste canal atp, igual en todos los frames para cuatificar*)
	atpsAdj = Map[ImageAdjust[#, {0.2,0.5,1}, {30,300}]&, atps];
	
	(*Extracellular mask*)
	extMaskTrans = ImageTransformation[extMask, # + coord &, DataRange -> Full];
	extMaskTrans = ImageTransformation[extMask, TranslationTransform[coord], DataRange->Full];
	
	(*Nuclear mask*)
	nucMask = Binarize[Blur[hos[[1]],3]//ImageAdjust]//FillingTransform//DeleteSmallComponents;
	
	(*Citoplasm mask*)
	citMask = 1-nucMask-extMask;
	
	(*Average intensity*)
	meanNuc = Mean[ImageMeasurements[Map[#*nucMask &, atpsAdj], "Mean"]];
	meanCit = Mean[ImageMeasurements[Map[#*citMask &, atpsAdj], "Mean"]];
	meanExt = Mean[ImageMeasurements[Map[#*extMaskTrans &, atpsAdj], "Mean"]];
	
	Return[{meanNuc, meanCit, meanExt}]
	
	];


(* ::Subsubsection::Closed:: *)
(*Iteraci\[OAcute]n sobre la lista*)


processFile[nBF_, natpHo_, extMaskList_, coord_, listas1_] := Module[
	{extMask, hos, atps, bf, meanNuc, meanCit, meanExt},
	
	(*Extracellular mask import*)
	extMask = Import[extMaskList,"BioImageFormat"];
	extMask = extMask[[1,1,1]];
	
	{hos, atps} = processImport[natpHo];
	{meanNuc, meanCit, meanExt} = averageInt[atps, hos, extMask, coord];
	
	listas = listas1;
	AppendTo[listas[[1]], meanNuc];
	AppendTo[listas[[2]], meanCit];
	AppendTo[listas[[3]], meanExt];
	
	Return[{listas}]
	
];


(*Lista m\[AAcute]scaras extracelulares*)
direct = "C:\\Users\\celia\\OneDrive\\Escritorio\\BioCompLab\\2024 - AECC\\nucleos_HeLa\\Nuevos resultados\\ATP_HO\\Confocal\\240729_Hela_ATP_Ho fijadas\\mask";
extMaskList = FileNames["*.jpeg",direct];

coord = {{-12,-10}, {-14,-8}, {-11,-8}, {-10,-8}, {-12,-7}, {-11,-9}};
listas = {{},{},{}};

(*N\[UAcute]mero de archivos*)
num = {1,2,3,4,5,6};

(*Iteraci\[OAcute]n sobre todo el directorio*)
Do[
	listas = processFile[nBF[[x]], natpHo[[x]], extMaskList[[x]], coord[[x]], listas][[1]],
	{x, num}	
	]


listas


(*Exportar lista medias*)
folder = "C:\\Users\\celia\\OneDrive\\Escritorio\\BioCompLab\\2024 - AECC\\nucleos_HeLa\\Nuevos resultados\\ATP_HO\\Confocal\\240729_Hela_ATP_Ho fijadas\\Resultados\\";
Export[folder<>"medias_videos.csv",listas];


(* ::Subsubsection::Closed:: *)
(*Representaci\[OAcute]n: boxplot*)


listas = Import["C:\\Users\\celia\\OneDrive\\Escritorio\\BioCompLab\\2024 - AECC\\nucleos_HeLa\\Nuevos resultados\\ATP_HO\\Confocal\\240729_Hela_ATP_Ho fijadas\\Resultados\\medias_videos.csv"]


(*Representaci\[OAcute]n*)
plot = BoxWhiskerChart[listas,
	"Outliers",
	ChartStyle->{Red, Green, Blue},
	Frame->True,
	FrameStyle->Directive[12,Black],
	FrameLabel->{"","\[LeftAngleBracket][ATP]\[RightAngleBracket]"},
	ImageSize->500,
	PlotLabel->"[ATP] fixed Hela",
	LabelStyle->Directive[15,Black],
	ChartLabels->{"Nuc","Cit","Ext"},
	PlotRange->{0,0.1}
	]


(*Exportar gr\[AAcute]fica*)
Export[folder<>"\\boxplot_videos.jpeg",plot];


(* ::Subsubsection::Closed:: *)
(*Test diferencia estad\[IAcute]stica*)


nuc = listas[[1]];
cit = listas[[2]];
ext = listas[[3]];

nucCitMWT = MannWhitneyTest[{nuc,cit}];
nucExtMWT = MannWhitneyTest[{nuc,ext}];
citExtMWT = MannWhitneyTest[{cit,ext}];


MWT = Dataset[{{"Nuc/Cit", "Nuc/Ext", "Cit/Ext"},{nucCitMWT, nucExtMWT, citExtMWT}}]
Export[folder<>"\\Mann_Whitney_Test.xlsx", MWT, "XLSX"];


(* ::Subsubsection::Closed:: *)
(*Test bootstraping*)


bootstrapping[values_, n_] := Module[
	{muestreo, medianas, med, desv},	
	
	muestreo = Table[RandomChoice[values, Length[values]],n];
	
	medianas = Map[Median, muestreo];
	med = Mean[medianas];
	desv = StandardDeviation[medianas];
	
	Return[{med,desv}]
	
	];


n = 1000;

nucMed = bootstrapping[nuc,n];
citMed = bootstrapping[cit,n];
extMed = bootstrapping[ext,n];


Boots = Dataset[{{"Nuclei", "Citoplasm", "Extracellular"},
	{ToString[nucMed[[1]]]<>" +/- "<>ToString[nucMed[[2]]],
	ToString[citMed[[1]]]<>" +/- "<>ToString[citMed[[2]]],
	ToString[extMed[[1]]]<>" +/- "<>ToString[extMed[[2]]]}}
	]


Export[folder<>"\\Bootstrapping.xlsx", Transpose[Boots], "XLSX"];
