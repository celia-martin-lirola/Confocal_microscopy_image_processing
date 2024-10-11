#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Procesado videos ATP + Ho (Particle tracking)*)


(* ::Text:: *)
(*-Operaciones: Algoritmo  de  procesamiento  en  funciones  para  los  v\[IAcute]deos  de  n\[UAcute]cleo  marcados  con  sonda  de  ATP  y  Hoesch . Realiza  los  test  de  actividad, estructura  y  difusividad, comaprando  			concentraciones  medias  de  ATP  y  cromatina, PIV  para  cromatina  y  WTV  para  ATP .*)
(*-Particle tracking: seguimiento de los puntos de inter\[EAcute]s sacados a partir de m\[AAcute]ximos locales. Representaci\[OAcute]n de difusividad y viscosidad.*)


(* ::Section::Closed:: *)
(*Input*)


(* ::Subsection::Closed:: *)
(*Importar v\[IAcute]deos*)


directory = "C:\\Users\\celia\\OneDrive\\Escritorio\\BioCompLab\\2024 - AECC\\nucleos_HeLa\\Nuevos resultados\\ATP_HO\\Confocal\\240729_Hela_ATP_Ho fijadas"; (*especificamos el directorio*)
fileList = FileNames["*.nd2",directory]; (*crea una lista con todos los archivos en formato .n2d*)


(* ::Text:: *)
(*Seleccionar los archivos con el nombre de ATP_Ho para no importar las im\[AAcute]genes de campo claro.*)


atpHoList = Select[fileList, (StringCases[FileNameTake[#],"ATP_Ho"]//Length)>0 &]


(*Carpeta de resultados*)
folder = directory <>"\\procesado_videos_OF_"<>FileNameTake[directory]<>"\\";
CreateDirectory[folder]


(* ::Section:: *)
(*Algoritmo en funciones*)


(* ::Text:: *)
(*Est\[AAcute]n cambiados los canales de sitio en Import video porque vinieron al rev\[EAcute]s en los \[UAcute]ltimos v\[IAcute]deos que mandaron.*)


(* ::Subsubsection::Closed:: *)
(*Import video*)


processImport[fullpath_] := Module[

	{file,newFolder,vid,cs,b,r},
		
	(* archivo *)
	file = FileNameTake[fullpath];
	(*Nueva carpeta de resultados para cada elemento de la lista*)
	newFolder = folder<>FileNameTake[file]<>"\\";
	
	(* Import Video *)
	Needs["BioFormatsLink`"];
	vid = Import[fullpath,"BioImageFormat"];
	
	(* referencia de niveles *)
	vid = vid[[1,All,1]]/2.^12;
	(*Separaci\[OAcute]n de canales*)
	cs = Map[ColorSeparate,vid];
	{r,b} = Transpose[cs];(*Canales cambiados de sitio*)
	
	Return[{file,newFolder,r,b}]
	
	];


(* ::Subsubsection::Closed:: *)
(*Ajuste y alineamiento*)


(* ::Text:: *)
(*Filtro gaussiano aplicado a un stack en 3D de todos los frames del v\[IAcute]deo para que todos los frames se ajusten por igual.*)
(*Los canales r y b para cuantificaci\[OAcute]n quedan con un ajuste de imagen igual para todas las c\[EAcute]lulas, para poder realizar comparaci\[OAcute]n (especificar el ajuste en la funci\[OAcute]n processVid). Se deja una variable llamada pAdj con el Image Adjust aplicado para utilizarlo con la m\[AAcute]scara y los dem\[AAcute]s procesos que no requieren cuantificaci\[OAcute]n.*)


processAdj[r_,b_,ajuste_] := Module[

	{rAdj,bAdj,pAdj,e,tf},
	
	(*Ajuste de imagen*)
	rAdj = Image3DSlices[GaussianFilter[Image3D[r],4] ajuste];
	bAdj = Image3DSlices[GaussianFilter[Image3D[b],4] ajuste];
	pAdj = Image3DSlices[ImageAdjust[GaussianFilter[Image3D[b],4]]];
	
	(*Alineamiento*)
	Table[
		{e,tf} = Transpose@FindGeometricTransform[pAdj[[i-1]],pAdj[[i]],Method->{"ImageAlign",{"Keypoints","SIFT"}},TransformationClass->"Rigid"];
		pAdj[[i]]=ImagePerspectiveTransformation[pAdj[[i]],tf,DataRange->Full];
		bAdj[[i]]=ImagePerspectiveTransformation[pAdj[[i]],tf,DataRange->Full];
		rAdj[[i]]=ImagePerspectiveTransformation[pAdj[[i]],tf,DataRange->Full];
		,{i,2,Length[pAdj],1}
	];
	
	Return[{bAdj,rAdj,pAdj}]
	];


(* ::Subsubsection::Closed:: *)
(*MinMaxDif*)


(* ::Text:: *)
(*Selecci\[OAcute]n de las posiciones de los puntos con intensidades m\[AAcute]ximas (es como aplicar un binarizado pero seleccionando solo los m\[AAcute]ximos).*)


processMinMax[bAdj_] := Module[
	{maxpos},	
		maxpos = ParallelMap[
			ImageAdjust[
				GaussianFilter[
					Dilation[MaxDetect[GaussianFilter[LocalAdaptiveBinarize[#,8],8]],DiskMatrix[2]] 
				,5]]&, bAdj];
	Return[maxpos]
	];


(* ::Subsubsection::Closed:: *)
(*M\[AAcute]scara de filtrado*)


(* ::Text:: *)
(*Selecci\[OAcute]n de la regi\[OAcute]n nuclear para descartar toda la actividad metab\[OAcute]lica citoplasm\[AAcute]tica. Se aplica una erosi\[OAcute]n para eliminar la fluctuaci\[OAcute]n de membrana.*)
(*Para evitar que la erosi\[OAcute]n no se aplique en los bordes de la imagen se a\[NTilde]ade un marco de p\[IAcute]xeles negros que luego se quita (funci\[OAcute]n Borders).*)


(*Erosi\[OAcute]n*)
Borders[img_,erosPix_] := Module[
		{img1,img2,img3},
			img1 = ImagePad[img,3];
			img2 = Erosion[img1,DiskMatrix[erosPix]];
			img3 = ImageCrop[img2,img//ImageDimensions];
			Return[img3];
			];

(*M\[AAcute]scara*)
processMask[bAdj_,erosPix_] := Module[{masc},
	masc = Map[FillingTransform[Binarize[GaussianFilter[#, 5],0.1]]&, bAdj];
		(*Erosi\[OAcute]n bordes*)
		masc = Map[Borders[#, erosPix]&, masc];
	Return[masc];
	];


(* ::Text:: *)
(*Para las im\[AAcute]genes nuevas para las que no funciona la m\[AAcute]scara nueva esto es cogiendo un par de coordenadas dentro de la imagen que est\[AAcute]n fuera del n\[UAcute]cleo para que seleccione en funci\[OAcute]n de eso.*)


(*1-RegionBinarize[ColorNegate[ImageAdjust[GaussianFilter[b[[1]],10]]],{{237.75, 227.908447265625}, {230.25, 25.408447265625}},0.25]*)


(* ::Subsubsection::Closed:: *)
(*PIV (Particle Image Velocimetry )*)


(* ::Text:: *)
(*Se divide cada frame de las posiciones m\[AAcute]ximas (MinMaxDif) en partes iguales (como en un grid), y se va comparando estas partes de la imagen con las del frame siguiente (para ver cu\[AAcute]nto se ha movido la part\[IAcute]cula). Se saca un vector de velocidad que refleja cu\[AAcute]nto se ha movido, del que luego se saca el m\[OAcute]dulo. Por \[UAcute]ltimo, se hace la media de los m\[OAcute]dulos de todos los frames. *)


(*Particle tracking*)	
	piv[pares_,initx_,iw_,ih_,\[Sigma]_]:= Module[
		{track,flow},
			(* Tracking *)
			track = ImageFeatureTrack[pares,initx,MaxIterations->20,MaxFeatureDisplacement->3];
			(* Flow *)
			flow = Reverse@Transpose[ArrayReshape[(track[[2]]-track[[1]])/. (x_+Missing[])->0.,{iw,ih,2}]];
			(* Filtrado espacial del flujo *)
			flow[[All,All,1]] = GaussianFilter[flow[[All,All,1]],\[Sigma]];
			flow[[All,All,2]] = GaussianFilter[flow[[All,All,2]],\[Sigma]];
			Return[flow]
		];

(*Funci\[OAcute]n PIV aplicado*)
processPIV[bBin_] := Module[
	{w,h,\[Delta],\[Sigma],\[Delta]t,initx,iw,ih,veli,normvelimean},
	
	(* Posiciones iniciales para el tracking *)
		{w,h} = ImageDimensions[bBin[[1]]];
		\[Delta] = 1; (* Resoluci\[OAcute]n espacial *)
		\[Sigma] = 20; (* Filtrado espacial del flujo *)
		\[Delta]t = 1; (* Espaciado temporal *)
		
	(* Posiciones iniciales para el tracking *)
		initx = Table[{x,y},{x,0.5,w-0.5,\[Delta]},{y,0.5,h-0.5,\[Delta]}];
		{iw,ih} = Dimensions[initx][[1;;2]];
		initx = Flatten[initx,1];
		
	(* Matriz de vectores para todos los frames *)
		veli = ParallelMap[piv[#,initx,iw,ih,\[Sigma]]&,Partition[bBin,2,1]];
				
	(* M\[OAcute]dulo de la media *)
		normvelimean = With[
			{normveli = Map[Norm,veli,{3}]},
			Mean[normveli]
			];
			
	Return[{veli,normvelimean}]
	];


(* ::Subsubsection::Closed:: *)
(*Waiting time variance (Tiempo de latencia)*)


(* ::Text:: *)
(*Se aplica al canal rojo del ATP ajustado. Se realiza el siguiente ajuste para obtener los tiempos de residencia:*)
(*	- Se obtiene una lista de datos de las intensidades de un mismo p\[IAcute]xel en el tiempo (en todos los frames del v\[IAcute]deo).*)
(*	- Se resta la media de intensidades a todos los datos (unos quedar\[AAcute]n negativos y otros positivos).*)
(*	- Los negativos se igualan a 0 y los positivos a 1.*)
(*	- Se cuentan los instantes de tiempo que cada p\[IAcute]xel permanece en 1 de forma continuada y se obtienen los tiempos de residencia. *)
(*	- Se aplica el c\[AAcute]lculo: varianza/media^2*)
(* A m\[AAcute]s varianza mayor movilidad en la imagen.*)


(*C\[AAcute]lculo del tiempo de latencia*)
	wtvar[signal_]:= Module[{s0=signal,split,bin,l0,l1},
		(* Referenciado a la media *)
		s0 = s0 - Mean[s0];
		(* Binarizando, eliminando 0.5 *)
		bin = DeleteCases[0.5 Sign[s0]+0.5,0.5];
		(* Recuento de tiempo pasado en estados 1 *)
		split = Split[bin];
		l0 = Length/@Select[split,(First@#)== 1 &];
		(* Output *)
		If[Length[l0]<=1,0,N[Variance[l0]/Mean[l0]^2]]
	];
	
(*Funci\[OAcute]n aplicaci\[OAcute]n*)
	processWT[rAdj_] := Module[{rRes,tstack,wtdata},
		tstack = Transpose[(ImageData/@rAdj),1<->3];
		wtdata = Map[wtvar,tstack,{2}];
		
		Return[wtdata]
		];


(* ::Section:: *)
(*Operaciones*)


(* ::Subsection::Closed:: *)
(*Procesado del v\[IAcute]deo y obtenci\[OAcute]n de par\[AAcute]metros*)


(* ::Text:: *)
(*La movilidad de la cromatina se calcula aplicando el PIV (particle tracking) a las regiones de m\[AAcute]xima intensidad (MinMaxDif) y calculando la media de los m\[OAcute]dulos de velocidad. Se trata de estructuras grandes que no presentan tanta movilidad como el ATP y a las que se puede seguir la trayectoria.*)
(*La movilidad del ATP se calcula con el tiempo de latencia (Waiting Time Variance) ya que se trata de mol\[EAcute]culas muy peque\[NTilde]as a las que no se puede seguir la trayectoria.*)


processVid[fullpath_] := Module[
	{file,newFolder,vidRaw,r,b,bAdj,rAdj,maxpos,erosPix,masc,bBin,rBin,bPIV,pAdj,pBin,ajuste,
	vidAdj,veli,normvelimean,vecImg,maskImg,mapaVec,rRes,wtdata,movATP,wtdataR,dataMovAtp,dataMovCrom,dataAtp,dataCrom},
	Print[fullpath];
	
	(*Video import*)
	Print["Video import"];
	{file,newFolder,r,b} = processImport[fullpath];
	CreateDirectory[newFolder];
	vidRaw = MapThread[ColorCombine[{#1, 0 #1, #2}]&,{r, b} ];
	Export[newFolder<>StringReplace[file,".nd2"->"_raw.gif"],Map[ImageAdjust,vidRaw]];
	
	(*Video adjust*)
	Print["Video adjust"];
	ajuste = 10;
	{bAdj,rAdj,pAdj} = processAdj[r,b,ajuste];
	
	(*MinMaxDif*)
	Print["MinMaxDif"];
	maxpos = processMinMax[pAdj];
	
	(*Mask*)
	Print["Mask"];
	erosPix = 4; (*P\[IAcute]xeles de erosi\[OAcute]n*)
	masc = processMask[pAdj,erosPix];
	Export[newFolder<>StringReplace[file,".nd2"->"_mask.gif"],masc];
		bBin = bAdj*masc;
		rBin = rAdj*masc;
		pBin = pAdj*masc;
		bPIV = maxpos*masc;
	vidAdj = MapThread[ColorCombine[{#1, 0 #1, #2}]&,{rBin, bBin} ];
	Export[newFolder<>StringReplace[file,".nd2"->"_adj.gif"],Map[ImageAdjust,vidAdj]];
	
	(*Particle tracking*)
	Print["PIV"];
	{veli,normvelimean} = processPIV[pBin];
			(*Mapa vectorial*)
			Print["Mapa vectorial"];
			vecImg = Show[Blur[Mean[b],2]//ImageAdjust,ListVectorPlot[Mean[veli],VectorPoints->40,VectorScaling->False,StreamPoints->10]];
			maskImg = ImageResize[Mean[masc], ImageDimensions[vecImg]];
			mapaVec = ImageMultiply[vecImg, maskImg];
			Export[newFolder<>StringReplace[file,".nd2"->"_mov_crom.jpeg"],mapaVec];
	
	(*Waiting time variance*)
	Print["Waiting time variance"];
	rRes = Map[ImageResize[#,Scaled[1/3]]&, rBin];
	wtdata = processWT[rRes]; (*Imagen reescalada*)
	movATP = Colorize[GaussianFilter[Image[wtdata],3],ColorFunction->"Rainbow"];
	Export[newFolder<>StringReplace[file,".nd2"->"_mov_atp.jpeg"],movATP];
	
	(*Promediado y reescalado*)
	Print["Promediado"];
	wtdataR = ImageResize[Image[wtdata],Scaled[1/4]];
	dataMovAtp =  Flatten[ImageData[ImageClip[wtdataR]]];
	dataMovCrom = Flatten[ImageData[ImageClip[ImageResize[Image[normvelimean],wtdataR//ImageDimensions]]]];
	dataAtp = Flatten[ImageData[ImageClip[ImageResize[Mean[rBin],wtdataR//ImageDimensions]]]];
	dataCrom = Flatten[ImageData[ImageClip[ImageResize[Mean[bBin],wtdataR//ImageDimensions]]]];
	
	Return[{newFolder,file,dataMovAtp,dataMovCrom,dataAtp,dataCrom}];
	Print["PROCESSED VIDEO"];

];


(* ::Subsection::Closed:: *)
(*Test y representaci\[OAcute]n*)


(*Funci\[OAcute]n media de puntos*)
media[lista_]:= Around[Mean[lista],StandardDeviation[lista]];


(* ::Subsubsection::Closed:: *)
(*Test actividad*)


(* ::Text:: *)
(*Concentraci\[OAcute]n de ATP en funci\[OAcute]n de la movilidad de la cromatina (PIV).*)


actTest[dataAtp_,dataMovCrom_] := Module[
	{data,xx,fit,ptoMedio,plot},
	
	data = Transpose[{dataAtp,dataMovCrom}];
	
	(*Ajuste*)
		fit = Fit[data,{xx,1},xx];
	(*Media de puntos*)
		ptoMedio = Map[media,Transpose[Select[data,#[[2]]>0.005 &]]];
		
	(*Representaci\[OAcute]n*)
	plot = Show[
		ListPlot[data,
			PlotStyle->Blue,
			Frame->True,
			FrameStyle->Directive[12,Black],
			FrameLabel->{"\[LeftAngleBracket]|\!\(\*SubscriptBox[\(v\), \(crom\)]\)|]\[RightAngleBracket]","\[LeftAngleBracket][ATP]\[RightAngleBracket]"},
			ImageSize->500,
			PlotLabel->"Test actividad",
			LabelStyle->Directive[15,Black]
			],
		Graphics[Point[ptoMedio[[All,1]]]],
		Plot[fit,
			{xx,0.0,Max[data[[All,1]]]},
			PlotStyle->Directive[Red,Thick,Dashed]	
			]
		];
		
	Return[{plot,fit,ptoMedio}]
	];


(* ::Subsubsection::Closed:: *)
(*Test estructural*)


(* ::Text:: *)
(*Concentraci\[OAcute]n de ATP en funci\[OAcute]n de la concentraci\[OAcute]n de cromatina.*)


strTest[dataAtp_,dataCrom_] := Module[
	{data,xx,fit,ptoMedio,plot},
	
	data = Transpose[{dataAtp,dataCrom}];
	
	(*Ajuste*)
		fit = Fit[data,{xx,1},xx];
	(*Media de puntos*)
		ptoMedio = Map[media,Transpose[Select[data,#[[2]]>0.005 &]]];
		
	(*Representaci\[OAcute]n*)
	plot = Show[
		ListPlot[data,
			PlotStyle->Blue,
			Frame->True,
			FrameStyle->Directive[12,Black],
			FrameLabel->{"\[LeftAngleBracket][Crom]\[RightAngleBracket]","\[LeftAngleBracket][ATP]\[RightAngleBracket]"},
			ImageSize->500,
			PlotLabel->"Test estructural",
			LabelStyle->Directive[15,Black]
			],
		Graphics[Point[ptoMedio[[All,1]]]],
		Plot[fit,
			{xx,0.0,Max[data[[All,1]]]},
			PlotStyle->Directive[Red,Thick,Dashed]	
			]
		];
		
	Return[{plot,fit,ptoMedio}]
	];


(* ::Subsubsection::Closed:: *)
(*Test difusividad*)


(* ::Text:: *)
(*Difusividad del ATP (tiempos de latencia) en funci\[OAcute]n de la movilidad de la cromatina (PIV)*)


difTest[dataMovAtp_,dataMovCrom_] := Module[
	{data,xx,fit,ptoMedio,plot},
	
	data = Transpose[{dataMovAtp,dataMovCrom}];
	
	(*Ajuste*)
		fit = Fit[data,{xx,1},xx];
	(*Media de puntos*)
		ptoMedio = Map[media,Transpose[Select[data,#[[2]]>0.005 &]]];
		
	(*Representaci\[OAcute]n*)
	plot = Show[
		ListPlot[data,
			PlotStyle->Blue,
			Frame->True,
			FrameStyle->Directive[12,Black],
			FrameLabel->{"\[LeftAngleBracket]|\!\(\*SubscriptBox[\(v\), \(crom\)]\)|]\[RightAngleBracket]","\[LeftAngleBracket]|\!\(\*SubscriptBox[\(v\), \(ATP\)]\)|]\[RightAngleBracket]"},
			ImageSize->500,
			PlotLabel->"Test difusividad",
			LabelStyle->Directive[15,Black]
			],
		Graphics[Point[ptoMedio[[All,1]]]],
		Plot[fit,
			{xx,0.0,Max[data[[All,1]]]},
			PlotStyle->Directive[Red,Thick,Dashed]	
			]
		];
		
	Return[{plot,fit,ptoMedio}]
	];


(* ::Subsubsection::Closed:: *)
(*Listas de pendiente y pto medio*)


(* ::Text:: *)
(*Se guardan las pendientes de todos los test y el punto medio en el eje x e y con su incertidumbre  para comparar la correlaci\[OAcute]n.*)


processTest[newFolder_,file_,dataMovAtp_,dataMovCrom_,dataAtp_,dataCrom_,listas_] := Module[
	{actPlot,actFit,actPtoMedio,actFitList,actMeanList,
	strPlot,strFit,strPtoMedio,strFitList,strMeanList,
	difPlot,difFit,difPtoMedio,difFitList,difMeanList
	},
	
	(*Test actividad: [ATP] vs movCrom*)
	Print["Test actividad"];
	{actPlot,actFit,actPtoMedio} = actTest[dataAtp,dataMovCrom];
		Export[newFolder<>StringReplace[file,".nd2"->"_act_test.jpeg"],actPlot];
		actFitList = Append[listas[[1,1]], {actFit[[1]],actFit[[2,1]]}];
		actMeanList = Append[listas[[2,1]], actPtoMedio];
	
	(*Test estructural: [ATP] vs [Crom]*)
	Print["Test estructural"];
	{strPlot,strFit,strPtoMedio} = strTest[dataAtp,dataCrom];
		Export[newFolder<>StringReplace[file,".nd2"->"_str_test.jpeg"],strPlot];
		strFitList = Append[listas[[1,2]], {strFit[[1]],strFit[[2,1]]}];
		strMeanList = Append[listas[[2,2]], strPtoMedio];
	
	(*Test difusividad: movATP vs movCrom*)
	Print["Test difusividad"];
	{difPlot,difFit,difPtoMedio} = difTest[dataMovAtp,dataMovCrom];
		Export[newFolder<>StringReplace[file,".nd2"->"_dif_test.jpeg"],difPlot];
		difFitList = Append[listas[[1,3]], {difFit[[1]],difFit[[2,1]]}];
		difMeanList = Append[listas[[2,3]], difPtoMedio];
	
	Print["Done "<>file];
	Return[{{actFitList,strFitList,difFitList},{actMeanList,strMeanList,difMeanList}}]
	];



(* ::Subsection::Closed:: *)
(*Iteraci\[OAcute]n sobre todos los archivos del directorio*)


(* ::Text:: *)
(*Se  realiza  el  mismo  proceso  para  todos  los archivos del directorio que se haya especificado.*)


(*Lista: {{actFitList,strFitList,difFitList},{actMeanList,strMeanList,difMeanList}} *)
listas = {{{},{},{}},{{},{},{}}};

Do[
	{newFolder,file,dataMovAtp,dataMovCrom,dataAtp,dataCrom} = processVid[fullpath];
	listas = processTest[newFolder,file,dataMovAtp,dataMovCrom,dataAtp,dataCrom,listas], 
	{fullpath,fileList}
	]
		
Print["DONE"]



(* ::Subsection::Closed:: *)
(*Ptos medios y pendientes*)


(* ::Text:: *)
(*Representaci\[OAcute]n de las listas de datos para observar la correlaci\[OAcute]n entre n\[UAcute]cleos celulares.*)


{{actFitList,strFitList,difFitList},{actMeanList,strMeanList,difMeanList}} = listas;


(* ::Subsubsection::Closed:: *)
(*Exportar datos*)


fits = Join[
	{{"actFit oo","actFit slope","strFit oo","strFit slope","difFit oo","difFit slope"}},
	Flatten/@Transpose[{actFitList,strFitList,difFitList}]
	];
	
ptoMed = Join[
	{{"actMean xpoint","actMean xerror","actMean ypoint","actMean yerror",
	"strMean xpoint","strMean xerror","strMean ypoint","strMean yerror",
	"difMean xpoint","difMean xerror","difMean ypoint","difMean yerror"}},
	Flatten/@Map[{First@#,Last@#}&,Flatten/@Transpose[{actMeanList,strMeanList,difMeanList}],{2}]
	];

Export[
	folder<>StringReplace[file,".nd2"->"_data_list.xlsx"],
	{
		{
			"Fits"-> fits,
			"Pto Medio"-> ptoMed
		}
	}
	,
	{{"Sheets"}}
]


(* ::Subsubsection::Closed:: *)
(*Representaci\[OAcute]n punto medio*)


(*Funci\[OAcute]n representaci\[OAcute]n para pendiente y ptomedio*)
dataRep[name_,xaxis_,yaxis_,data_] := Module[
	{fit,plot,xx},
		
	fit = Fit[data[[All,All,1]],{xx,1},xx];
	
	plot = Show[
		ListPlot[data,
			PlotStyle->Blue,
			Frame->True,
			FrameStyle->Directive[12,Black],
			FrameLabel->{xaxis,yaxis},
			ImageSize->500,
			PlotLabel->name,
			LabelStyle->Directive[15,Black]
			],
		Plot[fit,
			{xx,0.0,0.15},
			PlotStyle->Directive[Red,Thick,Dashed]	
			]
		];
		
	Return[{plot,fit}];
	];


(*Test de actividad*)	
	name = "Media test actividad";
	xaxis = "Std(\!\(\*SubscriptBox[\Crom\), \(B\)]\))";
	yaxis = "\[LeftAngleBracket][ATP]\[RightAngleBracket]";
	{actPlot,actFit} = dataRep[name,xaxis,yaxis,actMeanList];
	Export[folder<>StringReplace[file,".nd2"->"_act_test_ptomed.jpeg"],actPlot];
	
(*Test estructural*)	
	name = "Media test estructural";
	xaxis = "\[LeftAngleBracket][Crom]\[RightAngleBracket]";
	yaxis = "\[LeftAngleBracket][ATP]\[RightAngleBracket]";
	{strPlot,strFit} = dataRep[name,xaxis,yaxis,strMeanList];
	Export[folder<>StringReplace[file,".nd2"->"_str_test_ptomed.jpeg"],strPlot];
	
(*Test de difusi\[OAcute]n*)	
	name = "Media test difusi\[OAcute]n";
	xaxis = "Std(\!\(\*SubscriptBox[\(MinMaxDif Crom\), \(B\)]\))";
	yaxis = "Std(\!\(\*SubscriptBox[\(MinMaxDif ATP\), \(B\)]\))";
	{difPlot,difFit} = dataRep[name,xaxis,yaxis,difMeanList];
	Export[folder<>StringReplace[file,".nd2"->"_dif_test_ptomed.jpeg"],difPlot];


(* ::Subsubsection::Closed:: *)
(*Representaci\[OAcute]n pendiente*)


slopeRep = boxPlot = BoxWhiskerChart[{actFitList[[All,2]],strFitList[[All,2]],difFitList[[All,2]]}, 
		"Outliers", 
		PlotLabel->"ATP Crom correlation", 
		LabelStyle->{Black, 12},
		ChartLegends->{"Activity test", "Structural test", "Diffusion test"},
		ChartStyle->24
	]

Export[folder<>StringReplace[file,".nd2"->"_test_slopes.jpeg"],slopeRep];


(* ::Section::Closed:: *)
(*Particle Tracking*)


(* ::Subsection::Closed:: *)
(*Funci\[OAcute]n particle tracking*)


(* ::Text:: *)
(*Aplica una m\[AAcute]scara m\[AAcute]s sencilla al particle traking que se aplica solo al primer frame para sacar los interest key points (ikp). Los puntos de inter\[EAcute]s son los m\[AAcute]ximos locales en las intensidades de la imagen, que se ponen como punto dentro de la imagen. Una vez detectados estos puntos de inter\[EAcute]s, se van trackeando en cada frame, obteniendo una trayectoria para cada punto. En los puntos en los que no encuentra trayectoria a partir de un determinado frame se a\[NTilde]ade a la lista de trayectorias como Null. Estas trayectorias que se pierden a la mitad se eliminan de la lista de ikp (DeleteCases). *)
(**)
(*La funci\[OAcute]n devuelve los siguientes par\[AAcute]metros:*)
(*- imtrack: imagen reescalada utilizada para trackear*)
(*- mask: m\[AAcute]scara utilizada en el primer frame para obtener los ikp*)
(*- ikp: interest key points*)
(*-track:  lista de trayectorias de los ikp a lo largo de todos los frames*)
(*-traj: trayectorias limpias y colocadas a lo largo del tiempo y no de la imagen*)


processPTrack[pAdj_,scale_] := Module[
	{imtrack,mask,ikp,track,traj},
	imtrack = Map[ImageResize[#,Scaled[scale]]&,pAdj];
	mask = 1-DeleteSmallComponents[Binarize[ColorNegate[imtrack[[1]]],0.9]];
	
	(*interest key points*)
	ikp = ComponentMeasurements[MaxDetect[mask ImageAdjust[imtrack[[1]]],0.05],"Centroid"][[All,2]];
	
	track = ImageFeatureTrack[imtrack,ikp,Tolerance->0.1,MaxFeatureDisplacement->5];
	traj = Transpose[track];
	traj = Map[DeleteCases[#,Missing[]]&,traj];
	traj = Select[traj,Length[#]>3 &];
	
	Return[{imtrack,mask,ikp,traj,track}]
];


(* ::Subsubsection::Closed:: *)
(*Otras aproximaciones para particle tracking*)


(* ::Text:: *)
(*Utilizando funciones diferentes para detectar los puntos de inter\[EAcute]s (la que mejor funcionaba es la de arriba)*)


(*ikp1 = ImageKeypoints[ImageAdjust[bAdj[[1]]^2],Method->"KAZE",MaxFeatures->100,Masking->mask ];
HighlightImage[bAdj[[1]],ikp1]*)


(*ikp2 = ComponentMeasurements[MaxDetect[mask ImageAdjust[bAdj[[1]]],0.05],"Centroid"][[All,2]];
HighlightImage[bAdj[[1]],ikp2]*)


(* ::Subsection::Closed:: *)
(*Difusividad*)


(* ::Text:: *)
(*Se calcula:*)
(*	- \[CapitalDelta]x (trayectoria): desviaci\[OAcute]n est\[AAcute]ndar media en los ejes x e y*)
(*	- \[CapitalDelta]v (velocidad): calculada a partir de la desviaci\[OAcute]n est\[AAcute]ndar de la derivada de la trayectoria, promediado para ejes x e y*)
(*El producto de ambos par\[AAcute]metros permite obtener la difusividad de la part\[IAcute]cula, que es la inversa a la viscosidad (con un coeficiente): \[Eta] = 1/D*)


(* ::Text:: *)
(*Para compararlo en funci\[OAcute]n de la actividad metab\[OAcute]lica (ATP), se calcula la intensidad total del ATP en la regi\[OAcute]n de la part\[IAcute]cula de inter\[EAcute]s, siendo esta regi\[OAcute]n el \[AAcute]rea del cuadrado de 5 p\[IAcute]xeles alrededor (10x10). Para hacerlo m\[AAcute]s ajustado, se aplic\[OAcute] un filtro gaussiano al fragmento de imagen de inter\[EAcute]s, poniendo m\[AAcute]s peso a las intensidades centrales que a las de los extremos de la regi\[OAcute]n. Para que no haya problema con los puntos de inter\[EAcute]s cerca de los extremos, se aplica un marco de los mismos p\[IAcute]xeles.*)


processDif[trajx_,rRes_,scale_,mask_] := Module[
	{stdx,stdy,vx,vy,\[CapitalDelta]x,\[CapitalDelta]v,dif,pos,maskPad,meanAtp,width,coord,windowAtp,windowMask,totalI,atpI,atp},
	
	(*Difussion*)
	{stdx,stdy} = StandardDeviation[trajx];
	\[CapitalDelta]x = (stdx+stdy)/2;
	
	{vx,vy} = StandardDeviation[DerivativeFilter[trajx,{1},3]];
	\[CapitalDelta]v = (vx+vy)/2;
	
	dif = \[CapitalDelta]x \[CapitalDelta]v;
	
	(*ATP concentration*)
	width = 5; (*Tama\[NTilde]o cuadrado de la ventana, en p\[IAcute]xeles*)
	
	maskPad = ImagePad[mask,width+1];
	meanAtp = ImagePad[Mean[rRes],width+1];
	
	pos = Round[trajx[[1]]]+{width+1,width+1};
	coord = {{pos[[1]]-width+1, pos[[2]]-width+1},{pos[[1]]+width, pos[[2]]+width}};
	windowAtp = ImageTrim[meanAtp, coord];
	windowMask = ImageTrim[maskPad, coord];
		
	totalI = Total[Flatten[ImageData[windowMask]*GaussianMatrix[width]]];
	atpI = Total[Flatten[ImageData[windowAtp]*GaussianMatrix[width]]];
	atp = atpI/totalI;
	
	Return[{pos,\[CapitalDelta]x,\[CapitalDelta]v,dif,atp}]
];


(* ::Subsection::Closed:: *)
(*Representaci\[OAcute]n*)


(* ::Text:: *)
(*Representaci\[OAcute]n de la difusi\[OAcute]n y su inverso, la viscosidad, calculando tanto el fit como el punto medio (x,y).*)


(*Funci\[OAcute]n media de puntos*)
media[lista_]:= Around[Mean[lista],StandardDeviation[lista]];


(* ::Subsubsection::Closed:: *)
(*Difusi\[OAcute]n*)


difRes[difTable_] := Module[{dataDif,fit,ptoMedio,plot,xx},

	(*Data*)
	dataDif = Transpose[{difTable[[All,4]],difTable[[All,5]]}]; (*{Difusion,[ATP]}*)
		
	(*Ajuste*)
	fit = Fit[dataDif,{xx,1},xx];
	(*Punto medio*)
	ptoMedio = Map[media,Transpose[Select[dataDif,#[[2]]>0.005 &]]];
		
	(*Representaci\[OAcute]n*)
	plot = Show[
		ListLogLinearPlot[dataDif,
			PlotStyle->Blue,
			Frame->True,
			FrameStyle->Directive[12,Black],
			FrameLabel->{"\*SubscriptBox[\(D\), \(crom\)]","\[LeftAngleBracket][ATP]\[RightAngleBracket]"},
			ImageSize->500,
			PlotLabel->"Particle tracking",
			LabelStyle->Directive[15,Black]
			],
		Plot[fit,
			{xx,0.0,Max[dataDif[[All,1]]]},
			PlotStyle->Directive[Red,Thick,Dashed]	
			]
		];
		
	Return[{fit,ptoMedio,plot}]
];


(* ::Subsubsection::Closed:: *)
(*Viscosidad*)


difVis[difTable_] := Module[{dataVis,fit,ptoMedio,plot,xx},

	(*Data*)
	dataVis = Transpose[{1/difTable[[All,4]], difTable[[All,5]]}];(*{viscosity = 1/D,[ATP]}*)
		
	(*Ajuste*)
	fit = Fit[dataVis,{xx,1},xx];
	(*Punto medio*)
	ptoMedio = Map[media,Transpose[Select[dataVis,#[[2]]>0.005 &]]];
		
	(*Representaci\[OAcute]n*)
	plot = Show[
		ListPlot[dataVis,
			PlotStyle->Blue,
			Frame->True,
			FrameStyle->Directive[12,Black],
			FrameLabel->{"\*SubscriptBox[\(\[Eta]\), \(crom\)]","\[LeftAngleBracket][ATP]\[RightAngleBracket]"},
			ImageSize->500,
			PlotLabel->"Particle tracking",
			LabelStyle->Directive[15,Black]
			],
		Plot[fit,
			{xx,0.0,Max[dataVis[[All,1]]]},
			PlotStyle->Directive[Red,Thick,Dashed]	
			]
		];
		
	Return[{fit,ptoMedio,plot}]
];


(* ::Subsection:: *)
(*Operaciones*)


(* ::Subsubsection::Closed:: *)
(*Algoritmo*)


(* ::Text:: *)
(*Todas  las  operaciones  para  el particle tracking exportando:*)
(*- V\[IAcute]deo en crudo*)
(*- Ikp: interest key points*)
(*- Trayectoria de los puntos de inter\[EAcute]s*)
(*- V\[IAcute]deo ajustado*)
(*- Tabla en Excel de par\[AAcute]metros calculados*)
(*- Representaci\[OAcute]n:*)
(*	\[CenterDot] Difusi\[OAcute]n: ATP vs D*)
(*	\[CenterDot] Viscosidad: ATP vs \[Eta]*)
(*	*)
(*Devuelve las listas con la pendiente y el punto medio a\[NTilde]adido*)


processVidPT[fullpath_,listas_] := Module[
	{file,newFolder,b,r,vidRaw,bAdj,rAdj,ajuste,pAdj,
	imtrack,mask,ikp,traj,track,scale,ptoInt,ptoTraj,vidAdj,
	rRes,difTable,
	difFit,difPto,difPlot,difFitList,difPtoList,
	visFit,visPto,visPlot,visFitList,visPtoList},
	
	(*Video import*)
	Print["Video import and adjust"];
	{file,newFolder,r,b} = processImport[fullpath];
	CreateDirectory[newFolder];
	vidRaw = MapThread[ColorCombine[{#1, 0 #1, #2}]&,{r, b} ];
	Export[newFolder<>StringReplace[file,".nd2"->"_raw.gif"],Map[ImageAdjust,vidRaw]];
	
	(*Video adjust*)
	ajuste = 10;
	{bAdj,rAdj,pAdj} = processAdj[r,b,ajuste];
	
	(*Particle tracking*)
	Print["Particle tracking"];
	scale = 1/3;
	{imtrack,mask,ikp,traj,track} = processPTrack[pAdj,scale];
	ptoInt = HighlightImage[imtrack[[1]],ikp];
	Export[newFolder<>StringReplace[file,".nd2"->"_int_points.jpeg"],ptoInt];
	ptoTraj = Show[imtrack[[1]],ListPlot[traj,PlotStyle->Red,Joined->True]];
	Export[newFolder<>StringReplace[file,".nd2"->"_points_traj.jpeg"],ptoTraj];
	
	(*Adjusted video*)
	vidAdj = MapThread[ColorCombine[{#1, 0 #1, #2}]&,{Map[ImageResize[#,Scaled[scale]]&,rAdj]*mask, Map[ImageResize[#,Scaled[scale]]&,bAdj]*mask} ];
	Export[newFolder<>StringReplace[file,".nd2"->"_adj.gif"],Map[ImageAdjust,vidAdj]];
	
	(*Parameters*)
	Print["Parameters"];
	rRes = Map[ImageResize[#,Scaled[scale]]&,rAdj];
		(*pos,\[CapitalDelta]x,\[CapitalDelta]v,dif,atp*)
		difTable = Map[processDif[#,rRes,scale,mask]&,traj];
	Export[newFolder<>StringReplace[file,".nd2"->"_parameters.xlsx"],Grid[Join[{{"{x,y}","\[CapitalDelta]x","\[CapitalDelta]v","dif","atp"}},difTable]]];
		
	(*Representation*)
	{difFit, difPto, difPlot} = difRes[difTable];
		Export[newFolder<>StringReplace[file,".nd2"->"_dif_plot.jpeg"],difPlot];
		difFitList = Append[listas[[1]], {difFit[[1]],difFit[[2,1]]}];
		difPtoList = Append[listas[[3]], {difPto[[1]],difPto[[2,1]]}];
	{visFit, visPto, visPlot} = difVis[difTable];
		Export[newFolder<>StringReplace[file,".nd2"->"_vis_plot.jpeg"],visPlot];
		visFitList = Append[listas[[2]], {visFit[[1]],visFit[[2,1]]}];
		visPtoList = Append[listas[[4]], {visPto[[1]],visPto[[2,1]]}];
	
	Print["Done"<>file];	
	Return[{difFitList,visFitList,difPtoList,visPtoList}];
	
];


(* ::Subsection::Closed:: *)
(*Iteraci\[OAcute]n sobre todo el directorio*)


(*difFitList,visFitList,difPtoList,visPtoList*)
listas = {{},{},{},{}};


Do[
	listas = processVidPT[fullpath,listas],
	{fullpath,atpHoList[[1;;7]]}
]

Print["DONE"]


Export[folder<>StringReplace["240520_HeLa+ATP+Ho_.nd2",".nd2"->"_listas.dat"],listasAntiguas];
