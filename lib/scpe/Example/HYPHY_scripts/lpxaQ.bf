REPLACE_TREE_STRUCTURE = 1;
DataSet lpxa = ReadDataFile ("./lpxa70.phy");

tree ="((LPXA_ECO57,((LPXA_SALTI,LPXA_SALTY),((LPXA_YERPE,LPXA_YEREN), (((LPXA_WIGBR,CAD83354),LPXA_PROMI),((AAP96037,(LPXA_HAEIN,LPXA_PASMU)), ((LPXA_VIBPA,(LPXA_VIBCH,LPXA_VIBVU)),(Q9AQK4,(((((LPXA_XYLFT, LPXA_XYLFA),(LPXA_XANCP,LPXA_XANAC)),(AAQ20846,(LPXA_NEIMB,LPXA_NEIMA))), (LPXA_RALSO,(CAE41721,(CAE36841,CAE33110)))),((Q8EGG3,(Q88MG8, (Q886N1,LPXA_PSEAE))),((((((LPXA_CHLTR,LPXA_CHLMU),(LPXA_CHLCV, (AAP98605,LPXA_CHLPN))),Q8DGT9),((((LPXA_CYACA,Q85G74),CAD78188), ((AAQ00460,(CAE07073,CAE21585)),((Q8DI00,Q8YUR3),LPXA_SYNY3))), LPXA_CHLTE)),((Q8EZA6,LPXA_CHRVI),(Q8A016,LPXA_AQUAE))),(Q820W6, ((LPXA_FUSNN,EAA23397),(((LPXA_RICPR,(LPXA_RICCN,(EAA25923,LPXA_RICRI))), ((LPXA_RHILO,(Q8VQ21,(LPXA_BRUAB,LPXA_BRUME))),(((LPXA_AGRT5, LPXA_RHIME),Q89KQ5),LPXA_CAUCR))),(LPXA_CAMJE,(AAP77779,(LPXA_HELPJ, LPXA_HELPY)))))))))))))))),LPXA_ECOLI)";


/*Define and optimize H0 model= JTT*/
    #include "Jones1.mdl";
    lk=0;
    for(i=0;i<262;i=i+1){
        ExecuteCommands("Tree treeh0"+i+" = tree;");
        ExecuteCommands("DataSetFilter classh0"+i+" = CreateFilter (lpxa,1,\""+i+"\",\"\");");
        ExecuteCommands("LikelihoodFunction lkh0"+i+"= (classh0"+i+",treeh0"+i+");");
        ExecuteCommands("Optimize (paramsh0"+i+", lkh0"+i+");");
        ExecuteCommands("lk = paramsh0"+i+"[1][0] + lk;");
        ExecuteCommands("fprintf(lpxa-lrt-site.out,paramsh0"+i+"[1][0],\"\n\");");

   }

        fprintf (lpxa-lrt-site.out, "SCPE Log(L) H0 = ", lk,"\n");

/*Define and optimize H1 model= SCPE*/
    lk=0;
    for(i=0;i<262;i=i+1){
        ExecuteCommands("#include \"lpxa-"+i+".dat\";");
        ExecuteCommands("Tree treeh1"+i+" = tree;");
        ExecuteCommands("DataSetFilter classh1"+i+" = CreateFilter (lpxa,1,\""+i+"\",\"\");");
        ExecuteCommands("LikelihoodFunction lkh1"+i+"= (classh1"+i+",treeh1"+i+");");
        ExecuteCommands("Optimize (paramsh1"+i+", lkh1"+i+");");
        ExecuteCommands("lk = paramsh1"+i+"[1][0] + lk;");
        ExecuteCommands("fprintf(lpxa-lrt-site.out,paramsh1"+i+"[1][0],\"\n\");");

   }

        fprintf (lpxa-lrt-site.out, "SCPE Log(L) H1 = ", lk,"\n");

/*Parametric bootstraping*/

/*Simulation  H0=Monomeric   */

for (simCounter = 1; simCounter<=300; simCounter = simCounter+1)
{

  lkh0=0;
  lkh1=0;
  for (site = 0; site<262; site = site+1) {

    ExecuteCommands("DataSet simulated"+site+" = SimulateDataSet (lkh0"+site+");");
     
    ExecuteCommands("DataSetFilter classh0s"+site+" = CreateFilter (simulated"+site+",1);");
    ExecuteCommands("LikelihoodFunction lkh0"+site+"= (classh0s"+site+",treeh0"+site+");");
    ExecuteCommands("Optimize (paramsh0s"+site+", lkh0"+site+");");
    ExecuteCommands("fprintf(lpxa-lrt-datos-H0.out,paramsh0s"+site+"[1][0],\"\n\");");

    ExecuteCommands("DataSetFilter classh1s"+site+" = CreateFilter (simulated"+site+",1);");
    ExecuteCommands("LikelihoodFunction lkh1"+site+"= (classh1s"+site+",treeh1"+site+");");
    ExecuteCommands("Optimize (paramsh1s"+site+", lkh1"+site+");");
    ExecuteCommands("fprintf(lpxa-lrt-datos-H1.out,paramsh1s"+site+"[1][0],\"\n\");");


   }

    for(i=0;i<262;i=i+1) ExecuteCommands("lkh0 = paramsh0s"+i+"[1][0] + lkh0;");
    for(i=0;i<262;i=i+1) ExecuteCommands("lkh1 = paramsh1s"+i+"[1][0] + lkh1;");

    fprintf(lpxa-lrt-site.out,simCounter," ",lkh0," ",lkh1," ", 2*(lkh1-lkh0), "\n");

}
