LoadFunctionLibrary ("chooseGeneticCode");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);

ASSUME_REVERSIBLE_MODELS = 1;

settingsByClass = {};

ExecuteAFile ("QCAP.mdl");
settingsByClass + _qcap_settings;
LoadFunctionLibrary ("queryTree");

MULTIPLY_BY_FREQS = PopulateModelMatrix ("QCAP2", observedFreq, freqType, "class2.",0,0);
Model QCAPModel2 = (QCAP2,vectorOfFrequencies,MULTIPLY_BY_FREQS);
settingsByClass + _qcap_settings;

Tree tree2 = treeString;
ReplicateConstraint ("this1.?.?:=this2.?.?", tree2, givenTree);
global class_weight = 0.5; class_weight :< 1;
LikelihoodFunction lf_qcap = (filteredData, givenTree, filteredData, tree2, "Log(class_weight*SITE_LIKELIHOOD[0]+(1-class_weight)*SITE_LIKELIHOOD[1])");

MULTIPLY_BY_FREQS = PopulateModelMatrix ("QCAP3", observedFreq, freqType, "class3.",0,0);
Model QCAPModel3 = (QCAP3,vectorOfFrequencies,MULTIPLY_BY_FREQS);
settingsByClass + _qcap_settings;
Tree tree3 = treeString;
ReplicateConstraint ("this1.?.?:=this2.?.?", tree3, givenTree);
global class_weight2 = 0.5; class_weight2 :< 1;
LikelihoodFunction lf_qcap = (filteredData, givenTree, filteredData, tree2, filteredData, tree3, "Log(class_weight*SITE_LIKELIHOOD[0]+(1-class_weight)*class_weight2*SITE_LIKELIHOOD[1]+(1-class_weight)*(1-class_weight2)*SITE_LIKELIHOOD[2])");

MULTIPLY_BY_FREQS = PopulateModelMatrix ("QCAP4", observedFreq, freqType, "class4.",0,0);
Model QCAPModel4 = (QCAP4,vectorOfFrequencies,MULTIPLY_BY_FREQS);
settingsByClass + _qcap_settings;
Tree tree4 = treeString;
ReplicateConstraint ("this1.?.?:=this2.?.?", tree4, givenTree);
global class_weight3 = 0.5; class_weight3 :< 1;
LikelihoodFunction lf_qcap = (filteredData, givenTree, filteredData, tree2, filteredData, tree3, filteredData, tree4, 
            "Log(class_weight*SITE_LIKELIHOOD[0]+(1-class_weight)*class_weight2*SITE_LIKELIHOOD[1]+(1-class_weight)*(1-class_weight2)*class_weight3*SITE_LIKELIHOOD[2]+(1-class_weight2)*(1-class_weight3)*SITE_LIKELIHOOD[3])");

VERBOSITY_LEVEL = 1;
Optimize (res, lf_qcap);

fprintf (stdout, lf_qcap);

results = {};
results ["SETTINGS"]         = settingsByClass[0];

for (k = 0; k < Abs (settingsByClass); k+=1) {
    mles = {};
    ((settingsByClass[k])["Property Parameters"])["export_mles"][""];
    results [""+k] = {"Global model": {"LogL" : res[1][0], "DF":  res[1][1], "MLES": mles}};
}

SetDialogPrompt ("Save the results to:");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, results);

//------------------------------------------------------------------------------

function export_mles (key, value) {
    mles[value^{{"class[0-9]+\\.",""}}] = Eval (value);
    return 0;
}

