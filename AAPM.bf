ExecuteAFile        (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "chooseGeneticCode.def");
LoadFunctionLibrary ("GrabBag");

DataSet ds                 = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);

ASSUME_REVERSIBLE_MODELS = 1;

ExecuteAFile        ("QCAP.mdl");
ExecuteAFile        (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "queryTree.bf");

LikelihoodFunction lf_qcap = (filteredData, givenTree);
Optimize (res, lf_qcap);
fprintf (stdout, lf_qcap);

results = {};
results ["SETTINGS"]         = _qcap_settings;

mles = {};
(_qcap_settings["Property Parameters"])["export_mles"][""];
results ["0"]                = {"Global model": {"LogL" : res[1][0], "DF":  res[1][1], "MLES": mles}};

SetDialogPrompt ("Save the results to:");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, results);

//------------------------------------------------------------------------------

function export_mles (key, value) {
    mles[value] = Eval (value);
    return 0;
}
