LoadFunctionLibrary         ("chooseGeneticCode");
LoadFunctionLibrary         ("GrabBag");
SetDialogPrompt             ("A codon-alignment file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
sample_size                = Log(filteredData.species * filteredData.sites);

ChoiceList (_synRateVariation,"Synonymous rate variation",1,SKIP_NONE,
            "None","Synonymous substitution rates do not vary from site to site",
            "Yes","Model the variation in synonymous substitution rates using a general discrete distribution");

if (_synRateVariation) {
    _synRateVariation = prompt_for_a_value ("How many synonymous rate classes",3,1,8,1);
}

ASSUME_REVERSIBLE_MODELS = 1;

settingsByClass = {};
ExecuteAFile      ("QCAP.mdl");
LoadFunctionLibrary ("queryTree");

SetDialogPrompt ("Save results to: ");
fprintf         (PROMPT_FOR_FILE, CLEAR_FILE);
base_save_path = LAST_FILE_PATH;

LikelihoodFunction lf_qcap =  (filteredData, givenTree);

current_BIC = fitAndSpool (base_save_path);

class_count = 1;

frequency_multiplier = {"0":"1"};
lf_components = {};

lf_components + "filteredData";
lf_components + "givenTree";

do {
    BIC_to_beat = current_BIC;
    
    
    class_count += 1;
    fprintf (stdout, "\n[", class_count, " CLASSES]\n");
    
    MULTIPLY_BY_FREQS = PopulateModelMatrix ("QCAP" + class_count, observedFreq, freqType, "class"+class_count+".",3,0,0);
    
    class_str = "" + class_count;
    
    ExecuteCommands ("
    Model QCAPModel`class_str` = (QCAP`class_str`,vectorOfFrequencies,MULTIPLY_BY_FREQS);
    
    for (k = 0; k < Abs (_qcap_settings[\"Property Parameters\"]); k+=1) {
        ExecuteCommands ((_qcap_settings[\"Property Parameters\"])[k] + \"= 0.\");
    }
    
    settingsByClass + _qcap_settings;
    
    Tree tree`class_str` = treeString;
    ReplicateConstraint (\"this1.?.?:=this2.?.?\", tree`class_str`, givenTree);
    global class_weight`class_str` = 0.95; class_weight`class_str` :< 1;
    ");
    
    frequency_multiplier [class_count-1] = frequency_multiplier [class_count-2];
    frequency_multiplier [class_count-2] += "*class_weight`class_str`";
    frequency_multiplier [class_count-1] += "*(1-class_weight`class_str`)";
 
    lf_components + "filteredData";
    lf_components + "tree`class_str`";
    
    mixingString = {class_count,1};
    for (k = 0; k < class_count; k += 1) {
        mixingString[k] = frequency_multiplier[k] + "*SITE_LIKELIHOOD[" + k + "]";
    }
    
    ExecuteCommands ("LikelihoodFunction lf_qcap = (" + Join (",",lf_components) + ",\"Log(" + Join ("+", mixingString) + ")\")");
    
    current_BIC = fitAndSpool (base_save_path + "." + class_count);
    
}
while (BIC_to_beat > current_BIC);


//------------------------------------------------------------------------------

function fitAndSpool (save2) {
    Optimize (res, lf_qcap); 
    fprintf (stdout, lf_qcap);
    results = {};
    results ["SETTINGS"]         = settingsByClass[0];
    
    for (k = 0; k < Abs (settingsByClass); k+=1) {
        mles = {};
        ((settingsByClass[k])["Property Parameters"])["export_mles"][""];
        results [""+k] = {"Global model": {"LogL" : res[1][0], "DF":  res[1][1], "MLES": mles,
                            "Sample size" : Log(filteredData.species)+Log(filteredData.sites)}};
    }
    
    fprintf (save2, CLEAR_FILE, results);
    BIC = -2*res[1][0] + sample_size * res[1][1];
    
    
    fprintf (stdout, 
"
LogL = ", res[1][0],"
degF = ", res[1][1],"
samp = ", sample_size,"
BIC  = ", BIC, "\n");

    return BIC;
}

//------------------------------------------------------------------------------

function export_mles (key, value) {
    mles[value^{{"class[0-9]+\\.",""}}] = Eval (value);
    return 0;
}

