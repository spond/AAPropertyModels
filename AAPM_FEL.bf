MESSAGE_LOGGING = 0;

/******************************************************************************************************/
// A LOOONG SETUP PREAMBLE, STANDARD FOR FEL/SLAC etc       
/******************************************************************************************************/

    LoadFunctionLibrary ("GrabBag");
    ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "qndhelper1.ibf");
    modelDesc = ModelTitle; 
        // ModelTitle comes from qndhelper1; setting modelDesc here will stop QCAP.mdl from propmting for it
    ExecuteAFile        ("QCAP.mdl"); 
    
    ChoiceList (_fel_test_type,"Which FEL tests should be run?",1,SKIP_NONE,
            "Individual component testing","FEL will be run to test which of the property weights matter at a given site (i.e. are significantly different from 0)",
            "Positive selection testing","The unconstrained model will be tested against the model that limits omega <= 1");
            
    if (_fel_test_type < 0) {
        return -1;
    }

    fprintf (stdout, "Using the ", SELECTION_STRINGS , " metric\n");    
    
    REPLACE_TREE_STRUCTURE = 1;
    rOptions = 4;
    LoadFunctionLibrary("qndhelper2")
    
    saveNucs = {{AC__,AT__,CG__,CT__,GT__}};
                                    
    ClearConstraints (codonTree,AC,AT,CG,CT,GT); 
    LoadFunctionLibrary ("CodonTools.def");
    
    ASSUME_REVERSIBLE_MODELS = 1;
     
    OPTIMIZATION_METHOD     = 0;
    USE_LAST_RESULTS        = 1;
    
    fprintf (stdout, "\n[RETUNING BRANCH LENGTHS AND NUCLEOTIDE RATES UNDER THE CODON MODEL]\n");
    T0 = Time(1);
    Optimize (codonLF, lf);
    OPTIMIZATION_TIME_HARD_LIMIT = (Time(1)-T0)*4;
    
    fprintf (stdout, "IMPROVED Log(L) BY ", codonLF[1][0]-resC[1][0], " POINTS\n");
    
    OPTIMIZATION_METHOD = 4;
    SHORT_MPI_RETURN    = 1;
    USE_LAST_RESULTS    = 0;

/******************************************************************************************************/
// DEFINE ALL THE INDIVIDUAL PARAMETERS THAT WILL BE TESTED AT EACH SITE HERE       
/******************************************************************************************************/

SetDialogPrompt     ("Save raw results to");
fprintf             (PROMPT_FOR_FILE, CLEAR_FILE);
_felRawResultFile = LAST_FILE_PATH;

if (_fel_test_type == 0) {    
    _felTests = {
    "Full model"   : {"alpha_0" : 1, "alpha_1" : 1, "alpha_2": 1, "alpha_3" : 1, "alpha_4" : 1}, 
        // only these parameters (and the default syn rate scaler) will be fitted at a site; all others will be constrained to current values
    "Test property 1" : {"Parameter": "alpha_0", "Constraint" : "alpha_0 := 0", "Cleanup" : "alpha_0 = 1"},
    "Test property 2" : {"Parameter": "alpha_1", "Constraint" : "alpha_1 := 0", "Cleanup" : "alpha_1 = 1"},
    "Test property 3" : {"Parameter": "alpha_2", "Constraint" : "alpha_2 := 0", "Cleanup" : "alpha_2 = 1"},
    "Test property 4" : {"Parameter": "alpha_3", "Constraint" : "alpha_3 := 0", "Cleanup" : "alpha_3 = 1"},
    "Test property 5" : {"Parameter": "alpha_4", "Constraint" : "alpha_4 := 0", "Cleanup" : "alpha_4 = 1"}
    };
    
    results = runFELTestsOnPartition ("filteredData", "QCAPModel", treeString, "codonTree", "synRate", "synRate", _felTests, 1, _felRawResultFile);
    results["SETTINGS"] = _qcap_settings;
    
} else {

    fullSettings =  _qcap_settings;
    MULTIPLY_BY_FREQS = PopulateModelMatrix ("QCAPNeutNeg", observedFreq, freqType, "", 0, 1 ); // the last "1" argument forces the model to cap 'omega' at 1.

    Model QCAPModelNeutNeg = (QCAPNeutNeg,vectorOfFrequencies,MULTIPLY_BY_FREQS);

    _felTests = {
    "Full model"   : {"alpha_0" : 1, "alpha_1" : 1, "alpha_2": 1, "alpha_3" : 1, "alpha_4" : 1}, 
        // only these parameters (and the default syn rate scaler) will be fitted at a site; all others will be constrained to current values
    };
    
    nullModel  = runFELTestsOnPartition ("filteredData", "QCAPModelNeutNeg", treeString, "codonTree", "synRate", "synRate", _felTests, 1, _felRawResultFile);
    results    = runFELTestsOnPartition ("filteredData", "QCAPModel", treeString, "codonTree", "synRate", "synRate", _felTests, 1, _felRawResultFile);
    
    for (_copyVarValue = 0; _copyVarValue < Abs(nullModel); _copyVarValue += 1) {
        _thisParameter = nullModel ["INDEXORDER"][_copyVarValue];
        (results [_thisParameter])["Test positive selection"] = (nullModel[_thisParameter])["Full model"];
    }
    results["SETTINGS"]      = fullSettings;
    results["NULL SETTINGS"] = _qcap_settings;
}
fprintf (_felRawResultFile, CLEAR_FILE, results);

/******************************************************************************************************/

function runFELTestsOnPartition (filterName, modelName, treeString, copyLengthsFrom, blParameter1, blParameter2, testsToRun, verbose, resultFile) {

    FEL._optimizationTaskInformation = {};
 
    assert (Type (testsToRun["Full model"]) == "AssociativeList", "Incorrectly specified test parameters, missing 'Full model'");
    (testsToRun["Full model"])["_felScaler"] = 1;
    FEL._fullParameterList = testsToRun["Full model"];
    
    _testNames = Rows (testsToRun);

    ExecuteCommands ("UseModel (`modelName`);
                      Tree	 _felSiteTree = treeString;
                      global _felScaler   = 1;
                      ReplicateConstraint (\"this1.?.`blParameter1`:=_felScaler*this2.?.`blParameter2`__\",_felSiteTree,`copyLengthsFrom`);");
                      
    _site_count = Eval ("`filterName`.sites");   
    
    ExecuteCommands ("GetDataInfo (_mapSitesToUnuiquePatterns, `filterName`);
                      GetDataInfo (_filterParameters, `filterName`, \"PARAMETERS\");
                     ");
    
    //_filterParameters ["EXCLUSIONS"]

    FEL._felUniqueSitePatternResults = {};
    _alreadyProcessedPatterns    = {};

    if (verbose) {
        fprintf (stdout, "[Running FEL tests on ", _site_count, " sites]\n");
    }
    
    if (_site_count == 0) {
        return FEL._felUniqueSitePatternResults;
    }
    
    for (node_id = 0; node_id < Max (1, MPI_NODE_COUNT); node_id += 1) {
        FEL._optimizationTaskInformation [node_id] = {};
    }

    /* make the likelihood function here */
    
    COUNT_GAPS_IN_FREQUENCIES = 0;
    _hasLFBeenSetup = 0;
    
    for (_site_index = 0; _site_index < filteredData.sites; _site_index += 1) {
        _site_pattern_index = _mapSitesToUnuiquePatterns [_site_index];
        if (verbose) {
            fprintf (stdout, "...[Working on site ", _site_index, " (site pattern ", _site_pattern_index, ")]\n");
        }
        
        if (_alreadyProcessedPatterns[_site_pattern_index] > 0) { 
            continue; 
        } 
        
        FEL._felUniqueSitePatternResults [_site_pattern_index] = {};
        _alreadyProcessedPatterns    [_site_pattern_index] += 1;
        
    // check to see if the site is variable 
        filterString = "" + _site_index*3 + "-" + (_site_index*3+2);
        ExecuteCommands ("DataSetFilter  _felSiteFilter = CreateFilter (`filterName`,3,\"`filterString`\",\"\",_filterParameters [\"EXCLUSIONS\"])");
        HarvestFrequencies (_siteFreqs, _felSiteFilter, 3, 3, 0);
        
    // skip invariable sites
        if (+_siteFreqs["_MATRIX_ELEMENT_VALUE_>0"] <= 1) {
            if (verbose) {
                fprintf (stdout, ".....This site is invariable; models will not be fitted for this site\n");
            }
            (FEL._felUniqueSitePatternResults [_site_pattern_index]) ["CONSTANT"] = 1;
            continue;
        }  
        
    // only setup the likelihood function once
        if (_hasLFBeenSetup == 0) {
            if (verbose) {
                fprintf (stdout, "...[Setting up the FEL likelihood function]\n");
            }
            
            LikelihoodFunction FEL._felSiteLF = (_felSiteFilter, _felSiteTree);
            GetString (lfInfo, FEL._felSiteLF, -1);
            
            FEL._siteGlobalParameters    = lfInfo["Global Independent"];
            
            for (_global_var = 0; _global_var < Columns (FEL._siteGlobalParameters); _global_var += 1) {
                _init_value = (testsToRun["Full model"])[FEL._siteGlobalParameters[_global_var]];
                if (_init_value == 0) {
                    if (verbose) {
                        fprintf (stdout, "......Constraining ",FEL._siteGlobalParameters[_global_var],"\n");
                    }                    
                    ExecuteCommands (FEL._siteGlobalParameters[_global_var] + ":=" + FEL._siteGlobalParameters[_global_var] + "__");
                } else {
                    ExecuteCommands (FEL._siteGlobalParameters[_global_var] + "=" + _init_value);
                }
            }
            
            _siteGlobalParametersAVL = stringMatrixToAVL ("FEL._siteGlobalParameters");
            for (_testName = 0; _testName < Abs (testsToRun); _testName += 1) {
                _thisTest = _testNames [_testName];
                if (_thisTest != "Full model") {
                    assert (Type (testsToRun[_thisTest]) == "AssociativeList", "Incorrectly specified test parameters, `_thisTest` is not a properly formed test");
                    assert (_siteGlobalParametersAVL[(testsToRun[_thisTest])["Parameter"]] > 0, "Test parameter '" + (testsToRun[_thisTest])["Parameter"] + "' is not a part of the model/likelihood function");
                }
            }
                        
            _hasLFBeenSetup = 1;
        }
        
        for (_testName = 0; _testName < Abs (testsToRun); _testName += 1) {
            _thisTest = _testNames [_testName];
            _felSendAJob (_site_pattern_index, _thisTest, testsToRun[_thisTest], verbose);
        }
        
        fprintf (resultFile, CLEAR_FILE, FEL._felUniqueSitePatternResults);
    }    
    
    _stillToDo = 0;
    
    for (_leftOver = 0; _leftOver < Abs (FEL._optimizationTaskInformation); _leftOver += 1) {
        _stillToDo += Abs (FEL._optimizationTaskInformation[_leftOver]) > 0;
    }  
    
    for (; _stillToDo>0; _stillToDo += -1) {
        _felReceiveAJob (verbose);
    }
    
    
    
    FEL._felSiteResults = {};
    for (_site_index = 0; _site_index < filteredData.sites; _site_index += 1) {
        FEL._felSiteResults[_site_index] = FEL._felUniqueSitePatternResults[_mapSitesToUnuiquePatterns [_site_index]];        
    }
    
    return FEL._felSiteResults;                 
    
}

/******************************************************************************************************/
function _felSendAJob (_patternID, _testName, _testSpecification, verbose) {
    if (verbose) {
        fprintf (stdout, "...[Running ", _testName, " on site pattern ", _patternID, "]\n");
    }
    
    if (MPI_NODE_COUNT <= 1) {
        _nodeID = 0;
    } else {
        for (_nodeID = 0; _nodeID < MPI_NODE_COUNT-1; _nodeID += 1) {
            if (Abs(FEL._optimizationTaskInformation [_nodeID])==0) {
                break;
            }
        }
        if (_nodeID == MPI_NODE_COUNT - 1) {
            _nodeID = _felReceiveAJob (verbose);
        }
    }
    FEL._optimizationTaskInformation [_nodeID] = {"Pattern" : _patternID, "Name": _testName};
    if (Type (_testSpecification["Constraint"]) == "String") {
        ExecuteCommands (_testSpecification["Constraint"]);
    }
    if (MPI_NODE_COUNT <= 1) {
        Optimize (FEL._felSiteLF_MLES, FEL._felSiteLF);
        
    
        LIKELIHOOD_FUNCTION_OUTPUT = 7;
        fprintf ("" + _patternID + "." + _testName, CLEAR_FILE, FEL._felSiteLF);
        
        
        FEL._felSiteLF_MLE_VALUES = {};
        for (_copyVarValue = 0; _copyVarValue < Abs(FEL._fullParameterList); _copyVarValue += 1) {
            _thisParameter = FEL._fullParameterList ["INDEXORDER"][_copyVarValue];
            FEL._felSiteLF_MLE_VALUES [_thisParameter] = Eval (_thisParameter);
        }
        _felReceiveAJob (verbose);
    } else {
        if (verbose) {
            fprintf (stdout, "......[Sending job to node ", _nodeID+1,"]\n");
        }
        
        /*
        LIKELIHOOD_FUNCTION_OUTPUT = 7;
        fprintf ("" + _patternID + "." + _testName, CLEAR_FILE, FEL._felSiteLF);
        */
        
        MPISend (_nodeID+1, FEL._felSiteLF);
    }
    if (Type (_testSpecification["Cleanup"]) == "String") {
        ExecuteCommands (_testSpecification["Cleanup"]);
    }
    
    return None;
    
} 

/******************************************************************************************************/

function _felReceiveAJob (verbose) {
    if (MPI_NODE_COUNT <= 1) {
        _fromID = 0;
        _fromNode = "";
    } else {
        MPIReceive (-1,_fromID,_serializedResult);
        ExecuteCommands (_serializedResult);
        //fprintf (stdout, _serializedResult);
        _fromID = _fromID - 1;
        _fromNode = " from node " + _fromID;
    }
    
    finishedTestName        = (FEL._optimizationTaskInformation [_fromID])["Name"];
    finishedSitePattern     = 0 + (FEL._optimizationTaskInformation [_fromID])["Pattern"];
    
    if (verbose) {
        fprintf (stdout, "......[Received ", finishedTestName," on site pattern ", finishedSitePattern , _fromNode, "]\n");
    }
    
    FEL._optimizationTaskInformation [_fromID] = {};
    (FEL._felUniqueSitePatternResults[finishedSitePattern])[finishedTestName] = {};
    ((FEL._felUniqueSitePatternResults[finishedSitePattern])[finishedTestName])["LogL"] = FEL._felSiteLF_MLES [1][0];
    ((FEL._felUniqueSitePatternResults[finishedSitePattern])[finishedTestName])["DF"] = FEL._felSiteLF_MLES [1][1];
    _MLES = {};
    
    for (_copyVarValue = 0; _copyVarValue < Abs(FEL._fullParameterList); _copyVarValue += 1) {
        _thisParameter = FEL._fullParameterList ["INDEXORDER"][_copyVarValue];
        _MLES [_thisParameter] = FEL._felSiteLF_MLE_VALUES[_thisParameter];
    }
    ((FEL._felUniqueSitePatternResults[finishedSitePattern])[finishedTestName])["MLES"] = _MLES;
    
    return _fromID;
} 

