LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("PS_Plotters");
LoadFunctionLibrary ("WriteDelimitedFiles");
LoadFunctionLibrary ("chooseGeneticCode", {"0": "Universal"});

SetDialogPrompt ("Load a property-based FEL result file");
fscanf (PROMPT_FOR_FILE, "Raw", _fel_results);

_fel_results = Eval (_fel_results);


assert (Abs(_fel_results["SETTINGS"]) >= 3, "Missing analysis settings: " + _fel_results["SETTINGS"]);

fprintf (stdout, "\n[MODEL SETTINGS: '", (_fel_results["SETTINGS"])["Properties"], "' properties under the '", (_fel_results["SETTINGS"])["Metric"], "' weighting function]\n");

_is_positive_selection_test = Abs (_fel_results["NULL SETTINGS"]) >= 3;
if (_is_positive_selection_test) {
    fprintf (stdout, "[Positive selection testing results]\n"); 
} else {
    fprintf (stdout, "[Individual component testing results]\n"); 
}

_site_count = Abs (_fel_results) - 1 - _is_positive_selection_test;

fprintf (stdout, "[Loaded information on ", _site_count, " sites]\n");

assert (_site_count > 0, "No information about individual tested sites");

_test_names = {0,0};
_test_count = 0;

_multipliers = (_fel_results["SETTINGS"])["Multipliers"];
assert (Abs(_multipliers) > 0, "Missing rate matrix information");


/*SetDialogPrompt ("Write a detailed graphical report in PostScript to:");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN, 	_HYPSTextCommands(0), "\n/box {\n0 20 rlineto \n20 0 rlineto \n0 -20 rlineto \nclosepath } def\n");
_psFile = LAST_FILE_PATH;
*/

for (_site = 0; _site < _site_count; _site += 1) {
    _site_info = _fel_results[_site];
    

    if (_site_info ["CONSTANT"] == 1) {
        continue;
    }   
    
    if (Rows(_test_names) == 0) { // get test names
        _test_names = _sortStrings(Rows (_site_info));
        _test_count = Rows(_test_names);
        _variable_names = Rows(((_site_info["Full model"])["MLES"]));
        _resultColumns  = 3+Columns(_variable_names)+_test_count-1; 
            /* site 
               LogL
               DF
               MLE for each variable under the full model
               individual test 
            */
            
        _headerNames = {1, _resultColumns};
        _headerNames[0] = "Site";
        _headerNames[1] = "LogL";
        _headerNames[2] = "ParamCount"; _slider = 3;
        for (_column_count = 0; _column_count < Columns(_variable_names); _column_count += 1) {
            _headerNames[_slider] = _variable_names[_column_count];
            _slider += 1;
        }
        for (_column_count = 0; _column_count < _test_count; _column_count += 1) {
            if (_test_names[_column_count] == "Full model") {
                continue;
            }
            _headerNames[_slider] = "Corrected p-value for " + _test_names[_column_count];
            _slider += 1;
        }
        
        _resultMatrix = {_site_count, _resultColumns};
        for (_k = 0; _k  < _site_count; _k += 1) {
            _resultMatrix [_k][0] = _k + 1;
        }
        
    }
    
     _resultMatrix [_site][0] = _site + 1;

    ((_site_info["Full model"])["MLES"])["set_variable_values"][""];
    _numeric = {};
    _multipliers["compute_multiplers"][""];
    
    _fullModelLogL = (_site_info["Full model"])["LogL"];
    _fullModelDF   = (_site_info["Full model"])["DF"];
    
     _resultMatrix [_site][1] = _fullModelLogL;
     _resultMatrix [_site][2] = _fullModelDF;
     
     
    _test_pvalues = {_test_count, 2}["_MATRIX_ELEMENT_ROW_"];
    _slider = 3;
    for (_column_count = 0; _column_count < Columns(_variable_names); _column_count += 1) {
        _resultMatrix[_site][_slider] = ((_site_info["Full model"])["MLES"])[_variable_names[_column_count]];
        _slider += 1;
    }
    
    for (_aTest = 0; _aTest < _test_count; _aTest += 1) {
        if (_test_names[_aTest] != "Full model") {
            _testModelLogL = (_site_info[_test_names[_aTest]])["LogL"];
            _testModelDF = (_site_info[_test_names[_aTest]])["DF"];
            if (_is_positive_selection_test) {
                _test_pvalues [_aTest][0] = 1-CChi2 (2*(_fullModelLogL-_testModelLogL), 4);
            } else {
                _test_pvalues [_aTest][0] = 1-CChi2 (2*(_fullModelLogL-_testModelLogL), _fullModelDF-_testModelDF);            
            }
        } else {
            _test_pvalues [_aTest][0] = 1;
        }
    }
    
    _test_results_corrected = holm_multiple_testing_correction(_test_pvalues, _test_count-1) % 1;
    for (_aTest = 0; _aTest < _test_count; _aTest += 1) {
        if (_test_names[_aTest] != "Full model") {
            _resultMatrix[_site][_slider] = _test_results_corrected[_aTest][0];
            _slider += 1;
            if (_test_results_corrected[_aTest][0] <= 0.05) {
                fprintf (stdout, _test_names [_aTest], " is significant (corrected p = ", _test_results_corrected[_aTest][0], ") at site ", _site + 1, "\n");
            }
        }
    }
    //fprintf (_psFile,"\n", rateMatrixToPS (_numeric,"Site " + (1+_site),((_site_info["Full model"])["MLES"]),_test_names,_test_results_corrected));
}

fprintf (_psFile, CLOSE_FILE);
SetDialogPrompt ("Write a .csv result file to");
WriteSeparatedTable ("",_headerNames,_resultMatrix,0,",");


//------------------------------------------------------------------------------

function set_variable_values (key, value) {
    ExecuteCommands (key + " = " + value);
    return 0;
}

//------------------------------------------------------------------------------

function compute_multiplers (key, value) {
    _numeric [key] = Eval (value);
    
    return 0;
}

//------------------------------------------------------------------------------

function max_value (key, value) {
    _max = Max (value, _max);
    return 0;
}

//------------------------------------------------------------------------------

function holm_multiple_testing_correction (_test_pvalues, test_count) {
    _test_pvalues = _test_pvalues % 0;
    for (test_id = 0; test_id < test_count; test_id += 1) {
        _test_pvalues [test_id][0] = Min(1,_test_pvalues [test_id][0] * (test_count-test_id));
    }
    return _test_pvalues;
}

//------------------------------------------------------------------------------

function rateMatrixToPS (rate_list, title, MLES, test_names, test_p_values)
{
	
	psFigure = "";
	psFigure * 8192;
	
	table_rows = 2 + Abs (MLES) + Rows(test_p_values)-1;

	psFigure * _HYPSPageHeader (445,470 + 15*table_rows, "Protein Rate Matrix Plot");
	psFigure * "\n";
	psFigure * _HYPSSetFont ("Times-Roman", 12);
	psFigure * ("\n 0 " + (15*table_rows) +" translate\n");
	
	offset = 24;
	
	_max = 1e-10; rate_list ["max_value"][""];
	maxVal = _max;
	
	reordering = {}; // Stanfel Classification
	reordering["A"] = 0;
	reordering["C"] = 1;
	reordering["G"] = 2;
	reordering["I"] = 3;
	reordering["L"] = 4;
	reordering["M"] = 5;
	reordering["P"] = 6;
	reordering["S"] = 7;
	reordering["T"] = 8;
	reordering["V"] = 9;
	reordering["D"] = 10;
	reordering["E"] = 11;
	reordering["N"] = 12;
	reordering["Q"] = 13;
	reordering["F"] = 14;
	reordering["W"] = 15;
	reordering["Y"] = 16;
	reordering["H"] = 17;
	reordering["K"] = 18;
	reordering["R"] = 19;
	
	rk = Rows (reordering);
	
	_allPairs = Rows (rate_list);
	for (h=0; h<20; h+=1) {
		label = rk[h];
		psFigure * ("0 setgray\n10 "+(offset+(19-h)*20+6)+" \n("+label+") centertext\n");
		psFigure * ("0 setgray\n"+(offset+410)+" "+(offset+(19-h)*20+6)+" ("+label+") centertext\n");
		psFigure * ("0 setgray\n"+(offset+h*20+10)+" 10 ("+label+") centertext\n");
		psFigure * ("0 setgray\n"+(offset+h*20+10)+" "+(405+offset)+" ("+label+") centertext\n");
    	    
	}
	for (_pair = 0; _pair < Abs (rate_list); _pair += 1) {
	    _thisPair = _allPairs [_pair];
	    h2 = reordering[_thisPair[0]];
        v2 = reordering[_thisPair[1]];  
        psFigure * ("newpath\n"+(offset+v2*20)+" "+(offset+(19-h2)*20)+" moveto\n");
        greyColor = 1-rate_list[_thisPair]/maxVal;
		psFigure * (""+greyColor+" setgray\nbox fill\n");
        psFigure * ("newpath\n"+(offset+h2*20)+" "+(offset+(19-v2)*20)+" moveto\n");
		psFigure * (""+greyColor+" setgray\nbox fill\n");
	}
	
	psFigure * ("\n"+offset+" "+offset+" translate\n0 setgray\nnewpath\n0 0 moveto\n0 400 lineto\n400 400 lineto\n400 0 lineto\n0 0 lineto\nstroke\nclosepath");
	
	psFigure * ("\n\nnewpath\n0 200 moveto\n200 200 lineto\n200 400 lineto\nstroke");
	psFigure * ("\n\nnewpath\n200 120 moveto\n200 200 lineto\n280 200 lineto\n280 120 lineto\nclosepath stroke");
	psFigure * ("\n\nnewpath\n280 120 moveto\n340 120 lineto\n340 60 lineto\n280 60 lineto\nclosepath stroke");
	psFigure * ("\n\nnewpath\n340 60 moveto\n400 60 lineto\n400 0 lineto\n340 0 lineto\nclosepath stroke");
	psFigure * ("\n225.5 425 (Residue exchangeabilities for " + title + ") centertext\n 0 -" + (15*table_rows) +" translate\n");
	tableText = {2,1};
	tableText [0] = "Shading indicates relative substitution rates (black = max, white = 0)";
	tableText [1] = "Amino-acids are grouped into 4 Stanfel classification clusters";
	cellBorders = {{0}{0}};
	psFigure * _HYPSTextTable (400,30,12,tableText,cellBorders);
	_paramCount = Abs (MLES);
	_paramNames = Rows (MLES);
	
	psFigure * ("\n0 "+ (30) +" translate\n");
	
	tableText = {_paramCount,2};
	cellBorders = {_paramCount, 2};
	for (_pair = 0; _pair < _paramCount; _pair += 1) {
	    tableText[_pair][0] = _paramNames[_pair];
	    tableText[_pair][1] = Format(MLES[_paramNames[_pair]], 8, 2);
	}

	psFigure * _HYPSTextTable (200,15*_paramCount,12,tableText,cellBorders);
	psFigure * ("\nshowpage");
	psFigure * 0;
	return psFigure;
}
