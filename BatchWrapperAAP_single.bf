_stdinOverLoadOptions = {};

SetDialogPrompt ("Listing of directory paths");
fscanf		(PROMPT_FOR_FILE,"Lines", _inDirectoryPaths);

pathCount = Columns (_inDirectoryPaths);

fprintf 	(stdout, "[READ ", pathCount , " file path lines]\n");

modelCount	   =  Abs(modelTypes);

MPI_NODE_STATUS = {MPI_NODE_COUNT-1,2};
resultsByFile	= {};

fileToExec 		= PATH_TO_CURRENT_BF + "AAPM.bf";

doneM = 0;

using_propset       = "Atchley";
function_form       = "Interacting";

for (_fileLine = 0; _fileLine < Columns (_inDirectoryPaths); _fileLine += 1)
{
    IS_TREE_PRESENT_IN_DATA = 0;
	resultsByFile [_inDirectoryPaths[_fileLine]] = {};
	
	_stdinOverLoadOptions = {};
	_stdinOverLoadOptions["00"] = "Universal";
	_stdinOverLoadOptions["01"] = _inDirectoryPaths[_fileLine];
	_stdinOverLoadOptions["02"] = using_propset;
	_stdinOverLoadOptions["03"] = function_form + " properties";
	_stdinOverLoadOptions["04"] = "CF3x4";
	_stdinOverLoadOptions["05"] = "012345";
	_stdinOverLoadOptions["06"] = "y";
	_stdinOverLoadOptions["07"] = _inDirectoryPaths[_fileLine] + ".`using_propset`.`function_form`.out";
	
	modelDesc = ""; 
	
	if (!_stdinOverLoadOptions["07"]) {
	    fprintf (stdout, "Already done with " + _stdinOverLoadOptions["01"] + "\n");
	    continue;
	}
	
	ExecuteAFile (fileToExec, _stdinOverLoadOptions);
	
}
