_stdinOverLoadOptions = {};

SetDialogPrompt ("Listing of directory paths");
fscanf		(PROMPT_FOR_FILE,"Lines", _inDirectoryPaths);

pathCount = Columns (_inDirectoryPaths);

fprintf 	(stdout, "[READ ", pathCount , " file path lines]\n");

modelCount	   =  Abs(modelTypes);

MPI_NODE_STATUS = {MPI_NODE_COUNT-1,2};
resultsByFile	= {};

fileToExec 		= PATH_TO_CURRENT_BF + "AAPM_FEL.bf";

doneM = 0;

for (_fileLine = 0; _fileLine < Columns (_inDirectoryPaths); _fileLine += 1)
{
	resultsByFile [_inDirectoryPaths[_fileLine]] = {};
	
	_stdinOverLoadOptions = {};
	_stdinOverLoadOptions["00"] = "Universal";
	_stdinOverLoadOptions["01"] = "New Analysis";
	_stdinOverLoadOptions["02"] = _inDirectoryPaths[_fileLine];
	_stdinOverLoadOptions["03"] = "Custom";
	_stdinOverLoadOptions["04"] = "012345";
	_stdinOverLoadOptions["05"] = "y";
	_stdinOverLoadOptions["06"] = _inDirectoryPaths[_fileLine] + ".rev";
	_stdinOverLoadOptions["07"] = "Atchley";
	_stdinOverLoadOptions["08"] = "Independent properties";
	_stdinOverLoadOptions["09"] = "CF3x4";
	_stdinOverLoadOptions["10"] = "Positive selection testing";
	_stdinOverLoadOptions["11"] = _inDirectoryPaths[_fileLine] + ".aapm.fel";
	
	if (!_stdinOverLoadOptions["11"]) {
	    fprintf (stdout, "Already done with " + _stdinOverLoadOptions["02"] + "\n");
	    continue;
	}
	
	ExecuteAFile (fileToExec, _stdinOverLoadOptions);
}
