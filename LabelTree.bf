LoadFunctionLibrary ("TreePartitioning.ibf");

SetDialogPrompt ("Tree");
fscanf (PROMPT_FOR_FILE, "Tree", T);

fprintf (stdout, "Regexp:");
fscanf  (stdin, "String", regexp);

SetDialogPrompt ("Labeled branches written to:");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, Join("\n", labelBranchesByRegExp ("T", regexp, "tagged")));