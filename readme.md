This is the github repository of the bachelor thesis "Untersuchung tryptischer Peptide von Intrinsically Disordered-Proteinen in Proteom-Datenbanken und ihre Charakterisierung in Massenspektrometrie-Datensätzen".
The code and data used for the thesis can be found in this repository. <br>
There is a little help for the first use of the program in german in the appendix of the thesis. Other than that I tried to make the "help" option of the program as helpful as possible.
Therefore you should try using "DisorderAnalyzer.py --help" or "DisorderAnalyzer.py {analyze/update} --help" to get some extra information. <br>

The TextExporterConv.py program was written to create a tsv-file that can be used by the DisorderAnalyzer. It just converts the necessary columns in the required format and writes the result to a new file called "<old_filename>_mod.tsv" or "<old_filename>_mod.csv" depending on the input file.