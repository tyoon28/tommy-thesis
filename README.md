# tommy-thesis

Project file structure - orca and dynamic graphlet counting will have to be downloaded separately

thesis
|- tommy-thesis (this repo)
|- |- All scripts and other code
|- |- MD data (e.g. R1-30-closed, R1-15-closed)

|- orca
|- |- orca.exe
|- |- input
|- |- output
|- |- runorca.sh**

|- dynamic_graphlets
|- |- source code from https://www3.nd.edu/~cone/DG/ (the executable did not work on my machine)
|- |- input
|- |- output
|- |- run_count_tmp_graphlets.sh**

**included in this repo, you'll have to move it to the right directory for it to work


Many files in this repo are just tests or contain unused functions. The ones that ended up in the thesis are:

HELPER FUNCTIONS
dynamic_graphlets.py
graphlets.py
network_generation.py
network_statistics.py
visualization.py

VISUALIZATION
chimera_script.py
makebonds.tcl

RUNNING ANALYSIS
angles.py
analyze_orca_output.py (2)
contactprob.py
make_figures.py
make_orca_input.py (1)
run_dganalysis.py (2)
run_dyngraphlets.py (1)

These have to be run in order, but the rest should all run on their own:
(make_orca_input.py > runorca.sh > analyze_orca_output.py)
(run_dyngraphlets.py > run_count_tmp_graphlets.sh > run_dganalysis.py)

Apologies for the mess, any questions about anything feel free to contact tbyoon28@gmail.com