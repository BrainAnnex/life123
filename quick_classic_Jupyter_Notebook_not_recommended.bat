:: Quick start for Jupyter notebook.
:: ***NOT*** RECOMMENDED!  See "quicklab.bat" instead

:: This batch file can be run, for example, by typing its name (exclusive of the .bat) in
::     the PyCharm TERMINAL window tab at the bottom (not to be confused with the Python console!)

::IMPORTANT - FIRST CHANGE THE FILE NAMES BELOW TO WHAT YOU HAVE ON YOUR MACHINE!!

:: Add the root of the project files to the value of the sys.path seen inside the execution of the notebooks
set PYTHONPATH=\Docs\- MY CODE\BioSimulations\life123-Win7\

:: Start Jupyter notebook at a convenient startup folder;
:: warning: one will not be able to navigate to the parent of that startup folder!
:: Ditch the  --notebook-dir  argument if you want to see everything
jupyter notebook --notebook-dir="\Docs\- MY CODE\BioSimulations\life123-Win7\experiments\life_1D"