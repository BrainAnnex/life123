:: Quick start for Jupyter Lab (see IMPORTANT NOTE below: must change the folder name below!)

:: This batch file can be run, for example, by typing its name (exclusive of the .bat) in
::     the PyCharm TERMINAL window tab at the bottom (not to be confused with the Python console!)
::      On Win11, also need to prefix  ".\"          EXAMPLE:   .\quick

:: ****** FIRST CHANGE THE FOLDER NAME BELOW TO THE LOCATION ON YOUR MACHINE, and then uncomment
:: Add the root of the life123 library files to the value of the sys.path seen inside the execution of the notebooks
:: Note: only needed if you use a local copy of the life123 libraries;
::       NOT necessary if you pip install life123 in your venv
set PYTHONPATH=\Docs\- MY CODE (Win11)\Life123\Life123 development\

::@echo off
::echo *** Additional LOG FILES TO BE FOUND UNDER THE FOLDER: %PYTHONPATH%experiments
::echo on

:: Start Jupyter Lab (if a port other than the default 8888 is desired, use the option --port YOUR_PORT_NUMBER)
.\venv\Scripts\jupyter-lab