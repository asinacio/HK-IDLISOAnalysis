# TreeConverter

Script to take output from WCSim, and convert it to a text file for input to the main analysis.
Requires WCSim to be correctly setup and included. [Will add full instructions for that later]

To compile: (will change to Makefile later)

First source script to setup WCSim paths, then:
```
source buildTreeConverter.sh
```
_Will aim to replace this with a Makefile once I've figured out how to do that._


To run:
```
TreeConverter -i [inputfile.root] -o [outputfile.root] -m [macro.mac]
```
Possible arguments:
- `-i` : Input file
- `-o` : Output file
- `-m` : Macro file (the one used to run the input simulation)
- `-b` : Use B&L PMTs only
- `-v` : Print debug info
- `-h` : Display the help message

Only `-i` and `-m` are actually required for the converter to run, everything else will take default values.