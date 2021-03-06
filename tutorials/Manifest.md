
# Manifest

## Files
* `README.md`: contains a step by step docking tutorial to introduce the idea of containerizing docking methods using Docker.
* `BuildYourOwnDockingContainer.md`: details how to adapt the tutorial in `README.md` to write your own SAMPL challenge docking container.
* `BuildYourOwnLogDContainer.md`: details how to adapt the tutorial in `README.md` to write your own SAMPL challenge LogD container.
* `CondaEnvInstructions.pdf`: details how to build a conda environment that is compatible with SAMPL challenges. Please read this file as traditional ways of installing the conda environment may cause issues with SMILES string parsing.
* `ever_giveUsage.md`: details how to use the pip installable package `ever_given` to run your SAMPL containers.

## Directories
* `data`: contains sample input files for docking.
* `ever_given`: contains `ever_given` wrapper files, which will eventually be packaged into a pip installable package
* `templates`: contains the files to use as starting points to write your own docking container.
