
# Manifest

## Files
* `README.md`: contains a step by step docking tutorial to introduce the idea of containerizing docking methods using Docker.
* `BuildYourOwnDockingContainer.md`: details how to adapt the tutorial in `README.md` to write your own SAMPL challenge docking container.
* `BuildYourOwnLogDContainer.md`: details how to adapt the tutorial in `README.md` to write your own SAMPL challenge LogD container.
* `CondaEnvInstructions.pdf`: details how to build a conda environment that is compatible with SAMPL challenges. Please read this file as traditional ways of installing the conda environment may cause issues with SMILES string parsing.
* `ever_giveUsage.md`: details how to use the pip installable package `ever_given` to run your SAMPL containers.

## Directories
* `adv-base`: contains the files to build the AutoDock Vina base container, which has all the necessary environment package and command line program installations.
* `adv`: contains the files to build a run-able Autodock Vina docking container.
* `data`: contains sample input files for docking.
* `template`: contains the files to use as starting points to write your own docking container.
