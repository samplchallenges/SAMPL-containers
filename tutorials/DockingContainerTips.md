## Tips for using command line programs
* Some common command line programs (such as AutoDock Vina) already have docker containers made by other people or organizations. It may be worth it to search for pre-made docker containers to inherit from or build off of. (see [AutoDock Vina Docker](https://hub.docker.com/r/taccsciapps/autodock-vina))
* Some common command line programs may also have Python API's (see [AutoDock Vina API](https://pypi.org/project/vina/)) 
* If the above bullet doesn't work, you can install the command line program into your container by copying the files into the container and running the installation steps in the Dockerfile
    * Please see [`SAMPL-league/examples/adv-base/Dockerfile`](https://github.com/samplchallenges/SAMPL-league/blob/84ec83f00a637f9c79b6d2e3a1a336ea91837b7a/examples/adv-base/Dockerfile#L14)
* To run a command line program from within a Python module, consider using [os.system()](https://docs.python.org/3/library/os.html?highlight=os%20system#os.system) or similar from the Python3 library
    * Please see [`SAMPL-league/examples/adv/autodock.py`](https://github.com/samplchallenges/SAMPL-league/blob/84ec83f00a637f9c79b6d2e3a1a336ea91837b7a/examples/adv/autodock.py#L166)
