# Special Arguments Documentation

The Special Arguments allows you to add extra file arguments that will be passed to you container at run time, for example a license file. This implementation allows participants to keep information in file form separate from their container code. Non-executable files that are essential to run a container (such as a license file) can be stored elsewhere and only passed to the container when needed. This means the container on its own would be useless, and only those with the essential files can successfully run it. 

## Adding Special Arguments to your Container
* To use the Special Arguments option, your container must handle the required command line arguments for the challenge as well as your unique special argument(s).
* You may create as many special arguments as you want/need.
* Your argument must follow the conventions below:
    * only lower case letters, numbers, dashes/hyphens, or underscores are used in the naming convention of the argument
        * Examples: `--your-argument` or `--your_argument` NOT `--yourArgument`
    * Argument expects a path to a FILE as its value

```
import click
import os.path


# the following are decorators related to the MAIN FUNCTION, these '@click.command()' and '@click.option()' decorators
# MUST remain directly above your main function
@click.command()

@click.option("--required-arg1")
@click.option("--required-arg2")

@click.option("--output-dir",help="Output directory for receptor and docked_ligand files")

@click.option("--your-argument", type=click.Path(exists=True), help="Any special file arguments you")
def docking_main(required_arg1, required_arg2,  output_dir, your_argument) -> None:
    pass
```
## Other Important Information
* Please remember either copy the file(s) passed to the special arguments into your container's `/opt/app` folder at run time or set your license path to check for the license file at the path specified.
