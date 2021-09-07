import argparse
from pathlib import Path

from ever_given import wrapper
""" --file-key will get passed into file_kwargs as --key= in file_kwargs
into evergiven; container will get --key=filename

Example command:
python ../ever_given/run.py robbason/calc-coords:latest --file-molfile tests/data/ChEBI_16716.mdl --output-keys conformation

(assuming you are in app directory and you have built the container in  example_container/coords as robbason/calc-coords using its build.sh)
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("container_uri")
    parser.add_argument("--command", default="")
    parser.add_argument("--output-keys", help="comma separated list of output keys to be loaded from container output", default="")
    _parsed, unknown = parser.parse_known_args()
    #print("unknown", unknown)
    if len(unknown) % 2 != 0:
        raise ValueError("Must pass key/value pairs")
    for idx in range(len(unknown)//2):
        arg = unknown[idx*2]
        parser.add_argument(arg)
    args = parser.parse_args()
    kwargs = {}
    file_kwargs = {}
    FILE_PREFIX = "file_"
    for key, value in vars(args).items():
        if key.startswith(FILE_PREFIX):
            file_kwargs[key[len(FILE_PREFIX):]] = value
        else:
            if key not in ("command", "container_uri", "output_keys"):
                kwargs[key] = value
    print("command?", args.command)
    print("file kwargs", file_kwargs)
    print("kwargs", kwargs)
    output_file_keys = args.output_keys.split(",")
    output_dir = "evergiven_output"
    print("Putting output into", output_dir)
    results = {}
    for k, v in wrapper.run(args.container_uri, command=args.command, file_kwargs=file_kwargs, kwargs=kwargs, output_dir=str(output_dir), output_file_keys=output_file_keys):
        results[k] = v
    print("Results:", results)



if __name__ == "__main__":
    main()
