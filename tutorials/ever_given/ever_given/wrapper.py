import copy
import io
import shlex
import sys
import threading
from pathlib import Path

import docker

GUEST_INPUT_DIR = Path("/mnt") / "inputs"
GUEST_OUTPUT_DIR = Path("/mnt") / "outputs"


def _guest_input_path(filename):
    return GUEST_INPUT_DIR / Path(filename).name


def _prepare_commandline(command, args_dict):
    return (
        command
        + " "
        + " ".join(
            [f"--{key} {shlex.quote(str(value))}" for key, value in args_dict.items()]
        )
    )


def _parse_output(host_path, raw_text, output_file_keys):
    result = {}
    for line in raw_text.decode("utf-8").splitlines():
        lineparts = line.split(maxsplit=1)
        if len(lineparts) == 2:
            key, value = lineparts
            value = value.strip()
            if output_file_keys and key in output_file_keys:
                print(host_path, key, value)
                result[key] = _convert_guest_to_host_path(host_path, value)
            else:
                result[key] = value
        else:
            raise ValueError(f"Cannot parse output {line}, needs KEY VALUE format")
    return result


def _convert_guest_to_host_path(host_path, filepath):
    filepath = GUEST_OUTPUT_DIR / filepath
    relative_path = filepath.relative_to(GUEST_OUTPUT_DIR)
    # so we don't care whether container outputs relative or absolute paths
    return str(Path(host_path) / relative_path)


def _convert_file_kwargs(file_kwargs):
    """
    file_kwargs has key: path, where path is on the host
    generate a dict of host: guest path mappings that can be used to generate bind commands and return a file_kwargs pointing to the files with corresponding paths that are accessible to the container
    """
    dirpaths = {}
    final_file_kwargs = {}
    for idx, (key, pathname) in enumerate(file_kwargs.items()):
        path = Path(pathname).resolve()
        dirpath = path.parent
        basename = path.name
        if dirpath not in dirpaths:
            dirpaths[dirpath] = Path("/mnt") / f"inputs{idx}"

        key = str(key)
        final_file_kwargs[key] = str(dirpaths[dirpath] / basename)
    return dirpaths, final_file_kwargs


class PrintLogHandler:
    def handle_stdout(self, log):
        print("stdout", log)

    def handle_stderr(self, log):
        print("stderr", log)


def _read_stdout(container, log_handler, output_buffer):
    for log in container.logs(stdout=True, stderr=False, stream=True):
        log_handler.handle_stdout(log)
        output_buffer.write(log)


def _read_stderr(container, log_handler):
    for log in container.logs(stdout=False, stderr=True, stream=True):
        log_handler.handle_stderr(log)


def run(
    container_uri,
    command="",
    *,
    file_kwargs,
    kwargs,
    output_dir=None,
    output_file_keys=None,
    log_handler=None,
):
    """
    kwargs will be passed to container as --key=value
    where value will be shell escaped
    file_kwargs must be a dict of key to Python file-like object
    the underlying files will be made available to the container as bound mount points (readonly).
    Note that we must ensure the directories thus made available to the container only have challenge-wide data, no outputs from other submissions!
    output_dir, if not None, will be mounted r/w for container to write outputs
    output_file_keys is a set of keys for file output types. Their values will be mapped from paths on the guest to paths on the host.
    command is optional
    Iterate over key/values of results.
    """
    if log_handler is None:
        log_handler = PrintLogHandler()

    if not (
        hasattr(log_handler, "handle_stdout") and hasattr(log_handler, "handle_stderr")
    ):
        raise TypeError("log_handler must have handle_stdout and handle_stderr methods")

    input_dir_map, final_file_kwargs = _convert_file_kwargs(file_kwargs)
    final_kwargs = copy.deepcopy(kwargs)
    final_kwargs.update(final_file_kwargs)

    final_command = _prepare_commandline(command, final_kwargs)

    running_container = run_container(
        container_uri, final_command, input_dir_map, output_dir=output_dir
    )

    output_buffer = io.BytesIO()
    out_thread = threading.Thread(
        target=_read_stdout,
        name="stdout",
        args=(running_container, log_handler, output_buffer),
    )
    err_thread = threading.Thread(
        target=_read_stderr, name="stderr", args=(running_container, log_handler)
    )
    err_thread.start()
    out_thread.start()

    # TODO: timeouts here?
    err_thread.join()
    out_thread.join()
    result = output_buffer.getvalue()
    yield from _parse_output(output_dir, result, output_file_keys).items()

    running_container.reload()
    if running_container.status != "exited":
        print(f"Container status is {running_container.status}", file=sys.stderr)
    running_container.remove()


def run_container(container_uri, command, inputdir_map=None, output_dir=None):
    client = docker.from_env()
    volumes = {}
    for inputdir, guest_input_dir in inputdir_map.items():
        volumes[str(inputdir)] = {"bind": str(guest_input_dir), "mode": "ro"}
    if output_dir:
        output_dir = Path(output_dir).resolve()
        volumes[str(output_dir)] = {"bind": str(GUEST_OUTPUT_DIR), "mode": "rw"}
        command = f" --output-dir {GUEST_OUTPUT_DIR} {command}"

    running_container = client.containers.run(
        container_uri,
        command,
        volumes=volumes,
        network_disabled=True,
        network_mode="none",
        remove=False,
        detach=True,
    )
    return running_container
