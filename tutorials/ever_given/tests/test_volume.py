import os.path
import tempfile

import pytest

import ever_given.wrapper


def test_run_inputfile_only():
    test_mdlfile_rel = "data/ChEBI_16716.mdl"
    test_mdlfile_abs = os.path.join(os.path.dirname(__file__),
                                    test_mdlfile_rel)
    container_uri = "ghcr.io/robbason/calc-molwt:latest"
    file_kwargs = {"molfile": test_mdlfile_abs}
    results = {key: value
               for key, value in ever_given.wrapper.run(container_uri, kwargs={}, file_kwargs=file_kwargs)}

    assert set(results.keys()) == {"numAtoms", "numBonds", "molWeight"}
    molWeight = float(results["molWeight"].strip())
    assert pytest.approx(molWeight, 78.046950192)
