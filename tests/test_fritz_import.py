"""Offline regression test for https://github.com/MickaelRigault/ztfquery/issues/47

`from ztfquery import fritz` must not crash when the parent directory of the
local data dir (default ./Data) does not exist yet: it should create the
missing directory, including parents.
"""
import os
import subprocess
import sys

import pytest


def test_import_fritz_without_data_dir(tmp_path):
    code = "import ztfquery.fritz; print('IMPORT_OK')"
    env = {**os.environ, "ZTFDATA": "./Data"}
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"import failed:\n{result.stderr}"
    assert "IMPORT_OK" in result.stdout
    assert (tmp_path / "Data" / "fritz").is_dir(), "data dir was not created"
