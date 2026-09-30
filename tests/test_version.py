import re
import tomllib
from pathlib import Path

import artharness

ROOT = Path(__file__).resolve().parent.parent


def test_versions_agree_and_stay_in_the_0_0_x_series():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
    assert project == artharness.__version__
    assert re.fullmatch(r"0\.0\.\d+", project), "agents release only 0.0.x"
    lock = (ROOT / "uv.lock").read_text()
    locked = re.search(r'name = "artharness"\nversion = "([^"]+)"', lock).group(1)
    assert locked == project
