"""Guards on where the package and its data are resolved from.

Both checks cover failure modes that are silent rather than loud: a stale
editable install elsewhere on the machine shadows these sources without any
error, and a misresolved ``DATA_DIR`` makes the ``vlab_utils`` loaders fall
back to synthetic data for a notebook the reader believes is showing a
measurement.
"""
from pathlib import Path

import vlab_xray
from vlab_xray import vlab_utils as V

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_package_is_imported_from_this_repository():
    assert Path(vlab_xray.__file__).resolve().parents[2] == REPO_ROOT


def test_data_dir_is_the_repository_data_directory():
    assert Path(vlab_xray.DATA_DIR).resolve() == REPO_ROOT / "data"
    assert Path(vlab_xray.DATA_DIR).is_dir()


def test_vlab_utils_shares_the_package_data_dir():
    assert Path(V.DATA_DIR).resolve() == Path(vlab_xray.DATA_DIR).resolve()


def test_loader_prefers_a_real_file_over_the_generator():
    def _generator():
        raise AssertionError("generator called although a real file exists")

    x, y, path = V.load_1d("cxro_iron_00595u", gen=_generator)
    assert path is not None and Path(path).name == "cxro_iron_00595u.txt"
    assert len(x) == len(y) > 0


def test_loader_falls_back_to_the_generator_when_absent():
    x, y, path = V.load_1d("no_such_export", gen=lambda: ([0.0, 1.0], [2.0, 3.0]))
    assert path is None
    assert list(x) == [0.0, 1.0]
