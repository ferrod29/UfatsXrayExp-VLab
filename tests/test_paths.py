"""Guards on where the package and its data are resolved from.

Both checks cover failure modes that are silent rather than loud: a stale
editable install elsewhere on the machine shadows these sources without any
error, and a misresolved ``DATA_DIR`` makes the ``vlab_utils`` loaders fall
back to synthetic data for a notebook the reader believes is showing a
measurement.
"""
from pathlib import Path

import numpy as np
import pytest

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


def test_processed_dir_sits_under_the_data_dir():
    processed = Path(V.PROCESSED_DIR)
    root = Path(vlab_xray.DATA_DIR) / "processed"
    # Either data/processed/ itself, or one dated session inside it.
    assert processed == root or processed.parent == root


def test_sessions_sort_by_date_even_when_a_folder_name_is_yyyyddmm():
    # 20263108 is 31 August 2026 written back to front; plain string sorting
    # would place it after September and make it the default session.
    order = sorted(["20260901", "20263108", "20260226"], key=V._session_date)
    assert order == ["20260226", "20263108", "20260901"]


def test_use_session_rejects_a_session_that_was_never_built():
    with pytest.raises(FileNotFoundError, match="no built session"):
        V.use_session("19000101")


def _write_processed(tmp_path, monkeypatch, **arrays):
    """Point the loaders at a throwaway processed/ directory."""
    processed = tmp_path / "processed"
    processed.mkdir()
    for name, value in arrays.items():
        np.save(processed / f"{name}.npy", value)
    monkeypatch.setattr(V, "DATA_DIR", str(tmp_path))
    monkeypatch.setattr(V, "PROCESSED_DIR", str(processed))
    return processed


def test_loaders_read_npy_arrays_from_the_processed_directory(tmp_path, monkeypatch):
    # scripts/build_processed_data.py writes .npy, which the text-only loaders
    # used to skip straight past on their way to the synthetic fallback.
    _write_processed(
        tmp_path, monkeypatch,
        trace=np.column_stack([np.arange(4.0), np.arange(4.0) ** 2]),
        frame=np.ones((3, 3)),
        table=np.arange(8.0).reshape(4, 2),
    )

    def _fail():
        raise AssertionError("generator called although a processed array exists")

    x, y, path = V.load_1d("trace", gen=_fail)
    assert path is not None and path.endswith("trace.npy")
    assert list(y) == [0.0, 1.0, 4.0, 9.0]

    image, path = V.load_image("frame", gen=_fail)
    assert path is not None and image.shape == (3, 3)

    table, path = V.load_table("table", gen=_fail)
    assert path is not None and table.shape == (4, 2)


def test_load_delay_series_reads_the_processed_npy_bundle(tmp_path, monkeypatch):
    _write_processed(
        tmp_path, monkeypatch,
        transients_energy=np.array([7040.0, 7050.0, 7060.0]),
        # deliberately out of order: the loader must sort by delay
        transients_delays=np.array([100.0, -50.0, 0.0]),
        transients_matrix=np.array([[3.0, 3.0, 3.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]),
    )
    energy, delays, matrix, path = V.load_delay_series("transients")
    assert path is not None and path.endswith("transients_matrix.npy")
    assert list(delays) == [-50.0, 0.0, 100.0]
    assert [row[0] for row in matrix] == [1.0, 2.0, 3.0]
    assert len(energy) == matrix.shape[1]


def test_loaders_still_fall_back_when_the_processed_directory_is_absent(tmp_path, monkeypatch):
    monkeypatch.setattr(V, "DATA_DIR", str(tmp_path))
    monkeypatch.setattr(V, "PROCESSED_DIR", str(tmp_path / "processed"))
    x, y, path = V.load_1d("nothing_here", gen=lambda: ([0.0], [1.0]))
    assert path is None and list(y) == [1.0]
    with pytest.raises(FileNotFoundError):
        V.load_table("nothing_here")
