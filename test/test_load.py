import multiseq
from pathlib import Path
from scanpy import AnnData


def test_loading_str():
    x = multiseq.load_citeseq_count_matrix("test-data")
    assert isinstance(x, AnnData)


def test_loading_path():
    x = multiseq.load_citeseq_count_matrix(Path("test-data"))
    assert isinstance(x, AnnData)
