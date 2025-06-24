from cif_site_analyzer import cif_reader
import os


def test_cif_readers():
    CIF_DIR = os.path.join(os.path.dirname(__file__), "NaCl,cF8,225")
    for cif in os.listdir(CIF_DIR):
        if not cif.endswith("cif"):
            continue
        assert cif_reader.read_cif(
            f"{CIF_DIR}{os.sep}{cif}"
        ), f"Error reading {cif}"
