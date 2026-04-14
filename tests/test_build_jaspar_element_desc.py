"""Tests for scripts/build_jaspar_element_desc.py."""

from __future__ import annotations

from pathlib import Path

import pytest

from scripts.build_jaspar_element_desc import (
    DEFAULT_BUCKET,
    JasparMotif,
    build_argument_parser,
    build_description_rows,
    classify_family,
    main,
    parse_meme_motifs,
    write_xlsx,
)

pytestmark = pytest.mark.unit


MEME_SAMPLE = """\
MEME version 5.5.5

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF MA0020.1 Dof2
letter-probability matrix: alength= 4 w= 7 nsites= 20 E= 0
 0.25 0.25 0.25 0.25

MOTIF MA0044.1 HAHB4
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
 0.25 0.25 0.25 0.25

MOTIF MA0082.1 SQUA
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
 0.25 0.25 0.25 0.25

MOTIF MA0104.1
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
 0.25 0.25 0.25 0.25
"""


@pytest.fixture
def meme_bundle(tmp_path: Path) -> Path:
    path = tmp_path / "JASPAR_mini.meme"
    path.write_text(MEME_SAMPLE)
    return path


def test_parse_meme_motifs_extracts_motif_and_alt_ids(meme_bundle: Path) -> None:
    motifs = parse_meme_motifs(meme_bundle)

    assert motifs == [
        JasparMotif(motif_id="MA0020.1", alt_id="Dof2"),
        JasparMotif(motif_id="MA0044.1", alt_id="HAHB4"),
        JasparMotif(motif_id="MA0082.1", alt_id="SQUA"),
        JasparMotif(motif_id="MA0104.1", alt_id=""),
    ]


def test_parse_meme_motifs_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="MEME bundle not found"):
        parse_meme_motifs(tmp_path / "does_not_exist.meme")


def test_parse_meme_motifs_not_a_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="not a file"):
        parse_meme_motifs(tmp_path)


def test_display_name_falls_back_to_motif_id() -> None:
    assert JasparMotif(motif_id="MA0104.1", alt_id="").display_name == "MA0104.1"
    assert JasparMotif(motif_id="MA0020.1", alt_id="Dof2").display_name == "Dof2"


@pytest.mark.parametrize(
    "name,expected_fragment",
    [
        ("MYB15", "MYB"),
        ("myb3r1", "MYB"),
        ("WRKY40", "WRKY"),
        ("NAC002", "NAC"),
        ("ATAF1", "NAC"),
        ("ERF1", "AP2/ERF"),
        ("DREB2A", "AP2/ERF"),
        ("RAP2.3", "AP2/ERF"),
        ("bHLH112", "bHLH"),
        ("PIF3", "bHLH"),
        ("bZIP11", "bZIP"),
        ("ABF2", "bZIP"),
        ("HY5", "bZIP"),
        ("TGA1", "bZIP"),
        ("Dof3", "DOF"),
        ("GATA1", "GATA"),
        ("BBX19", "GATA"),
        ("TCP4", "TCP"),
        ("SPL9", "SBP"),
        ("HSF1", "HSF"),
        ("ARF5", "ARF"),
        ("ARR1", "ARR"),
        ("AGL3", "MADS"),
        ("SEP3", "MADS"),
        ("AP1", "MADS"),
        ("ATHB2", "Homeobox"),
        ("WOX5", "Homeobox"),
        ("LHY", "circadian"),
        ("CCA1", "circadian"),
        ("IDD5", "C2H2"),
        ("ZAT6", "C2H2"),
        ("LFY", "LFY"),
        ("GRF1", "GRF"),
        ("PHL1", "GARP"),
        ("GLK1", "GARP"),
        ("CAMTA3", "CAMTA"),
        ("PBF", "PBF"),
        ("AHL20", "AT-hook"),
        ("BPC1", "BBR/BPC"),
        ("E2F1", "E2F"),
        ("FAR1", "FAR1"),
        ("AT1G12345", "Species-specific"),
        ("Zm00001d012345", "Species-specific"),
    ],
)
def test_classify_family_known_patterns(name: str, expected_fragment: str) -> None:
    assert expected_fragment in classify_family(name)


def test_classify_family_default_bucket() -> None:
    assert classify_family("QQQ123") == DEFAULT_BUCKET
    assert classify_family("") == DEFAULT_BUCKET
    assert classify_family("   ") == DEFAULT_BUCKET


def test_build_description_rows_deduplicates_and_sorts() -> None:
    motifs = [
        JasparMotif(motif_id="M1", alt_id="WRKY40"),
        JasparMotif(motif_id="M2", alt_id="MYB15"),
        JasparMotif(motif_id="M3", alt_id="WRKY40"),  # duplicate
        JasparMotif(motif_id="M4", alt_id=""),  # falls back to M4 (classified via motif_id)
    ]

    rows = build_description_rows(motifs)

    elements = [element for element, _ in rows]
    assert elements == sorted(elements)
    assert elements.count("WRKY40") == 1
    assert "MYB15" in elements
    assert "M4" in elements

    row_dict = dict(rows)
    assert "MYB" in row_dict["MYB15"]
    assert "WRKY" in row_dict["WRKY40"]


def test_build_description_rows_skips_entirely_empty_names() -> None:
    motifs = [JasparMotif(motif_id="", alt_id="")]
    assert build_description_rows(motifs) == []


def test_write_xlsx_round_trip(tmp_path: Path) -> None:
    openpyxl = pytest.importorskip("openpyxl")

    rows = [
        ("Dof3", "DOF zinc finger family"),
        ("MYB15", "MYB transcription factor family"),
    ]
    output_path = tmp_path / "nested" / "jaspar_element.desc.xlsx"

    returned = write_xlsx(rows, output_path)

    assert returned == output_path
    assert output_path.exists()

    workbook = openpyxl.load_workbook(output_path)
    sheet = workbook.active
    assert sheet.title == "jaspar_elements"
    data = list(sheet.iter_rows(values_only=True))
    assert data[0] == ("element", "description")
    assert data[1] == ("Dof3", "DOF zinc finger family")
    assert data[2] == ("MYB15", "MYB transcription factor family")


def test_build_argument_parser_requires_meme_and_output() -> None:
    parser = build_argument_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([])


def test_main_end_to_end(tmp_path: Path, meme_bundle: Path, capsys: pytest.CaptureFixture[str]) -> None:
    pytest.importorskip("openpyxl")
    output_path = tmp_path / "jaspar_element.desc.xlsx"

    exit_code = main(["--meme", str(meme_bundle), "--output", str(output_path)])

    assert exit_code == 0
    assert output_path.exists()
    captured = capsys.readouterr()
    assert "element descriptions" in captured.out


def test_main_missing_meme(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    exit_code = main(
        [
            "--meme",
            str(tmp_path / "does_not_exist.meme"),
            "--output",
            str(tmp_path / "out.xlsx"),
        ]
    )
    assert exit_code == 1
    assert "Error" in capsys.readouterr().err


def test_main_empty_meme(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    empty = tmp_path / "empty.meme"
    empty.write_text("MEME version 5.5.5\n\nALPHABET= ACGT\n")

    exit_code = main(
        [
            "--meme",
            str(empty),
            "--output",
            str(tmp_path / "out.xlsx"),
        ]
    )
    assert exit_code == 1
    assert "no MOTIF entries" in capsys.readouterr().err
