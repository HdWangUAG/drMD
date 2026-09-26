"""
Checks that drGeometry refuses a topology PDB whose elements would be guessed wrong.

A tleap-written `*_solvated.pdb` leaves PDB columns 73-80 empty AND writes every atom name
left-justified from column 13. The PDB convention is that a two-letter element symbol starts in
column 13 and a one-letter symbol in column 14, so with both signals gone there is nothing to tell
`CA` the alpha carbon from `CA` calcium, and mdtraj takes the two-letter reading. Measured on a
162,983-atom UcFatB1 build (HFA_spec/WT_C12), the tleap PDB and the Amber prmtop of the *same* atoms
in the *same* order disagree on 2,266 elements, of which 820 are hydrogens promoted to heavy atoms:

    HG/HG1/HG2/HG3 -> mercury (200.59 amu)   546 atoms
    HE/HE1/HE2/HE3 -> helium  (4.00 amu)     274 atoms
    CA -> calcium, CD -> cadmium, CE -> cerium, NE -> neon, ND -> neodymium, SG -> seaborgium

Nothing in drGeometry can notice: the mis-elemented hydrogens pass the `element.symbol == "H"` test
in select_atoms and heavy_atom_indices, so they enter `sidechain` and `heavy` selections and drag
every minDistance short (up to ~0.9 A on an Arg-phosphate anchor, with fracIntact inflated), while
the wrong masses capture every mass-weighted comDistance centroid (the Val144/His140/Ala243 floor
group weighs 497 amu instead of 124). No warning, no traceback, just wrong published numbers.

So such a topology must be an error, an element-complete topology (prmtop/parm7/PSF/mmCIF) must be
loadable directly, and a conventionally written PDB whose blank element column happens to be
unambiguous must still be accepted - a check that cried wolf would just be deleted.

No MD and no real trajectory: the fixtures are a handful of hand-written PDB lines.

Run with:  python tests/test_geometry_topology_elements.py   (or pytest tests/)
"""
import sys
import tempfile
from os import path as p

ROOT = p.dirname(p.dirname(p.abspath(__file__)))
sys.path.insert(0, p.join(ROOT, "src"))
from ExaminationRoom import drGeometry

## One arginine, the residue that made the defect visible: HG2/HG3 become mercury and HE becomes
## helium when the element column is blank, so a `sidechain` selection gains three fake heavy atoms.
ARG_ATOMS = [
    ("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C"), ("HB2", "H"), ("HB3", "H"),
    ("CG", "C"), ("HG2", "H"), ("HG3", "H"), ("CD", "C"), ("HD2", "H"), ("HD3", "H"),
    ("NE", "N"), ("HE", "H"), ("CZ", "C"), ("NH1", "N"), ("NH2", "N"),
]


def write_pdb(path, withElements, leftJustified=True):
    """A one-residue PDB, built column by column, with or without columns 77-78 filled.

    tleap writes the atom name left-justified from column 13, which is what the real
    `*_solvated.pdb` looks like, so that is the default. `leftJustified=False` writes the
    conventional layout instead, where a one-letter element symbol is indented to column 14 - the
    layout that makes the element guessable even with the element field blank.
    """
    with open(path, "w") as fh:
        for serial, (name, element) in enumerate(ARG_ATOMS, start=1):
            x = 0.9 * serial
            field = f"{name:<4s}" if leftJustified or len(name) == 4 else f" {name:<3s}"
            line = (f"ATOM  "                  # 1-6
                    f"{serial:5d}"             # 7-11
                    f" "                       # 12
                    f"{field}"                 # 13-16 atom name
                    f" "                       # 17 altLoc
                    f"ARG"                     # 18-20 resName
                    f" "                       # 21
                    f"A"                       # 22 chainId
                    f" 167"                    # 23-26 resSeq
                    f"    "                    # 27-30 iCode + pad
                    f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}"   # 31-54
                    f"{1.0:6.2f}{0.0:6.2f}"    # 55-66
                    f"{'':10s}")               # 67-76
            line += f"{element:>2s}" if withElements else "  "   # 77-78
            assert len(line) == 78, (len(line), line)
            fh.write(line + "\n")
        fh.write("TER\nEND\n")
    return path


def test_blank_element_column_is_counted():
    with tempfile.TemporaryDirectory() as tmp:
        blank = write_pdb(p.join(tmp, "blank.pdb"), withElements=False)
        filled = write_pdb(p.join(tmp, "filled.pdb"), withElements=True)
        nBlank, total = drGeometry.count_blank_element_records(blank)
        print(f"blank.pdb: {nBlank} of {total} records have no element")
        assert total == len(ARG_ATOMS) and nBlank == len(ARG_ATOMS)
        nBlank, total = drGeometry.count_blank_element_records(filled)
        print(f"filled.pdb: {nBlank} of {total} records have no element")
        assert total == len(ARG_ATOMS) and nBlank == 0


def test_ambiguous_element_names_are_refused():
    with tempfile.TemporaryDirectory() as tmp:
        blank = write_pdb(p.join(tmp, "blank.pdb"), withElements=False)
        try:
            drGeometry.load_topology(blank)
        except ValueError as error:
            print(f"blank element column raises: {str(error)[:140]}...")
            message = str(error)
            assert "element column" in message, message
            ## the message must name the atoms it would have got wrong, and the way out
            assert "HG2" in message and "HE" in message, message
            assert "prmtop" in message, message
        else:
            raise AssertionError("a topology PDB with a blank element column must not be accepted")


def test_blank_element_column_can_be_forced_but_warns():
    with tempfile.TemporaryDirectory() as tmp:
        blank = write_pdb(p.join(tmp, "blank.pdb"), withElements=False)
        topology, source = drGeometry.load_topology(blank, allowInferredElements=True)
        print(f"forced: elements from {source!r}")
        assert "guessed" in source.lower(), source
        ## and this is what the guess costs on one arginine: three hydrogens become heavy metals,
        ## and three carbons/nitrogens get a mass 3-9x too large
        wrong = {atom.name: atom.element.symbol for atom in topology.atoms
                 if atom.element is not None and atom.element.symbol not in ("C", "N", "O", "H")}
        print(f"mis-elemented atoms under inference: {wrong}")
        assert wrong == {"CA": "Ca", "CD": "Cd", "NE": "Ne",
                         "HG2": "Hg", "HG3": "Hg", "HE": "He"}, wrong
        ## the three hydrogens-as-metals are the damaging ones: they are not element "H" any more
        for name in ("HG2", "HG3", "HE"):
            atom = next(a for a in topology.atoms if a.name == name)
            assert atom.element.symbol != "H", name
            assert atom.element.mass > 3.0, (name, atom.element.mass)


def test_filled_element_column_is_accepted_and_gives_right_elements():
    with tempfile.TemporaryDirectory() as tmp:
        filled = write_pdb(p.join(tmp, "filled.pdb"), withElements=True)
        topology, source = drGeometry.load_topology(filled)
        print(f"filled.pdb: elements from {source!r}, {topology.n_atoms} atoms")
        assert source == "PDB element column", source
        elements = {atom.name: atom.element.symbol for atom in topology.atoms}
        for name, element in ARG_ATOMS:
            assert elements[name] == element, (name, elements[name], element)


def test_hydrogens_are_excluded_from_sidechain_only_with_right_elements():
    """The bug in one assertion: a blank element column puts HG2/HG3/HE into a sidechain selection."""
    with tempfile.TemporaryDirectory() as tmp:
        blank = write_pdb(p.join(tmp, "blank.pdb"), withElements=False)
        filled = write_pdb(p.join(tmp, "filled.pdb"), withElements=True)
        selection = {"chain": "A", "resId": 167, "atoms": "sidechain"}

        good, _ = drGeometry.load_topology(filled)
        indices = drGeometry.select_atoms(good, drGeometry.residue_labels(good, "A:1-1:166"), selection)
        rightNames = sorted(good.atom(i).name for i in indices)
        print(f"element-complete sidechain: {rightNames}")
        assert rightNames == ["CB", "CD", "CG", "CZ", "NE", "NH1", "NH2"], rightNames

        bad, _ = drGeometry.load_topology(blank, allowInferredElements=True)
        indices = drGeometry.select_atoms(bad, drGeometry.residue_labels(bad, "A:1-1:166"), selection)
        wrongNames = sorted(bad.atom(i).name for i in indices)
        print(f"guessed-element sidechain:   {wrongNames}")
        assert set(wrongNames) - set(rightNames) == {"HE", "HG2", "HG3"}, wrongNames
        ## those three also survive the heavy-atom filter that minDistance relies on
        heavy = drGeometry.heavy_atom_indices(bad, indices, selection, "Arg167__phosphate")
        heavyNames = sorted(bad.atom(i).name for i in heavy)
        print(f"guessed-element heavy atoms: {heavyNames}")
        assert "HG2" in heavyNames and "HE" in heavyNames, heavyNames


def test_a_conventional_pdb_without_elements_is_still_accepted():
    """The check must not cry wolf: an indented one-letter name is unambiguous even with no element."""
    with tempfile.TemporaryDirectory() as tmp:
        conventional = write_pdb(p.join(tmp, "conventional.pdb"), withElements=False, leftJustified=False)
        nBlank, total = drGeometry.count_blank_element_records(conventional)
        assert nBlank == total == len(ARG_ATOMS)
        topology, source = drGeometry.load_topology(conventional)
        print(f"conventional blank-element PDB: elements from {source!r}")
        assert "unambiguous" in source, source
        assert drGeometry.ambiguous_element_atoms(topology) == {}
        wrong = {atom.name: atom.element.symbol for atom in topology.atoms
                 if atom.element is not None and atom.element.symbol not in ("C", "N", "O", "H")}
        print(f"mis-elemented atoms: {wrong}")
        assert wrong == {}, wrong


def test_unambiguous_two_letter_elements_are_not_flagged():
    """ZN -> zinc is a correct two-letter read; Z is no element, so it must not be called ambiguous."""
    with tempfile.TemporaryDirectory() as tmp:
        path = p.join(tmp, "zinc.pdb")
        with open(path, "w") as fh:
            fh.write("ATOM      1 ZN   ZN  A 999       0.000   0.000   0.000  1.00  0.00\n")
            fh.write("TER\nEND\n")
        topology, source = drGeometry.load_topology(path)
        print(f"zinc.pdb: elements from {source!r}, "
              f"{[(a.name, a.element.symbol) for a in topology.atoms]}")
        assert drGeometry.ambiguous_element_atoms(topology) == {}


def test_element_complete_suffixes_include_amber_topologies():
    for suffix in (".prmtop", ".parm7"):
        assert suffix in drGeometry.ELEMENT_COMPLETE_SUFFIXES, drGeometry.ELEMENT_COMPLETE_SUFFIXES
    print(f"element-complete suffixes: {drGeometry.ELEMENT_COMPLETE_SUFFIXES}")
    assert set("HCNOPS") <= drGeometry.AMBIGUOUS_LEADING_ELEMENTS, drGeometry.AMBIGUOUS_LEADING_ELEMENTS


if __name__ == "__main__":
    test_blank_element_column_is_counted()
    test_ambiguous_element_names_are_refused()
    test_blank_element_column_can_be_forced_but_warns()
    test_filled_element_column_is_accepted_and_gives_right_elements()
    test_hydrogens_are_excluded_from_sidechain_only_with_right_elements()
    test_a_conventional_pdb_without_elements_is_still_accepted()
    test_unambiguous_two_letter_elements_are_not_flagged()
    test_element_complete_suffixes_include_amber_topologies()
    print("ALL GEOMETRY TOPOLOGY ELEMENT TESTS PASSED")
