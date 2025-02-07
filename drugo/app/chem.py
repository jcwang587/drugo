from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdDepictor


def draw_molecule(mol: Chem.Mol, name="labeled"):
    mol = Chem.RemoveHs(mol)

    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)

    Chem.RemoveStereochemistry(mol)

    num_atoms = mol.GetNumHeavyAtoms()
    size = min(1000, num_atoms * 30)

    mol_draw = Draw.rdMolDraw2D.PrepareMolForDrawing(mol, addChiralHs=False)
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(size, size)

    opts = drawer.drawOptions()

    for i in range(mol_draw.GetNumAtoms()):
        opts.atomLabels[i] = f"{mol_draw.GetAtomWithIdx(i).GetSymbol()}<sub>{i+1}</sub>"

    drawer.drawOptions().prepareMolsBeforeDrawing = False
    drawer.drawOptions().maxFontSize = 30

    drawer.drawOptions().clearBackground = False
    drawer.DrawMolecule(mol_draw)

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    with open(f"{name}.svg", "w") as f:
        f.write(svg)

    return svg
