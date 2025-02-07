from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdDepictor
from PIL import Image
import io
import cairosvg
import base64


def prepare_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.RemoveHs(mol)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    Chem.RemoveStereochemistry(mol)
    return Draw.rdMolDraw2D.PrepareMolForDrawing(mol, addChiralHs=False)


def draw_molecule_size(smiles, size=1000, rotate=False):
    mol_draw = prepare_molecule(smiles)

    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(size, size)

    opts = drawer.drawOptions()

    for i in range(mol_draw.GetNumAtoms()):
        opts.atomLabels[i] = f"{mol_draw.GetAtomWithIdx(i).GetSymbol()}<sub>{i+1}</sub>"

    if rotate:
        opts.rotate = 90

    drawer.drawOptions().prepareMolsBeforeDrawing = False
    drawer.drawOptions().maxFontSize = 30
    drawer.drawOptions().clearBackground = False

    drawer.DrawMolecule(mol_draw)
    drawer.FinishDrawing()

    svg = drawer.GetDrawingText()
    png_data = cairosvg.svg2png(bytestring=svg.encode("utf-8"))
    img = Image.open(io.BytesIO(png_data)).convert("RGBA")

    # Crop the image to the bounding box
    bbox = img.getbbox()
    img_cropped = img.crop(bbox)
    width, height = img_cropped.size

    # Add a white background
    background = Image.new("RGBA", (width, height), (255, 255, 255, 255))
    background.paste(img_cropped, (0, 0), img_cropped)

    # Add a white margin
    margin_size = 25
    new_size = (background.width + 2 * margin_size, background.height + 2 * margin_size)
    img_margin = Image.new("RGBA", new_size, (255, 255, 255, 255))
    img_margin.paste(background, (margin_size, margin_size))

    # Convert the image to a bytes object
    img_bytes = io.BytesIO()
    img_margin.save(img_bytes, format="PNG")
    img_bytes.seek(0)

    png_content = base64.b64encode(img_bytes.read()).decode("utf-8")

    return png_content, width, height


def draw_smiles(smiles):
    mol_draw = prepare_molecule(smiles)

    # Calculate the size based on the number of atoms
    num_atoms = mol_draw.GetNumHeavyAtoms()

    if num_atoms > 50:
        size = min(900, num_atoms * 25)
    elif num_atoms > 40:
        size = min(800, num_atoms * 25)
    elif num_atoms > 30:
        size = min(700, num_atoms * 25)
    elif num_atoms > 20:
        size = min(600, num_atoms * 25)
    else:
        size = min(500, num_atoms * 25)

    png_content, width, height = draw_molecule_size(smiles, size, rotate=False)

    if height > 1.1 * width:
        png_content, width, height = draw_molecule_size(smiles, size=size, rotate=True)

    if num_atoms > 50 and width > 1.3 * height:
        png_content, width, height = draw_molecule_size(smiles, 1000, rotate=False)
    elif num_atoms > 40 and width > 1.3 * height:
        png_content, width, height = draw_molecule_size(smiles, 900, rotate=False)
    elif num_atoms > 30 and width > 1.3 * height:
        png_content, width, height = draw_molecule_size(smiles, 800, rotate=False)
    elif num_atoms > 20 and width > 1.3 * height:
        png_content, width, height = draw_molecule_size(smiles, 700, rotate=False)
    else:
        pass

    return png_content
