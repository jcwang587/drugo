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

def draw_molecule_size(smiles, size=1000, rotate=False, highlight_atoms=None):
    """
    Draws the molecule as an SVG, converts it to PNG, and returns the PNG
    content (base64 encoded) along with the image width and height.
    
    Parameters:
      smiles (str): The SMILES string for the molecule.
      size (int): The size (width and height in pixels) for the drawing.
      rotate (bool): Whether to rotate the molecule 90 degrees.
      highlight_atoms (list or None): A list of atom indices (0-indexed) to highlight.
    """
    mol_draw = prepare_molecule(smiles)
    highlight_atoms = [atom - 1 for atom in highlight_atoms]

    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(size, size)
    opts = drawer.drawOptions()

    # Label atoms with their element symbol and 1-indexed number
    for i in range(mol_draw.GetNumAtoms()):
        opts.atomLabels[i] = f"{mol_draw.GetAtomWithIdx(i).GetSymbol()}<sub>{i+1}</sub>"

    if rotate:
        opts.rotate = 90

    opts.prepareMolsBeforeDrawing = False
    opts.maxFontSize = 30
    opts.clearBackground = False

    if highlight_atoms is not None:
        light_blue = (173/255, 216/255, 230/255)
        highlight_atom_colors = {atom: light_blue for atom in highlight_atoms}
    else:
        highlight_atom_colors = {}

    # Draw the molecule, passing the highlight information if any
    drawer.DrawMolecule(mol_draw,
                          highlightAtoms=highlight_atoms,
                          highlightAtomColors=highlight_atom_colors,
                          highlightBonds=None)
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

    # Convert the image to a bytes object and then base64 encode it
    img_bytes = io.BytesIO()
    img_margin.save(img_bytes, format="PNG")
    img_bytes.seek(0)
    png_content = base64.b64encode(img_bytes.read()).decode("utf-8")

    return png_content, width, height

def draw_smiles(smiles, som_list):
    """
    Draws a molecule from a SMILES string, automatically chooses the size,
    and highlights atoms given in som_list (assumed to be 0-indexed).
    
    Parameters:
      smiles (str): The SMILES string.
      som_list (list): A list of atom indices to highlight.
    """
    mol_draw = prepare_molecule(smiles)
    num_atoms = mol_draw.GetNumHeavyAtoms()

    if num_atoms > 50:
        size = min(900, num_atoms * 25)
    elif num_atoms > 40:
        size = min(800, num_atoms * 25)
    elif num_atoms > 30:
        size = min(700, num_atoms * 25)
    elif num_atoms > 20:
        size = min(600, num_atoms * 25)
    elif num_atoms > 10:
        size = min(500, num_atoms * 30)
    else:
        size = min(400, num_atoms * 30)
    size = max(size, 300)

    # Initial drawing with the given size and highlights
    png_content, width, height = draw_molecule_size(smiles, size, rotate=False, highlight_atoms=som_list)

    # Rotate if the height is substantially larger than the width
    if height > width:
        rotate_flag = True
        png_content, width, height = draw_molecule_size(smiles, size, rotate=rotate_flag, highlight_atoms=som_list)
    else:
        rotate_flag = False

    # Adjust the size if the molecule is wide (and many atoms are present)
    if num_atoms > 50 and width > 1.3 * height:
        png_content, width, height = draw_molecule_size(smiles, 1000, rotate=rotate_flag, highlight_atoms=som_list)
    elif num_atoms > 40 and width > 1.3 * height:
        png_content, width, height = draw_molecule_size(smiles, 900, rotate=rotate_flag, highlight_atoms=som_list)
    elif num_atoms > 30 and width > 1.3 * height:
        png_content, width, height = draw_molecule_size(smiles, 800, rotate=rotate_flag, highlight_atoms=som_list)
    elif num_atoms > 20 and width > 1.3 * height:
        png_content, width, height = draw_molecule_size(smiles, 700, rotate=rotate_flag, highlight_atoms=som_list)
    elif num_atoms > 10 and width > 1.3 * height:
        png_content, width, height = draw_molecule_size(smiles, 600, rotate=rotate_flag, highlight_atoms=som_list)

    return png_content
