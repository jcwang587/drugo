from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdDepictor
from PIL import Image
import io
import cairosvg
import base64

def draw_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.RemoveHs(mol)

    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    Chem.RemoveStereochemistry(mol)
    mol_draw = Draw.rdMolDraw2D.PrepareMolForDrawing(mol, addChiralHs=False)

    num_atoms = mol.GetNumHeavyAtoms()
    size = min(1000, num_atoms * 30)
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
    png_data = cairosvg.svg2png(bytestring=svg.encode('utf-8'))
    img = Image.open(io.BytesIO(png_data)).convert("RGBA")

    # Crop the image to the bounding box
    bbox = img.getbbox()
    img_cropped = img.crop(bbox)

    # Add a white background
    background = Image.new("RGBA", img_cropped.size, (255, 255, 255, 255))
    background.paste(img_cropped, (0, 0), img_cropped)

    # Add a white margin
    margin_size = 25  
    new_size = (background.width + 2 * margin_size, background.height + 2 * margin_size)
    img_margin = Image.new("RGBA", new_size, (255, 255, 255, 255))
    img_margin.paste(background, (margin_size, margin_size))

    # Convert the image to a bytes object
    img_bytes = io.BytesIO()
    img_margin.save(img_bytes, format='PNG')
    img_bytes.seek(0)

    png_content = base64.b64encode(img_bytes.read()).decode("utf-8")

    return png_content
