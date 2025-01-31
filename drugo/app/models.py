from .server import db


class Drugs(db.Model):
    __tablename__ = "drugs"
    id = db.Column(db.Integer, primary_key=True)
    drug_id = db.Column(db.String)
    drug_title = db.Column(db.String)
    smiles = db.Column(db.String)

class Molecules(db.Model):
    __tablename__ = "molecules"
    id = db.Column(db.Integer, primary_key=True)
    som = db.Column(db.Integer)
    som_level = db.Column(db.Integer)
    som_element = db.Column(db.String)
    drug_id = db.Column(db.String)

class References(db.Model):
    __tablename__ = "references"
    id = db.Column(db.Integer, primary_key=True)
    reference = db.Column(db.String)
    doi = db.Column(db.String)
    drug_id = db.Column(db.String)