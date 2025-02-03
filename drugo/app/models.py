from .server import db


class Drugs(db.Model):
    __tablename__ = "drugs"
    drug_id = db.Column(db.String, primary_key=True)
    drug_title = db.Column(db.String)
    smiles = db.Column(db.String)


class Molecules(db.Model):
    __tablename__ = "molecules"
    drug_id = db.Column(db.String, primary_key=True)
    som = db.Column(db.Integer)
    som_element = db.Column(db.String)
    som_level = db.Column(db.Integer)


class References(db.Model):
    __tablename__ = "references"
    drug_id = db.Column(db.String, primary_key=True)
    reference = db.Column(db.String)
    doi = db.Column(db.String)
