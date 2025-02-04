from .server import db


class Drugs(db.Model):
    __tablename__ = "drugs"
    drug_id = db.Column(db.String, primary_key=True)
    drug_title = db.Column(db.String)
    smiles = db.Column(db.String)

    # Relationships
    molecules = db.relationship("Molecules", backref="drug", lazy=True)
    references = db.relationship("References", backref="drug", lazy=True)


class Molecules(db.Model):
    __tablename__ = "molecules"
    # Add a new primary key column
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    # drug_id is now a foreign key, and not the primary key
    drug_id = db.Column(db.String, db.ForeignKey("drugs.drug_id"), nullable=False)
    som = db.Column(db.Integer)
    som_element = db.Column(db.String)
    som_level = db.Column(db.Integer)


class References(db.Model):
    __tablename__ = "references"
    # Add a new primary key column
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    # drug_id is now a foreign key, and not the primary key
    drug_id = db.Column(db.String, db.ForeignKey("drugs.drug_id"), nullable=False)
    reference = db.Column(db.String)
    doi = db.Column(db.String)
