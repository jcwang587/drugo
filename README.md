# drugo 💊🗄️

[![Release](https://img.shields.io/github/v/release/jcwang587/drugo)](https://github.com/jcwang587/drugo/releases/latest)
[![Heroku Status](https://img.shields.io/badge/Heroku-5A1BA9?logo=heroku)](https://drugo-a54338d8b0d8.herokuapp.com/)
[![Powered by SQLAlchemy](https://img.shields.io/badge/powered%20by-SQLAlchemy-CB2222.svg)](https://github.com/sqlalchemy/sqlalchemy)
[![Powered by Zotero](https://img.shields.io/badge/powered%20by-Zotero-CD2533.svg)](https://github.com/zotero/zotero)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Drugo is a database primarily composed of drug molecules, designed to serve as a benchmark dataset for predicting sites of metabolism. The current goal of this project is to provide a database containing molecular structural information using SMILES strings and structurally assigned literature references. Currently, the priority for the construction of the database is focusing on the CYP450 3A4 substrate database.

## References

If you are using the Drugo database for academic work, please cite the following original papers:

```bibtex

@article{JChemInfModel2013,
  author = {Zaretzki, Jed and Matlock, Matthew and Swamidass, S. Joshua},
  title = {XenoSite: Accurately Predicting CYP-Mediated Sites of Metabolism with Neural Networks},
  journal = {Journal of Chemical Information and Modeling},
  volume = {53},
  number = {12},
  pages = {3373-3383},
  year = {2013},
  doi = {10.1021/ci400518g}
}
```
