from flask import Flask
from flask_sqlalchemy import SQLAlchemy
import os


class Config(object):
    # Construct the absolute path to the database file
    BASE_DIR = os.path.abspath(os.path.dirname(__file__))
    SQLALCHEMY_DATABASE_URI = (
        f"sqlite:///{os.path.join(BASE_DIR, '../database/v2025.1/drugo_3a4.db')}"
    )
    SQLALCHEMY_TRACK_MODIFICATIONS = False


db = SQLAlchemy()


def create_app():
    server = Flask(__name__)
    server.config.from_object(Config)
    db.init_app(server)

    from .dashboard import create_dashapp

    dash_app = create_dashapp(server)
    return server
