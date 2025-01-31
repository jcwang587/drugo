from flask import Flask
from flask_sqlalchemy import SQLAlchemy


class Config(object):
    SQLALCHEMY_DATABASE_URI = 'sqlite:///../database/v2025.1/drugo_3a4.db'
    SQLALCHEMY_TRACK_MODIFICATIONS = False

db = SQLAlchemy()


def create_app():
    server = Flask(__name__)
    server.config.from_object(Config)
    db.init_app(server)
    
    from .dashboard import create_dashapp
    dash_app = create_dashapp(server)
    return server 