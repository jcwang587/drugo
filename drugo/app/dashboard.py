import sqlalchemy as sa
import plotly.graph_objects as go
import dash_bootstrap_components as dbc
from dash import dash_table

from .server import db
from .models import Drugs, Molecules, References
from dash import Dash, html, dcc, Input, Output


def create_dashapp(server):
    # Initialize the app with a Bootstrap theme
    app = Dash(__name__, server=server, external_stylesheets=[dbc.themes.BOOTSTRAP])    

    # App layout
    app.layout = dbc.Container([
        html.H1('Drug Oxidation Database', className='text-center my-3'), 
        dbc.Row([
            dbc.Col([
                html.H2('Tables', className='text-center'),
                html.Hr(),
                dcc.RadioItems(
                    options=[
                        {'label': 'Drugs', 'value': 'drugs'},
                        {'label': 'Molecules', 'value': 'molecules'},
                        {'label': 'References', 'value': 'references'}
                    ],
                    value='drugs',
                    id='controls',
                    inline=False
                )
            ], width=2, className='bg-light p-3'), 

            dbc.Col([
                html.Div(id='first-graph')
            ], width=10) 
        ])
    ], fluid=True)

    @app.callback(
        Output(component_id='first-graph', component_property='children'),
        Input(component_id='controls', component_property='value')
    )
    def update_table(table_chosen):
        # Map the chosen value to the actual model class
        table_map = {
            'drugs': Drugs,
            'molecules': Molecules,
            'references': References
        }
        
        # Get the model class based on the chosen table
        model_class = table_map.get(table_chosen)
        
        if model_class is None:
            return html.Div("No data available")  # Return a message if no valid model is found
        
        # Query the selected table
        query = db.session.query(model_class).all()
        
        # Convert query results to a list of dictionaries
        data = [{column.name: getattr(row, column.name) for column in model_class.__table__.columns} for row in query]
        
        # Define column widths for each table
        column_widths = {
            'drugs': {'name': '200px', 'type': '150px', 'description': '300px'},
            'molecules': {'name': '250px', 'formula': '200px', 'weight': '150px'},
            'references': {'title': '300px', 'author': '200px', 'year': '100px'}
        }

        # Get the column widths for the selected table
        widths = column_widths.get(table_chosen, {})

        # Create a Dash DataTable
        table = dash_table.DataTable(
            columns=[{"name": i, "id": i} for i in data[0].keys()],
            data=data,
            style_table={'overflowX': 'hidden'},  # Hide horizontal scrollbar
            style_cell={'textAlign': 'left'},
            style_data_conditional=[
                {
                    'if': {'column_id': column_id},
                    'minWidth': width,
                    'maxWidth': width,
                    'width': width
                } for column_id, width in widths.items()
            ],
            style_header={
                'backgroundColor': 'rgb(230, 230, 230)',
                'fontWeight': 'bold'
            },
            page_size=20  
        )
        
        return table

    return app