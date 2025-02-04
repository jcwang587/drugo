import dash_bootstrap_components as dbc
from dash import dash_table
from dash import Dash, html, dcc, Input, Output, callback_context

from .server import db
from .models import Drugs, Molecules, References


def create_dashapp(server):
    # Initialize the app with a Bootstrap theme
    app = Dash(
        __name__,
        server=server,
        external_stylesheets=[dbc.themes.BOOTSTRAP],
        title="Drug Oxidation Database",
    )

    with app.server.app_context():
        drug_options = [
            {
                "label": (
                    (drug.drug_title[:25] + "...")
                    if len(drug.drug_title) > 25
                    else drug.drug_title
                ),
                "value": drug.drug_id,
                "title": drug.drug_title,
            }
            for drug in db.session.query(Drugs).all()
        ]

    # App layout
    app.layout = dbc.Container(
        [
            html.H1(
                "Drug Oxidation Database",
                className="text-center my-3",
                style={
                    "fontSize": "24px",
                    "fontWeight": "bold",
                    "color": "#343a40",
                    "fontFamily": "Roboto, sans-serif",
                },
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            html.H2(
                                "Tables",
                                className="text-center",
                                style={
                                    "fontSize": "20px",
                                    "fontWeight": "bold",
                                    "color": "#343a40",
                                    "fontFamily": "Roboto, sans-serif",
                                },
                            ),
                            html.Hr(),
                            dcc.RadioItems(
                                options=[
                                    {"label": " Drugs", "value": "drugs"},
                                    {"label": " Molecules", "value": "molecules"},
                                    {"label": " References", "value": "references"},
                                ],
                                value="drugs",
                                id="controls",
                                inline=False,
                                style={"fontFamily": "Roboto, sans-serif"},
                            ),
                            html.Br(),
                            dbc.Row(
                                [
                                    dbc.Col(
                                        html.Label(
                                            "Items per page",
                                            style={"fontFamily": "Roboto, sans-serif"},
                                        ),
                                        width="auto",
                                    ),
                                    dbc.Col(
                                        dcc.Dropdown(
                                            id="items-per-page-dropdown",
                                            options=[
                                                {"label": str(i), "value": i}
                                                for i in [10, 15, 25, 50, 100]
                                            ],
                                            value=15,  # Default value
                                            clearable=False,
                                            style={
                                                "fontFamily": "Roboto, sans-serif",
                                                "fontSize": "14px",
                                            },
                                        ),
                                        width=True,
                                    ),
                                ],
                                align="center",
                            ),
                            html.Br(),
                            html.Br(),
                            html.H2(
                                "Molecular View",
                                className="text-center",
                                style={
                                    "fontSize": "20px",
                                    "fontWeight": "bold",
                                    "color": "#343a40",
                                    "fontFamily": "Roboto, sans-serif",
                                },
                            ),
                            html.Hr(),
                            dcc.Dropdown(
                                id="dropdown",
                                options=drug_options,
                                value=None,
                                style={
                                    "fontFamily": "Roboto, sans-serif",
                                    "fontSize": "14px",
                                },
                            ),
                            html.Br(),
                        ],
                        width=2,
                        className="bg-light p-3",
                    ),
                    dbc.Col([html.Div(id="first-graph")], width=10),
                ]
            ),
        ],
        fluid=True,
    )

    # This callback returns two outputs:
    # 1. The table (or an empty Div if the molecular view is selected)
    # 2. The value of the molecular view dropdown
    # When the user clicks on the radio buttons, we clear the molecular view selection,
    # which causes the table to be shown again.
    @app.callback(
        Output("first-graph", "children"),
        Output("dropdown", "value"),
        Input("controls", "value"),
        Input("items-per-page-dropdown", "value"),
        Input("dropdown", "value"),
    )
    def update_table(table_chosen, items_per_page, selected_molecule):
        ctx = callback_context
        triggered_id = (
            ctx.triggered[0]["prop_id"].split(".")[0] if ctx.triggered else None
        )

        # When switching table type, clear the molecular view selection
        if triggered_id == "controls":
            selected_molecule = None

        # If a molecule has been chosen, hide the table
        if selected_molecule is not None:
            return html.Div(), selected_molecule

        table_map = {"drugs": Drugs, "molecules": Molecules, "references": References}
        model_class = table_map.get(table_chosen)
        if model_class is None:
            return html.Div("No data available"), selected_molecule

        query = db.session.query(model_class).all()
        data = [
            {
                column.name: getattr(row, column.name)
                for column in model_class.__table__.columns
            }
            for row in query
        ]

        # Define column widths for each table
        column_widths = {
            "drugs": {"drug_id": "150px", "drug_title": "400px", "smiles": "1600px"},
            "molecules": {
                "drug_id": "150px",
                "drug_title": "150px",
                "som": "150px",
                "som_element": "150px",
                "som_level": "150px",
            },
            "references": {"drug_id": "150px", "reference": "800px", "doi": "1200px"},
        }
        widths = column_widths.get(table_chosen, {})

        # Redefine column names
        column_names = {
            "drug_id": "drug id",
            "drug_title": "drug title",
            "smiles": "smiles",
            "som": "som",
            "som_level": "som level",
            "som_element": "som element",
            "reference": "reference",
            "doi": "doi",
        }

        # If the table is molecules or references, remove the "id" column
        if data:
            if table_chosen in ("molecules", "references"):
                columns = [
                    {"name": column_names.get(i, i), "id": i}
                    for i in data[0].keys()
                    if i != "id"
                ]
                # Also remove the "id" key from each row
                data = [{k: v for k, v in row.items() if k != "id"} for row in data]
            else:
                columns = [
                    {"name": column_names.get(i, i), "id": i} for i in data[0].keys()
                ]
        else:
            columns = []

        table = dash_table.DataTable(
            columns=columns,
            data=data,
            sort_action="native",
            style_table={"overflowX": "hidden"},
            style_cell={
                "textAlign": "left",
                "fontFamily": "Roboto, sans-serif",
                "fontSize": "14px",
            },
            style_data_conditional=[
                {
                    "if": {"column_id": column_id},
                    "minWidth": width,
                    "maxWidth": width,
                    "width": width,
                }
                for column_id, width in widths.items()
            ],
            style_header={
                "backgroundColor": "rgb(230, 230, 230)",
                "fontWeight": "bold",
                "fontFamily": "Roboto, sans-serif",
                "fontSize": "16px",
            },
            page_size=items_per_page,
            page_current=0,
            page_action="native",
        )

        return table, selected_molecule
