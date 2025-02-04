import dash_bootstrap_components as dbc
from dash import dash_table
from dash import Dash, html, dcc, Input, Output, callback_context

from .server import db
from .models import Drugs, Molecules, References


def create_dashapp(server):
    # Initialize the app with a dark Bootstrap theme
    app = Dash(
        __name__,
        server=server,
        external_stylesheets=[dbc.themes.DARKLY],
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
            dbc.Navbar(
                dbc.Container(
                    [
                        dbc.NavbarBrand("Drug Oxidation Database", className="ms-2"),
                        dbc.Nav(
                            [
                                dbc.NavItem(dbc.NavLink("Home", href="#")),
                                dbc.NavItem(dbc.NavLink("About", href="#")),
                            ],
                            className="ms-auto",
                        ),
                    ],
                    fluid=True,
                ),
                color="dark",
                dark=True,
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Tabs(
                                [
                                    dbc.Tab(label="Tables", tab_id="tab-tables"),
                                    dbc.Tab(label="Molecular View", tab_id="tab-molecular"),
                                ],
                                id="tabs",
                                active_tab="tab-tables",
                            ),
                            html.Div(id="tab-content"),
                        ],
                        width=12,
                    ),
                ]
            ),
        ],
        fluid=True,
    )

    @app.callback(
        Output("tab-content", "children"),
        Input("tabs", "active_tab"),
    )
    def render_tab_content(active_tab):
        if active_tab == "tab-tables":
            return dbc.Row(
                [
                    dbc.Col(
                        [
                            html.H2(
                                "Available Tables",
                                className="text-left",
                                style={
                                    "fontSize": "16px",
                                    "fontWeight": "bold",
                                    "color": "#343a40",
                                    "fontFamily": "Roboto, sans-serif",
                                },
                            ),
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
                                            style={
                                                "fontFamily": "Roboto, sans-serif",
                                                "fontSize": "16px",
                                                "fontWeight": "bold",
                                                "color": "#343a40",
                                            },
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
                                                "color": "#343a40",
                                            },
                                        ),
                                        width=True,
                                    ),
                                ],
                                align="center",
                            ),

                            html.Hr(),
                            dcc.Dropdown(
                                id="dropdown",
                                options=drug_options,
                                value=None,
                                style={
                                    "fontFamily": "Roboto, sans-serif",
                                    "fontSize": "14px",
                                    "color": "#343a40",  # Darker font color for dropdown text
                                },
                            ),
                            html.Br(),
                        ],
                        width=2,
                        className="bg-light p-3",
                    ),
                    dbc.Col([html.Div(id="first-graph")], width=10),
                ]
            )
        elif active_tab == "tab-molecular":
            return dbc.Row(
                [
                    dbc.Col(
                        [
                            html.H2(
                                "Molecular Controls",
                                className="text-left",
                                style={
                                    "fontSize": "16px",
                                    "fontWeight": "bold",
                                    "color": "#343a40",
                                    "fontFamily": "Roboto, sans-serif",
                                },
                            ),
                            dcc.Dropdown(
                                id="molecular-dropdown",
                                options=[
                                    {"label": "Option 1", "value": "option1"},
                                    {"label": "Option 2", "value": "option2"},
                                ],
                                value="option1",
                                style={
                                    "fontFamily": "Roboto, sans-serif",
                                    "fontSize": "14px",
                                    "color": "#343a40",
                                },
                            ),
                            html.Br(),
                            dbc.Button(
                                "Apply",
                                id="apply-button",
                                color="primary",
                                className="me-1",
                            ),
                        ],
                        width=2,
                        className="bg-light p-3",
                    ),
                    dbc.Col(
                        html.Div("Molecular View Content"),
                        width=10,
                    ),
                ]
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
            "drugs": {"drug_id": "150px", "drug_title": "300px", "smiles": "1000px"},
            "molecules": {
                "drug_id": "150px",
                "drug_title": "150px",
                "som": "150px",
                "som_element": "150px",
                "som_level": "150px",
            },
            "references": {"drug_id": "150px", "reference": "500px", "doi": "1000px"},
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
                "color": "#343a40",
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
                "color": "#343a40",
            },
            page_size=items_per_page,
            page_current=0,
            page_action="native",
            cell_selectable=False,
        )

        return table, selected_molecule
