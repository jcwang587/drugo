import dash_bootstrap_components as dbc
from dash import dash_table
from dash import Dash, html, dcc, Input, Output
from .server import db
from .models import Drugs, Molecules, References
from .chem import draw_molecule


def create_dashapp(server, db_version):
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
                    (drug.drug_title[:20] + "...")
                    if len(drug.drug_title) > 20
                    else drug.drug_title
                ),
                "value": drug.drug_id,
                "title": drug.drug_title,
            }
            for drug in db.session.query(Drugs).all()
        ]

    # App layout with tabs: Tables, Molecular View, Resources
    app.layout = dbc.Container(
        [
            dbc.Navbar(
                dbc.Container(
                    [
                        dbc.NavbarBrand(
                            [
                                "Drug Oxidation Database ",
                                html.Span(
                                    db_version,
                                    style={
                                        "fontSize": "14px",
                                        "color": "lightgrey",
                                    },
                                ),
                            ],
                            className="ms-2",
                        ),
                        dbc.Nav(
                            [
                                dbc.NavItem(
                                    dbc.NavLink(
                                        "Home",
                                        href="https://drugo-a54338d8b0d8.herokuapp.com/",
                                    )
                                ),
                                dbc.NavItem(
                                    dbc.NavLink(
                                        "GitHub Repository",
                                        href="https://github.com/jcwang587/drugo",
                                    )
                                ),
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
                                    dbc.Tab(
                                        label="Molecular View", tab_id="tab-molecular"
                                    ),
                                    dbc.Tab(label="Resources", tab_id="tab-resources"),
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
                                className="text-white",
                                style={
                                    "fontSize": "16px",
                                    "fontWeight": "bold",
                                    "fontFamily": "Roboto, sans-serif",
                                },
                            ),
                            dcc.RadioItems(
                                options=[
                                    {"label": " Drugs", "value": "drugs"},
                                    {
                                        "label": " Sites of Metabolism",
                                        "value": "molecules",
                                    },
                                    {"label": " References", "value": "references"},
                                ],
                                value="drugs",
                                id="controls",
                                inline=False,
                                style={
                                    "fontFamily": "Roboto, sans-serif",
                                    "color": "white",
                                },
                            ),
                            html.Br(),
                            dbc.Row(
                                [
                                    dbc.Col(
                                        html.Label(
                                            "Items per page",
                                            className="text-white",
                                            style={
                                                "fontFamily": "Roboto, sans-serif",
                                                "fontSize": "16px",
                                                "fontWeight": "bold",
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
                                            value=15,
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
                            html.Br(),
                        ],
                        width=2,
                        className="bg-dark p-3",
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
                                "Molecules",
                                className="text-white",
                                style={
                                    "fontSize": "16px",
                                    "fontWeight": "bold",
                                    "fontFamily": "Roboto, sans-serif",
                                },
                            ),
                            dcc.Dropdown(
                                id="molecule-dropdown",
                                options=drug_options,
                                value=None,
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
                        className="bg-dark p-3",
                    ),
                    dbc.Col(
                        html.Div(id="molecular-svg", className="text-white"),
                        width=9,
                    ),
                ],
                align="start",
            )

    @app.callback(
        Output("molecular-svg", "children"),
        Input("molecule-dropdown", "value"),
    )
    def update_molecular_svg(selected_molecule_id):
        if selected_molecule_id is None:
            return "Select a molecule to view its details."

        # Fetch the selected molecule's name from the database
        selected_molecule = (
            db.session.query(Drugs).filter_by(drug_id=selected_molecule_id).first()
        )
        if not selected_molecule:
            return "Molecule not found."

        # Create an SVG with the molecule's name
        # svg_content = f"""
        # <svg width="200" height="100" xmlns="http://www.w3.org/2000/svg">
        #     <rect width="100%" height="100%" fill="lightgrey" />
        #     <text x="10" y="50" font-family="Roboto, sans-serif" font-size="20" fill="black">
        #         {selected_molecule.drug_title}
        #     </text>
        # </svg>
        # """

        png_content = draw_molecule(selected_molecule.smiles)

        return html.Div(
            html.Img(
                src=f"data:image/png;base64,{png_content}",
                style={"border": "none", "height": "100%"},
            )
        )

    @app.callback(
        Output("first-graph", "children"),
        Input("controls", "value"),
        Input("items-per-page-dropdown", "value"),
    )
    def update_table(table_chosen, items_per_page):
        table_map = {"drugs": Drugs, "molecules": Molecules, "references": References}
        model_class = table_map.get(table_chosen)

        if model_class is None:
            return html.Div("No data available", className="text-white")

        query = db.session.query(model_class).all()
        data = [
            {
                column.name: getattr(row, column.name)
                for column in model_class.__table__.columns
            }
            for row in query
        ]

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

        if data:
            if table_chosen in ("molecules", "references"):
                if table_chosen == "references":
                    for row in data:
                        doi_value = row.get("doi")
                        if doi_value:
                            # Convert the DOI value (full URL) into a clickable Markdown link
                            row["doi"] = f"[{doi_value}]({doi_value})"
                # Build the columns while skipping the internal "id" column
                columns = []
                for key in data[0].keys():
                    if key == "id":
                        continue
                    col = {"name": column_names.get(key, key), "id": key}
                    # Specify Markdown presentation for the DOI column
                    if table_chosen == "references" and key == "doi":
                        col["presentation"] = "markdown"
                    columns.append(col)
                # Remove the "id" key from each row
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
            style_table={
                "overflowX": "auto",
                "tableLayout": "fixed",
                "width": "100%",
            },
            # Data cells: light background and dark text
            style_cell={
                "textAlign": "left",
                "fontFamily": "Roboto, sans-serif",
                "fontSize": "14px",
                "backgroundColor": "white",
                "color": "#343a40",
                "padding": "5px",
                "overflow": "hidden",
                "textOverflow": "ellipsis",
                "maxWidth": 0,
            },
            # Header: dark background with white text
            style_header={
                "backgroundColor": "#343a40",
                "fontWeight": "bold",
                "fontFamily": "Roboto, sans-serif",
                "fontSize": "16px",
                "color": "white",
            },
            page_size=items_per_page,
            page_current=0,
            page_action="native",
            cell_selectable=False,
        )

        return table

    return app
