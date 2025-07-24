import geopandas as gpd
from shapely.geometry import Point

# Percorso file di input e output
INPUT_PATH = "../Data/Raw/Punti_di_interesse/farmacie_raw.geojson"
OUTPUT_PATH = "../Data/Clean/Punti_di_interesse/farmacie_clean.geojson"

# Carica il file GeoJSON e Rinomina le colonne
gdf = gpd.read_file(INPUT_PATH)
column = {
    "DESCRIZIONE_FARMACIA": "Descrizione_farmacia",
    "CAP":"CAP",
    "MUNICIPIO": "Municipio",
    "NIL": "Quartiere",
    "LONGITUDINE": "Longitudine",
    "LATITUDINE": "Latitudine",
    "geometry": "geometry"
}
gdf = gdf[column.keys()]
gdf["icon"] = "control-on"
gdf = gdf.rename(columns=column)
gdf.to_file(OUTPUT_PATH, driver="GeoJSON")