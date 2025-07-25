import pandas as pd
import geopandas as gpd

# Percorso file di input e output
INPUT_PATH = "../Data/Raw/Punti_di_interesse/aree_gioco_raw.geojson"
OUTPUT_PATH = "../Data/Clean/Punti_di_interesse/aree_gioco_clean.geojson"

gdf = gpd.read_file(INPUT_PATH)

# Modifico i valori mancanti di "data_ini"
valori_mancanti = gdf["data_ini"].isna() | (gdf["data_ini"].astype(str).str.strip() == "")

def get_fill_date(row, reference_gdf):
    area = row["area"]
    localita = row["località"]
    
    match = reference_gdf[
        (reference_gdf["area"] == area) &
        (reference_gdf["località"] == localita) &
        (~reference_gdf["data_ini"].isna()) &
        (reference_gdf["data_ini"].astype(str).str.strip() != "")
    ]
    
    if not match.empty:
        return match.iloc[0]["data_ini"]
    else:
        return "2025-01-01"

gdf.loc[valori_mancanti, "data_ini"] = gdf[valori_mancanti].apply(lambda row: get_fill_date(row, gdf), axis=1)
gdf["data_ini"] = pd.to_datetime(gdf["data_ini"]).dt.year

# Filtro gli anni che non mi interessano
gdf = gdf[~gdf["data_ini"].isin([2024, 2025])]

# Elimino colonne superflue
gdf = gdf.drop(columns=["id_area", "obj_id", "codice", "descrizione_codice"])

# Salvo GeoJSON pulito
gdf.to_file(OUTPUT_PATH, driver="GeoJSON")