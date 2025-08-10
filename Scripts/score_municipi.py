import geopandas as gpd
import pandas as pd
from my_paths import *

# Carica dati
municipi = gpd.read_file(PATH_MUNICIPI_CLEAN)
parchi = gpd.read_file(PATH_PARCHI_CLEAN)
fontane = gpd.read_file(PATH_FONTANE_CLEAN)
impianti = gpd.read_file(PATH_IMPIANTI_SPORTIVI_CLEAN)
scuole = gpd.read_file(PATH_SCUOLE_CLEAN)
aree_gioco = gpd.read_file(PATH_AREE_GIOCO_CLEAN)
biblioteche = gpd.read_file(PATH_BIBLIOTECHE_CLEAN)
stazioni_bikemi = gpd.read_file(PATH_BIKEMI_CLEAN)
case_acqua = gpd.read_file(PATH_CASE_ACQUA_CLEAN)
cinema = gpd.read_file(PATH_CINEMA_CLEAN)
farmacie = gpd.read_file(PATH_FARMACIE_CLEAN)
musei = gpd.read_file(PATH_MUSEI_CLEAN)
ospedali = gpd.read_file(PATH_OSPEDALI_CLEAN)
teatri = gpd.read_file(PATH_TEATRI_CLEAN)

# Aggiungiamo i POI 
parchi["tipo"] = "parchi"
fontane["tipo"] = "fontane"
impianti["tipo"] = "impianti_sportivi"
scuole["tipo"] = "scuole"
aree_gioco["tipo"] = "aree_gioco"
biblioteche["tipo"] = "biblioteche"
stazioni_bikemi["tipo"] = "stazioni_bikemi"
case_acqua["tipo"] = "case_acqua"
cinema["tipo"] = "cinema"
farmacie["tipo"] = "farmacie"
musei["tipo"] = "musei"
ospedali["tipo"] = "ospedali"
teatri["tipo"] = "teatri"

# Unisco POI in un unico GeoDataFrame
pois = pd.concat([
    parchi, fontane, impianti, scuole, aree_gioco, biblioteche,
    stazioni_bikemi, case_acqua, cinema, farmacie, musei, ospedali, teatri
], ignore_index=True)

# Join spaziale per associare POI al municipio
pois_con_municipio = gpd.sjoin(pois, municipi, how="left", predicate="intersects")

# Raggruppa per MUNICIPIO e tipo POI
conteggi = pois_con_municipio.groupby(['MUNICIPIO_right', 'tipo']).size().unstack(fill_value=0).reset_index()

# Assicura che MUNICIPIO nei municipi e MUNICIPIO_right in conteggi siano dello stesso tipo
municipi['MUNICIPIO'] = municipi['MUNICIPIO'].astype(int)
conteggi['MUNICIPIO_right'] = conteggi['MUNICIPIO_right'].astype(int)

# Unisci i conteggi ai municipi con il giusto join
municipi = municipi.merge(conteggi, left_on='MUNICIPIO', right_on='MUNICIPIO_right', how='left')

# Riempi NaN con zero e converti in int per tutte le colonne dei POI
poi_columns = conteggi.columns.drop('MUNICIPIO_right')
for col in poi_columns:
    municipi[col] = municipi[col].fillna(0).astype(int)

# Rimuovi colonna MUNICIPIO_right perché non serve più
municipi = municipi.drop(columns=['MUNICIPIO_right'])


# Salva GeoJSON pronto per Kepler
municipi.to_file("score_municipi.geojson", driver="GeoJSON")
