import geopandas as gpd

# Caricamento dati
gdf_sport = gpd.read_file("C:/Users/dalla/OneDrive/Documenti/Project/Progetto-Bike/Data/Clean/Analisi/sport_tempo_libero_clean.geojson")
gdf_municipi = gpd.read_file("C:/Users/dalla/OneDrive/Documenti/Project/Progetto-Bike/Data/Clean/Punti_di_interesse/municipi_clean.geojson")

# Uniformare CRS
gdf_sport = gdf_sport.to_crs(gdf_municipi.crs)

segmenti = []

for _, sport_row in gdf_sport.iterrows():
    geom = sport_row["geometry"]
    municipi_intersectati = []

    # Trova municipi intersecati
    for _, muni_row in gdf_municipi.iterrows():
        muni_geom = muni_row["geometry"]
        if muni_geom.intersects(geom):
            municipi_intersectati.append((muni_row["MUNICIPIO"], muni_row.get("NOME", None), muni_geom))

    if not municipi_intersectati:
        nuova_riga = sport_row.copy()
        nuova_riga["MUNICIPIO"] = None
        nuova_riga["NOME_MUNICIPIO"] = None
        segmenti.append(nuova_riga)
        continue

    # Spezza geometria per ogni municipio
    for municipio, nome_municipio, muni_geom in municipi_intersectati:
        parte = muni_geom.intersection(geom)
        if parte.is_empty:
            continue

        if parte.geom_type == "GeometryCollection":
            for sub_geom in parte.geoms:
                if sub_geom.is_empty or (sub_geom.length == 0 and sub_geom.area == 0):
                    continue
                nuova_riga = sport_row.copy()
                nuova_riga["geometry"] = sub_geom
                nuova_riga["MUNICIPIO"] = municipio
                nuova_riga["NOME_MUNICIPIO"] = nome_municipio
                segmenti.append(nuova_riga)
        else:
            nuova_riga = sport_row.copy()
            nuova_riga["geometry"] = parte
            nuova_riga["MUNICIPIO"] = municipio
            nuova_riga["NOME_MUNICIPIO"] = nome_municipio
            segmenti.append(nuova_riga)

gdf_segmenti = gpd.GeoDataFrame(segmenti, crs=gdf_sport.crs, geometry="geometry")

# Reproiezione a CRS metrico per calcolare lunghezza corretta (esempio EPSG:32632 per Milano)
gdf_segmenti_metrico = gdf_segmenti.to_crs(epsg=32632)
gdf_segmenti["lunghezza_km"] = gdf_segmenti_metrico.length / 1000

#Torna al CRS originale se vuoi mantenere lat/lon
gdf_segmenti = gdf_segmenti.to_crs(gdf_sport.crs)

# Salvataggio GeoJSON
output_path = "C:/Users/dalla/OneDrive/Desktop/nuovi_punti_di_interesse/sport_segmentati_per_municipio.geojson"
gdf_segmenti.to_file(output_path, driver="GeoJSON")

