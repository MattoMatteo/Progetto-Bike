import json

with open("Impianti sportivi\impianti_sportivi_raw.geojson", "r", encoding="utf-8") as f:
    data = json.load(f)

# Conversione da MultiPolygon a Polygon, ha senso farlo perchè ogni MultiPolygon conteneva un solo poligono
for feature in data["features"]:
    geom = feature["geometry"]
    if geom["type"] == "MultiPolygon":
        first_polygon = geom["coordinates"][0]
        feature["geometry"] = {
            "type": "Polygon",
            "coordinates": first_polygon
        }

#Pulizia Dataset
cleaned_features = []

for feature in data["features"]:
    props = feature["properties"]
    geometry = feature["geometry"]

    # Salta se 'data_ini' è mancante o è nel 2024
    data_ini = props.get("data_ini")
    if not data_ini or data_ini.startswith("2024"):
        continue

    new_props = {
        "municipio": props.get("municipio"),
        "località": props.get("località"),
        "id_impianto": props.get("obj_id"), 
        "data_ini": data_ini,
        "area_mq": props.get("area_mq"),
        "perim_m": props.get("perim_m")
    }

    # Crea la nuova feature con proprietà filtrate e geometria
    cleaned_features.append({
        "type": "Feature",
        "properties": new_props,
        "geometry": geometry
    })

cleaned_geojson = {
    "type": "FeatureCollection",
    "features": cleaned_features
}

# Salvataggio file in .geojson
with open("impianti_sportivi_staging.geojson", "w", encoding="utf-8") as f:
    json.dump(cleaned_geojson, f, ensure_ascii=False, indent=2)
