import json

# Percorso file di input e output
INPUT_PATH = "Data/Raw/Punti_di_interesse/teatri_raw.json"
OUTPUT_PATH = "Data/Clean/Punti_di_interesse/teatri_clean.geojson"

# Carico il file JSON originale
with open(INPUT_PATH, "r", encoding="utf-8") as f:
    data = json.load(f)

# Filtro per sottotipo "TEATRO" o "AUDITORIUM" ed elimino l'unica riga con "zona":"0" (Chiuso definitivamente)
filtrati = [record for record in data if record["sottotipo"].upper() in ["TEATRO", "AUDITORIUM"] and str(record.get("zona", "0")).strip() not in ["0", ""]]

# Creo le feature per il file GeoJSON
features = []
for r in filtrati:
    feature = {
        "type": "Feature",
        "geometry": {
            "type": "Point",
            "coordinates": [
                float(r["glongitude"]),
                float(r["glatitude"])
            ]
        },
        "properties": {
            "nome_teatro": r["denominaz"],
            "sottotipo": r["sottotipo"],
            "id_teatro": r["id"],
            "municipio": r["zona"],
            "longitudine": r["glongitude"],
            "latitudine": r["glatitude"]
        }
    }
    features.append(feature)

# Converto in GeoJSON e salvo il file
geojson_data = {
    "type": "FeatureCollection",
    "features": features
}
with open(OUTPUT_PATH, "w", encoding="utf-8") as f_out:
    json.dump(geojson_data, f_out, ensure_ascii=False, indent=2)
