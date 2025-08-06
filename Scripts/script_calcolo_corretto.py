import pickle
import osmnx as ox
from collections import defaultdict

with open("C:/Users/dalla/OneDrive/Documenti/Project/Progetto-Bike/Data/Raw/Info_bici/strade_ciclabili.pickle", 'rb') as x:    
    G = (pickle.load(x))
gdf = ox.graph_to_gdfs(G, nodes=False)

# Trasformazione righe combinate in liste
import ast

def safe_parse_highway(hwy):
    if isinstance(hwy, list):
        return hwy
    if isinstance(hwy, str) and hwy.startswith("["):
        try:
            return ast.literal_eval(hwy)
        except:
            return [hwy]
    return [hwy]

gdf["highway"] = gdf["highway"].apply(safe_parse_highway)

#gdf = gdf.explode("highway").reset_index(drop=True)


# 1. Raw data filtrato (solo tipi validi)
raw_data = {
    "primary": 2298,
    "tertiary": 9358,
    "residential": 25570,
    "unclassified": 14870,
    "secondary": 4055,
    "footway": 178667,
    "motorway_link": 230,
    "motorway": 81,
    "construction": 173,
    "steps": 1072,
    "pedestrian": 9802,
    "primary_link": 137,
    "service": 36717,
    "cycleway": 7543,
    "tertiary_link": 131,
    "secondary_link": 104,
    "living_street": 169,
    "path": 2110,
    "proposed": 87,
    "trunk_link": 90,
    "bridleway": 50,
    "track": 1212,
    "trunk": 16,
    "busway": 28,
    "services": 38,
    "demolished": 4,
    "corridor": 6,
    "platform": 78,
    "bus_stop": 4,
    "emergency_bay": 2
}

# 2. Mappa macrocategorie in base a linee guida
macro_map = {
    "road": [
        "motorway", "trunk", "primary", "secondary", "tertiary",
        "residential", "unclassified", "living_street"
    ],
    "link": [
        "motorway_link", "trunk_link", "primary_link",
        "secondary_link", "tertiary_link"
    ],
    "pedestrian": [
        "footway", "pedestrian", "steps", "platform", "corridor"
    ],
    "cycleway": [
        "cycleway"
    ],
    "mixed": [
        "service", "track", "path", "bridleway"
    ],
    "construction": [
        "construction", "proposed", "demolished"
    ],
    "public_transport": [
        "busway", "bus_stop", "platform"
    ],
    "other": [
        "services", "emergency_bay"
    ]
}

# 3. Pesi personalizzati
custom_weights = {
    "motorway": 99,
    "trunk": 99,
    "primary": 1.2,
    "secondary": 1.1,
    "tertiary": 1.25,
    "residential": 1.3,
    "unclassified": 1.3,
    "living_street": 1.3,
    "motorway_link": 99,
    "trunk_link": 99,
    "primary_link": 1.2,
    "secondary_link": 1.1,
    "tertiary_link": 1.25,
    "service": 99,
    "footway": 1.2,
    "pedestrian": 1.5,
    "steps": 99,
    "platform": 99,
    "path": 1.1,
    "cycleway": 1,
    "track": 1,
    "bridleway": 1.6,
    "construction": 99,
    "proposed": 99,
    "demolished": 99,
    "busway": 99,
    "bus_stop": 99,
    "corridor": 99,
    "services": 99,
    "emergency_bay": 99
}


# 4. Costruzione dizionario e ordinamento

category_by_type = {}
for macro, types in macro_map.items():
    for t in types:
        category_by_type[t] = macro

result = defaultdict(dict)

for k in raw_data:
    if k.startswith("["):  # ignora chiavi ambigue
        continue
    weight = custom_weights.get(k)
    if weight is None:
        continue
    category = category_by_type.get(k, "unknown")
    result[category][k] = weight

# Ordina per peso decrescente in ciascuna categoria
for macro in result:
    sorted_items = sorted(result[macro].items(), key=lambda x: -x[1])
    result[macro] = dict(sorted_items)

def get_weight(hwy):
    if isinstance(hwy, list):
        # Trova il valore più vicino a 1
        return min([custom_weights.get(i, 0) for i in hwy], key=lambda x: abs(x - 1))
    else:
        return custom_weights.get(hwy, 0)
    
gdf["weight"] = gdf["highway"].apply(get_weight)

# Stessa funzione di prima ma qui prendiamo il max
#def get_weight(hwy):
    #if isinstance(hwy, list):
        #return max([custom_weights.get(i, 0) for i in hwy])  # usa il valore massimo tra i tipi
    #else:
        #return custom_weights.get(hwy, 0)  # default 0 se non trovato

#6. Esporta in GeoJSON
output_path = "Aree gioco/strade_peso_corretto.geojson"
gdf.to_file(output_path, driver='GeoJSON')