import pickle
import json
import math
from collections import defaultdict

import geopandas as gpd
from shapely import Point, LineString, MultiLineString, MultiPolygon, Polygon, GeometryCollection
from shapely.prepared import prep
import networkx as nx
from networkx.algorithms.approximation import steiner_tree
import osmnx as ox
from scipy.spatial import cKDTree
import numpy as np

from my_paths import *

# Funzioni Private

def _explode_gdf_with_MultiPolygon_to_Polygon(gdf:gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Converte un GeoDataFrame che contiene dei MultiPolygon in Polygon, creando nuove righe per i nuovi Polygon separati,
    mantenendo le medesime informazioni delle altre colonne su tutte le righe.
    Tutte le geometrie che non sono MultiPolygon rimarranno inalterate.
    """
    # Esplodiamo i Multypoligon se presenti
    exploded_rows = []
    for idx, row in gdf.iterrows():
        geometry = row["geometry"]
        if isinstance(geometry, MultiPolygon):
            for polygon in geometry.geoms:
                new_row = row.copy()
                new_row["geometry"] = polygon
                exploded_rows.append(new_row)
        else:
            exploded_rows.append(row)
    return gpd.GeoDataFrame(exploded_rows, crs=gdf.crs).reset_index(drop=True)

def _get_points_of_gdf(gdf:gpd.GeoDataFrame) -> dict:
    """
    Dato un GeoDataFrame con geometria di tipo Point o Polygon, restituisce una collezione dei
    suoi Point, in formato dizionario che ha come chiave l'indice di riga del GeoDataFrame e come
    valore una tupla di long/lat: {idx: (long, lat), ...}
    """
    poi_nodes = {}
    for idx, row in gdf.iterrows():
        geometry = row["geometry"]
        if isinstance(geometry, Polygon):
            poi_nodes[idx] = (geometry.centroid.x, geometry.centroid.y)
        elif isinstance(geometry, Point):
            poi_nodes[idx] = (geometry.x, geometry.y)
        else:
            print(f"[WARNING] Geometria non valida (id={idx}): {type(geometry)} — ignorata.")
    return poi_nodes

def _change_edge_weight_in_gdf(G:nx.MultiDiGraph, gdf:gpd.GeoDataFrame, weight:float=0) -> nx.MultiDiGraph:
    geometry_collection = prep(gdf.union_all())
    for _, __, data in G.edges(data=True):
        geom = data.get("geometry", None)
        if geom is None:
            continue

        if geometry_collection.intersects(geom):
            data["weight"] = weight

    return G

def _custom_graph_to_gdf(G:nx.MultiDiGraph) -> gpd.GeoDataFrame:
    gdf_edges = ox.graph_to_gdfs(G, edges=True, nodes=False).reset_index().drop(["u", "v", "key"], axis=1)
    # formato geojson non supporta "liste" come valori possibili nelle colonne quindi li convertiamo in stringhe
    colonne_da_esplodere = ["name", "highway", "access", "tunnel", "service"]
    for colonna in colonne_da_esplodere:    
        gdf_edges[colonna] = gdf_edges[colonna].apply(lambda x: " ".join(x) if isinstance(x, list) else x)
    
    def convert_maxspeed(x):
        if x is None or (isinstance(x, float) and math.isnan(x)):
            return None
        elif isinstance(x, list):
            max_speed_list = []
            for y in x:
                if isinstance(y, str):
                    if y == "IT:urban":
                        y = 50
                max_speed_list.append(int(y))
            return max(max_speed_list)
        elif isinstance(x, str):
            if x == "IT:urban":
                return 50
            else:
                return int(x)
        return int(x)

    gdf_edges["maxspeed"] = gdf_edges["maxspeed"].apply(convert_maxspeed)

    gdf_municipi = gpd.read_file(PATH_MUNICIPI_CLEAN)

    segmenti = []
    for _, sport_row in gdf_edges.iterrows():
        geom = sport_row["geometry"]
        municipi_intersectati = []

        # Trova municipi intersecati
        for _, muni_row in gdf_municipi.iterrows():
            muni_geom = muni_row["geometry"]
            if muni_geom.intersects(geom):
                municipi_intersectati.append((muni_row["MUNICIPIO"], muni_geom))
        # Spezza geometria per ogni municipio
        for municipio, muni_geom in municipi_intersectati:
            parte = muni_geom.intersection(geom)

        if isinstance(parte, LineString):
            nuova_riga = sport_row.copy()
            nuova_riga["geometry"] = parte
            nuova_riga["MUNICIPIO"] = municipio
            segmenti.append(nuova_riga)

        elif isinstance(parte, MultiLineString):
            for subgeom in parte.geoms:
                nuova_riga = sport_row.copy()
                nuova_riga["geometry"] = subgeom
                nuova_riga["MUNICIPIO"] = municipio
                segmenti.append(nuova_riga)

    gdf_segmenti = gpd.GeoDataFrame(segmenti, crs=CRS_GRAD, geometry="geometry")
    gdf_segmenti["length"] = gdf_segmenti.to_crs(CRS_METR).length
    return gdf_segmenti

def _get_proportional_reduce_n_pois(poi_group:list[str],
                                    max_pois:int,
                                    peso_inquinamento_vs_incidenti:float = 0.5,
                                    peso_fattori_esterni_vs_distribuzione_poi:float = 0.5):
    """
    Dato un gruppo di "poi" e un numero massimo di poi da "tenere", 
    restituisce un dataframe con i nuovi valori, ridotti in base ai pesi assegnati, per ogni categoria e municipio.
    """
    # Carichiamo da file il gdf dei conteggi di tutti i poi che abbiamo, ragruppati per Municipio
    gdf_conteggi = gpd.read_file(PATH_SCORE_X_MUNICIPI)
    # Calcoliamo i totali sommando per riga tutte le colonne dei "poi"
    gdf_conteggi_gruppo = gdf_conteggi[poi_group].copy()
    gdf_conteggi_gruppo["Totale"] = gdf_conteggi_gruppo.sum(axis=1)
    gdf_conteggi_gruppo["geometry"] = gdf_conteggi["geometry"].to_numpy()
    
    def converti_totale(row:gpd.GeoSeries, tot_poi:int, max_poi:int):
        """
        Funzione da usare in ".apply" per ottenere i valori risultanti dalla
        proporzione, avendo il nuovo max_poi:
        x = numero_colonna * numero_massimo_poi / numero_totale
        """
        # Da inserire nella formula anche le variabili di:
        # Inquinamento e incidenti
        if max_poi > tot_poi:
            max_poi = tot_poi
        risultato = int(round((row["Totale"]*max_poi/tot_poi)))
        # Almeno 1 per tipo
        if risultato == 0:
            return 1
        return risultato

    # Riduciamo proporzionalmente i totali di ogni municipio, creando "nuovo_totale"
    totale_poi_gdf = int(gdf_conteggi_gruppo["Totale"].sum())   # Da dare alla funzione
    gdf_conteggi_gruppo["nuovo_totale"] = gdf_conteggi_gruppo.apply(lambda row: converti_totale(row, totale_poi_gdf, max_pois), axis=1)
    gdf_conteggi_gruppo["MUNICIPIO"] = gdf_conteggi_gruppo.index +1

    # -------------------------------------------------------------#
    # ---       Con i nuovi totali possiamo ridistribuire       ----
    # ---   in base a inquinamento e incidenti i nuovi totali   ----

    # Uniamo le colonne dei dataframe che tanto hanno lo stesso numero di righe ordinate per municipio:
    # conteggi + incidenti + inquinamento
    gdf_conteggi_gruppo["pm25"] = gpd.read_file(PATH_INQUINAMENTO_CLEAN)["pm25"]
    gdf_conteggi_gruppo["Incidenti"] = gpd.read_file(PATH_INCIDENTI_CLEAN)["Incidenti"]
    # Normalizziamo i valori di incidenti e inquinamento -> in % rispetto al valore massimo in tabella
    gdf_conteggi_gruppo['pm25_norm'] = gdf_conteggi_gruppo['pm25'] / gdf_conteggi_gruppo['pm25'].max()
    gdf_conteggi_gruppo['incidenti_norm'] = gdf_conteggi_gruppo['Incidenti'] / gdf_conteggi_gruppo['Incidenti'].max()
    # Unifichiamo i 2 fattori in un solo "indice_priorità" utilizzando i rispettivi pesi
    w_pm = peso_inquinamento_vs_incidenti # peso dell’inquinamento
    w_inc = 1 - peso_inquinamento_vs_incidenti  # peso degli incidenti
    gdf_conteggi_gruppo['indice_priorita'] = w_pm * gdf_conteggi_gruppo['pm25_norm'] + w_inc * gdf_conteggi_gruppo['incidenti_norm']
    
    # Rispetto ai propri totali, calcoliamo le percentuali per ogni municipio
    #  sia per il fattore "territorio" (incidenti + inquinamenti) che per quello "base" (distribuzione poi)
    gdf_conteggi_gruppo['peso_territorio'] = gdf_conteggi_gruppo['indice_priorita'] / gdf_conteggi_gruppo['indice_priorita'].sum()
    gdf_conteggi_gruppo['peso_base'] = gdf_conteggi_gruppo['nuovo_totale'] / gdf_conteggi_gruppo['nuovo_totale'].sum()

    # Uniamo i 2 pesi dei 2 fattori in uno unico "peso finale"
    peso_fattori_esterni = peso_fattori_esterni_vs_distribuzione_poi
    peso_distribuzione_poi = 1 - peso_fattori_esterni_vs_distribuzione_poi
    gdf_conteggi_gruppo['peso_finale'] = peso_distribuzione_poi * gdf_conteggi_gruppo['peso_base'] + peso_fattori_esterni * gdf_conteggi_gruppo['peso_territorio']
    
    # Calcoliamo grazie al peso finale il nuovo totale per ogni municipio
    totale_poi = gdf_conteggi_gruppo["nuovo_totale"].sum()
    gdf_conteggi_gruppo['poi_teorico'] = (totale_poi * gdf_conteggi_gruppo['peso_finale']).round().astype(int)
    gdf_conteggi_gruppo.sort_values("poi_teorico", ascending=False, inplace=True)

    # Sistema per non eccede con i valori di attribuzione per ogni municipio
    gdf_conteggi_gruppo['poi_finali'] = 0
    accumulo = 0
    for idx, row in gdf_conteggi_gruppo.iterrows():
        if row["Totale"] - row["poi_teorico"] > 0:
            gdf_conteggi_gruppo.loc[idx, "poi_finali"] = row["poi_teorico"]
        else:
            gdf_conteggi_gruppo.loc[idx, "poi_finali"] = row["Totale"]
            accumulo += row["poi_teorico"] - row["Totale"]
    for idx, row in gdf_conteggi_gruppo.iterrows():
        if row["poi_finali"] < row["Totale"]:
            if row["Totale"] - row["poi_finali"] > accumulo:
                gdf_conteggi_gruppo.loc[idx, "poi_finali"] += accumulo
                break
            else:
                accumulo -= row["Totale"] - row["poi_finali"]
                gdf_conteggi_gruppo.loc[idx, "poi_finali"] = row["Totale"]

    # --------------------------------------------------------------------------------#
    # ----      Basandoci sui nuovi totali ridistributi per municipio,      -----------
    # ---- settiamo i singoli valori per ogni "tipo" di poi per municipio   -----------
    for idx, row in gdf_conteggi_gruppo.iterrows():
        for colonna in poi_group:
            gdf_conteggi_gruppo.loc[idx, colonna] = round(row[colonna] / row["Totale"] * row["poi_finali"])
    return gdf_conteggi_gruppo

# Funzioni Pubbliche

def geojson_to_graph(gdf: gpd.GeoDataFrame, **attr) -> nx.MultiDiGraph:
    """
    Dato un geojson caricato come GeoDataFrame avente geometry di tipo LineString, lo converte in un MultiDiGraph, prendendo come
    nodi, ogni Point delle LineString e come archi, ogni Point con quello successivo (p[x], p[x+1]) di ogni LineString.
    """
    G = nx.MultiDiGraph()
    node_id = 0
    coord_to_node = {}  # Mappa coord (x,y) → ID nodo

    # Aggiunge edge artificiali
    for _, row in gdf.iterrows():
        geom = row.geometry
        coords = list(geom.coords)
        for i in range(len(coords) - 1):
            u_coord = coords[i]
            v_coord = coords[i + 1]

            for coord in [u_coord, v_coord]:
                if coord not in coord_to_node:
                    coord_to_node[coord] = node_id
                    x, y = coord
                    G.add_node(node_id, x=x, y=y, **attr)
                    node_id += 1

            u = coord_to_node[u_coord]
            v = coord_to_node[v_coord]

            length = gpd.GeoDataFrame([LineString([u_coord, v_coord])], columns=["geometry"], crs=CRS_GRAD).to_crs(CRS_METR)["geometry"][0].length
            weight_multipler = row["weight_multipler"]
            column = row.drop(["geometry", "weight_multipler", "weight", "length"]).to_dict()
            G.add_edge(u, v, length=length, weight_multipler = weight_multipler, weight = length * weight_multipler, geometry=LineString([u_coord, v_coord]), **column, **attr)
    
    # Aggiunge crs
    G.graph['crs'] = gdf.crs
    return G

def add_edge_near_nodes(G:nx.MultiDiGraph,
                        distance:int = 5,
                        weight_moltiplicator:float = 1.0,
                        tipo_da_connettere: str = None,
                        tipo_su_cui_connettere: str = None,
                        **attr) -> nx.MultiDiGraph:
    """
    Aggiunge Archi ad un MultiDiGraph a partire dalla Geometry di tipo LineString di un GeoDataFrame.
    distance è la distanza in metri dei nodi. Se non specificato è 5 metri
    """
    gdf_nodes = ox.graph_to_gdfs(G, nodes=True, edges=False)   # Per farlo ci serve ottenere una lista

    gdf_nodes_METR = gdf_nodes.to_crs(CRS_METR)
    gdf_nodes_GRAD = gdf_nodes.to_crs(CRS_GRAD)

    coords = np.array([[geom.x, geom.y] for geom in gdf_nodes_METR.geometry])
    node_ids = gdf_nodes.index.to_list()

    tree = cKDTree(coords)
    # .query_pairs -> Trova le coppie di nodi vicine alla distanza data
    pairs = tree.query_pairs(r=distance)  # distanza massima in metri

    for i, j in pairs:
        node_u = node_ids[i]
        node_v = node_ids[j]

        tipo_u = G.nodes[node_u].get("tipo")
        tipo_v = G.nodes[node_v].get("tipo")

        # Connetti solo se i tipi sono diversi nel modo richiesto
        if tipo_da_connettere is not None and tipo_su_cui_connettere is not None:
            if (not ((tipo_u == tipo_da_connettere and tipo_v == tipo_su_cui_connettere) or
                    (tipo_v == tipo_da_connettere and tipo_u == tipo_su_cui_connettere))):
                continue


        if not G.has_edge(node_u, node_v):
            # Calcoliamo la lunghezza in metri
            point_u_metr = gdf_nodes_METR.iloc[i].geometry
            point_v_metr = gdf_nodes_METR.iloc[j].geometry
            length = point_u_metr.distance(point_v_metr)

            # Prendiamo le coordinate lat/long
            point_u_grad = gdf_nodes_GRAD.iloc[i].geometry
            point_v_grad = gdf_nodes_GRAD.iloc[j].geometry
            geometry = LineString([point_u_grad, point_v_grad])

            G.add_edge(node_u, node_v, geometry=geometry, length=length,
                       weight_moltiplicator = weight_moltiplicator,
                       weight = length * weight_moltiplicator, **attr)
    return G

def connect_poi_nodes_to_graph(G: nx.MultiDiGraph, poi_gdf: gpd.GeoDataFrame, **attr) -> nx.MultiDiGraph:
    """
    Connette al grafo un GeoDataFrame, creando come nodi ogni Point presente nelle geometrie.
    Accetta sia geometria di tipo: Point che di tipo: Polygon e MultiPolygon. Su quest'ultimi verrà
    presa in considerazione il loro "centroid".
    """
    # Convertiamo eventuali MultiPolygon
    poi_gdf = _explode_gdf_with_MultiPolygon_to_Polygon(poi_gdf)

    # Raccogliamo i Point
    poi_nodes = _get_points_of_gdf(poi_gdf)

    # Estraiamo i nodi dal grafo, prendiamo le coordinate in metri e creiamo "tree" per poter fare le query
    gdf_nodes = ox.graph_to_gdfs(G, nodes=True, edges=False)
    coords = np.array([[geom.x, geom.y] for geom in gdf_nodes.to_crs(CRS_METR).geometry])
    tree = cKDTree(coords)

    # Teniamo anche gli indici del gdf dei vari nodi
    node_ids = gdf_nodes.index.to_list()

    # ID di partenza per i nuovi nodi
    poi_node_start_id = max(G.nodes) + 1  

    # Inseriamo i poi (nodi e archi)
    for poi_idx, (lon, lat) in poi_nodes.items():

        # Conversione dei punti in metri per eseguire la query
        point = gpd.GeoSeries([Point(lon, lat)], crs=poi_gdf.crs).to_crs(CRS_METR).iloc[0]
        poi_coord_metr = np.array([[point.x, point.y]])

        #.query -> trova i punti più vicini (k=1) alle coordinate date all'albero generato con i nodi del grafo.
        # Restituisce la distanza in metri (perché abbiamo dato tutto in metri) e l'INDICE del nodo
        dist, idx = tree.query(poi_coord_metr, k=1)

        # Aggiungiamo il nodo POI
        G.add_node(poi_node_start_id, x=lon, y=lat, **attr, **poi_gdf.iloc[poi_idx].drop("geometry").to_dict())

        # LineString geometria in EPSG:4326
        nearest_node = node_ids[idx[0]]
        geom_line = LineString([Point(lon, lat), gdf_nodes.loc[nearest_node].geometry])

        # Aggiungiamo arco di collegamento (può essere bidirezionale)
        G.add_edge(poi_node_start_id,
                nearest_node,
                length=dist[0],
                weight = 0,
                geometry=geom_line,
                **attr)
        G.add_edge(nearest_node,
                poi_node_start_id,
                length=dist[0],
                weight = 0,
                geometry=geom_line,
                **attr)

        poi_node_start_id += 1

    G = _change_edge_weight_in_gdf(G, poi_gdf, 0)

    return G

def connetti_due_grafi(G1:nx.MultiDiGraph, G2:nx.MultiDiGraph, weight_moltiplicator:float=0, **attr) -> nx.MultiDiGraph:
    # Connettiamo il 2 al 1

    gdf1_nodes = ox.graph_to_gdfs(G1, nodes=True, edges=False)   # Per farlo ci serve ottenere una lista
    idx1 = gdf1_nodes.index.to_list()
    gdf1_nodes_METR = gdf1_nodes.to_crs(CRS_METR)
    gdf1_nodes_GRAD = gdf1_nodes.to_crs(CRS_GRAD)

    gdf2_nodes = ox.graph_to_gdfs(G2, nodes=True, edges=False)   # Per farlo ci serve ottenere una lista
    idx2 = gdf2_nodes.index.to_list()
    gdf2_nodes_METR = gdf2_nodes.to_crs(CRS_METR)
    gdf2_nodes_GRAD = gdf2_nodes.to_crs(CRS_GRAD)

    m_coords_1 = np.array([[geom.x, geom.y] for geom in gdf1_nodes_METR.geometry])
    m_coords_2 = np.array([[geom.x, geom.y] for geom in gdf2_nodes_METR.geometry])

    tree_1 = cKDTree(m_coords_1)
    dists, idxs = tree_1.query(m_coords_2, k=1)

    G = nx.compose(G1, G2)
    for i, (idx1_match, dist) in enumerate(zip(idxs, dists)):
        if dist > 50:
            continue
        node_G2 = idx2[i]
        node_G1 = idx1[idx1_match]

        # Coordinate in lat/lon per costruire LineString
        point_G2 = gdf2_nodes_GRAD.loc[node_G2].geometry
        point_G1 = gdf1_nodes_GRAD.loc[node_G1].geometry
        geometry = LineString([point_G2, point_G1])

        # Lunghezza in metri
        point_G2_m = gdf2_nodes_METR.loc[node_G2].geometry
        point_G1_m = gdf1_nodes_METR.loc[node_G1].geometry
        length = point_G2_m.distance(point_G1_m)

        G.add_edge(node_G2, node_G1,
                   geometry=geometry,
                   length=length,
                   weight=length * weight_moltiplicator,
                   weight_moltiplicator=weight_moltiplicator,
                   **attr)
    return G

def etl_strade_ciclabili(custom_weights:dict = None) -> nx.MultiDiGraph:
    """
    Processo di ETL sul Grafo Raw scaricato da openstreet maps.
    Il proceso comprende:
    - Selezione attributi da mantenere: name, highway, lenght, maxspeed, tunnel, access, service, geometry
    - Aggiunta di attributi custom:
        - weight_multipler: Un moltiplicatore che verrà applicato a weigth in base a "highway" tramite
        il parametro "custom_weights", oppure da file se non specificato.
        N.B: di base cycleway dovrebbe assumere valore: 0, per dare massima priorità nel calcolo dei percorsi.
        - weight: Sarà di base: length * weight_multipler.
        - tipo: Un tag che viene assegnato a tutti i grafi del progetto. In questo caso assumerà valore "Strade_ciclabili".
        - artificial: Un tag che viene assegnato a tutit i grafi del progetto. Essendo questi archi reali: False
        - poi: Un tag che viene assegnato a tutti i grafi del progetto. Essendo queste le strade: False
    - Se Geometria non presente, aggiungiamo una LineString "dritta" tra i due nodi
    (può capitare che non ci sia tra nodi molto vicini).
    - Infine il grafo viene reso bidirezionale senza ridondanza di archi tra due punti.

    Args:
        custom_weights (dict): del tipo {"nome_strada": 1.3, "default": 1, ...} è possibile assegnare
            sempre un "default" che verrà assegnato in caso di mancanza di un nome_strada. Se non viene
            specificato un "default" verrà assegnato il valore "99", cioè massimo peso.
    """
    # Carichiamo grafo raw
    with open(PATH_STRADE_CICLABILI_PICKLE_RAW, "rb") as f:
        G_strade = pickle.load(f)

    # Se non specificato, carica da file i pesi delle highway
    if not custom_weights:
        with open(PATH_CUSTOM_WEIGHTS_STRADE_RAW, 'r') as f:
            custom_weights = json.load(f)

    # Iniziamo manipolazione dei "data" degli archi del grafo
    for u, v, k, data in G_strade.edges(data=True, keys=True):
        # Copio i data ed elimino tutto il dizionario.
        data = data.copy()
        G_strade[u][v][k].clear()
        # Scelgo gli attributi da mantenere
        G_strade[u][v][k]['length'] = data['length']
        G_strade[u][v][k]['name'] = data.get('name', None)
        G_strade[u][v][k]['highway'] = data["highway"]
        # ------------ Da capire se tenere ----------------
        G_strade[u][v][k]["maxspeed"] = data.get("maxspeed", None)
        G_strade[u][v][k]["tunnel"] = data.get("tunnel", None)
        G_strade[u][v][k]["access"] = data.get("access", None)
        G_strade[u][v][k]["service"] = data.get("service", None)
        # ------------------------------------------------
        # Calcolo pesi
        if isinstance(data["highway"], str):
            highway = [data["highway"]]
        max_weight = 99
        for h in highway:
            new_weight = custom_weights.get(h, custom_weights.get("default", 99))
            if new_weight < max_weight:
                max_weight = new_weight
        weight_multipler = max_weight
        # Aggiungo gli attributi custom
        if "cycleway" in data["highway"]:
            G_strade[u][v][k]["ciclabile"] = True
        else:
            G_strade[u][v][k]["ciclabile"] = False
        G_strade[u][v][k]["weight_multipler"] = weight_multipler
        G_strade[u][v][k]['weight'] = data['length'] * weight_multipler
        G_strade[u][v][k]["tipo"] = "Strade_ciclabili"
        G_strade[u][v][k]["artificial"] = False
        G_strade[u][v][k]["poi"] = False
        # Conversione geometry in Linestring
        G_strade[u][v][k]["geometry"] = data.get("geometry", LineString([(G_strade.nodes[u]["x"], G_strade.nodes[u]["y"]),
                                            (G_strade.nodes[v]["x"], G_strade.nodes[v]["y"])]))
    nx.set_node_attributes(G_strade, "Strade_ciclabili", "tipo")
    return G_strade.to_undirected()

# Funzione generazione percorsi automatica

def auto_analysis_poi(list_gdfs_poi:list[dict],
                      custom_weights:dict = None,
                      PATH_GEOJSON:str = None,
                      PATH_PICKLE:str = None):
    """
    Esegue in automatico: Creazione grafo combinato di strade e ciclabili con i pesi dati.
    Successivamente aggiunge tutti i "poi" dei gdf dati e crea il grafo dei percorsi minimi,
    salvando nei files. Leggere la descrizione degli attributi per maggiori info.
    Args:
        dict_gdfs_poi (dict): del tipo: [{"gdf": GeoDataFrame, "tipo": str, "attr": dict}, ...]
        custom_weights (dict): del tipo {"nome_strada": 1.3, "default": 1, ...} è possibile assegnare
            sempre un "default" che verrà assegnato in caso di mancanza di un nome_strada. Se non viene
            specificato un "default" verrà assegnato il valore "99", cioè massimo peso.
            Se non viene assegnato nessun custom_weights, verrà usato il grafo staging salvato su pickle
            (Operazione molto più veloce).
        PATH_GEOJSON (str): Il path dove verrà salvato il geoDataFrame convertito dal grafo.
        PATH_PICKLE (str): Il path dove verrà salvato il Grafo in formato pickle (se dovesse servire usarlo).
            E' opzionale se non si desidera salvare il pickle.
    """
    # Caricamento Grafo strade
    if custom_weights:   # Se necessario ricreare il grafo strade con i nuovi pesi
        G_strade = etl_strade_ciclabili(custom_weights)
    else:               # Se è possibile riutilizzare quello standard salvato su file
        with open(PATH_STRADE_CICLABILI_PICKLE_STAGING, "rb") as f:
            G_strade = pickle.load(f)

    # Aggiunta dei "poi" al grafo delle strade
    for dict_gdf in list_gdfs_poi:
        gdf = dict_gdf["gdf"]
        tipo = dict_gdf["tipo"]
        attr = dict_gdf.get("attr", None)
        kwargs = {"tipo": tipo, "poi": True, "artificial": True}
        if attr:
            kwargs.update(attr)
        G_compose = connect_poi_nodes_to_graph(G_strade, gdf, **kwargs)
    
    # uso di algoritmo per trovare i percorsi minimi tra i poi
    poi_nodes = [n for n, d in G_compose.nodes(data=True) if d.get("poi") == True]
    G_sport_tempo_libero = nx.MultiDiGraph(steiner_tree(G_compose, terminal_nodes=poi_nodes, weight='weight')).to_undirected()

    # salvataggio risultati
    if PATH_PICKLE: # Il grafo in pickle se richiesto
        with open(PATH_PICKLE, "wb") as f:
            pickle.dump(G_sport_tempo_libero, f)

    if PATH_GEOJSON: # In geodataframe in geojson se richiesto
        _custom_graph_to_gdf(G_sport_tempo_libero).to_file(PATH_GEOJSON, driver="GeoJSON")

def auto_filter_poi(G_percorsi:nx.MultiDiGraph, max_pois:int,
                    peso_inquinamento_vs_incidenti:float = 0.5,
                    peso_fattori_esterni_vs_distribuzione_poi:float = 0.5,
                    custom_weights:dict = None,
                      PATH_GEOJSON:str = None,
                      PATH_PICKLE:str = None):
    """
    Dato un grafo ottenuto da "auto_analysis_poi" (recuperare da file pickle in staging), crea
    un nuovo grafo (come auto_analysis_poi) filtratro con un massimo di pois dato in input.
    I pois verranno ridotti proporzionalmente per ogni municipio e mantenuti un minimo di 1 per "tipo".
    Verranno poi scelti i più vicini ad una ciclabile per ogni categoria.
    Scelti i punti verrà generato un nuovo grafo di percorsi utilizzando solo quei punti.
    Per maggiori info sui parametri, guardare auto_analysis_poi.
    Args:
        G_percorsi (nx.MultiDiGraph): Il grafo risultante da "auto_analysis_poi" che si vuole filtrare.
            Eseguire "auto_analysis_poi" specificando: tipo e priorità (valori crescenti hanno meno priorità).
        max_pois (int): Il numero massimo di pois in totale tra tutti i municipi di cui si vogliono ottenere
            i percorsi.
        peso_inquinamento_vs_incidenti (float): Un numero da 0 a 1 che andrà ad influenzare la quantità di "poi"
            che verranno assegnati per ogni Municipio.
            Si riferisce al peso dato all'inquinamento, rispetto a quello degli incidenti (quindi incidenti sarà
            calcolato da: 1 - inquinamento). Nell'insieme questi fattori saranno identificati come: "pesi fattori esterni".
        peso_fattori_esterni_vs_distribuzione_poi (float): Un numero da 0 a 1 che influenzare la quantità di "poi" che
            verranno assegnati per ogni Municipio.
            Questi fattori sono l'insieme di fattori esterni (incidenti e inquinamento) che andrà a contrapporsi
            con la normale distribuzione dei "poi" nel territorio (proporzionalmente al max pois).
        custom_weights (dict): Parametro di "auto_analysis_poi", dizionario dei pesi delle strade.
        PATH_GEOJSON (str): Parametro di "auto_analysis_poi", percorso in cui salvare il gdf.
        PATH_PICKLE (str): Parametro di "auto_analysis_poi", percorso in cui salvare il grafo.
    """
    # Otteniamo il numero di poi ridotto proporzionalmente per ogni gruppo in base
    # al massimo ricevuto come parametro, tramite la funzione: _get_proportional_reduce_n_pois
    colonne = {}
    for _, data in G_percorsi.nodes.items():
        colonne[data["tipo"]] = data.get("priorita", 99)
    colonne.pop("Strade_ciclabili")
    colonne_sorted = sorted(colonne, key=lambda x: colonne[x])
    gdf_n_pois = _get_proportional_reduce_n_pois(poi_group=colonne_sorted,
                                                 max_pois=max_pois,
                                                 peso_inquinamento_vs_incidenti=peso_inquinamento_vs_incidenti,
                                                 peso_fattori_esterni_vs_distribuzione_poi=peso_fattori_esterni_vs_distribuzione_poi)

    # Facciamo un municipio alla volta creando il sotto grafo
    selected_poi = []
    for idx, row in gdf_n_pois.iterrows():
        # Creiamo il sub grafo del municipio e prendiamo i nodi dei poi
        poligono_municipio:Polygon = row["geometry"]

        poi_nodes_municipio_categoria = defaultdict(list)
        nodes_in_municipio = []
        for n, data in G_percorsi.nodes(data=True):
            if poligono_municipio.contains(Point(data["x"], data["y"])):
                nodes_in_municipio.append(n)
                if data.get("poi", False):
                    poi_nodes_municipio_categoria[data["tipo"]].append(n)
        G_sub_grafo_municipio = G_percorsi.subgraph(nodes_in_municipio).copy()

        # Troviamo i nodi che hanno un arco "ciclabile" -> Quindi nodi di ciclabile
        ciclabile_nodes = set() 
        for u, v, k, data in G_sub_grafo_municipio.edges(keys=True, data=True):
            if data.get("ciclabile", False):
                ciclabile_nodes.add(u)
                ciclabile_nodes.add(v)
        ciclabile_nodes = list(ciclabile_nodes)
        distances = nx.multi_source_dijkstra_path_length(G_sub_grafo_municipio, sources=ciclabile_nodes, weight="length")

        for tipo, nodes in poi_nodes_municipio_categoria.items():
            nodes_sorted = sorted(nodes, key=lambda n: distances.get(n, float("inf")))
            max_num = int(gdf_n_pois.loc[idx, tipo])
            selected_poi.extend(nodes_sorted[:max_num])

    new_gdf_filtrato = defaultdict(list)
    for node in selected_poi:
        row = G_percorsi.nodes[node].copy()
        priorita = G_percorsi.nodes[node]["priorita"]
        tipo = G_percorsi.nodes[node]["tipo"]
        geometria = Point(G_percorsi.nodes[node]["x"], G_percorsi.nodes[node]["y"])
        row.pop("priorita")
        row.pop("tipo")
        row.pop("x")
        row.pop("y")
        row.pop("poi")
        row.pop("artificial")
        row["geometry"] = geometria
        new_gdf_filtrato[(tipo, priorita)].append(row)

    final_gdf = []
    for (tipo, priorita), v in new_gdf_filtrato.items():
        final_gdf.append({"gdf": gpd.GeoDataFrame(v, crs=CRS_GRAD),
                         "tipo": tipo,
                         "priorita": priorita
                         })
    auto_analysis_poi(list_gdfs_poi=final_gdf,
                      custom_weights=custom_weights,
                      PATH_GEOJSON=PATH_GEOJSON,
                      PATH_PICKLE=PATH_PICKLE)
