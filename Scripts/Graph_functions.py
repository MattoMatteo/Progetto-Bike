import geopandas as gpd
from shapely import Point, LineString, MultiPolygon, Polygon, GeometryCollection
from shapely.prepared import prep
import networkx as nx
import osmnx as ox
from scipy.spatial import cKDTree
import numpy as np

from my_paths import CRS_GRAD, CRS_METR

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
    for u, v, k, data in G.edges(keys=True, data=True):
        geom = data.get("geometry", None)
        if geom is None:
            continue

        if geometry_collection.intersects(geom):
            data["weight"] = 0

    return G

# Funzioni Pubbliche

def calcola_moltiplicatore_peso_strade(highway: str = None) -> float:
    return 2 # al momento

def geojson_to_graph(gdf: gpd.GeoDataFrame, weight_moltiplicator:float = 1.0, **attr) -> nx.MultiDiGraph:
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
                    G.add_node(node_id, x=x, y=y)
                    node_id += 1

            u = coord_to_node[u_coord]
            v = coord_to_node[v_coord]

            length = gpd.GeoDataFrame([LineString([u_coord, v_coord])], columns=["geometry"], crs=CRS_GRAD).to_crs(CRS_METR)["geometry"][0].length
            column = row.drop("geometry").to_dict()
            G.add_edge(u, v, distance=length, weight = length * weight_moltiplicator, geometry=LineString([u_coord, v_coord]), **column, **attr)
    
    # Aggiunge crs
    G.graph['crs'] = gdf.crs
    return G

def add_edge_near_nodes(G:nx.MultiDiGraph, distance:int = 5, weight_moltiplicator:float = 1.0, **attr) -> nx.MultiDiGraph:
    """
    Aggiunge Archi ad un MultiDiGraph a partire dalla Geometry di tipo LineString di un GeoDataFrame.
    distance è la distanza in metri dei nodi. Se non specificato è 5 metri
    """
    gdf_nodes, _ = ox.graph_to_gdfs(G)   # Per farlo ci serve ottenere una lista

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

        if not G.has_edge(node_u, node_v):
            # Calcoliamo la lunghezza in metri
            point_u_metr = gdf_nodes_METR.iloc[i].geometry
            point_v_metr = gdf_nodes_METR.iloc[j].geometry
            length = point_u_metr.distance(point_v_metr)

            # Prendiamo le coordinate lat/long
            point_u_grad = gdf_nodes_GRAD.iloc[i].geometry
            point_v_grad = gdf_nodes_GRAD.iloc[j].geometry
            geometry = LineString([point_u_grad, point_v_grad])

            G.add_edge(node_u, node_v, geometry=geometry, distance=length, weight = length * weight_moltiplicator, **attr)
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
                distance=dist[0],
                weight = 0,
                geometry=geom_line,
                **attr)
        G.add_edge(nearest_node,
                poi_node_start_id,
                distance=dist[0],
                weight = 0,
                geometry=geom_line,
                **attr)

        poi_node_start_id += 1

    G = _change_edge_weight_in_gdf(G, poi_gdf, 0)

    return G
