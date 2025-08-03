import geopandas as gpd
from shapely import LineString
import networkx as nx
import osmnx as ox
from scipy.spatial import cKDTree
import numpy as np

from my_paths import CRS_GRAD, CRS_METR

def geojson_to_graph(gdf: gpd.GeoDataFrame) -> nx.MultiDiGraph:
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
            G.add_edge(u, v, weight=length, geometry=LineString([u_coord, v_coord]))
    
    # Aggiunge crs
    G.graph['crs'] = gdf.crs
    return G

def add_edge_near_nodes(G:nx.MultiDiGraph, distance:int = 5) -> nx.MultiDiGraph:
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

            G.add_edge(node_u, node_v, geometry=geometry, artificial=True, weight=length)
    return G
