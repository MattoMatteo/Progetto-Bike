# Sezione addetti ai lavori

## Struttura repository
```
📦 Progetto
├── 📁 Data
│   ├── 📁 Clean
│   │   ├── 📁 Analisi
│   │   ├── 📁 Incidenti
│   │   ├── 📁 Inquinamento
│   │   ├── 📁 Punti_di_interesse
│   │   ├── 📁 Rete_ciclabile_stradale
│   │   └── 📁 Visualization
│   │
│   ├── 📁 Raw
│   │   ├── 📁 Incidenti
│   │   ├── 📁 Inquinamento
│   │   ├── 📁 Istruzione
│   │   ├── 📁 Punti_di_interesse
│   │   └── 📁 Rete_ciclabile_stradale
│   │
│   └── 📁 Staging
│       ├── 📁 Analisi
│       ├── 📁 Inquinamento
│       └── 📁 Rete_ciclabile_stradale
│
├── 📁 Scripts
│   ├── 📓 Analisi.ipynb
│   ├── 📓 Git.ipynb
│   ├── 🐍 Graph_functions.py
│   ├── 📓 Incidenti.ipynb
│   ├── 📓 Inquinamento.ipynb
│   ├── 🐍 my_paths.py
│   ├── 📓 Punti_di_interesse.ipynb
│   ├── 📓 Rete_ciclabile_stradale.ipynb
│   └── 📓 Visualization.ipynb
│
└── 📁 Tableau
    └── 📊 project.twb
```
### Consiglio di utilizzare questo schema per i nomi dei file:
**{**<span style="color:yellow">categoria</span>**}`_`{**<span style="color:yellow">nome_file</span>}`_`{<span style="color:yellow">cartella destinazione</span>}`.`{<span style="color:yellow">estensione</span>}  
esempio:    **inquinamento_stazioni_raw.json**

### Raw, Staging, Clean
**Raw**     →   File senza nessuna modifica. Scaricati così come sono oppure ottenuti da scraping (cambiamo solo il nome)  
**Staging** →   File Raw che hanno subito delle modifiche ETL ma che non sono ancora pronti (fase intermedia se serve)  
**Clean**   →   File pronti all'uso su kepler  

---

### Ho creato il file my_path.py
Inseriremo tutti i percorsi dei vari notebbok (raw, staging e clean) e importeremo direttamente il file  
in modo che questi percorsi risultino "centralizzati" all'interno del progetto e quindi si possano usare  
da più file contemporanemente (esempio i file di ingestion/etl e da quelli di analisi).  
Esempio di utilizzo:

```python

import geopandas as gpd

from my_paths import *

gpd.read_file(PATH_FILE_RAW)

```
"PATH_FILE_RAW" verrà automaticamente trovato all'interno del file my_paths e importato in quello attuale  
consentendo l'utilizzo di tutti gli oggetti presenti al suo interno.  
N.b. usando **"import *"** si sta dicendo di importare tutto il suo contenuto

## Sintesi librerie Grafi

### 📦 OSMnx – lavoro su grafi geografici da OpenStreetMap

```python
import osmnx as ox

ox.graph_from_place()	#Scarica un grafo (strade, piste, ecc.) da OSM per una città o zona.
ox.graph_from_polygon()	#Scarica il grafo OSM all'interno di un poligono (es. confine Milano).
ox.distance.nearest_nodes(G, lon, lat)	#Trova il nodo più vicino a una coordinata geografica.
ox.distance.shortest_path(G, source, target, weight='length')	#Trova il percorso più corto tra due nodi, pesato per lunghezza o altro.
ox.utils_graph.get_undirected()	#Converte un grafo MultiDiGraph in grafo non orientato. Utile per MST.
ox.plot_graph() / ox.plot_graph_route()	#Visualizza il grafo o un percorso sopra di esso.
```

### 📦 NetworkX – manipolazione di grafi generici

```python
import networkx as nx

nx.Graph() / nx.MultiDiGraph()	#Crea un grafo (semplice o con direzioni e archi multipli).
nx.compose(G1, G2)	#Unisce due grafi (nodi e archi) in uno solo.
nx.has_path(G, u, v)	#Verifica se esiste un percorso tra due nodi.
nx.shortest_path(G, source, target, weight='length')	#Trova il percorso più corto (come OSMnx).
nx.shortest_path_length()	#Dà la lunghezza del percorso più corto.
nx.minimum_spanning_tree(G)	#Restituisce l’albero di copertura minimo (MST).
nx.connected_components(G)	#Trova i componenti connessi in un grafo non orientato.
nx.degree(G, node)	#Dà il grado di un nodo (quanti archi ha).
nx.add_edge(u, v, attr_dict)	#Aggiunge un arco tra due nodi.
nx.add_node(node, attr_dict)	#Aggiunge un nodo al grafo.
```

### 📦 Scipy-spatial: cKDTree

```python
from scipy.spatial import cKDTree

cKDTree(data)	#Costruisce l’albero da una lista di punti [[x1, y1], [x2, y2], ...]
.query(point, k)	#Trova i k vicini più vicini di point
.query_ball_point(point, r)	#Trova tutti i punti entro un raggio r da point
.query_pairs(r)	#Trova tutte le coppie di punti entro distanza r (come nel tuo esempio)
```

### 📦 Scikit-learn – clustering

```python
from scipy.cluster import DBSCAN, KMeans

DBSCAN(eps, min_samples)	#Algoritmo per clusterizzare punti vicini. Non richiede sapere in anticipo il numero di cluster.
KMeans(n_clusters)	        #Clusterizza in n gruppi definiti.
.fit_predict(X)	            #Assegna un'etichetta (cluster) a ogni punto.
```

--- 

# 🚴‍♂️ Progetto Bike – Analisi per la Pianificazione di Piste Ciclabili a Milano

> Progetto del corso **Big Data Processing and Data Engineering**  
> Master in *Artificial Intelligence & Data Analytics for Business* – A.A. 2025/2026

---

## 📌 Obiettivo del Progetto

Il progetto si propone di supportare la pianificazione urbana nella città di **Milano**, con l’obiettivo di individuare **nuove aree strategiche per la creazione di piste ciclabili**.  
L'analisi si basa sull'integrazione di **dati geospaziali** provenienti da fonti eterogenee, al fine di migliorare l'infrastruttura ciclabile, promuovere la mobilità sostenibile e aumentare la sicurezza dei ciclisti.

---

## 🔍 Dataset Analizzati

Il gruppo sta raccogliendo, pulendo e analizzando una serie di dataset geolocalizzati, tra cui:

- 📍 **Stazioni di bike sharing** (es. BikeMi)
- 🌳 **Aree verdi e parchi**
- 🚰 **Fontanelle e punti d’acqua pubblici**
- 🏫 **Luoghi di interesse** (scuole, università, attrazioni culturali)
- 🛣️ **Infrastrutture esistenti** (piste ciclabili attuali, zone pedonali)
- 📊 Altri dati open source forniti dal Comune di Milano e portali pubblici

---

## ⚙️ Tecnologie Utilizzate

- **Python** (pandas, geopandas, shapely, etc.)
- **Jupyter Notebook** per l’esplorazione e la documentazione dell’analisi
- **Git/GitHub** per la collaborazione e il versionamento
- **Kepler** per la visualizzazione dei dati geospaziali

---

## 👥 Collaborazione

Questa repository è pensata per la **collaborazione tra i membri del gruppo**, con una suddivisione delle attività e delle analisi.  
Sono presenti notebook, script e materiali utili per lo sviluppo coordinato.

---

## Studenti coinvolti:
- Matteo Calautti
- Arianna Peralti 
- Lorenzo Dalla Rosa

---