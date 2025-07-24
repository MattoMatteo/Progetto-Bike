# Sezione addetti ai lavori

## Struttura repository
```
📦 Progetto  
├──── 📁 Data  
│      ├──── 📁 Clean  
│      │     ├──── 📁 Incidenti  
│      │     ├──── 📁 Info_bici  
│      │     ├──── 📁 Inquinamento  
│      │     └──── 📁 Punti_di_interesse  
│      │  
│      ├──── 📁 Raw  
│      │     ├──── 📁 Incidenti  
│      │     ├──── 📁 Info_bici  
│      │     ├──── 📁 Inquinamento  
│      │     └──── 📁 Punti_di_interesse  
│      │           ├──── 📁 Municipi  
│      │           └──── 📁 Parchi  
│      │  
│      └──── 📁 Staging  
│            └──── 📁 Inquinamento  
│  
└──── 📁 Scripts  
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