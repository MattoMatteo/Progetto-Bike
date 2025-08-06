##--------------- Info bici --------------------------------------------#

# RAW
PATH_CICLABILI_RAW = "../Data/Raw/Rete_ciclabile_stradale/bike_ciclabili_raw.geojson"
PATH_STRADE_CICLABILI_PICKLE_RAW = "../Data/Raw/Rete_ciclabile_stradale/strade_ciclabili.pickle"

# STAGING
        #Ciclabili
PATH_CICLABILI_PICKLE_STAGING = "../Data/Staging/Rete_ciclabile_stradale/ciclabile.pickle"
PATH_CICLABILI_GEOJSON_STAGING = "../Data/Staging/Rete_ciclabile_stradale/ciclabile.geojson"
        #Strade
PATH_STRADE_CICLABILI_PICKLE_STAGING = "../Data/Staging/Rete_ciclabile_stradale/strade_ciclabili.pickle"
PATH_STRADE_CICLABILI_GEOJSON_STAGING = "../Data/Staging/Rete_ciclabile_stradale/strade_ciclabili.geojson"

# CLEAN
PATH_RETE_CICLABILE_COMPLETA_PICKLE_CLEAN = "../Data/Clean/Rete_ciclabile_stradale/rete_ciclabile_completa.pickle"
PATH_RETE_CICLABILE_COMPLETA_GEOJSON_CLEAN = "../Data/Clean/Rete_ciclabile_stradale/rete_ciclabile_completa.geojson"

#------------------------------------------------------------------------

#----------------------------- Punti di interesse ------------------------#

#PATH RAW
PATH_BIKEMI_RAW = '../Data/Raw/Punti_di_interesse/bikemi_stazioni_raw.csv'
PATH_MUNICIPI_RAW = "../Data/Raw/Punti_di_interesse/Municipi/Municipi.shp"
PATH_PARCHI_RAW = "../Data/Raw/Punti_di_interesse/Parchi/Parchi_WGS84_Milano_1.shp"
PATH_FONTANE_RAW = "../Data/Raw/Punti_di_interesse/fontane_raw.csv"
PATH_BIBLIOTECHE_RAW = "../Data/Raw/Punti_di_interesse/biblioteche_raw.geojson"
PATH_IMPIANTI_SPORTIVI_RAW = "../Data/Raw/Punti_di_interesse/impianti_sportivi_raw.geojson"
PATH_FARMACIE_RAW = "../Data/Raw/Punti_di_interesse/farmacie_raw.geojson"
PATH_AREE_GIOCO_RAW = "../Data/Raw/Punti_di_interesse/aree_gioco_raw.geojson"
PATH_TEATRI_RAW = "../Data/Raw/Punti_di_interesse/teatri_raw.json"
PATH_CINEMA_RAW = "../Data/Raw/Punti_di_interesse/cinema_raw.json"
PATH_MUSEI_RAW = "../Data/Raw/Punti_di_interesse/musei_raw.json"
PATH_OSPEDALI_RAW = '../Data/Raw/Punti_di_interesse/ospedali_raw.geojson'
PATH_CASE_ACQUA_RAW = "../Data/Raw/Punti_di_interesse/case_acqua_raw.geojson"
#   Scuole
PATH_SCUOLE_PRIMARIE_RAW = "../Data/Raw/Istruzione/Scuole/scuole-primarie_raw.geojson"
PATH_SCUOLE_SECONDARIE_PRIMO_RAW = "../Data/Raw/Istruzione/Scuole/scuole-secondarie-primo-grado_raw.geojson"
PATH_SCUOLE_SECONDARIE_SECONDO_RAW = "../Data/Raw/Istruzione/Scuole/scuole-secondarie-secondo-grado_raw.geojson"
PATH_UNIVERSITA_RAW = "../Data/Raw/Istruzione/Scuole/universita_raw.geojson"

#PATH CLEAN
PATH_BIKEMI_CLEAN = "../Data/Clean/Punti_di_interesse/bikemi_stazioni_clean.geojson"
PATH_MUNICIPI_CLEAN = "../Data/Clean/Punti_di_interesse/municipi_clean.geojson"
PATH_PARCHI_CLEAN = "../Data/Clean/Punti_di_interesse/parchi_clean.geojson"
PATH_FONTANE_CLEAN = "../Data/Clean/Punti_di_interesse/fontane_clean.geojson"
PATH_BIBLIOTECHE_CLEAN = "../Data/Clean/Punti_di_interesse/biblioteche_clean.geojson"
PATH_IMPIANTI_SPORTIVI_CLEAN = "../Data/Clean/Punti_di_interesse/impianti_sportivi_clean.geojson"
PATH_SCUOLE_CLEAN = "../Data/Clean/Istruzione/Scuole/scuole_clean.geojson"
PATH_FARMACIE_CLEAN = "../Data/Clean/Punti_di_interesse/farmacie_clean.geojson"
PATH_AREE_GIOCO_CLEAN = "../Data/Clean/Punti_di_interesse/aree_gioco_clean.geojson"
PATH_TEATRI_CLEAN = "../Data/Clean/Punti_di_interesse/teatri_clean.geojson"
PATH_CINEMA_CLEAN = "../Data/Clean/Punti_di_interesse/cinema_milano.geojson"
PATH_MUSEI_CLEAN = "../Data/Clean/Punti_di_interesse/musei_milano.geojson"
PATH_OSPEDALI_CLEAN = '../Data/Clean/Punti_di_interesse/ospedali_clean.geojson'
PATH_CASE_ACQUA_CLEAN = "../Data/Clean/Punti_di_interesse/case_acqua_clean.geojson"

#-----------------------------------------------------------------------------------------------#

#--------------------------------------- Inquinamento ------------------------------------------#

# RAW
PATH_INQUINAMENTO_INGESTION_RAW = "../Data/Raw/Inquinamento/Airnet/inquinamento-airnet_media_raw.geojson"

# CLEAN
PATH_INQUINAMENTO_INGESTION_CLEAN = "../Data/Clean/Inquinamento/Airnet/inquinamento-airnet_ingestion_clean.geojson"

#-----------------------------------------------------------------------------------------------#

#-------------------------------------- Analisi ----------------------------------------------#
PATH_SPORT_TEMPO_LIBERO_ANALISI_PICKLE_STAGING = "../Data/Staging/Analisi/sport_tempo_libero_clean.pickle"
PATH_SPORT_TEMPO_LIBERO_ANALISI_CLEAN = "../Data/Clean/Analisi/sport_tempo_libero_clean.geojson"

#-----------------------------------------------------------------------------------------------#

#----------------------------- CRS (lo so non sono dei path) ----------------------------------------#


CRS_GRAD = "EPSG:4326"  # Lat-Long
CRS_METR = "EPSG:32632" # Metri