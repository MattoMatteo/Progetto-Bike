##--------------- Info bici --------------------------------------------#

# RAW
PATH_STRADE_CICLABILI_PICKLE_RAW = "../Data/Raw/Rete_ciclabile_stradale/strade_ciclabili_raw.pickle"

# STAGING
PATH_STRADE_CICLABILI_PICKLE_STAGING = "../Data/Staging/Rete_ciclabile_stradale/strade_ciclabili_staging.pickle"

# CLEAN
PATH_CUSTOM_WEIGHTS_STRADE_CLEAN = "../Data/Clean/Rete_ciclabile_stradale/custom_weight_strade_clean.json"
PATH_STRADE_CICLABILI_GEOJSON_CLEAN = "../Data/Clean/Rete_ciclabile_stradale/strade_ciclabili_clean.geojson"

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
PATH_SCUOLE_CLEAN = "../Data/Clean/Punti_di_interesse/scuole_clean.geojson"
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
PATH_INQUINAMENTO_RAW = "../Data/Raw/Inquinamento/inquinamento_raw.geojson"

# STAGING
PATH_INQUINAMENTO_STAGING = "../Data/Staging/Inquinamento/inquinamento_staging.geojson"

# CLEAN
PATH_INQUINAMENTO_CLEAN = "../Data/Clean/Inquinamento/inquinamento_clean.geojson"

#-------------------------------------- Incidenti -------------------------------------------------#
# RAW
PATH_INCIDENTI_RAW = "../Data/Raw/Incidenti/incidenti_mese_municipio.csv"

# CLEAN
PATH_INCIDENTI_CLEAN = "../Data/Clean/Incidenti/incidenti_clean.geojson"

#-------------------------------------- Analisi ----------------------------------------------#
PATH_SCORE_X_MUNICIPI = "../Data/Clean/Analisi/score_x_municipi_clean.geojson"

# STAING (pickle)
PATH_EXTENDED_CICLABILI_COMPLETO_PICKLE_STAGING = "../Data/Staging/Analisi/Extended_ciclabili/extended_ciclabili_completo_staging.pickle"

PATH_SPORT_TEMPO_LIBERO_COMPLETO_PICKLE_STAGING = "../Data/Staging/Analisi/Sport_tempo_libero/sport_tempo_libero_completo_staging.pickle"
PATH_SPORT_TEMPO_LIBERO_INQUINAMENTO_PICKLE_STAGING = "../Data/Staging/Analisi/Sport_tempo_libero/sport_tempo_libero_inquinamento_staging.pickle"
PATH_SPORT_TEMPO_LIBERO_INCIDENTI_PICKLE_STAGING = "../Data/Staging/Analisi/Sport_tempo_libero/sport_tempo_libero_incidenti_staging.pickle"

PATH_ISTRUZIONE_COMPLETO_PICKLE_STAGING = "../Data/Staging/Analisi/Istruzione/istruzione_completo_staging.pickle"
PATH_ISTRUZIONE_INQUINAMENTO_PICKLE_STAGING = "../Data/Staging/Analisi/Istruzione/istruzione_inquinamento_staging.pickle"
PATH_ISTRUZIONE_INCIDENTI_PICKLE_STAGING = "../Data/Staging/Analisi/Istruzione/istruzione_incidenti_staging.pickle"

PATH_CULTURA_SPETTACOLO_COMPLETO_PICKLE_STAGING = "../Data/Staging/Analisi/Cultura_spettacolo/cultura_spettacolo_completo_staging.pickle"
PATH_CULTURA_SPETTACOLO_INQUINAMENTO_PICKLE_STAGING = "../Data/Staging/Analisi/Cultura_spettacolo/cultura_spettacolo_inquinamento_staging.pickle"
PATH_CULTURA_SPETTACOLO_INCIDENTI_PICKLE_STAGING = "../Data/Staging/Analisi/Cultura_spettacolo/cultura_spettacolo_incidenti_staging.pickle"

PATH_SANITA_COMPLETO_PICKLE_STAGING = "../Data/Staging/Analisi/Sanita/sanita_completo_staging.pickle"
PATH_SANITA_INQUINAMENTO_PICKLE_STAGING = "../Data/Staging/Analisi/Sanita/sanita_inquinamento_staging.pickle"
PATH_SANITA_INCIDENTI_PICKLE_STAGING = "../Data/Staging/Analisi/Sanita/sanita_incidenti_staging.pickle"

# CLEAN "completi" (geojson)
PATH_SPORT_TEMPO_LIBERO_COMPLETO_CLEAN = "../Data/Clean/Analisi/Sport_tempo_libero/sport_tempo_libero_completo_clean.geojson"
PATH_EXTENDED_CICLABILI_COMPLETO_CLEAN = "../Data/Clean/Analisi/Extended_ciclabili/extended_ciclabili_completo_clean.geojson"
PATH_ISTRUZIONE_COMPLETO_CLEAN = "../Data/Clean/Analisi/Istruzione/istruzione_completo_clean.geojson"
PATH_CULTURA_SPETTACOLO_COMPLETO_CLEAN = "../Data/Clean/Analisi/Cultura_spettacolo/cultura_spettacolo_completo_clean.geojson"
PATH_SANITA_COMPLETO_CLEAN = "../Data/Clean/Analisi/Sanita/sanita_completo_clean.geojson"

# CLEAN "filtrati" (geojson)
PATH_EXTENDED_CICLABILI_INCIDENTI_CLEAN = "../Data/Clean/Analisi/Extended_ciclabili/extended_ciclabili_filtrato_incidenti_clean.geojson"
PATH_EXTENDED_CICLABILI_INQUINAMENTO_CLEAN = "../Data/Clean/Analisi/Extended_ciclabili/extended_ciclabili_filtrato_inquinamento_clean.geojson"

PATH_CULTURA_SPETTACOLO_INCIDENTI_CLEAN = "../Data/Clean/Analisi/Cultura_spettacolo/cultura_spettacolo_filtrato_incidenti_clean.geojson"
PATH_CULTURA_SPETTACOLO_INQUINAMENTO_CLEAN = "../Data/Clean/Analisi/Cultura_spettacolo/cultura_spettacolo_filtrato_inquinamento_clean.geojson"

PATH_ISTRUZIONE_INCIDENTI_CLEAN = "../Data/Clean/Analisi/Istruzione/istruzione_filtrato_incidenti_clean.geojson"
PATH_ISTRUZIONE_INQUINAMENTO_CLEAN = "../Data/Clean/Analisi/Istruzione/istruzione_filtrato_inquinamento_clean.geojson"

PATH_SANITA_INCIDENTI_CLEAN = "../Data/Clean/Analisi/Sanita/sanita_filtrato_incidenti_clean.geojson"
PATH_SANITA_INQUINAMENTO_CLEAN = "../Data/Clean/Analisi/Sanita/sanita_filtrato_inquinamento_clean.geojson"

PATH_SPORT_TEMPO_LIBERO_INCIDENTI_CLEAN = "../Data/Clean/Analisi/Sport_tempo_libero/sport_tempo_libero_filtrato_incidenti_clean.geojson"
PATH_SPORT_TEMPO_LIBERO_INQUINAMENTO_CLEAN = "../Data/Clean/Analisi/Sport_tempo_libero/sport_tempo_libero_filtrato_inquinamento_clean.geojson"

# CLEAN "uniti"
PATH_EXTENDED_CICLABILI_UNITI_CLEAN = "../Data/Clean/Analisi/Extended_ciclabili/extended_ciclabili_unito_clean.geojson"

PATH_CULTURA_SPETTACOLO_UNITI_CLEAN = "../Data/Clean/Analisi/Cultura_spettacolo/cultura_spettacolo_unito_clean.geojson"
PATH_CULTURA_SPETTACOLO_POI_UNITI_CLEAN = "../Data/Clean/Analisi/Cultura_spettacolo/cultura_spettacolo_POI_unito_clean.geojson"

PATH_ISTRUZIONE_UNITI_CLEAN = "../Data/Clean/Analisi/Istruzione/istruzione_unito_clean.geojson"
PATH_ISTRUZIONE_POI_UNITI_CLEAN = "../Data/Clean/Analisi/Istruzione/istruzione_POI_unito_clean.geojson"

PATH_SANITA_UNITI_CLEAN = "../Data/Clean/Analisi/Sanita/sanita_unito_clean.geojson"
PATH_SANITA_POI_UNITI_CLEAN = "../Data/Clean/Analisi/Sanita/sanit_POI_unito_clean.geojson"

PATH_SPORT_TEMPO_LIBERO_UNITI_CLEAN = "../Data/Clean/Analisi/Sport_tempo_libero/sport_tempo_libero_unito_clean.geojson"
PATH_SPORT_TEMPO_LIBERO_POI_UNITI_CLEAN = "../Data/Clean/Analisi/Sport_tempo_libero/sport_tempo_libero_POI_unito_clean.geojson"
#-----------------------------------------------------------------------------------------------#

#----------------------------- Data Visualization ----------------------------------------#

PATH_GRUPPI_ANALISI_UNITI_CLEAN = "../Data/Clean/Visualization/gruppi_analisi_uniti_clean.geojson"
PATH_POI_GRUPPI_ANALISI_UNITI_CLEAN = "../Data/Clean/Visualization/poi_gruppi_analisi_uniti_clean.geojson"

#-----------------------------------------------------------------------------------------------#

#----------------------------- CRS (lo so non sono dei path) ----------------------------------------#
CRS_GRAD = "EPSG:4326"  # Lat-Long
CRS_METR = "EPSG:32632" # Metri