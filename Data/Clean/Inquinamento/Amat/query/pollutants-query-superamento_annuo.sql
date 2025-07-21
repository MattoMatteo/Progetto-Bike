SELECT strftime('%Y', data) AS anno, stazione.nome, inquinante.simbolo, AVG(misurazione.valore) as Media_annua,
		soglia.valore, soglia.periodo_di_riferimento, soglia.max_superamenti_anno
FROM misurazione
JOIN stazione ON stazione.id_stazione == misurazione.id_stazione
JOIN bollettino ON bollettino.id_bollettino == misurazione.id_bollettino
JOIN inquinante ON inquinante.id_inquinante == misurazione.id_inquinante
JOIN soglia ON soglia.id_inquinante == inquinante.id_inquinante
WHERE soglia.periodo_di_riferimento == "Annuale"
GROUP BY anno, stazione.nome, inquinante.simbolo