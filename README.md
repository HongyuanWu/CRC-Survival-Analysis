# CRC-Survival-Analysis
Survival Analysis on CRC samples subsetted by different oncology progression.

RICHIEDE : 
- Ruby
- R
- PiCnIc (o dati ottenuti da esso)

TESTATO SU : Rstudio 1.0.143 con R 3.4 / R 3.4.1 e Ruby 2.3

CARTELLE E FILE :
- TCGA-DATA : dati precedentemente utilizzati, ottenuti e restituiti da PiCnIc-COADREAD
- PICNICGRAPH : cartella d'esportazione dei grafi ottenuti da PiCnIc in formato xml
- CSVGRAPH : cartella di salvataggio dei dataset in csv ottenuti dalla versione xml dei grafi
	(ottenuti tramite uno script in Ruby: XmlGraphParser.rb) 
- CSVLIFETABLES : cartella per il salvataggio delle life tables create dal programma
- SCRIPTS : cartella contenente gli script utilizzati dal programma

CODICE:
- SurvivalAnalysisCRC.Rproj
- SurvivalAnalysisMain.r = script principale, carica librerie e altri scripts
Nella cartella Scripts : 
- GettingData : improta i dati da TCGA-DATA
- CreateSurvivalObject : sfrutta la libreria survival per creare gli oggetti di sopravvivenza
- XmlGraphParser.r e XmlGraphParser.rb : il primo script invoca il secondo che ritorna nodi e
	archi dei grafi ottenuti da PiCnIc come dataset per R.
	Successivamente scorre i grafi per ottenere dei raggruppamenti.
- SubsettingData.r : esclude determianti pazienti dal dataset iniziale e crea dei subsets
	in accordo con i gruppi ottenuti dai grafi dallo script precedente
- MSI_StatusDrivenSurvivalAnalysis.r : utilizza il campione intero per eseguire una 
	analisi di sopravvivenza generica, guidata da MSI_status (MSS or MSI-HIGH)
	E' inserita più a livello di esempio per illustrare l'utilizzo di molte funzioni
	relative all'analisi di sopravvivenza.
- MutationsDrivenSurvivalAnalysis.r : esegue un'analisi di sopravvivenza sui dataset 
	dei precedenti raggruppamenti guidata dalla presenza di determinate mutazioni o
	dalla presenza del paziente in determinati cammini/gruppi/livelli ottenuti dai grafi.
- functions.r : contiene le principali funzioni utilizzate dal resto degli scripts.

ESECUZIONE : 
Aprire 'SurvivalAnalysisCRC.Rproj' e eseguire 'SurvivalAnalysisMain.R'.

RISULTATI : 
(OSS: per ogni cox model usare survfit(coxmodel) per ottenere la sopravvivenza e poterla 
      visualizzare graficamente)
- Life Tables in 'CSVLifeTables' relative a:
	- campione totale
	- risultati per kaplan-meier sulla parte di campione con 'MSI-status = MSS' raggruppati
	  secondo 'levelGrouping' (guardare 'SurvivalAnalysisForTCGA-COADREAD.pdf')
	- risultati per kaplan-meier sulla parte di campione con 'MSI-status = MSS' raggruppati
	  secondo 'pathLevelGrouping' (guardare 'SurvivalAnalysisForTCGA-COADREAD.pdf')

- Risultati in R sull'intero campione basandosi su 'MSI-status': 
	- km_tot : kaplan meier sull'intero campione
	- km_mss_msi : kaplan meier sul campione diviso per MSI_status
	- log_rank : log rank test su km_mss_msi (MSS e MSI-HIGH le due parti confrontate)
	- hazard_tot : hazard function per tutto il campione
	- hazard_mss : hazard function per quella parte di campione con MSI_status = MSS
	- hazard_msi : hazard function per quella parte di campione con MSI_status = MSI-HIGH
	- cumulative_hazard_tot : cumulative hazard function per tutto il campione
	- cox_mss_msi : modello di cox con MSI_status come covariata 
	- baseline_cumulative_hazard : baseline cumulative hazard function per cox_mss_msi
	- strata_cox_mss_msi : modello di cox stratificato per MSI_status
	- baseline_cumulative_hazard_strata : baseline cumulative hazard function per 
					      il modello strata_cox_mss_msi
	- baseline_cumulative_hazard_mss : baseline cumulative hazard function per quella parte
					   di campione con MSI_status = MSS del modello
				           strata_cox_mss_msi
	- baseline_cumulative_hazard_msi : baseline cumulative hazard function per quella parte
					   di campione con MSI_status = MSI-HIGH del modello
				           strata_cox_mss_msi
	- haz_mss : hazard function per quella parte del campione con MSI_status = MSS
	- haz_msi : hazard function per quella parte del campione con MSI_status = MSI-HIGH
	- strata_cox_mss_msi_fake : modello di cox stratificato per MSI_status con due covariate
				    casuali generate a livello di esempio
	- res_surv : funzione dei residui da confrontare con la retta y=x per controllare
		     l'overall fit del modello strata_cox_mss_msi_fake
	- ssr_check_ph : Scaled Schoenfeld Residuals per controllare la proporzionalità
			 delle covariate di strata_cox_mss_msi_fake

- Risultati in R sui raggruppamenti ottenuti dal grafo
	- MSSnodes : lista di nodi per il grafo di MSS ottenuto da PiCnIc
	- MSInodes : lista di nodi per il grafo di MSI ottenuto da PiCnIc
	- MSSedges : lista di archi per il grafo di MSS ottenuto da PiCnIc
	- MSIedges : lista di archi per il grafo di MSI-HIGH ottenuto da PiCnIc
	- MSSgroups_grouping1 : mutazioni disposte per gruppo secondo raggruppamento1 per MSS
				(guardare 'SurvivalAnalysisForTCGA-COADREAD.pdf')
		      		struttura : [[group1] ... [groupN]], accesso : MSSlevels[[i]] 
				con i = 1, ..., N
	- MSIgroups_grouping1 : mutazioni disposte per gruppo secondo raggruppamento1 per MSI-HIGH
				(guardare 'SurvivalAnalysisForTCGA-COADREAD.pdf')
		      		struttura : [[group1] ... [groupN]], accesso : MSSlevels[[i]] 
				con i = 1, ..., N
	- MSSlevels : mutazioni disposte per livello (guardare 'SurvivalAnalysisForTCGA-COADREAD.pdf')
		      struttura : [[level1] [level2] [level3]], accesso : MSSlevels[[i]] con i = 1,2,3
	- OSS: manca MSSpathLevels perchè il dataset sarà ottenuto da MSSlevels rimuovendo 
	  i pazienti che non seguono un cammino (guardare 'SurvivalAnalysisForTCGA-COADREAD.pdf');
          mentre le mutazioni per livello sarebbero le stesse di MSSlevels

- Risultati in R dell'analisi di sopravvivenza sui precedenti raggruppamenti
	- MSSgroup_cox : lista di liste relativa all'analisi di sopravvivenza tramite modello di Cox
			 per quella parte di campione con MSI_status = MSS per ogni gruppo in 
			 MSSgroups_grouping1 a cui sono precedentemente stati assegnati i relativi
			 pazienti.
			 struttura : [[group1] ... [groupN]], accesso: MSSgroups_cox[[i]] con i = 1..N
			 struttura interna --> per ogni gruppo saranno salvati diversi campi : 
				- MSSgroups_cox$models : modello di cox con le mutazioni come covariate
				- MSSgroups_cox$survival : sopravvivenza di quel gruppo
				- MSSgroups_cox$baseline_cumulative_hazard : baseline cumulative 
					hazard function per il modello di cox di quel gruppo
				- MSSgroups_cox$cox_snell_residuals_surv_fun : funzione dei residui per
					il modello di cox di quel gruppo da confrontare con y=x per 
					controllare l'overall fit del modello
				- MSSgroups_cox$scaled_schoenfeld_residuals : scaled schoenfeld residuals
					per controllare la proporzionalità di ogni covariata
	- MSIgroup_cox : lista di liste relativa all'analisi di sopravvivenza tramite modello di Cox
			 per quella parte di campione con MSI_status = MSI-HIGH per ogni gruppo in 
			 MSIgroups_grouping1 a cui sono precedentemente stati assegnati i relativi
			 pazienti.
			 struttura : [[group1] ... [groupN]], accesso: MSSgroups_cox[[i]] con i = 1..N
			 struttura interna --> per ogni gruppo saranno salvati diversi campi : 
				- MSSgroups_cox$models : modello di cox con le mutazioni come covariate
				- MSSgroups_cox$survival : sopravvivenza di quel gruppo
				- MSSgroups_cox$baseline_cumulative_hazard : baseline cumulative 
					hazard function per il modello di cox di quel gruppo
				- MSSgroups_cox$cox_snell_residuals_surv_fun : funzione dei residui per
					il modello di cox di quel gruppo da confrontare con y=x per 
					controllare l'overall fit del modello
				- MSSgroups_cox$scaled_schoenfeld_residuals : scaled schoenfeld residuals
					per controllare la proporzionalità di ogni covariata
	- MSScox_strata : modello di Cox stratificato per MSSgroups_grouping1 con tutte le mutazioni dei 
			  livelli come covariate
	- MSScox_strata2 : modello di Cox stratificato per MSSgroups_grouping1 senza covariate
	- cox_g2g4 : cox model stratificato per gruppo solo per MSSgroups_groupin1[[2]] o 
		     MSSgroups_groupin1[[4]]
	- log_rank_test_g2g4 : log rank test fra MSSgroup_cox[[2]] e MSSgroup_cox[[4]]
	- km_level1 : kaplan-meier per i pazienti nel 'level1' in 'levelGrouping', ovvero
		      per quelli in MSSlevel[[1]]
	- km_level2 : kaplan-meier per i pazienti nel 'level2' in 'levelGrouping', ovvero
		      per quelli in MSSlevel[[2]]
	- km_level3 : kaplan-meier per i pazienti nel 'level3' in 'levelGrouping', ovvero
		      per quelli in MSSlevel[[3]]
	- log_rank_levels : log rank test fra km_level1, km_level2 e km_level3
	- log_rank_levels_1_2 : log rank test fra km_level1 e km_level2
	- log_rank_levels_2_3 : log rank test fra km_level2 e km_level3
	- log_rank_levels_1_3 : log rank test fra km_level1 e km_level3
	- cox_levels : cox model per i pazienti in ogni MSSlevels con l'appartenenza 
		       ad un determinato livello come covariata
	- l_res_surv : funzione dei residui da confrontare con y=x per l'overall fit del modello
		       cox_levels
	- ssr_cox_levels : scaled schoenfeld residuals per cox_levels per controllare la
			   proporzionalità di ogni covariata nel modello
	- cox_levels_2 : cox model per i pazienti in ogni MSSlevels con i vari gruppi come 
			 diverse covariate
	- l2_res_surv : funzione dei residui da confrontare con y=x per l'overall fit del modello
			cox_levels_2
	- ssr_cox_levels_2 : scaled schoenfeld residuals per cox_levels_2 per controllare la 
			     proporzionalità di ogni covariata del modello
	- cox_levels_strata : cox model stratificato per livello per tutti i pazienti appartenenti
			      ad un MSSlevels
	- l_strata_res_surv : funzione dei residui da confrontare con y=x per l'overall fit del modello
			      cox_levels_strata
	- km_pathlevel1 : kaplan-meier per i pazienti nel 'pathLevel1' in 'pathLevelGrouping', ovvero
		          per quelli in MSSpathLevel[[1]]
	- km_pathlevel2 : kaplan-meier per i pazienti nel 'pathLevel2' in 'pathLevelGrouping', ovvero
		          per quelli in MSSpathLevel[[2]]
	- km_pathlevel3 : kaplan-meier per i pazienti nel 'pathLevel3' in 'pathLevelGrouping', ovvero
		          per quelli in MSSpathLevel[[3]]
	- log_rank_pathlevels : log rank test fra km_pathlevel1, km_pathlevel2 e km_pathlevel3
	- log_rank_pathlevels_1_2 : log rank test fra km_pathlevel1 e km_pathlevel2
	- log_rank_pathlevels_2_3 : log rank test fra km_pathlevel2 e km_pathlevel3
	- log_rank_pathlevels_1_3 : log rank test fra km_pathlevel1 e km_pathlevel3
	- cox_pathlevels : cox model per i pazienti in ogni MSSpathLevels con l'appartenenza 
		           ad un determinato livello come covariata
	- pl_res_surv : funzione dei residui da confrontare con y=x per l'overall fit del modello
		        cox_pathlevels
	- ssr_cox_pathlevels : scaled schoenfeld residuals per cox_pathlevels per controllare la 
			       proporzionalità di ogni covariata
	- cox_pathlevels_2 : cox model per i pazienti in ogni MSSpathLevels con i vari gruppi come 
			     diverse covariate
	- pl2_res_surv : funzione dei residui da confrontare con y=x per l'overall fit del modello
			 cox_pathlevels_2
	- ssr_cox_pathlevels_2 : scaled schoenfeld residuals per cox_pathlevels_2 per controllare la 
			         proporzionalità di ogni covariata del modello
	- cox_pathlevels_strata : cox model stratificato per livello per tutti i pazienti appartenenti
			      ad un MSSpathLevels
	- pl_strata_res_surv : funzione dei residui da confrontare con y=x per l'overall fit del modello
			       cox_pathlevels_strata
