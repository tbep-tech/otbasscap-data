# otbasscap-data

Raw data:

./data-raw/_WIN_WAVES_MEDINA_M_3_20240227101624_74965.txt
=========================================================
Raw nitrogen and phosphorus concentration data retrieved from FDEP WIN database (https://prodenv.dep.state.fl.us/DearWin/public/welcomeGeneralPublic).
Search query:
  COUNTY = HILLSBOROUGH, FLORIDA, PINELLAS, FLORIDA
  PRIMARY TYPE = Surface Water
  DEP ANALYTE GROUP = Nutrients
  Report Run on FEBRUARY 27, 2024

./data-raw/EPC_WQ.csv
=====================
Raw water quality data retrieved from Hillsborough County Environmental Protection Commmission (https://epcbocc.sharepoint.com/sites/Share/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FShare%2FShared%20Documents%2FOutbound%2FERM%2FWQM%5FReports&p=true&ga=1).

./data-raw/normalization.xlsx
=============================
Excel file visualizes hydrologic normalization of nitrogen loads to Old Tampa Bay, for Task 1 Draft Report (3/29/2024):
  Sheet 1 demonstrates the current linear adjustment procedure and the alternative logistic procedure (Figure 10).
  Sheet 2 visualizes the compliance load adjustment factor over time (Figure 9).

./data-raw/OTB_loads_seagrass.xlsx
==================================
Excel file visualizes nitrogen loads and provides proof-of-concept visualization for hysteresis analysis, for Task 1 Draft Report (3/29/2024).
- Sheet 1 contains TBEP's monthly nitrogen load data.
- Sheet 2 contains the OTB subset of data on Sheet 1, a pivot table of load types by year, and calculated percentages for each load type between 2012-2021 (Figure 3).
- Sheet 3 contains TBEP's annual seagrass transect data and N load data, and shows hysteresis-related plots (Figure 7).

./data-raw/TB_hydro_monthly.csv
===============================
Tidy dataset for monthly hydrologic loads to Tampa Bay segments, 1985-2021. Retrieved from TBEP (https://tbep.org/tampa-bay-nitrogen-loads/).
Columns:
  year: 1985-2021
  month: 1-12
  bay_segment: All Segments, Hillsborough Bay, Lower Tampa Bay, Middle Tampa Bay, Old Tampa Bay, Remainder Lower Tampa Bay.
  hy_load_106_m3_mo: Hydrologic load (10e6 m3)

./data-raw/TB_loads_monthly.csv
==============================
Tidy dataset for monthly TN, TP, TSS & BOD load for Tampa Bay segments, 1995-2021. Retrieved from TBEP (https://tbep.org/tampa-bay-nitrogen-loads/).
Columns:
  year: 1995-2021
  month: 1-12
  bay_segment: All Segments, Hillsborough Bay, Lower Tampa Bay, Middle Tampa Bay, Old Tampa Bay, Remainder Lower Tampa Bay.
  source: AD (atmospheric deposition), DPS (domestic points source), GWS (groundwater), IPS (industrial point source), NPS (nonpoint source).
  tn_load: Total nitrogen load (tons)
  tp_load: Total phosphorus load (tons)
  tss_load: Total suspended solids (tons)
  bod_load: Biological oxygen demand (tons)
