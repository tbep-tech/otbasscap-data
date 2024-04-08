# otbasscap-data

## ./data/winwq.RData

Water quality data retrieved from FDEP WIN database (https://prodenv.dep.state.fl.us/DearWin/public/welcomeGeneralPublic). Raw data files stored on Google Drive were processed by ./R/dat_proc.R.

Search query:

  * Requested report - WIN WAVES
  * PRIMARY TYPE = Surface Water
  * ACTIVITY TYPE = Field, Sample
  * DEP ANALYTE GROUP = General Physical-Chemical, Oxygen Demand, Field Observation, Biological, Nutrients
  * COUNTY = HILLSBOROUGH, FLORIDA, PASCO, FLORIDA, PINELLAS, FLORIDA
  * Report Run on APRIL 5, 2024

## ./data/epcwq.RData

Water quality data from Hillsborough County Environmental Protection Commission(https://epcbocc.sharepoint.com/sites/Share/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FShare%2FShared%20Documents%2FOutbound%2FERM%2FWQM%5FReports&p=true&ga=1). Raw data files stored on Google Drive were processed by ./R/dat_proc.R.

## ./data/sg.RData
Seagrass transect data for Tampa Bay, obtained from the Tampa Bay Estuary Program via the `tbeptools` API:

  * `sg_dat`: transect data (returned by `tbeptools::read_transect()`)
  * `sg_lines`: transect locations (returned by `tbeptools::trnlns`)
  * `sg_pts`: transect points (returned by `tbeptools::trnpts`)

## ./data-raw/normalization.xlsx

Excel file visualizes hydrologic normalization of nitrogen loads to Old Tampa Bay, for Task 1 Draft Report (3/29/2024):

  * Sheet 1 demonstrates the current linear adjustment procedure and the alternative logistic procedure (Figure 10).
  * Sheet 2 visualizes the compliance load adjustment factor over time (Figure 9).

## ./data-raw/OTB_loads_seagrass.xlsx

Excel file visualizes nitrogen loads and provides proof-of-concept visualization for hysteresis analysis, for Task 1 Draft Report (3/29/2024).

  * Sheet 1 contains TBEP's monthly nitrogen load data.
  * Sheet 2 contains the OTB subset of data on Sheet 1, a pivot table of load types by year, and calculated percentages for each load type between 2012-2021 (Figure 3).
  * Sheet 3 contains TBEP's annual seagrass transect data and N load data, and shows hysteresis-related plots (Figure 7).

## ./data-raw/TB_hydro_monthly.csv

Tidy dataset for monthly hydrologic loads to Tampa Bay segments, 1985-2021. Retrieved from TBEP (https://tbep.org/tampa-bay-nitrogen-loads/).

Columns:
  * `year`: 1985-2021
  * `month`: 1-12
  * `bay_segment`: All Segments, Hillsborough Bay, Lower Tampa Bay, Middle Tampa Bay, Old Tampa Bay, Remainder Lower Tampa Bay.
  * `hy_load_106_m3_mo`: Hydrologic load (10e6 m3)

## ./data-raw/TB_loads_monthly.csv

Tidy dataset for monthly TN, TP, TSS & BOD load for Tampa Bay segments, 1995-2021. Retrieved from TBEP (https://tbep.org/tampa-bay-nitrogen-loads/).

Columns:
  * `year`: 1995-2021
  * `month`: 1-12
  * `bay_segment`: All Segments, Hillsborough Bay, Lower Tampa Bay, Middle Tampa Bay, Old Tampa Bay, Remainder Lower Tampa Bay.
  * `source`: AD (atmospheric deposition), DPS (domestic points source), GWS (groundwater), IPS (industrial point source), NPS (nonpoint source).
  * `tn_load`: Total nitrogen load (tons)
  * `tp_load`: Total phosphorus load (tons)
  * `tss_load`: Total suspended solids (tons)
  * `bod_load`: Biological oxygen demand (tons)

## ./R/dat_proc.R

Queries raw EPC and WIN data from Google Drive and exports RData to ./data/