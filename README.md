# otbasscap-data

## Summary HTML docs

* Initial summary of existing paradigm [link](https://tbep-tech.github.io/otbasscap-data/eval_paradigm.html)

* Assessment of attainment with hypothetical chlorophyll upper limits [link](https://tbep-tech.github.io/otbasscap-data/upperchla.html)

* GLM of monthly chlorophyll attainment for sug-segments [link](https://tbep-tech.github.io/otbasscap-data/chlmoattain.html)

* Explanation of monthly load data assignment to OTB subsegments [link](https://tbep-tech.github.io/otbasscap-data/otb_subbasin_loads.html)

* Load increase projection [link](https://tbep-tech.github.io/otbasscap-data/load_increase.html)

----------------------

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

Water quality data from Hillsborough County Environmental Protection Commission (https://epcbocc.sharepoint.com/sites/Share/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FShare%2FShared%20Documents%2FOutbound%2FERM%2FWQM%5FReports&p=true&ga=1). Raw data files stored on Google Drive were processed by ./R/dat_proc.R.

## ./data/sg.RData
Seagrass transect data for Tampa Bay, obtained from the Tampa Bay Estuary Program via the `tbeptools` API:

  * `sg_dat`: transect data (returned by `tbeptools::read_transect()`)
  * `sg_lines`: transect locations (returned by `tbeptools::trnlns`)
  * `sg_pts`: transect points (returned by `tbeptools::trnpts`)

## ./data-raw/DataDownload_2999342_row.csv

Water quality data for Old Tampa Bay monitoring strata, obtained from Pinellas County Water Atlas (https://pinellas.wateratlas.usf.edu/datadownload).

Search query:
  * Data Type: Surface water quality
  * Water Atlas: Pinellas County Water Atlas
  * Date Range: 2000-01-01 to 2024-04-08
  * Parameter: BOD, Biochemical oxygen demand; Chlorophyll (a+b+c); Chlorophyll a (probe relative fluorescence); Chlorophyll a, uncorrected for pheophytin; Chlorophyll a, corrected for pheophytin; Chlorophyll b; Chlorophyll c; Chlorophyll/Pheophytin ratio; Dissolved oxygen (DO); Dissolved oxygen saturation; Nitrogen, inorganic; Nitrogen, ammonia as N; Nitrogen, ammonia as N; Nitrogen, ammonia as N; Nitrogen, ammonia (NH3) as NH3; Nitrogen, ammonia (NH3) as NH3; Nitrogen, ammonia (NH3) + ammonium (NH4); Nitrogen, ammonium (NH4) as N; Nitrogen, ammonia (NH3) + ammonium (NH4); Nitrogen ion (N); Nitrogen, Nitrite (NO2) as N; Nitrogen, Nitrite (NO2) as N; Nitrogen, Nitrite (NO2) as NO2; Nitrogen, Nitrite (NO2) as NO2; Nitrogen, Nitrite (NO2) as N; Nitrogen, Nitrate (NO3) as N; Nitrogen, Nitrate (NO3) as N; Nitrogen, Nitrate (NO3) as N; Nitrogen, Nitrate (NO3) as NO3; Nitrogen, Nitrate (NO3) as NO3; Nitrogen, Nitrate (NO3) as N; Nitrogen, organic; Nitrogen, organic; Nitrogen, Nitrite (NO2) + Nitrate (NO3) as N; Nitrogen, Nitrite (NO2) + Nitrate (NO3) as N; Salinity; Salinity; Temperature, water; Temperature, water; Nitrogen, Kjeldahl; Nitrogen, Kjeldahl; Nitrogen, mixed forms (NH3)+(NH4)+organic+(NO2)+(NO3); Nitrogen; Nitrogen; Total Suspended Solids (TSS); Turbidity; Turbidity; Turbidity; Velocity - stream.

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

## ./data-raw/fwcpyro20112023.csv

Monthly *Pyrodinium bahamense* sampling in Old Tampa Bay, 2011-2023. From https://myfwc.com/research/redtide/monitoring/database/, C. Lopez.

Columns: 
  * `yr`: Year, 2011-2023
  * `date`: Date of sampling, YYYY-MM-DD
  * `Latitude`: Latitude of sampling location, CRS 4326
  * `Longitude`: Longitude of sampling location, CRS 4326
  * `pyro`: Pyrodinium cell count, cells/L, `NA` is sampled but not present/detected
  * `doy`: Day of year, 1-365
  * `pyrocat`: Bloom intensity category from cells/L, <= 1e4 No bloom, > 1e4 & <= 1e5 Low, > 1e5 & <= 1e6 Medium, > 1e6 High

## ./data/otbsub.RData

Approximate boundaries for six sub-segments in OTB. Attribute column `subseg` defines the sub-segment name in cardinal directions, eg., NW: northwest, etc. CRS WGS 84.

## ./R/dat_proc.R

Queries raw EPC and WIN data from Google Drive and exports RData to ./data/