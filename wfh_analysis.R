library(readr)
library(tidyr)
library(dplyr)
library(zoo)
library(tidycensus)
library(stringr)
library(ipumsr)
library(httr)
### FUNCTIONS ###
get_bls_data <- function(url, email) {
  bls_res <- GET(url = url, user_agent(email))
  stop_for_status(bls_res)
  
  bls_content <- content(bls_res, 
                         as = "parsed",
                         type = "text/tab-separated-values",
                         encoding = "UTF-8",
                         col_names = T,
                         col_types = cols(.default = col_character()),
                         trim_ws = T
  )
  return(bls_content)
  
}

### ANALYSIS ###

kastle_back_to_work_barometer <- read_csv(
  "./data/raw_data/kastle_20_24_long.csv",
  col_names = T,
  col_types = "Dcd"
)

# Looking office occupancy rates in 10 metro areas tracked in Kastle data
# via office keycard swipes. Taking 4-week moving average to smooth out
# weekly data

kastle_office_occupancy <- kastle_back_to_work_barometer %>% 
  arrange(desc(data_week)) %>% 
  pivot_wider(names_from = city, values_from = occupancy) %>% 
  mutate(across(-data_week, ~rollmean(.x, k = 4, fill = NA, align = "left"))) %>% 
  rowwise() %>% 
  mutate(`10 Metro Avg.` = mean(c_across(Washington:Austin))) %>% 
  ungroup() %>% 
  filter(data_week >= base::as.Date("2020-04-15"))

new_order <- order(kastle_office_occupancy[1,2:12], decreasing = T)

kastle_office_occupancy_cities <- kastle_office_occupancy[,2:12]

kastle_reordered <- kastle_office_occupancy_cities[, new_order]

kastle_office_occupancy_ordered <- bind_cols(
  tibble(data_week = kastle_office_occupancy$data_week), 
  kastle_reordered)


write_csv(kastle_office_occupancy, 
          "./data/clean_data/kastle_office_occ_roll_four_dw_line.csv")


# Getting WFH rates from 2022 1-year ACS for the same 10 metro areas as the Kastle data:
# Reading in Census api key
con <- file(description = "~/.api_keys/census.txt", open = "rt", blocking = F)
CENSUS_API_KEY <- readLines(con, n = 1)
close(con)

msa_wfh_22 <- get_acs(
  geography = "cbsa",
  variables = c("wfh" = "B08006_017"),
  year = 2022,
  survey = "acs1",
  geometry = F,
  summary_var = "B08006_001", # My population here is _workers_ 16 years and older.
  key = CENSUS_API_KEY
)

# Now the national numbers
us_wfh_22 <- get_acs(
  geography = "us",
  variables = c("wfh" = "B08006_017"),
  year = 2022,
  survey = "acs1",
  geometry = F,
  summary_var = "B08006_001", # My population here is _workers_ 16 years and older.
  key = CENSUS_API_KEY
)

# Writing out data to ./raw_data/ folder
bind_rows(msa_wfh_22, us_wfh_22) %>% 
  write_csv("./data/raw_data/msa_us_wfh_22_acs.csv")

# Filter to 10 Kastle metros & calculate WFH percentage:
# FIPS code reference: https://www2.census.gov/programs-surveys/cps/methodology/2015%20Geography%20Cover.pdf
kastle_fips <- c("12420",
                 "16980",
                 "19100",
                 "26420",
                 "31080",
                 "35620",
                 "37980",
                 "41860",
                 "41940",
                 "47900",
                 "1")

msa_wfh_22_prop <- msa_wfh_22 %>% 
  bind_rows(us_wfh_22) %>% 
  filter(GEOID %in% kastle_fips) %>% 
  mutate(wfh_prop = (estimate / summary_est) * 100, # getting derived margins of error
         wfh_prop_moe = moe_prop(estimate, summary_est, moe, summary_moe) * 100,
         wfh_90_high = wfh_prop + wfh_prop_moe, # Calculating these for the Datawrapper: https://blog.datawrapper.de/confidence-intervals-value-markers-bar-charts/
         wfh_90_low = wfh_prop - wfh_prop_moe,
         msa_short_name = str_remove(str_remove(NAME, "-.*$"), ",.*$")) %>% 
  select(GEOID, NAME, msa_short_name, wfh_prop, wfh_90_low, wfh_90_high) %>% 
  arrange(desc(wfh_prop))

write_csv(
  msa_wfh_22_prop,
  "./data/clean_data/census_acs_22_1_wfh_dw_bar.csv"
)

# Getting IPUMS CPS microdata for telework estimates:
con <- file(description = "~/.api_keys/ipums.txt", open = "rt", blocking = F)
IPUMS_API_KEY <- readLines(con, n = 1)
close(con)

cps_samples_df <- get_sample_info("cps", api_key = IPUMS_API_KEY)

# Oct. '22 was when telework questions were added to the CPS
cps_samples <- cps_samples_df$name[which(cps_samples_df$description == "IPUMS-CPS, October 2022"):nrow(cps_samples_df)]

cps_wfh_msa_def <- define_extract_micro(
  collection = "cps",
  description = "CPS extract for Telework Analysis of Kastle Metros",
  samples = cps_samples,
  variables = c("WTFINL", "HWTFINL", "TELWRKPAY", "TELWRKHR", 
                "TELWRKBFCVD", "TELWRKDIFFCVD", "IND", "OCC", "METFIPS")
)

cps_wfh_msa_submitted <- submit_extract(cps_wfh_msa_def, api_key = IPUMS_API_KEY)
cps_wfh_msa_complete <- wait_for_extract(cps_wfh_msa_submitted, api_key = IPUMS_API_KEY)

cps_wfh_msa_complete$status
cps_wfh_msa_num <- cps_wfh_msa_complete$number
names(cps_wfh_msa_complete$download_links)

download_extract(cps_wfh_msa_complete, 
                 download_dir = "./microdata/",
                 api_key = IPUMS_API_KEY)

ddi <- read_ipums_ddi(paste0("./microdata/cps_000", cps_wfh_msa_num, ".xml"))

cps_wfh_data <- read_ipums_micro(ddi)

# Calculating weighted average of teleworking rates
cps_kastle_wfh_long <- cps_wfh_data %>% 
  filter(METFIPS %in% kastle_fips, TELWRKPAY != 0) %>% # Need to filter out to only those in the universe and in our MSAs of interest
  group_by(YEAR, MONTH, METFIPS) %>% 
  summarize(wfh_pct = weighted.mean(TELWRKPAY == 1, WTFINL)) %>% 
  ungroup() %>% 
  filter(!is.na(wfh_pct))

cps_kastle_wfh_long_us <- cps_wfh_data %>% 
  filter(TELWRKPAY != 0) %>% 
  group_by(YEAR, MONTH) %>% 
  summarize(wfh_pct = weighted.mean(TELWRKPAY == 1, WTFINL)) %>% 
  ungroup() %>% 
  filter(!is.na(wfh_pct)) %>% 
  mutate(METFIPS = 1)

# Writing out data to ./raw_data/ folder
bind_rows(cps_kastle_wfh_long, cps_kastle_wfh_long_us) %>% 
  write_csv("./data/raw_data/kastle_msa_us_wfh_22_24_cps.csv")
  
kastle_fips_df <- msa_wfh_22_prop %>% 
  select(GEOID, NAME, msa_short_name) %>% 
  mutate(GEOID = as.character(GEOID))

msa_wfh_22_24_prop <- bind_rows(cps_kastle_wfh_long,
          cps_kastle_wfh_long_us) %>% 
  mutate(METFIPS = as.character(METFIPS), 
         date = base::as.Date(paste0(YEAR, "-", MONTH, "-01")),
         wfh_pct = wfh_pct * 100) %>% 
  inner_join(kastle_fips_df, by = c("METFIPS" = "GEOID")) %>% 
  select(date, msa_short_name, wfh_pct) %>% 
  pivot_wider(names_from = msa_short_name, values_from = wfh_pct)

write_csv(
  msa_wfh_22_24_prop,
  "./data/clean_data/census_bls_cps_22_24_wfh_dw_line.csv"
)

# Looking at teleworking rates for workers by industry for Philly and the country
# as a whole.

# I'll need to use the Census industry to NAICS industry crosswalk to get 
# fewer NAICS codes.
# https://www.census.gov/topics/employment/industry-occupation/guidance/code-lists.html
# https://www.bls.gov/sae/additional-resources/naics-supersectors-for-ces-program.htm

cps_wfh_data <- cps_wfh_data %>% 
  mutate(supersector = case_when(
    between(IND, 170, 490) ~ "10|Natural Resources and Mining",
    IND == 770 ~ "20|Construction",
    between(IND, 1070, 3990) ~ "30|Manufacturing",
    (between(IND, 4070, 5790) | between(IND, 6070, 6390) | between(IND, 570, 690)) ~ "40|Trade, Transportation, and Utilities",
    between(IND, 6470, 6780) ~ "50|Information",
    between(IND, 6870, 7190) ~ "55|Financial Activities",
    between(IND, 7270, 7790) ~ "60|Professional and Business Services",
    between(IND, 7860, 8470) ~ "65|Education and Health Services",
    between(IND, 8561, 8690) ~ "70|Leisure and Hospitality",
    between(IND, 8770, 9290) ~ "80|Other Services",
    between(IND, 9370, 9890) ~ "90|Government",
    IND == 0 ~ "00|Unknown",
    T ~ NA
  )) %>% 
  separate_wider_delim(cols = supersector, names = c("supersector_code", "supersector_name"), delim = "|")
  
philly_msa_wfh_sprsctr <- cps_wfh_data %>% 
  filter(METFIPS == 37980, TELWRKPAY != 0) %>% 
  group_by(supersector_code, supersector_name) %>% 
  summarize(philly_wfh_pct = weighted.mean(TELWRKPAY == 1, WTFINL) * 100) %>% 
  ungroup()
  
us_wfh_sprsctr <- cps_wfh_data %>% 
  filter(TELWRKPAY != 0) %>% 
  group_by(supersector_code, supersector_name) %>% 
  summarize(us_wfh_pct = weighted.mean(TELWRKPAY == 1, WTFINL) * 100) %>% 
  ungroup()

us_philly_wfh_sprsctr <- inner_join(philly_msa_wfh_sprsctr, 
           us_wfh_sprsctr, 
           by = c("supersector_code", "supersector_name")) %>% 
  filter(supersector_code != "00") %>% 
  arrange(desc(us_wfh_pct))

write_csv(us_philly_wfh_sprsctr,
          "./data/clean_data/census_bls_cps_22_24_philly_us_dw_dot.csv"
          )

# Looking at WFH rates "office workers" aka Information, Financial Activities,
# and Professional and Business Services that have the three highest WFH percentages.

msa_cps_office_workers_wfh <- cps_wfh_data %>% 
  filter(METFIPS %in% kastle_fips, 
         TELWRKPAY != 0,
         supersector_code %in% c("50", "55", "60")) %>% # Need to filter out to only those in the universe and in our MSAs of interest and are "office workers"
  group_by(METFIPS) %>% 
  summarize(wfh_pct = weighted.mean(TELWRKPAY == 1, WTFINL) * 100) %>% 
  ungroup() %>% 
  mutate(METFIPS = as.character(METFIPS)) %>% 
  inner_join(kastle_fips_df, by = c("METFIPS" = "GEOID")) 

us_cps_office_workers_wfh <- cps_wfh_data %>% 
  filter(TELWRKPAY != 0,
         supersector_code %in% c("50", "55", "60")) %>% # Need to filter out to only those in the universe and in our MSAs of interest and are "office workers"
  summarize(wfh_pct = weighted.mean(TELWRKPAY == 1, WTFINL) * 100,
            METFIPS = "1") %>% 
  ungroup() %>% 
  inner_join(kastle_fips_df, by = c("METFIPS" = "GEOID"))

cps_office_workers_wfh <- bind_rows(msa_cps_office_workers_wfh,
                                    us_cps_office_workers_wfh) %>% 
  arrange(desc(wfh_pct))

write_csv(
  cps_office_workers_wfh,
  "./data/clean_data/census_bls_cps_22_24_msa_us_office_workers_dw_bar.csv"
  )

# Recreating monthly time series above but now only with office worker's teleworking rate.
cps_kastle_office_workers_wfh_long <- cps_wfh_data %>% 
  filter(METFIPS %in% kastle_fips, TELWRKPAY != 0, supersector_code %in% c("50", "55", "60")) %>% # Need to filter out to only those in the universe and in our MSAs of interest and for office workers
  group_by(YEAR, MONTH, METFIPS) %>% 
  summarize(wfh_pct = weighted.mean(TELWRKPAY == 1, WTFINL)) %>% 
  ungroup() %>% 
  filter(!is.na(wfh_pct))

cps_kastle_office_workers_wfh_long_us <- cps_wfh_data %>% 
  filter(TELWRKPAY != 0, supersector_code %in% c("50", "55", "60")) %>% 
  group_by(YEAR, MONTH) %>% 
  summarize(wfh_pct = weighted.mean(TELWRKPAY == 1, WTFINL)) %>% 
  ungroup() %>% 
  filter(!is.na(wfh_pct)) %>% 
  mutate(METFIPS = 1)

# Writing out raw data
bind_rows(cps_kastle_office_workers_wfh_long, cps_kastle_office_workers_wfh_long_us) %>% 
  write_csv("./data/raw_data/kastle_msa_us_office_workers_wfh_22_24_cps.csv")

msa_wfh_22_24_prop <- bind_rows(cps_kastle_office_workers_wfh_long,
                                cps_kastle_office_workers_wfh_long_us) %>% 
  mutate(METFIPS = as.character(METFIPS), 
         date = base::as.Date(paste0(YEAR, "-", MONTH, "-01")),
         wfh_pct = wfh_pct * 100) %>% 
  inner_join(kastle_fips_df, by = c("METFIPS" = "GEOID")) %>% 
  select(date, msa_short_name, wfh_pct) %>% 
  pivot_wider(names_from = msa_short_name, values_from = wfh_pct)

# Getting BLS SAE data to get supersector proportion
user_email <- "adriannesta@gmail.com"
sae_supersector <- get_bls_data("https://download.bls.gov/pub/time.series/sm/sm.supersector",
                             email = user_email)

sae_all <- get_bls_data("https://download.bls.gov/pub/time.series/sm/sm.data.1.AllData",
                        email = user_email)

sae_parsed <- sae_all %>%
  mutate(seas_adj = str_sub(series_id, 3, 3),
         state_code = str_sub(series_id, 4, 5),
         area_code = str_sub(series_id, 6, 10),
         supersector_code = str_sub(series_id, 11, 12),
         naics_code = str_sub(series_id, 13, 18),
         industry_code = str_sub(series_id, 11, 18),
         measure_code = str_sub(series_id, 19, 20),
         date = base::as.Date(paste0(year, "-", str_sub(period, 2, 3), "-01")),
         value = as.numeric(value)) %>% 
  filter(period == "M13", 
         year == "2023",
         naics_code == "000000",
         seas_adj == "U",
         area_code %in% kastle_fips,
         measure_code == "01")  # getting latest year annual averages for msas and measure I want

write_csv(
  sae_parsed,
  "./data/raw_data/kastle_msas_23_emp_supersector.csv"
)

sae_supersector_prop_msa <- sae_parsed %>%   
  inner_join(kastle_fips_df, by = c("area_code" = "GEOID")) %>% 
  inner_join(sae_supersector, by = "supersector_code") %>% 
  filter(supersector_code %in% c("15", "30", "40", "50", "55", "60", "65", "70", "80", "90"),
         area_code %in% c("16980", "19100", "35620", "37980")
         ) %>%  # filter to only supersectors needed for viz and four cities I want to visualize
  group_by(area_code, msa_short_name) %>% 
  mutate(sector_prop = prop.table(value) * 100) %>% 
  ungroup() %>% 
  select(msa_short_name, supersector_name, sector_prop) %>% 
  pivot_wider(names_from = msa_short_name, values_from = sector_prop) %>% 
  arrange(desc(Philadelphia))

write_csv(sae_supersector_prop_msa,
          "./data/clean_data/bls_sae_23_msa_supersector_prop_dw_dot.csv")

# Making Oct. '22 to May '24 WFH rates for office worker industries vs. office attendance by metro.

kastle_avg_occ_1022_0524 <- kastle_back_to_work_barometer %>% 
  filter(between(data_week, base::as.Date("2022-10-01"), base::as.Date("2024-05-31"))) %>% 
  group_by(city) %>% 
  summarize(kastle_avg_occ = mean(occupancy))

kastle_office_occ_cps_wfh_1022_0524 <- kastle_avg_occ_1022_0524 %>% 
  inner_join(msa_cps_office_workers_wfh, by = c("city" = "msa_short_name")) %>% 
  rename(`Office Occupancy` = kastle_avg_occ, `Telework Rate` = wfh_pct)

write_csv(
  kastle_office_occ_cps_wfh_1022_0524,
  "./data/clean_data/kastle_cps_office_occ_wfh_1022_0524_dw_scatter.csv"
)

# Making county and PUMS maps of WFH percentage for the Philly Area
philly_msa_counties <- c(
  "34005",
  "34007",
  "34015",
  "34033",
  "42017",
  "42029",
  "42045",
  "42091",
  "42101",
  "10003",
  "24015"
)

county_wfh_22 <- get_acs(
  geography = "county",
  variables = c("wfh" = "B08006_017"),
  year = 2022,
  survey = "acs1",
  geometry = F,
  summary_var = "B08006_001", # My population here is _workers_ 16 years and older.
  key = CENSUS_API_KEY
)

philly_msa_county_wfh_22 <- county_wfh_22 %>% 
  filter(GEOID %in% philly_msa_counties) %>% 
  mutate(wfh_prop = (estimate / summary_est) * 100, # getting derived margins of error
         msa_short_name = str_remove(str_remove(NAME, "-.*$"), ",.*$")) %>% 
  select(GEOID, NAME, msa_short_name, wfh_prop) %>% 
  arrange(desc(wfh_prop))

write_csv(philly_msa_county_wfh_22,
          "./data/clean_data/census_acs_22_1_wfh_philly_msa_county_dw_map.csv"
          )
