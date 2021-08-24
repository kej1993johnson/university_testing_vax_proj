# university_testing_vax_proj

There are three Rmd files in this repository that do the following:
1. Estimate initial infection prevalence on a university campus using number of students per county in the US and NYT case count data
2. Estimate vaccination rate on a university campusing using either TX DSHS 18-49 year old county-level vaccination data or CDC 18+ county-level data
3. Run the transmission model with the cost analysis to project cases at each vaccination rate and testing frequency, and find for each vaccination rate the
frequency and number of tests needed to keep cases under a certain threshold, and the associated cost of those different scenarios. Also performs projections of
isolation facility needs for on campus housing assuming symptomatic testing only is performed. All parameters (epidemiological and scenario-specific) can be modified (and should be as we learn about VE). 

The names of these files and their required input data sets if any are:
1. UT_case_introduction_estimates.Rmd 
      student_origins: zip or county of students
      uszips.csv: zip to county mapping
      county_pops.csv: population size in each county
      cases_by_county: https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties-recent.csv (NYT case counts)
2. UT_vaccine_coverage_estimates.Rmd
      student_origins: zip or county of students
      uszips.csv: zip to county mapping
      county_vaccinations_age_date_2021-07-26.csv: vaccinations by age group and county in TX from DSHS (from spencerwoody's github)
      COVID-19_Vaccinations_in_the_United_States_County.csv : vaccinations by county in US from CDC
3. UT_projections_report.Rmd
      no inputs needed, but can start in the middle by inputing par_table, cost_table, and par_bounds or any of the dataframes needed to make the figures
