  #' Get Descendant Concept IDs
  #' 
  #' @param cdm A CDM reference object from CDMConnector
  #' @param concept_ids A numeric vector of ancestor concept IDs
  #' 
  #' @return A numeric vector of descendant concept IDs
get_concept_descendants <- function(cdm, concept_ids) {
  if (length(concept_ids) == 0) {
    warning("No concept IDs provided")
    return(numeric(0))
  }
  
  descendants <- cdm$concept_ancestor |>
    dplyr::filter(ancestor_concept_id %in% concept_ids) |>
    dplyr::select(descendant_concept_id) |>
    dplyr::distinct() |>
    dplyr::pull(descendant_concept_id)
  
  return(descendants)
}

dbms <- Sys.getenv("pgh_dbms")
server <- Sys.getenv("pgh_server")
user <- Sys.getenv("pgh_user")
password <- Sys.getenv("pgh_password")
port <- Sys.getenv("pgh_port")
host <- Sys.getenv("pgh_host")
dbname <- Sys.getenv("pgh_dbname")
driver <- Sys.getenv("pgh_driver")
cdmSchema <- Sys.getenv("pgh_cdm_schema")
writeSchema <- Sys.getenv("pgh_scratch_schema")

con <- DBI::dbConnect(
  drv = RPostgres::Postgres(),
  dbname = dbname,
  host = host,
  port = as.integer(port),
  user = user,
  password = password,
  bigint = c("numeric")
)

message("Connected to database")

cdm <- CDMConnector::cdmFromCon(
  con,
  cdmName = "papageorgiou",
  cdmSchema = cdmSchema,
  writeSchema = writeSchema,
  cohortTables = c("outcomes", "cohorts")
)

ckd_concept_ids <- get_concept_descendants(cdm, 46271022)
all_ckd_concepts <- c(ckd_concept_ids, 46271022)

diabetes_concept_ids <- get_concept_descendants(cdm, 201820)
all_diabetes_concepts <- c(diabetes_concept_ids, 201820)

liver_concept_ids <- get_concept_descendants(cdm, 194984)
all_liver_concepts <- c(liver_concept_ids, 194984)

defaultSpecification <- FeatureExtraction::getDefaultTable1Specifications()

specification <- data.frame(
  label = c(
    "Age",
    "Gender: female",
    "Medical history: Liver disease",
    "Medical history: Diabetes",
    "Medical history: Kidney disease",
    "Medical history: Cardiovascular disease",
    "Medication use"
  ),
  analysisId = c(2, 1, 210, 210, 210, 210, 410),
  covariateIds = c(
    NA,
    "8532001", # Gender
    paste(paste0(all_liver_concepts, 210), collapse = ","),
    paste(paste0(all_diabetes_concepts, 210), collapse = ","),
    paste(paste0(all_ckd_concepts, 210), collapse = ","),
    paste(defaultSpecification[6, ]$covariateIds, collapse = ","),
    paste(paste(defaultSpecification[8, ]$covariateIds, collapse = ","), "956874410", sep = ",")
  )
)

connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = dbms,
  user = user,
  password = password,
  server = file.path(server, dbname),
  pathToDriver = driver
)

covariateSettings <- FeatureExtraction::createCovariateSettings(
  useDemographicsGender = TRUE,
  useDemographicsAge = TRUE,
  useDemographicsAgeGroup = TRUE,
  useConditionGroupEraLongTerm  = TRUE,
  useDrugGroupEraLongTerm = TRUE,
  useConditionEraAnyTimePrior = TRUE,
  useDrugExposureAnyTimePrior = TRUE,
  addDescendantsToInclude = TRUE
)

table1CovariateSettings <- FeatureExtraction::createTable1CovariateSettings(
  specifications = specification,
  covariateSettings = covariateSettings,
  addDescendantsToInclude = TRUE,
  excludedCovariateConceptIds = c(),
  addDescendantsToExclude = FALSE
)

covariateData <- FeatureExtraction::getDbCovariateData(
  connectionDetails = connectionDetails,
  tempEmulationSchema = NULL,
  cdmDatabaseSchema = cdmSchema,
  cdmVersion = "5",
  cohortTable = "cohorts",
  cohortDatabaseSchema = writeSchema,
  cohortTableIsTemp = FALSE,
  cohortIds = 1,
  rowIdField = "subject_id",
  covariateSettings = covariateSettings,
  aggregated = TRUE
)

FeatureExtraction::createTable1(
  covariateData1 = covariateData,
  covariateData2 = NULL,
  cohortId1 = NULL,
  cohortId2 = NULL,
  specifications = specification,
  output = "two columns",
  showCounts = FALSE,
  showPercent = TRUE,
  percentDigits = 1,
  valueDigits = 1,
  stdDiffDigits = 2
)
