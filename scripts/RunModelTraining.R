source("code/Helper.R")

library(PatientLevelPrediction)

dbms <- Sys.getenv("pgh_dbms")
user <- Sys.getenv("pgh_user")
password <- Sys.getenv("pgh_password")
server <- Sys.getenv("pgh_server")
port <- Sys.getenv("pgh_port")
driver <- Sys.getenv("pgh_driver")
host <- Sys.getenv("pgh_host")
dbname <- Sys.getenv("pgh_dbname")

cdmDatabaseSchema <- "omop"
workDatabaseSchema <- "analysis"
cdmDatabaseName <- "pgh"
cohortTable <- "cohorts"
outcomeTable <- "outcomes"

connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = dbms,
  server = glue::glue("{ host }/{ dbname }"),
  user = user,
  password = password,
  port = port,
  pathToDriver = driver
)

databaseDetails <- createDatabaseDetails(
  connectionDetails = connectionDetails, 
  cdmDatabaseSchema = cdmDatabaseSchema, 
  cdmDatabaseName = cdmDatabaseName, 
  cohortDatabaseSchema = workDatabaseSchema, 
  cohortTable = cohortTable, 
  outcomeDatabaseSchema = workDatabaseSchema, 
  outcomeTable = outcomeTable,
  cdmVersion = 5
)

restrictPlpDataSettings <- createRestrictPlpDataSettings()

# predict outcome within 1 to 180 days after index
# remove people with outcome prior and with < 365 days observation
populationSettings <- createStudyPopulationSettings(
  binary = TRUE, 
  firstExposureOnly = FALSE, 
  washoutPeriod = 1, 
  removeSubjectsWithPriorOutcome = FALSE,
  priorOutcomeLookback = 9999,
  requireTimeAtRisk = FALSE, 
  riskWindowStart = 1, 
  riskWindowEnd = 30
)

# use age/gender in groups and condition groups as features
covariateSettings <- FeatureExtraction::createCovariateSettings(
  useDemographicsGender = TRUE, 
  useDemographicsAge = TRUE,
  useConditionOccurrenceAnyTimePrior = TRUE,
  useDrugExposureAnyTimePrior = TRUE,
  includedCovariateConceptIds = c(
    956874, 19010482, 950637,
    19081224, 45775324, 19080985,
    316139, 197320, 46271022,
    434610, 3006906, 4064161,
    201820, 201620, 132797,
    81902, 132797
  ),
  addDescendantsToInclude = TRUE
)

modelDesign1 <- createModelDesign(
  targetId = 1, 
  outcomeId = 1, 
  restrictPlpDataSettings = restrictPlpDataSettings, 
  populationSettings = populationSettings, 
  covariateSettings = covariateSettings, 
  featureEngineeringSettings = createFeatureEngineeringSettings(),
  sampleSettings = createSampleSettings(), 
  splitSettings = createDefaultSplitSetting(), 
  preprocessSettings = createPreprocessSettings(), 
  modelSettings = setLassoLogisticRegression()
)

modelDesign2 <- createModelDesign(
  targetId = 1, 
  outcomeId = 1, 
  restrictPlpDataSettings = restrictPlpDataSettings, 
  populationSettings = populationSettings, 
  covariateSettings = covariateSettings, 
  featureEngineeringSettings = createRandomForestFeatureSelection(
    ntrees = 500,
    maxDepth = 7
  ),
  sampleSettings = createSampleSettings(), 
  splitSettings = createDefaultSplitSetting(), 
  preprocessSettings = createPreprocessSettings(), 
  modelSettings = setRandomForest()
)

modelDesign3 <- createModelDesign(
  targetId = 1, 
  outcomeId = 1, 
  restrictPlpDataSettings = restrictPlpDataSettings, 
  populationSettings = populationSettings, 
  covariateSettings = covariateSettings, 
  featureEngineeringSettings = createRandomForestFeatureSelection(),
  sampleSettings = createSampleSettings(), 
  splitSettings = createDefaultSplitSetting(), 
  preprocessSettings = createPreprocessSettings(), 
  modelSettings = setGradientBoostingMachine()
)

results <- runMultiplePlp(
  databaseDetails = databaseDetails,
  modelDesignList = list(
    modelDesign1,
    modelDesign2,
    modelDesign3
  ),
  onlyFetchData = FALSE,
  logSettings = createLogSettings(),
  saveDirectory = "./PlpMultiOutput"
)
