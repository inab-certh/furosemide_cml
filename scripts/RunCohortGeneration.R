dbms <- Sys.getenv("pgh_dbms")
server <- Sys.getenv("pgh_server")
user <- Sys.getenv("pgh_user")
password <- Sys.getenv("pgh_password")
port <- Sys.getenv("pgh_port")
host <- Sys.getenv("pgh_host")
dbname <- Sys.getenv("pgh_dbname")
driver <- Sys.getenv("pgh_driver")

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
  cdmSchema = "omop",
  writeSchema = "analysis"
)

ancestors <- list(
  "essential_hypertension" = 320128
)

allDescendants <- purrr::map(ancestors, function(concept_id) {
  CodelistGenerator::getDescendants(
    cdm = cdm,
    conceptId = concept_id,
    withAncestor = TRUE
  ) |>
    dplyr::pull(concept_id)
})

message("Defined condition concepts")

cdm$cohorts <- CohortConstructor::conceptCohort(
  cdm = cdm,
  conceptSet = allDescendants,
  name = "cohorts",
  table = "condition_occurrence"
  # subsetCohort = "exposures"
)
message("Generated table: cohorts")

chronic_kidney_disease_concept_ids <- c(46271022) 
ckd_descendants <- cdm$concept_ancestor |>
  dplyr::filter(ancestor_concept_id %in% chronic_kidney_disease_concept_ids) |>
  dplyr::select(descendant_concept_id) |>
  dplyr::distinct() |>
  dplyr::pull(descendant_concept_id)

# Combine parent concept(s) with all descendants
all_ckd_concepts <- c(chronic_kidney_disease_concept_ids, ckd_descendants)

message(paste("Total CKD concepts (including descendants):", length(all_ckd_concepts)))

# Step 3: Exclude patients with chronic kidney disease (including all descendants)
# --------------------------------------------------------------------------------
# requireConceptIntersect with intersections = 0 excludes patients
cdm$cohorts <- cdm$cohorts |>
  CohortConstructor::requireConceptIntersect(
    conceptSet = list("chronic_kidney_disease" = all_ckd_concepts),
    indexDate = "cohort_start_date",
    window = list(c(-Inf, Inf)),  # Check entire patient history
    intersections = c(0, 0),  # Require 0 occurrences (i.e., exclude)
    name = "cohorts"
  )

message("Excluded chronic kidney disease patients")


cdm$outcomes <- CohortConstructor::measurementCohort(
  cdm = cdm,
  name = "outcomes",
  conceptSet = list("creatinine" = c(3016723)),
  valueAsNumber = list(
    "creatinine_increased" = list("8840" = c(1.2, Inf))
  ),
  subsetCohort = "cohorts"
)

message("Generated table: outcome")
message("Run completed successfully")
