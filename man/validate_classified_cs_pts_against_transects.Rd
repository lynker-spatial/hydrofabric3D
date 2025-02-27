% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validators.R
\name{validate_classified_cs_pts_against_transects}
\alias{validate_classified_cs_pts_against_transects}
\title{Validate Classified Cross Section Points Against Transects
Ensure all cross section points are valid relative to a set of transects.
This validates the points in the same manner as validate_cs_pts_against_transects() but also checks that
classification columns ('class', 'point_type', 'valid_banks', 'has_relief') exist.}
\usage{
validate_classified_cs_pts_against_transects(
  cs_pts,
  transects,
  crosswalk_id = NULL
)
}
\arguments{
\item{cs_pts}{sf object, cross section points}

\item{transects}{sf object, transects}

\item{crosswalk_id}{character, column name of the crosswalk id}
}
\value{
logical, TRUE if all validations pass, FALSE otherwise
}
\description{
Validate Classified Cross Section Points Against Transects
Ensure all cross section points are valid relative to a set of transects.
This validates the points in the same manner as validate_cs_pts_against_transects() but also checks that
classification columns ('class', 'point_type', 'valid_banks', 'has_relief') exist.
}
