utils::globalVariables(
  c(".", "hy_id", "cs_id", "pt_id", "Z", "middle_index", "point_type", "minZ", 
    "maxZ", "minZ_bottom", "maxZ_left_bank", "maxZ_right_bank", "valid_left_bank", 
    "valid_right_bank", "bottom", "left_bank", "right_bank", "valid_banks", 
    "relative_distance", "cs_lengthm", "default_middle", "has_relief", 
    "max_relief", "braid_id", "geometry",
    
    "comid", "fromnode", "tonode", 
    "tocomid", "divergence", "cycle_id", "node", "braid_vector", "totdasqkm", 
    "changed", "relative_position", "head_distance", "tail_distance", 
    "component_id", "cs_measure", "ds_distance", "along_channel", "euclid_dist", 
    "sinuosity", "points_per_cs", "Z_at_bottom", "lower_bound", "upper_bound", 
    "ge_bottom", "is_near_bottom", "pts_near_bottom", "total_valid_pts", 
    "pct_near_bottom", 
    "member_braids",  "braid_members", "diff_pts", "is_extended", 
    "new_cs_id", "split_braid_ids",
    
    "braid_length", 
    "crosswalk_id", 
    "lengthm", 
    "check_z_values", 
    "geom", 
    "is_same_Z", 
    "is_multibraid", 
    "channel", "unique_count",
    "left_bank_count", "right_bank_count", "channel_count", "bottom_count", 
    "terminalID",
    "tmp_id",
    "make_geoms_to_cut_plot",
    "Y", "improved", "length_vector_col", "median", "min_ch", "new_validity_score",
    "old_validity_score", "transects", "validity_score", "x",
    "A", "DEPTH", "DINGMAN_R", "TW", "X", "X_end", "X_start", "Y_end", "Y_start",
    "ahg_a", "ahg_index", "ahg_x", "ahg_y", 
    "bottom_end", "bottom_length", "bottom_midpoint", 
    "bottom_start", "cs_partition", "distance_interval", "fixed_TW", 
    "has_new_DEPTH", "has_new_TW", "ind", "is_dem_point", "left_max", 
    "left_start", "max_right_position", "new_DEPTH", "new_TW", "next_X_is_missing", "next_Y_is_missing",
    "parabola", "partition", "prev_X_is_missing", 
    "prev_Y_is_missing", "right_start", "right_start_max", "start_or_end", "start_pt_id",
    "cs_source", 
    "partition_lengthm", "left_fema_index", "right_fema_index", 
    "left_is_within_fema", "right_is_within_fema", "left_distance", "right_distance",
    "new_cs_lengthm", "polygon_index",
    "crosswalk_id",  
    "anchors", "deriv_type", "edge", "extension_distance", 
    "left_is_extended", "right_is_extended", "to_node", "verbose", 
    "toindid", "indid", "toid", "is", "internal_is_braided2"
  )
)

#' Add a "middle_index" column denoting the middle index of a specific point_type 
#' The middle index is relative to the rest of the cross section points in each hy_id, cs_id
#'
#' @param cs_pts cross section points dataframe with a "point_type" column and "hy_id", "cs_id" columns
#' @param crosswalk_id character, unique ID column 
#' @param point_type character, which point type to get the middle index for. Must be one of  "left_bank", "bottom", "right_bank", or "channel". Default is "bottom"
#' @param default_col_name logical, whether the output column should be named "middle_index" or if 
#' the new column should take the point_type string and use that in the column name (i.e. "left_bank_middle_index" instead of "middle_index"). 
#' Default is TRUE and adds a column named "middle_index"
#'
#' @return dataframe of the input cs_pts with an added middle index column
#' @importFrom dplyr group_by mutate n ungroup select across any_of
#' @noRd
#' @keywords internal
add_middle_index_by_point_type <- function(
    cs_pts, 
    crosswalk_id = NULL, 
    point_type = "bottom",
    default_col_name = TRUE
) {
  
  # Throw an error if an invalid "point_type" value is given
  if(!point_type %in% c("left_bank", "bottom", "right_bank", "channel")) {{
    stop("Invalid 'point_type' value, 'point_type' must be one of: 'bottom', 'channel', 'left_bank', 'right_bank'")
  }}
  
  # Add a middle index column for the given point type, 
  # if a cross section does NOT have the given point type, then 
  # the middle index of the entire cross section is used as a default value
  cs_pts <- 
    cs_pts %>% 
    dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id, "cs_id")))) %>% 
    dplyr::mutate(
      default_middle = dplyr::n() %/% 2,
      middle_index   =  ifelse(
        identical(which(point_type == !!point_type)[ceiling(length(which(point_type == !!point_type)) / 2)], integer(0)), 
        default_middle,
        which(point_type == !!point_type)[ceiling(length(which(point_type == !!point_type)) / 2)]
      )
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-default_middle)
  
  # Make the custom column name to replace the default "middle_index" column name
  if(!default_col_name) {
    new_column_name <- paste0(point_type, "_middle_index")
    names(cs_pts)[names(cs_pts) == "middle_index"] <- new_column_name
  }
  
  return(cs_pts)
  
}

#' Add the degree angle between the middle index of a specific point type and the maximum XY points to the left and right of the given middle index
#' Uses Law of Cosines to determine angle from a given point given a set of 3 points that can form a triangle the rest of the cross section points in each hy_id, cs_id
#' @param cs_pts cross section points dataframe with a "point_type" column and "hy_id", "cs_id" columns
#' @param crosswalk_id character, unique ID column 
#' @param angle_at character, which point type to get the degree angle for. Must be one of  "left_bank", "bottom", "right_bank", or "channel". Default is "bottom"
#' @param default_col_name logical, whether the output column should be named "angle_at" or 
#' if the new column should take the "point_type" string and use. Default is TRUE and adds a column named "angle_at"
#'
#' @return dataframe, the cs_pts dataframe with an added "angle_at" column
#' @importFrom dplyr group_by mutate ungroup select across any_of
#' @noRd
#' @keywords internal
add_angle_at_point_type <- function(cs_pts, 
                                    crosswalk_id = NULL,
                                    # point_type = "bottom",
                                    angle_at = "bottom",
                                    default_col_name = TRUE
) {

  # Throw an error if an invalid "angle_at" value is given
  if(!angle_at %in% c("left_bank", "bottom", "right_bank", "channel")) {{
    stop("Invalid 'angle_at' value, 'angle_at' must be one of: 'bottom', 'channel', 'left_bank', 'right_bank'")
  }}
  
  cs_pts <- 
    cs_pts %>% 
    add_middle_index_by_point_type(crosswalk_id = crosswalk_id, 
                                   point_type = angle_at) %>% 
    dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id, "cs_id")))) %>% 
    dplyr::mutate(
      angle_at = angle_at_index(pt_id, Z, middle_index[1])
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-middle_index)
  
  # Make the custom column name to replace the default "middle_index" column name
  if(!default_col_name) {
    new_column_name <- paste0("angle_at_", point_type)
    names(cs_pts)[names(cs_pts) == "angle_at"] <- new_column_name
  }
  
  return(cs_pts)
}

#' Function to calculate the angle using the Law of Cosines at a given index of X, Y, points
#'
#' @param x numeric vector of size n
#' @param y numeric vector of size n 
#' @param middle_index numeric value, indicating middle index X, Y point to calculate the angle at (can be obtained from add_middle_index_by_point_type())
#'
#' @return numeric angle in degrees between the middle_index point and the maximum Y value XY points to the left and right of middle_index point
#' @noRd
#' @keywords internal
angle_at_index <- function(x, y, middle_index = NULL) {
  
  # get the number of points
  n <- length(x)
  
  # if no index is given, just use the minimum Y value index
  if(is.null(middle_index)) {
    middle_index <- which.min(y)
  }
  
  # Check if the index is valid
  if (middle_index < 1 || middle_index > n) {
    stop("Index out of range")
  }
  
  # Find the maximum Y value to the left of the given index (Y1 for triangle)
  Y_left <- ifelse(middle_index == 1, 
                   NA,
                   max(y[1:(middle_index-1)], na.rm = TRUE)
  )
  
  # find the Y value at the given index (Y2 for triangle)
  Y_middle <- y[middle_index]
  
  # Find the maximum Y value to the right of the given index (Y3 for triangle)
  Y_right <- ifelse(middle_index == n, 
                    NA,
                    max(y[(middle_index+1):n], na.rm = TRUE)
  )
  
  # Find corresponding X values for the maximum Y value on the LEFT side of line (X1 for triangle)
  X_left <- ifelse(middle_index == 1, 
                   NA,
                   x[1:(middle_index-1)][
                     which.max(y[1:(middle_index-1)])
                   ]
  )
  
  # Find corresponding X values at the MIDDLE index (X2 for triangle)
  X_middle <- x[middle_index]
  
  # Find corresponding X values for the maximum Y value on the RIGHT side of line (X3 for triangle)
  X_right <- ifelse(middle_index == n, 
                    NA,
                    x[(middle_index+1):n][
                      which.max(y[(middle_index+1):n])
                    ]
  )
  
  # TODO: come up with a better way of dealing/outputting invalid scenario (NA might be the best option but this needs more work/thought)
  # if invalid (missing a proper left, right, middle point(s)), return an NA value
  if (any(is.na(c(X_left, X_middle, X_right, Y_left, Y_middle, Y_right)))) {
    return(NA)
  }
  
  # Calculate the Euclidean distances for all 3 points of a triangle: 
  # ---> triangle points are the middle (specified index) and the 2 maximum points to the LEFT and RIGHT of the middle
  
  # distance from max on LEFT side to middle (bottom)
  ab <- sqrt(((X_left - X_middle)**2) + ((Y_left - Y_middle)**2))
  
  # distance from middle (bottom) to max on RIGHT side
  bc <- sqrt(((X_middle - X_right)**2) + ((Y_middle - Y_right)**2))
  
  # distance from max on LEFT side to max on RIGHT side
  ac <- sqrt(((X_left - X_right)**2) + ((Y_left - Y_right)**2))
  
  # TODO: come up with a better way of dealing/outputting invalid scenario (NA might be the best option but this needs more work/thought)
  # if missing values, return NA
  if (any(is.na(c(ab, bc, ac)))) {
    return(NA)
  }
  
  # check triangle inequality is satisfied
  if (!(ab + bc <= ac || ab + ac <= bc || bc + ac <= ab)) {
    angle_radians <- acos((ac^2 - ab^2 - bc^2) / (-2 * ab * bc))
    angle <- angle_radians * 180 / pi
  } else {
    angle <- NA
  }
  
  return(angle)
}
