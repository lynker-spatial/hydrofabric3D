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
    "toindid", "indid", "toid", "is", "internal_is_braided2",
    "index"
  )
)

#' Give a set of transecct linestrings and poylgons and get the minimum distance to extend each transect line (from both directions, to try and reach the edge of a "polygons")
#' Superseces old version of function (now named extend_transects_to_polygons2())
#' @param transect_lines Set of Sf linestrings to extend (only if the transect lines are ENTIRELLY within a polygons)
#' @param polygons set of sf polygons that transect lines should be exteneded 
#' @param flowlines set of Sf linestrings
#' @param crosswalk_id character, flowline ID that matches flowlines with transect lines. This crosswalk_id must appear are a column in both flowlines and transect_lines.
#' @param grouping_id character, name of a column in flowlines that should be used to group each transect with 1 or more flowlines. 
#' That is, when transects are checked to make sure they don't intersect 
#' other transects or other flowlines, this group ID will distinguise which flowlines a transect should be checked against.
#' The intersect_group_id must appear as a column in both flowlines and transect_lines dataframes
#' @param max_extension_distance numeric, maximum distance (meters) to extend a transect line in either direction to try and intersect one of the "polygons". Default is 3000m
#' @param tolerance A minimum distance to use for simplification on polygons. Use a higher value for more simplification on the polygons. Default is NULL which will apply no simplification to polygons.
#' @param keep_lengths logical whether to keep a record of the original transect lengths or not, default is FALSE, original lengths are not kept 
#' @param reindex_cs_ids logical, whether to reindex the cs_ids to ensure each crosswalk_id has cs_ids of 1-number of transects. Default is TRUE, which makes sure if any cross sections were removed from a crosswalk_id, 
#' then the cs_ids are renumbered so there are no gaps between cs_ids within a crosswalk_id. 
#' Setting this to FALSE will make sure crosswalk_id/cs_ids remain untouched as they were given in the input data.
#' @param verbose logical, whether to output messages or not. Default is TRUE, and messages will output 
#' 
#' @return sf linestring, with extended transect lines
#' @importFrom rmapshaper ms_simplify
#' @importFrom geos as_geos_geometry geos_intersects_matrix geos_simplify_preserve_topology geos_within_matrix geos_empty geos_point_start geos_point_end
#' @importFrom sf st_as_sf st_cast st_segmentize st_length st_drop_geometry st_geometry
#' @importFrom dplyr mutate case_when select left_join relocate n any_of group_by ungroup arrange across
#' @importFrom lwgeom st_linesubstring
#' @importFrom wk wk_crs 
#' @importFrom vctrs vec_c
#' @export
extend_transects_to_polygons <- function(
    transect_lines, 
    polygons, 
    flowlines, 
    crosswalk_id = NULL,  
    grouping_id = 'mainstem',
    max_extension_distance = 3000,
    tolerance = NULL,
    keep_lengths = FALSE,
    reindex_cs_ids = TRUE,
    verbose = TRUE
) {
  
  # ----------------------------------------------------------------------------------
  # ----------- Input checking ------
  # ----------------------------------------------------------------------------------
 
  # standardize geometry name
  transect_lines <- hydroloom::rename_geometry(transect_lines, "geometry")
  flowlines      <- hydroloom::rename_geometry(flowlines, "geometry")
  
  is_valid_flowlines <- validate_df(flowlines, 
                             c(crosswalk_id, grouping_id, "geometry"), 
                             "flowlines")
  
  is_valid_transects <- validate_df(transect_lines, 
                             c(crosswalk_id, "cs_id", grouping_id, "geometry"), 
                             "transect_lines")
  
  # ----------------------------------------------------------------------------------
  
  # stash input UIDs, used in a check at the end of this function to make sure all unique hy_id/cs_id in the INPUT are in the OUTPUT,
  # and raise an error if they're are missing hy_id/cs_ids
  input_uids    <- get_unique_tmp_ids(transect_lines, x = crosswalk_id)
  
  # stash the original starting transect line lengths
  starting_lengths <- 
    transect_lines %>%  
    add_length_col(length_col = "initial_length") %>% 
    sf::st_drop_geometry() %>% 
    dplyr::select(dplyr::any_of(crosswalk_id), cs_id, initial_length)
  
  # # # preserve the initial ordering
  # starting_order <- get_transect_initial_order(transect_lines, crosswalk_id)
  
  # get only the relevent polygons/transects
  transect_subset   <- subset_transects_in_polygons(transect_lines, polygons)
  polygons_subset   <- subset_polygons_in_transects(transect_lines, polygons)

  # get a dataframe that tells you how far to extend each line in either direction
  extensions_by_id  <- get_transect_extension_distances_to_polygons(transects = transect_subset, 
                                            polygons = polygons_subset, 
                                            crosswalk_id = crosswalk_id, 
                                            max_extension_distance = max_extension_distance,
                                            tolerance = tolerance,
                                            verbose = verbose
                                            )
  
  # TODO: Add left/right extension distancces to transect data
  # TODO: this can ultimately just be the "transects" variable, dont need to make new "transects_with_distances" variable
  transect_lines <-
    transect_lines %>% 
    dplyr::left_join(
      dplyr::select(extensions_by_id, 
                    dplyr::any_of(crosswalk_id), cs_id, left_distance, right_distance),
      by = c(crosswalk_id, "cs_id")
    ) %>% 
    # TODO: I think i want to keep the NAs and NOT fill w/ 0 
    dplyr::mutate(
      left_distance = dplyr::case_when(
        is.na(left_distance) ~ 0,
        TRUE                 ~ left_distance
      ),
      right_distance = dplyr::case_when(
        is.na(right_distance) ~ 0,
        TRUE                  ~ right_distance
      )
    ) %>% 
    hydrofabric3D::add_tmp_id(x = crosswalk_id, y = "cs_id") 
  
  transect_lines <- extend_transects_sides(
    transects    = transect_lines,
    flowlines    = flowlines,
    crosswalk_id = crosswalk_id,
    cs_id        = "cs_id",
    grouping_id  = grouping_id,
    direction    = "any_by_specific_distances"
  )
  
  # Set the is_extended flag based on if either the left OR the right side were extended
  transect_lines <-
    transect_lines %>%
    hydroloom::rename_geometry("geometry") %>%
    dplyr::mutate(
      is_extended = dplyr::case_when(
        left_is_extended | right_is_extended ~ TRUE,
        TRUE                                 ~ FALSE
      )
    ) 
  
  # shorten any transects that intersect multiple transects back to their original lengths
  transect_lines  <- shorten_multi_transect_intersecting_extended_transects(x = transect_lines, 
                                                          crosswalk_id = crosswalk_id)
  
  # shorten any transects that intersect multiple flowlines (or a flowline more than once) back to their original lengths
  transect_lines  <- shorten_multi_flowline_intersecting_extended_transects(x = transect_lines, 
                                                                   flowlines = flowlines,
                                                                   crosswalk_id = crosswalk_id)
  
  # remove transects that intersect with OTHER TRANSECTS
  transect_lines <-
    transect_lines %>%
    rm_multi_intersects()
  
  # remove transects that intersect multiple flowlines
  transect_lines <- 
    transect_lines %>% 
    rm_multiflowline_intersections(flowlines = flowlines) 
  
  # # make sure correct length column
  # transect_lines <- 
  #   transect_lines %>% 
  #   add_length_col("cs_lengthm")
  
  # check to make sure all unique hy_id/cs_id in the INPUT are in the OUTPUT,
  # and raise an error if they're are missing hy_id/cs_ids
  output_uids   <- get_unique_tmp_ids(transect_lines, x = crosswalk_id)
  
  # check all of the output_uids exist in the input UIDs
  has_all_uids           <- all(output_uids %in% input_uids)
  
  # throw an error if NOT all crosswalk_id/cs_ids are the same in the input and output data
  if (!has_all_uids) {
    stop("Missing unique crosswalk_id/cs_id UIDs from input transects in the output transects")
  }
  
  if (keep_lengths) {
    transect_lines <- 
      transect_lines %>% 
      dplyr::left_join(
        starting_lengths %>% 
          dplyr::select(dplyr::any_of(crosswalk_id), cs_id, initial_length),
        by = c(crosswalk_id, "cs_id")
      ) %>% 
      dplyr::select(
        dplyr::any_of(c(crosswalk_id, "cs_id", "cs_source", "initial_length")),
        cs_lengthm, cs_measure,
        left_distance, right_distance,
        geometry
      )
  }
  
  # re-index the cs_ids to make sure there are 1-number of transects for each crosswalk_id and that there are NO gaps between cs_ids
  if (reindex_cs_ids) {
    warning("Re-indexing cs_ids may result in a mismatch between unique crosswalk_id/cs_ids in the input 'transect_lines' and the output data's unique crosswalk_id/cs_ids")
    transect_lines <- renumber_cs_ids(transect_lines, crosswalk_id = crosswalk_id)
  }
  
  # ----------------------------------
  # ---- Final reorder ----
  # ----------------------------------
  
  # transect_lines <- 
  #   transect_lines %>% 
  #   dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id, "cs_id")))) %>% 
  #   dplyr::arrange(cs_id, .by_group = TRUE) %>% 
  #   dplyr::ungroup()
  
  # transect_lines <- 
  #   transect_lines %>% 
  #   dplyr::left_join(
  #     starting_order,
  #     by = crosswalk_id
  #   ) %>% 
  #   dplyr::arrange(initial_order, cs_id) %>% 
  #   dplyr::select(-initial_order)
  
  # ----------------------------------
  # ---- Final column select ----
  # ----------------------------------
  
  # remove added tmp_id column
  transect_lines <-
    transect_lines %>% 
      dplyr::select(
        dplyr::any_of(
          c(crosswalk_id, 
            "cs_id", 
            "cs_lengthm", 
            "cs_measure", 
            "ds_distance",
            "sinuosity", 
            
            "initial_length", 
            "left_distance", 
            "right_distance", 
            
            "cs_source",
            "geometry"
            )
        )
        # cs_lengthm, cs_measure,
        # left_distance, right_distance,
        # geometry
      )
      
  return(transect_lines)

}



#' Calculate the minimum distance a line would need to extend to reach the boundary of the polygon/line that the input geometries are entirely within 
#' VERSION 2
#' @param geos_geoms list of geos_geometrys
#' @param ids character vector
#' @param lines_to_cut geos_linestrings
#' @param lines_to_cut_indices numeric vector
#' @param direction character, either "head", "tail" or "both"
#' @param max_extension_distance numeric
#'
#' @return numeric vector, distance to extend each geos_geoms
#' @importFrom vctrs vec_c
#' @noRd
#' @keywords internal
calc_extension_distances <- function(
    geos_geoms, 
    ids, 
    lines_to_cut, 
    lines_to_cut_indices, 
    direction = "head", 
    max_extension_distance = 2500
) {
  
  if (!direction %in% c("head", "tail")) {
    stop("Invalid 'direction' value, must be one of 'head' or 'tail'")
  }
  
  distance_range               <- 1:max_extension_distance
  
  # preallocate vector that stores the extension. distances
  extension_dists              <- vctrs::vec_c(rep(0, length(ids)))
  
  # number of geometries that will be iterated over, keeping this variable to reference in message block  
  total                        <- length(ids)
  
  # output a message every ~10% intervals
  message_interval             <- total %/% 20
  
  make_progress <- make_progress_bar(TRUE, length(ids))
  
  for (i in seq_along(ids)) {
    # log percent complete
    make_progress()
    
    # # log percent complete
    # if (message_interval != 0 && i %% message_interval == 0) {
    #   # get the percent complete
    #   percent_done <- round(i/total, 2) * 100
    #   message(i, " > ", percent_done, "% ") 
    # }
    
    index_vect <- unlist(lines_to_cut_indices[[i]])
    
    distance_to_extend <- geos_bs_distance(
      distances    = distance_range,
      line         = geos_geoms[[i]],
      geoms_to_cut = lines_to_cut[index_vect],
      direction    = direction
    )
    
    extension_dists[i] <- distance_to_extend
    
  }
  
  return(extension_dists)
}

#' Given 2 geos_geometry point geometries, create a line between the 2 points
#'
#' @param start geos_geoemtry, point
#' @param end geos_geoemtry, point
#' @param line_crs crs
#' @importFrom geos geos_y geos_x geos_make_linestring
#' @return geos_geometry linestring
#' @noRd
#' @keywords internal
make_line_from_start_and_end_pts <- function(start, end, line_crs) {
  
  # Y_start <- geos::geos_y(start)
  # X_start <- geos::geos_x(start)
  # Y_end   <- geos::geos_y(end)
  # X_end   <- geos::geos_x(end)
  # 
  # # make the new transect line from the start and points 
  # geos_ls <- geos::geos_make_linestring(x = c(X_start, X_end),
  #                                       y = c(Y_start, Y_end), 
  #                                       crs = line_crs)
  # return(geos_ls)  
  
  # make the new transect line from the start and points 
  return(
    geos::geos_make_linestring(
      x   = c(geos::geos_x(start), geos::geos_x(end)),
      y   = c(geos::geos_y(start), geos::geos_y(end)), 
      crs = line_crs
    )
  )
  
}

#' Check if an updated transect line is valid relative to the other transects and flowlines in the network
#' The 'transect_to_check' should be 'used' (i.e. function returns TRUE) if 
#' the 'transect_to_check' does NOT interesect any other transects ('transect_lines') AND it only intersects a single flowline ONCE.
#' If the 'transect_to_check' intersects ANY other transects OR intersects a flowline more
#' than once (OR more than one flowline in the network) then the function returns FALSE.
#' @param transect_to_check geos_geometry, linestring
#' @param trans geos_geometry, linestring
#' @param flines geos_geometry, linestring
#'
#' @return TRUE if the extension should be used, FALSE if it shouldn't be used
#' @importFrom geos geos_intersection geos_type geos_intersects
#' @noRd
#' @keywords internal
is_valid_transect_line <- function(transect_to_check, trans, flines) {
  
  # ###   ##   ##   ##   ##   ##   ##   ##   ##   ##  
  # Define conditions to decide which version of the transect to use
  
  # 1. Use transect with extension in BOTH directions
  # 2. Use transect with LEFT extension only
  # 3. Use transect with RIGHT extension only
  
  # check for NULL / empty geometries
  if (is.null(transect_to_check) || is.null(trans) || is.null(flines)) {
    stop("invalid geometry: geometries cannot be NULL")
  }
  
  # Check that the extended transect lines only intersect a single flowline in the network only ONCE
  intersects_with_flowlines <- geos::geos_intersection(
    transect_to_check,
    flines
  )
  
  intersects_flowline_only_once <- sum(geos::geos_type(intersects_with_flowlines) == "point") == 1 && 
    sum(geos::geos_type(intersects_with_flowlines) == "multipoint") == 0 
  
  # NOTE: return early if if the transects does NOT intersect the flowline ONLY once
  # little optimization, avoids extra geos_intersects() calls 
  if (!intersects_flowline_only_once) {
    return(FALSE)
  }
  
  # check that the extended transect line does NOT intersect other transect lines (other than SELF)
  intersects_other_transects <- lengths(geos::geos_intersects_matrix(transect_to_check,  trans)) > 1
  # intersects_other_transects <- sum(geos::geos_intersects(transect_to_check, trans)) > 1
  
  # TRUE == Only one flowline is intersected a single time AND no other transect lines are intersected
  use_transect <- intersects_flowline_only_once  && !intersects_other_transects
  
  return(use_transect)
}

#' Check if an updated transect line is valid relative to the other transects and flowlines in the network (v2)
#' The 'transect_to_check' should be 'used' (i.e. function returns TRUE) if 
#' the 'transect_to_check' does NOT interesect any other transects ('transect_lines') AND it only intersects a single flowline ONCE.
#' If the 'transect_to_check' intersects ANY other transects OR intersects a flowline more
#' than once (OR more than one flowline in the network) then the function returns FALSE.
#' @param transect_to_check geos_geometry, linestring
#' @param trans geos_geometry, linestring
#' @param flines geos_geometry, linestring
#'
#' @return TRUE if the extension should be used, FALSE if it shouldn't be used
#' @importFrom geos geos_intersection geos_type geos_intersects
#' @noRd
#' @keywords internal
is_valid_transect_line2 <- function(transect_to_check, trans, flines) {
  
  # ###   ##   ##   ##   ##   ##   ##   ##   ##   ##  
  # transect_to_check = extended_trans
  # # trans = transects_geos[transect_group_id_array == transect_group_id_array[i]]
  # trans <-    transects_geos[transect_group_id_array == transect_group_id_array[i] & transect_uid_array != current_uid]
  # flines <- flowlines_geos[fline_group_id_array == transect_group_id_array[i]]

  # transects_geos[transect_group_id_array == transect_group_id_array[i] & transect_uid_array != current_uid]
  # transects_geos[transect_group_id_array == transect_group_id_array[i]]
  # flowlines_geos[fline_group_id_array == transect_group_id_array[i]]
  
  # Define conditions to decide which version of the transect to use
  
  # 1. Use transect with extension in BOTH directions
  # 2. Use transect with LEFT extension only
  # 3. Use transect with RIGHT extension only
  
  # check for NULL / empty geometries
  if (is.null(transect_to_check) || is.null(trans) || is.null(flines)) {
    stop("invalid geometry: geometries cannot be NULL")
  }
  
  # Check that the extended transect lines only intersect a single flowline in the network only ONCE
  intersects_with_flowlines <- geos::geos_intersection(
    transect_to_check,
    flines
  )
  
  intersects_flowline_only_once <- sum(geos::geos_type(intersects_with_flowlines) == "point") == 1 && 
    sum(geos::geos_type(intersects_with_flowlines) == "multipoint") == 0 
  
  # NOTE: return early if if the transects does NOT intersect the flowline ONLY once
  # little optimization, avoids extra geos_intersects() calls 
  if (!intersects_flowline_only_once) {
    return(FALSE)
  }
  
  # check that the extended transect line does NOT intersect other transect lines (other than SELF)
  intersects_other_transects <- lengths(
    geos::geos_intersects_matrix(transect_to_check,  trans)
    ) > 0
  # intersects_other_transects <- sum(geos::geos_intersects(transect_to_check, trans)) > 1
  
  # TRUE == Only one flowline is intersected a single time AND no other transect lines are intersected
  use_transect <- intersects_flowline_only_once  && !intersects_other_transects
  
  return(use_transect)
}

#' Get transects that intersect with the polygons 
#' @param transect_lines Set of Sf linestrigns to extend (only if the transect lines are ENTIRELLY within a polygons)
#' @param polygons set of sf polygons that transect lines should be exteneded 
#' @return sf linestring, with extended transect lines
#' @importFrom geos as_geos_geometry geos_intersects_matrix 
#' @noRd
#' @keywords internal
subset_transects_in_polygons <- function(transect_lines, polygons)  {
  
  transects_polygons_matrix <- geos::geos_intersects_matrix(geos::as_geos_geometry(transect_lines), geos::as_geos_geometry(polygons))
  
  # return(transect_lines[lengths(transects_polygons_matrix) != 0, ])
  return(transect_lines[lengths(transects_polygons_matrix) != 0, ])
  
}

#' Get polygons that intersect with the transects 
#' @param transect_lines Set of Sf linestrigns to extend (only if the transect lines are ENTIRELLY within a polygons)
#' @param polygons set of sf polygons that transect lines should be exteneded 
#' @return sf polygon dataframe
#' @importFrom geos as_geos_geometry geos_intersects_matrix 
#' @noRd
#' @keywords internal
subset_polygons_in_transects <- function(transect_lines, polygons)  {
  
  # TODO: this should be a function argument OR removed, shouldn't probably forcibly and silently simplify the input polygons without user knowing..
  # keep 10% of the original points for speed
  # polygons <- rmapshaper::ms_simplify(polygons, keep_shapes = T, keep = 0.01)
  
  # polygons_transects_matrix <- geos::geos_intersects_matrix(polygons_geos, transects_geos)
  polygons_transects_matrix <- geos::geos_intersects_matrix(geos::as_geos_geometry(polygons), geos::as_geos_geometry(transect_lines))
  
  # return(polygons_geos[lengths(polygons_transects_matrix) != 0])
  return(polygons[lengths(polygons_transects_matrix) != 0, ])
  
}

sf_polygons_to_geos_multilinestrings <- function(polygons, tolerance = 200) {
  
  mls <- 
    polygons %>% 
    geos::as_geos_geometry() %>% 
    geos::geos_node()
    
  if (!is.null(tolerance)) {
    mls <- 
      mls %>% 
      geos::geos_simplify_preserve_topology(tolerance)
  }
  
  return(mls)
  
}

#' Helper function for cleaning up transects and adding some meta data after going through partition_transects_for_extension
#'
#' @param partition dtaframe or sf dataframe
#' @param dir character, 'left' or 'right'. Default is 'left'
#' @param crosswalk_id character, unique column ID name
#'
#' @return dataframe or sf dataframe
#' @importFrom dplyr mutate relocate any_of
#' @importFrom sf st_length
#' @noRd
#' @keywords internal
wrangle_paritioned_transects <- function(partition, 
                                         dir = "left", 
                                         crosswalk_id = "hy_id"
) {
  
  partition <- 
    partition %>% 
    dplyr::mutate(
      partition         = dir,
      partition_lengthm = as.numeric(sf::st_length(.))
    ) %>%
    add_tmp_id(x = crosswalk_id, y = "cs_id") %>%
    dplyr::relocate(tmp_id,
                    dplyr::any_of(crosswalk_id),
                    # cs_source,
                    cs_id
                    # cs_measure,
                    # cs_lengthm,
                    # is_extended,
                    # partition, partition_lengthm, geometry
    ) 
  # dplyr::select(tmp_id, dplyr::any_of(crosswalk_id), cs_source, cs_id, cs_measure, cs_lengthm,
  # partition, partition_lengthm, geometry)  
  
  return(partition)
}

# From a set of transects and a set of polygons, subset the transects based on whether their start / end point is fully within a polygon
# This function will return the points transect lines that need to be extended in a given direction
# i.e. When dir = "left" then the returned transects have their starting points entirely within a polygon

#' From a set of transects and a set of polygons, subset the transects based on whether their start / end point is fully within a polygon
#'  This function will return the points transect lines that need to be extended in a given direction. ( i.e. When dir = "left" then the returned transects have their starting points entirely within a polygon)
#' @param transects sf dataframe w/ linestrings
#' @param polygons_subset sf dataframe of polygons 
#' @param dir character, direction to partition
#'
#' @return sf dataframe subset of 'transects' that are fully within polygons
#' @importFrom geos geos_point_start geos_point_end geos_within_matrix
#' @importFrom dplyr filter 
#' @noRd
#' @keywords internal
partition_transects_for_extension <- function(transects, polygons_subset, dir = "left") {
  
  # POINT START = LEFT
  # POINT END   = RIGHT
  
  if (!dir %in% c("left", "right")) {
    stop("Invalid 'dir' value '", dir, "', 'dir' must be either 'left' or 'right'")  
  }
  
  # get the function needed for a given direction, a transects starting point is the "left" and the ending point is the right  
  dir_function      <- ifelse(dir == "left", geos::geos_point_start, geos::geos_point_end)
  
  # determine transects whose starting/end points are within a polygon
  is_within_matrix  <- geos::geos_within_matrix(dir_function(transects), polygons_subset)
  
  is_within_vect    <- lapply(is_within_matrix, function(i) { if(length(i) > 0) { c(i) } else { c(NA) } })
  
  transects$polygon_index <- is_within_vect
  
  # return only the transects whose start / end point are within a polygon
  return(dplyr::filter(transects, !is.na(polygon_index)))
  
}

#' Get the left and right extension distances for a set of transects out to a set of polygons
#'
#' @param transects sf linestring dataframe
#' @param polygons sf polygon dataframe
#' @param crosswalk_id character
#' @param max_extension_distance numeric 
#' @param tolerance A minimum distance to use for simplification on polygons. Use a higher value for more simplification on the polygons. Default is NULL which will apply no simplification to polygons.
#' @param verbose logical, whether to output messages or not. Default is TRUE, and messages will output 
#' 
#' @return data.frame or tibble
#' @export
get_transect_extension_distances_to_polygons <- function(
    transects, 
    polygons, 
    crosswalk_id, 
    max_extension_distance,
    tolerance = NULL,
    verbose = TRUE
) {
  
  # split transects into left and right partitions
  left_partition <- partition_transects_for_extension(
    transects, 
    polygons, 
    dir = "left"
  ) %>% 
    wrangle_paritioned_transects(
      dir          = "left", 
      crosswalk_id = crosswalk_id
    )
  
  right_partition <- partition_transects_for_extension(
    transects, 
    polygons, 
    dir = "right"
  ) %>% 
    wrangle_paritioned_transects(
      dir          = "right", 
      crosswalk_id = crosswalk_id
    )
  
  # Convert the polygon to a MULTILINESTRING geometry for checking extension distances
  mls <- sf_polygons_to_geos_multilinestrings(polygons, tolerance)
  # mls <- sf_polygons_to_geos_multilinestrings(polygons, 100)
  
  message_if_verbose("Generating left side transects distances to polygons... ", verbose = verbose)

  left_distances <- calc_extension_distances(
    geos_geoms             = geos::as_geos_geometry(left_partition),
    ids                    = left_partition$tmp_id,
    lines_to_cut           = mls,
    lines_to_cut_indices   = left_partition$polygon_index,
    direction              = "head",
    max_extension_distance = max_extension_distance
  )
  
  message_if_verbose("Generating right side transects distances to polygons... ", verbose = verbose)
  
  right_distances <- calc_extension_distances(
    geos_geoms             = geos::as_geos_geometry(right_partition),
    ids                    = right_partition$tmp_id,
    lines_to_cut           = mls,
    lines_to_cut_indices   = right_partition$polygon_index,
    direction              = "tail",
    max_extension_distance = max_extension_distance
  )
  
  left_partition$left_distance    <- left_distances
  right_partition$right_distance  <- right_distances
  
  # Distance to extend LEFT and/or RIGHT for each hy_id/cs_id
  extensions_by_id <- dplyr::left_join(
    sf::st_drop_geometry(
      dplyr::select(left_partition, 
                    dplyr::any_of(crosswalk_id),
                    cs_id,
                    left_distance
      )
    ),
    sf::st_drop_geometry(
      dplyr::select(right_partition, 
                    dplyr::any_of(crosswalk_id),
                    cs_id, 
                    right_distance
      )
    ),
    by = c(crosswalk_id, "cs_id")
  )
  
  # add any missing crosswalk_id/cs_id that DID NOT have any extension distance w/ values of 0
  extensions_by_id <- dplyr::bind_rows(
    extensions_by_id, 
    transects %>% 
      sf::st_drop_geometry() %>% 
      hydrofabric3D::add_tmp_id(x = crosswalk_id) %>% 
      dplyr::filter(!tmp_id %in% hydrofabric3D::add_tmp_id(extensions_by_id, x = crosswalk_id)$tmp_id) %>% 
      dplyr::mutate(
        left_distance  = 0,
        right_distance = 0
      ) %>% 
      dplyr::select(
        dplyr::any_of(crosswalk_id), cs_id, 
        left_distance, right_distance
      ) 
  ) 
  # # TODO: i want to make any NA left/right distances set to 0 instead of having both NAs AND 0 values, need to fix tests for this function first
  # extensions_by_id <- 
  #   extensions_by_id %>% 
  #   dplyr::mutate(
  #     left_distance = dplyr::case_when(
  #       is.na(left_distance) ~ 0,
  #       TRUE                 ~ left_distance
  #     ),
  #     right_distance = dplyr::case_when(
  #       is.na(right_distance) ~ 0,
  #       TRUE                  ~ right_distance
  #     )
  #   ) 
  
  return(extensions_by_id)
}

#' Get the left and right extension distances for a set of transects out to a set of polygons
#' Deprecated, new name is get_transect_extension_distances_to_polygons()
#' @param transects sf linestring dataframe
#' @param polygons sf polygon dataframe
#' @param crosswalk_id character
#' @param max_extension_distance numeric 
#' @param tolerance A minimum distance to use for simplification on polygons. Use a higher value for more simplification on the polygons. Default is NULL which will apply no simplification to polygons.
#'
#' @return data.frame or tibble
#' @noRd
#' @keywords internal
get_extensions_by_id <- function(
    transects, 
    polygons, 
    crosswalk_id, 
    max_extension_distance,
    tolerance = NULL
    ) {
  
  # split transects into left and right partitions
  left_partition <- partition_transects_for_extension(
    transects, 
    polygons, 
    dir = "left"
  ) %>% 
    wrangle_paritioned_transects(
      dir          = "left", 
      crosswalk_id = crosswalk_id
    )
  
  right_partition <- partition_transects_for_extension(
    transects, 
    polygons, 
    dir = "right"
  ) %>% 
    wrangle_paritioned_transects(
      dir          = "right", 
      crosswalk_id = crosswalk_id
    )
  
  # Convert the polygon to a MULTILINESTRING geometry for checking extension distances
  mls <- sf_polygons_to_geos_multilinestrings(polygons, tolerance)
  # mls <- sf_polygons_to_geos_multilinestrings(polygons, 100)
  
  message("Generating left side distances....") 
  left_distances <- calc_extension_distances(
    geos_geoms             = geos::as_geos_geometry(left_partition),
    ids                    = left_partition$tmp_id,
    lines_to_cut           = mls,
    lines_to_cut_indices   = left_partition$polygon_index,
    direction              = "head",
    max_extension_distance = max_extension_distance
  )
  
  message("Generating right side distances...")
  right_distances <- calc_extension_distances(
    geos_geoms             = geos::as_geos_geometry(right_partition),
    ids                    = right_partition$tmp_id,
    lines_to_cut           = mls,
    lines_to_cut_indices   = right_partition$polygon_index,
    direction              = "tail",
    max_extension_distance = max_extension_distance
  )
  
  left_partition$left_distance    <- left_distances
  right_partition$right_distance  <- right_distances
  
  # Distance to extend LEFT and/or RIGHT for each hy_id/cs_id
  extensions_by_id <- dplyr::left_join(
                        sf::st_drop_geometry(
                          dplyr::select(left_partition, 
                                        dplyr::any_of(crosswalk_id),
                                        cs_id,
                                        left_distance
                          )
                        ),
                        sf::st_drop_geometry(
                          dplyr::select(right_partition, 
                                        dplyr::any_of(crosswalk_id),
                                        cs_id, 
                                        right_distance
                          )
                        ),
                        by = c(crosswalk_id, "cs_id")
                      )
  
  # add any missing crosswalk_id/cs_id that DID NOT have any extension distance w/ values of 0
  extensions_by_id <- dplyr::bind_rows(
                          extensions_by_id, 
                          transects %>% 
                            sf::st_drop_geometry() %>% 
                            hydrofabric3D::add_tmp_id(x = crosswalk_id) %>% 
                            dplyr::filter(!tmp_id %in% hydrofabric3D::add_tmp_id(extensions_by_id, x = crosswalk_id)$tmp_id) %>% 
                            dplyr::mutate(
                              left_distance  = 0,
                              right_distance = 0
                            ) %>% 
                            dplyr::select(
                              dplyr::any_of(crosswalk_id), cs_id, 
                              left_distance, right_distance
                            ) 
                        ) 
  # # TODO: i want to make any NA left/right distances set to 0 instead of having both NAs AND 0 values, need to fix tests for this function first
  # extensions_by_id <- 
  #   extensions_by_id %>% 
  #   dplyr::mutate(
  #     left_distance = dplyr::case_when(
  #       is.na(left_distance) ~ 0,
  #       TRUE                 ~ left_distance
  #     ),
  #     right_distance = dplyr::case_when(
  #       is.na(right_distance) ~ 0,
  #       TRUE                  ~ right_distance
  #     )
  #   ) 
  
  return(extensions_by_id)
}

#' Decide the start and end points for the final transect line given two extended versions of the same transect
#' Requires two logicals indicating what to do with the extensions (these are decided by checking for intersections with the rest of the network)
#' Internal helper function
#' @param left_extension geos_geometry linestring
#' @param right_extension geos_geometry linestring 
#' @param use_left logical, do we use the left extension
#' @param use_right logical, do we use the right extension
#' @importFrom geos geos_point_start geos_point_end
#' @return geos_geometry points, the start and end point of the final extension line
#' @noRd
#' @keywords internal
pick_extension_pts <- function(
    left_extension, 
    right_extension, 
    use_left, 
    use_right
) {
  
  use_both <- use_left && use_right
  
  # Get the start and end of both extended tranects
  left_start  <- geos::geos_point_start(left_extension)
  left_end    <- geos::geos_point_end(left_extension)
  right_start <- geos::geos_point_start(right_extension)
  right_end   <- geos::geos_point_end(right_extension)
  
  # Extend in BOTH directions
  if(use_both) {
    # message("Extend direction: BOTH")
    start  <- left_start
    end    <- right_end
    
    # extend ONLY the left side
  } else if(use_left && !use_right) {
    # message("Extend direction: LEFT")       
    start  <- left_start
    end    <- left_end
    
    # Extend ONLY the right side
  } else if(!use_left && use_right) {
    # message("Extend direction: RIGHT")       
    start  <- right_start
    end    <- right_end
    
    # DO NOT extend either direction
  } else {
    # message("No extension")   
    # TODO: Really dont need to do anything 
    # TODO: in this scenario because we just use the original transect line
    start  <- left_end
    end    <- right_start
  }
  
  return( c(start, end) )
  
}

#' Get the starting and ending points of geos_linestring 
#' Internal helper function
#' @param extension geos_geometry linestring
#' @param use_extension logical, do we use the extension
#' @importFrom geos geos_point_start geos_point_end
#' @return geos_geometry points, the start and end point of the final extension line
#' @noRd
#' @keywords internal
get_line_node_pts <- function(
    line
) {
  # Get the start and end of the geos_linestring (extended transect)
  return( c(geos::geos_point_start(line),  geos::geos_point_end(line)) )
}


#' Trim a set of transects to the bounds of polygons
#'
#' @param transect_lines sf dataframe
#' @param flowlines sf dataframe
#' @param polygons sf dataframe 
#' @param crosswalk_id character unique ID
#' @param dissolve logical, whether to dissolve polygon internal boundaries or not. Default is FALSE. 
#' @importFrom dplyr filter select any_of mutate bind_rows n slice_max group_by ungroup across distinct
#' @importFrom rmapshaper ms_explode
#' @importFrom sf st_intersection st_geometry_type
#' @importFrom hydroloom rename_geometry 
#' @return sf dataframe
#' @export
trim_transects_to_polygons <- function(transect_lines, 
                                       flowlines,
                                       polygons,
                                       crosswalk_id = NULL,
                                       dissolve = FALSE
                                       ) {
  
  # TODO: this is a hacky way of doing this, can definitely be improved so user wont get an error if there transects/flowlines have 'polygon_id' as the crosswalk_id
  # temporary ID to use for keeping track of polygon IDs
  POLYGON_ID <- "polygon_id"
  
  if(crosswalk_id == POLYGON_ID) {
    stop("Invalid crosswalk_id value of '", crosswalk_id, "'. 'crosswalk_id' can NOT be 'polygon_id'. 
         \nTry renaming your crosswalk_id column in the 'transect_lines' and 'flowlines' input datasets to a name that is not 'polygon_id'"
         )
  }
  
  # standardize geometry name and drop polygon_id if its a column in any of the data
  transect_lines <-
    transect_lines %>% 
    hydroloom::rename_geometry("geometry") %>% 
    dplyr::select(!dplyr::any_of(POLYGON_ID))
  
  flowlines <- 
    flowlines %>% 
    hydroloom::rename_geometry("geometry") %>% 
    dplyr::select(!dplyr::any_of(POLYGON_ID))
  
  
  if (dissolve) {
    polygons <- 
      polygons %>% 
      # # sf::st_union() %>% 
      # sf::st_difference() %>%
      # rmapshaper::ms_explode() %>%
      # sf::st_as_sf() %>%
      dissolve_overlaps()
      # cluster_dissolve()
  }
  
  polygons <- 
    polygons %>% 
    hydroloom::rename_geometry("geometry") %>% 
    dplyr::mutate(
      polygon_id = 1:dplyr::n()
    ) %>% 
    dplyr::select(polygon_id, geometry)
  
  # Validate the transects / flowlines input datasets
  is_valid_transects <- validate_df(transect_lines, 
                                    c(crosswalk_id, "cs_id", "geometry"), 
                                    "transect_lines")
  
  is_valid_flowlines <- validate_df(flowlines, 
                                    c(crosswalk_id, "geometry"), 
                                    "flowlines")
  
  # Figure out which transects intersect the polygons and which do NOT
  split_transects <- 
    add_intersects_ids(
      transect_lines, 
      polygons, 
      POLYGON_ID
    )
  
  suppressWarnings({
    
    # trim any of the transects that DO hit the polygons 
    trimmed_trans <-
      split_transects %>% 
      dplyr::filter(!is.na(.data[[POLYGON_ID]])) %>% 
      dplyr::select(-dplyr::any_of(POLYGON_ID)) %>% 
      sf::st_intersection(polygons) %>% 
      dplyr::select(!dplyr::any_of(c(POLYGON_ID, paste0(POLYGON_ID, ".1")))) %>% 
      dplyr::filter(sf::st_geometry_type(geometry) %in% c("LINESTRING", "MULTILINESTRING")) %>% 
      rmapshaper::ms_explode() %>% 
      add_intersects_ids(
        flowlines %>% 
          dplyr::mutate(
            new_id = .data[[crosswalk_id]]
          ), 
        "new_id"
      ) %>% 
      dplyr::filter(!is.na(new_id), .data[[crosswalk_id]] %in% strsplit(new_id, ", ")) %>% 
      dplyr::select(-new_id)
    
    # make sure all of the intersection transects are distinct and then select the longest of the transects if there are any duplicates
    trimmed_trans <- 
      trimmed_trans %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id, "cs_id")))) %>% 
      add_length_col("length_check") %>% 
      dplyr::slice_max(length_check, with_ties = FALSE) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(-length_check)
  
    
  })
  
  merged_transects <- 
    split_transects %>% 
    dplyr::filter(is.na(.data[[POLYGON_ID]])) %>% 
    dplyr::bind_rows(trimmed_trans) %>% 
    dplyr::select(!dplyr::any_of(c(POLYGON_ID)))
  
  # has_same_uids(merged_transects, transect_lines, crosswalk_id = crosswalk_id)
  missing_trans <- 
    transect_lines %>% 
    add_tmp_id(crosswalk_id) %>% 
    dplyr::filter(!tmp_id %in% get_unique_tmp_ids(merged_transects, crosswalk_id)) %>% 
    dplyr::select(-tmp_id)
  
  # TODO:
  # add any missing transects BACK to the merged_transects set, 
  # these transects COULDN'T be properly dealt with so they will just be thrown back in as is
  merged_transects <- 
    merged_transects %>% 
    dplyr::bind_rows(missing_trans)
  
  has_all_same_uids <- has_same_uids(merged_transects, transect_lines, crosswalk_id = crosswalk_id)
  
  if(!has_all_same_uids) {
    warning("Not all unique crosswalk_id/cs_ids were retained from input...")     
  }
  
  # put transects back in proper order
  merged_transects <- 
    merged_transects %>% 
    cs_arrange(crosswalk_id = crosswalk_id, order_by = "cs_id")
  
  # recalculate lengths
  merged_transects <- 
    merged_transects %>% 
    add_length_col("cs_lengthm")
  
  return(merged_transects)
}

#' Add an ID column from 'y' if it intersects with 'x'
#'
#' @param x sf dataframe
#' @param y sf dataframe
#' @param id_col character, unique ID column name in 'y'
#'
#' @return sf dataframe of x with id_col added if it intersects with y
#' @importFrom sf st_transform st_crs st_intersects
#' @export
add_intersects_ids <- function(x, y, id_col) {
  
  # TODO: Remove this function as an 'export' 
  
  # make sure the crs are tjhe same
  y <- sf::st_transform(y, sf::st_crs(x))
  
  # Perform the intersection
  intersections <- sf::st_intersects(x, y)
  
  # add the intersected values to the first dataframe
  x[[id_col]] <- unlist(lapply(intersections, function(idx) {
    if (length(idx) > 0) {
      paste0(unlist(y[[id_col]][idx]), collapse = ", ")
    } else {
      NA
    }
  }))
  
  return(x)
}

#' Create a 'group_id' for sf polygons based on an sf spatial predicate function
#' An internal method for doing a spatial predicate based 'group_by' 
#' 
#' @param polys sf dataframe of POLYGONS or MULTIPOLYGONS 
#' @param predicate sf geometry binary predicate function (i.e. 'st_intersects', 'st_within', etc.)
#' @importFrom fastmap fastmap
#' @importFrom dplyr bind_rows arrange group_by count ungroup left_join slice_max select
#' @return polys sf dataframe of POLYGONS or MULTIPOLYGONS with added 'group_id' column 
#' @noRd
#' @keywords internal
add_predicate_group_id <- function(polys, predicate) {
  
  relations <- predicate(polys)
  
  relations <- lapply(seq_along(relations), function(i) { as.character(sort(unique(c(relations[i][[1]], i)))) })
  
  group_ids_map <- fastmap::fastmap()
  ids_to_groups <- fastmap::fastmap()
  
  group_id <- 0
  
  for (i in seq_along(relations)) {
    
    predicate_ids <- relations[i][[1]]
    
    id_group_check <- ids_to_groups$has(predicate_ids)
    
    if(any(id_group_check)) {
      
      known_groups  <- ids_to_groups$mget(predicate_ids)
      known_group   <- known_groups[unname(sapply(known_groups , function(kg) {
        !is.null(kg)
      }))][[1]]
      
      past_group_ids     <- group_ids_map$get(known_group)[[1]]
      updated_group_ids  <- as.character(
        sort(as.numeric(unique(c(past_group_ids, predicate_ids))))
      )
      
      group_ids_map$set(known_group, list(updated_group_ids))
      
      new_ids <- predicate_ids[!predicate_ids %in% past_group_ids]
      
      # add any newly added IDs to the seen map
      for (seen_id in new_ids) {
        # message(seen_id)
        ids_to_groups$set(as.character(seen_id), as.character(group_id))
      }
      
    } else {
      # get a new group ID number
      group_id <- group_id + 1    
      
      # create a new key in the map with the predicate IDs list as the value
      group_ids_map$set(as.character(group_id), list(predicate_ids))
      
      # add each predicate ID to the map storing the seen indexes and their respective group IDs 
      for (seen_id in predicate_ids) {
        ids_to_groups$set(as.character(seen_id), as.character(group_id))
      }
    }
  }
  
  group_ids   <- group_ids_map$as_list() 
  
  grouping_df <- lapply(seq_along(group_ids), function(i) {
    grouping  <- group_ids[i] 
    group_id  <- names(grouping)
    indices   <- grouping[[1]][[1]]
    
    data.frame(
      index      = as.numeric(indices),
      group_id   = rep(group_id, length(indices))   
    )
    
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::arrange(i) 
  
  # count up the number of IDs for each group, well use this to determine which group 
  # to put any indices that had MULTIPLE groups they were apart of (use the group with the most other members)
  group_id_counts <- 
    grouping_df %>% 
    dplyr::group_by(group_id) %>% 
    dplyr::count() %>% 
    dplyr::ungroup()
  
  # select the IDs with the most other members
  grouping_df <- 
    grouping_df %>% 
    dplyr::left_join(
      group_id_counts, 
      by = 'group_id'
    ) %>% 
    dplyr::group_by(index) %>% 
    dplyr::slice_max(n, with_ties = FALSE) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-n) %>% 
    dplyr::arrange(-index) 
  
  polys$group_id <- grouping_df$group_id
  
  return(polys)
  
}