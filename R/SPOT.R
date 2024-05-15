#' build_K_matrix
#'
#' Construct a matrix of dimension \eqn{n \times P} where \eqn{n} is the number
#' of samples and \eqn{P} is the number of radii under consideration of the
#' values of Ripley's K or Besag's L evaluated at each radius for each sample.
#'
#' This function is primarily used within the \code{SPOT} function.
#'
#' @param data The dataset of single cell spatial data. Must include column names
#' PID, id, x, y, and type if \code{marked}. PID corresponds to the patient IDs
#' which should be replicated across ROIs for each patient. id corresponds to the
#' image ID within a patient. x and y give the coordinates of each cell. type (if given)
#' is a cell type label.
#' @param radii Vector of length \eqn{P} giving the radii at which to evaluate the
#' spatial summary measure. The first value must be zero.
#' @param use.K Boolean. If \code{TRUE}, calculates Ripley's K for each radius and
#' sample. If \code{FALSE}, calculates Besag's L for each radius and sample.
#' @param homogeneous Boolean. If \code{TRUE}, calculates the homogeneous versions of
#' Ripley's K or Besag's L. If \code{FALSE}, calculates the inhomogeneous versions of
#' these summaries.
#' @param marked Boolean. If \code{TRUE}, subsets the data to consider only cells with the
#' labels given in the \code{cell.type} parameter. If \code{FALSE}, considers all detected cells
#' included in \code{data}.
#' @param cell.type The cell type labels to consider. Can be one string, e.g. \code{cell.type = "CD4 T cell"}, or
#' a vector of two cell types, e.g. \code{cell.type = c("CD4 T cell", "CD8 T cell")}. If \code{!marked}, leave
#' this parameter as \code{NULL}.
#' @param K.diff Boolean. If \code{TRUE}, returns the difference between the empirical and theoretical values of
#' Ripley's K or Besag's L. If \code{FALSE}, returns the empirical value of the spatial summary.
#' @param nlarge Parameter for the \code{Kest} and \code{Lest} function for efficiency. If there are more cells
#' than \code{nlarge}, will apply a border correction for edge effects instead of the isotropic edge correction.
#' Default is 10000.
#' @param print.progress Boolean. Should progress in calculating the spatial summary across samples be printed?
#' If \code{TRUE}, prints a progress bar.
#'
#' @return Returns two objects: \code{Kr} and \code{Kr.df}. \code{Kr} is a list of length \eqn{n} containing
#' the \eqn{P} spatial summary values at each radius. \code{Kr.df} is the data.frame containing the spatial summaries
#' in \code{Kr} but includes the PIDs and image ids.
#'
#' @importFrom magrittr %>%
#' @importFrom tidyselect all_of
#' @importFrom dplyr across
#'
#' @export
build_K_matrix <- function(data, radii, use.K, homogeneous, marked, cell.type, K.diff, nlarge = 10000, print.progress = TRUE) {

  # Summarize the images within each PID
  pids.images <- data %>%
    dplyr::select(PID, id) %>%
    distinct()

  # Save the unique image ids
  image.ids <- pids.images$id

  # Initialize list for each image
  Kr <- lapply(1:length(image.ids), function(i) list())
  names(Kr) <- image.ids

  # Iterate through the images
  for (i in 1:length(image.ids)) {
    # Print out the progress
    if (print.progress) svMisc::progress(i/(length(image.ids)/100))

    # Save the current image ID
    id.i <- image.ids[i]

    # Should the point process marked?
    if (!marked) {

      # Save the data
      data.xy <- data %>%
        dplyr::filter(id == id.i) %>%
        dplyr::select(x, y)

      # Convex hull for image
      w <- spatstat.geom::convexhull.xy(data.xy$x, data.xy$y)

      # Convert to a point process
      data.ppp <- spatstat.geom::as.ppp(data.xy, W=w)

      # Apply Ripley's K/Besag's L
      if (homogeneous) {

        # Ripley's K
        if (use.K) {
          Ki <- Ki.hom <- spatstat.explore::Kest(data.ppp, r = radii, nlarge = nlarge, correction = "isotropic")
        }

        # Besag's L
        if (!use.K) {
          Ki <- Ki.hom <- spatstat.explore::Lest(data.ppp, r = radii, nlarge = nlarge, correction = "isotropic")
        }
      }

      if (!homogeneous) {

        # Ripley's K
        if (use.K) {
          Ki <- Ki.inhom <- spatstat.explore::Kinhom(data.ppp, r = radii, nlarge = nlarge, correction = "isotropic")
        }

        # Besag's L
        if (!use.K) {
          Ki <- Li.inhom <- spatstat.explore::Linhom(data.ppp, r = radii, nlarge = nlarge, correction = "isotropic")
        }
      }
    }

    if (marked) {

      # Filter to this id
      data.types <- data %>% dplyr::filter(id == id.i) %>% dplyr::select(x, y, type)

      # Convert the types to a factor
      data.types$type <- as.factor(data.types$type)

      # Convert to a ppp
      w <- spatstat.geom::convexhull.xy(data.types$x, data.types$y)
      data.types.ppp <- spatstat.geom::as.ppp(data.types, W = w, marks = data.types$type)

      # Subset to specific cell type
      data.types.ppp.subset <- subset(data.types.ppp, marks %in% cell.type)

      # Are there 2 cell types?
      bivariate <- length(cell.type) == 2

      # Ripely's K/Besag's L for a specific cell type
      if (!bivariate) {

        # Ripley's K
        if (use.K) {

          if (homogeneous) {
            Ki <- Ki.type <- spatstat::Kest(data.types.ppp.subset, r = radii, nlarge = nlarge, correction = "isotropic")
          }

          if (!homogeneous) {
            Ki <- Ki.type <-  spatstat::Kinhom(data.types.ppp.subset, r = radii, nlarge = nlarge, correction = "isotropic")
          }
        }

        # Besag's L
        if (!use.K) {

          if (homogeneous) {
            Ki <- Ki.type <-  spatstat::Lest(data.types.ppp.subset, r = radii, nlarge = nlarge, correction = "isotropic")
          }

          if (!homogeneous) {
            Ki <- Ki.type <-  spatstat::Linhom(data.types.ppp.subset, r = radii, nlarge = nlarge, correction = "isotropic")
          }

        }
      }

      if (bivariate) {

        # Check if both cell types are available
        cell.type.counts <- table(data.types.ppp.subset$marks)[table(data.types.ppp.subset$marks) != 0]

        # If the number of non-zero cell type counts is equal to the number of cell types
        if (length(cell.type.counts) == length(cell.type)) {

          if (use.K) {

            if (homogeneous) {
              Ki <- Ki.type <-  spatstat.explore::Kcross(data.types.ppp.subset, i = cell.type[1], j = cell.type[2],
                                       r = radii, nlarge = nlarge, correction = "isotropic")
            }

            if (!homogeneous) {
              Ki <- Ki.type <-  spatstat.explore::Kcross.inhom(data.types.ppp.subset, i = cell.type[1], j = cell.type[2],
                                             r = radii, nlarge = nlarge, correction = "isotropic")
            }
          }

          if (!use.K) {

            if (homogeneous) {
              Ki <- Ki.type <-  spatstat.explore::Lcross(data.types.ppp.subset, i = cell.type[1], j = cell.type[2],
                                       r = radii, nlarge = nlarge, correction = "isotropic")
            }

            if (!homogeneous) {
              Ki <- Ki.type <-  spatstat.explore::Lcross.inhom(data.types.ppp.subset, i = cell.type[1], j = cell.type[2],
                                             r = radii, nlarge = nlarge, correction = "isotropic")
            }
          }

        }

        # If both cell types are not available
        if (length(cell.type.counts) != length(cell.type)) {
          Ki <- Ki.type <- list(iso = rep(NA, length(radii)), theo = rep(NA, length(radii)))
        }
      }
    }

    # Save the result
    if (!K.diff) {
      Kr[[i]] <- Ki$iso
    }

    if (K.diff) {
      Kr[[i]] <- Ki$iso - Ki$theo
    }
  }

  # Create a data.frame
  Kr.df <- do.call(rbind, Kr)

  # Add column names
  colnames(Kr.df) <- paste0("radius.", radii)

  # Combine with PIDs.timepoint
  Kr.df <- cbind.data.frame(pids.images, Kr.df)

  # Return the results
  list(Kr = Kr, Kr.df = Kr.df)

}


#' SPatial Omnibus Test (SPOT)
#'
#' Perform a SPatial Omnibus Test (SPOT) to relate a spatial summary of each
#' image to patient outcomes.
#'
#' @param data The dataset of single cell spatial data. Must include column names
#' PID, id, x, y, and type if \code{marked}. PID corresponds to the patient IDs
#' which should be replicated across ROIs for each patient. id corresponds to the
#' image ID within a patient. x and y give the coordinates of each cell. type (if given)
#' is a cell type label. Should also contain an outcome column which will be specified in
#' the \code{outcome} parameter.
#' @param radii Vector of length \eqn{P} giving the radii at which to evaluate the
#' spatial summary measure. The first value must be zero.
#' @param outcome String specifying the column name in \code{data} containing the outcome.
#' Must be a survival outcome, a binary outcome, or a continuous outcome.
#' @param censor String specifying the column name in \code{data} containing the censoring indicator.
#' Leave this parameter \code{NULL} if the outcome is not survival.
#' @param model.type String specifying which predictive model to use. The options are: "survival" for a
#' Cox proportional hazards model, or "logistic" for a logistic regression model, or "linear" for a
#' standard linear model.
#' @param use.K Boolean. If \code{TRUE}, calculates Ripley's K for each radius and
#' sample. If \code{FALSE}, calculates Besag's L for each radius and sample.
#' @param homogeneous Boolean. If \code{TRUE}, calculates the homogeneous versions of
#' Ripley's K or Besag's L. If \code{FALSE}, calculates the inhomogeneous versions of
#' these summaries.
#' @param adjustments String(s) specifying the column names in \code{data} of covariates
#'  to adjust for in the predictive model. May be left \code{NULL}.
#' @param marked Boolean. If \code{TRUE}, subsets the data to consider only cells with the
#' labels given in the \code{cell.type} parameter. If \code{FALSE}, considers all detected cells
#' included in \code{data}.
#' @param cell.type The cell type labels to consider. Can be one string, e.g. \code{cell.type = "CD4 T cell"}, or
#' a vector of two cell types, e.g. \code{cell.type = c("CD4 T cell", "CD8 T cell")}. If \code{!marked}, leave
#' this parameter as \code{NULL}.
#' @param K.diff Boolean. If \code{TRUE}, returns the difference between the empirical and theoretical values of
#' Ripley's K or Besag's L. If \code{FALSE}, returns the empirical value of the spatial summary.
#' @param pick.roi String specifying how to select among ROIs if multiple are available per sample. Default is \code{all}
#' which should be used if only one image is available per person. If multiple are available, specify this parameter
#' as \code{"average"} to average the spatial summaries across ROIs within each radius or \code{random} to randomly select
#' an ROI for each sample.
#' @param seed Specify a seed to use if \code{pick.roi = "random"}.
#' @param nlarge Parameter for the \code{Kest} and \code{Lest} function for efficiency. If there are more cells
#' than \code{nlarge}, will apply a border correction for edge effects instead of the isotropic edge correction.
#' Default is 10000.
#' @param prop.nonzero How many spatial summaries should be non-zero in order for the corresponding radius to be
#' included in the predictive model? Specify as a proportion of the total number of samples. Default is 0.2.
#' @param print.progress Boolean. Should progress in calculating the spatial summary across samples be printed?
#' If \code{TRUE}, prints a progress bar.
#'
#' @return Returns a list with the following elements:
#' \item{overall.pval}{The overall p-value resulting from the Cauchy combination test characterizing the
#' strength of association between the spatial summary and the outcome across radii.}
#' \item{Kr.df}{The matrix of spatial summaries across samples and radii. If multiple ROIs were person and \code{pick.roi = "average" or "random"},
#' returns only the averaged spatial summaries or the summaries for the randomly-selected ROIs.}
#' \item{pval.df}{P-values for the association between the spatial summary at each radius and the outcome. }
#' \item{pids.images}{The PIDs and the corresponding image IDs within each patient. }
#' \item{images.missing.cells}{How many images were missing cells altogether?}
#' \item{n.images.association}{How many images were used in the predictive model?}
#' \item{arguments}{Returns the arguments specified in the function call.}
#'
#' @importFrom dplyr n
#' @importFrom tidyselect all_of
#' @importFrom dplyr across
#' @importFrom survival coxph.control
#' @importFrom magrittr %>%
#'
#' @export
spot <- function(data, radii, outcome, censor = NULL, model.type, use.K = TRUE,
                                 homogeneous = TRUE, adjustments = NULL, marked = FALSE,
                                 cell.type = NULL, K.diff = FALSE, pick.roi = "all", seed = 12345,
                                 nlarge = 10000, prop.nonzero = 0.2, print.progress = TRUE) {

  # Save the arguments
  arguments <- match.call()

  # -------------------------------------
  # Calculate the K(r)s or L(r)s
  # -------------------------------------

  # Summarize the images within each PID
  pids.images <- data %>%
    dplyr::select(PID, id) %>%
    dplyr::distinct()

  # Save the unique image ids
  image.ids <- pids.images$id

  # Output the Kr results
  Kr.results <- build_K_matrix(data = data, radii = radii, use.K = use.K, homogeneous = homogeneous,
                               marked = marked, cell.type = cell.type, K.diff = K.diff, nlarge = nlarge,
                               print.progress = print.progress)

  # Save
  Kr.df <- Kr.results$Kr.df

  # -------------------------------------
  # If there are multiple ROIs per person, choose one to continue with
  # -------------------------------------

  # Calculate the column indices for the K(r)s
  col.inds <- which(grepl("radius.", colnames(Kr.df)))

  # Remove NA rows due to missingness
  rows.missing.cells <- apply(Kr.df[,col.inds], 1, function(row) all(is.na(row)))
  Kr.df <- Kr.df[!rows.missing.cells,]

  # How many images were lost due to lack of cells?
  images.missing.cells <- sum(rows.missing.cells)

  # Filter to the number of images remaining
  pids.images <- pids.images %>% dplyr::filter(id %in% Kr.df$id)

  if (pick.roi != "all") {

    # Pick a random ROI for each individual
    if (pick.roi == "random") {

      # Set a seed
      if (is.null(seed)) seed <- 12345
      set.seed(seed)

      # For each subject ID, check whether there are multiple ROIs and count how many
      Kr.df <-  Kr.df %>%
        dplyr::group_by(PID) %>%
        dplyr::mutate(multiple.ROI = any(duplicated(PID)),
                      number.ROI = n(),
                      selected.ROI = as.numeric(id == sample(id, size = 1)))

      selected.ROI <- Kr.df$selected.ROI

      # Select just the selected ROIs again
      Kr.df <- Kr.df %>% dplyr::filter(selected.ROI == 1)

      # Filter to the number of images remaining
      pids.images <- pids.images %>% dplyr::filter(id %in% Kr.df$id)

    }

    # Take the mean of the L(r)s within each patient
    if (pick.roi == "average") {

      Kr.df <-  Kr.df %>%
        dplyr::select(-id) %>%
        dplyr::group_by(PID) %>%
        dplyr::summarise(across(everything(), mean))

    }
  }

  # Update col.inds if need be
  if (!identical(col.inds, which(grepl("radius.", colnames(Kr.df))))) {
    col.inds <- which(grepl("radius.", colnames(Kr.df)))
  }

  # -------------------------------------
  # Some radii will be mostly 0. Remove these here.
  # -------------------------------------

  # How many images were ultimately used in the association testing?
  n.images.association <- nrow(Kr.df)

  # How many non-zero K(r) or L(r)s are there for each radius?
  num.nonzero.Kr <- apply(Kr.df[,col.inds], 2, function(col) sum(col != 0, na.rm = TRUE))

  # Save the number of non-zero K(r)s or L(r)s needed
  required.nonzero.Kr <- round(prop.nonzero * n.images.association)

  # Which columns (radii) meet this threshold?
  radii.to.save <- names(num.nonzero.Kr)[num.nonzero.Kr >= required.nonzero.Kr]

  # Save the radii with the corresponding number
  if ("id" %in% colnames(Kr.df)) {
    Kr.df <- Kr.df[, c("PID", "id", radii.to.save)]
  }

  if (!("id" %in% colnames(Kr.df))) {
    Kr.df <- Kr.df[, c("PID", radii.to.save)]
  }

  # -------------------------------------
  # Test association with given outcome
  # if enough samples are available
  # -------------------------------------

  # Test association for each radius (multiple ROIs per person)
  radii.to.save.numeric <- as.numeric(sapply(strsplit(radii.to.save, "radius."), function(i) i[2]))
  pval.df <- data.frame(radius = radii.to.save.numeric, coef = NA, pval = NA)

  # If we have enough images
  if (n.images.association > 5) {

    # Iterate through radii
    for (i in 1:length(radii.to.save)) {

      # Select the radius
      radius.i <- radii.to.save[i]

      # Subset the K(r)s for this radius
      if ("id" %in% colnames(Kr.df)) {
        Kr.subset <- Kr.df %>% dplyr::select(id, PID)
      }

      if (!("id" %in% colnames(Kr.df))) {
        Kr.subset <- Kr.df[, c("PID")]
      }

      Kr.subset <- cbind.data.frame(Kr.subset, ripley = Kr.df[,radius.i,drop=TRUE])

      # Add in the outcome data
      outcome.data <- dplyr::right_join(data %>%
                                          dplyr::select(all_of(c("PID", outcome, censor, adjustments))) %>%
                                          dplyr::distinct(),
                                        Kr.subset,
                                        by = "PID")

      if (model.type == "linear") {
        # Save the model equation
        model.equation <- reformulate(c("ripley", adjustments), "get(outcome)")

        # Run the model
        res <- lm(model.equation, data = outcome.data)
      }

      if (model.type == "survival") {
        # Save the model equation
        model.equation <- reformulate(c("ripley", adjustments), "Surv(get(outcome), get(censor))")

        # Run the model
        res <- survival::coxph(model.equation, data = outcome.data, control = coxph.control(iter.max = 100))
      }

      if (model.type == "logistic") {
        # Save the model equation
        model.equation <- reformulate(c("ripley", adjustments), "get(outcome)")

        # Run the model
        res <- glm(model.equation, family = binomial(), data = outcome.data)
      }

      # Save the results
      coef <- summary(res)$coefficients
      row.ind <- rownames(coef) == "ripley"
      col.inds <- c(1, ncol(coef))
      pval.df$coef[i] <- coef[row.ind, col.inds[1]]
      pval.df$pval[i] <- coef[row.ind, col.inds[2]]
    }
  }

  # Remove any NA rows due to all the same value (not an issue after removing all 0 K(r) or L(r)s)
  pval.df <- pval.df[!is.na(pval.df$pval),]

  # Initialize the overall p-values at NA
  overall.pval <- NA

  # Combine results with and without rounding
  if (!is.logical(pval.df$pval)) {
    pval.df$pval.round <- pval.df$pval

    # Combine results (no rounding)
    overall.pval <- ACAT::ACAT(pval.df$pval)
  }

  # Return
  list(overall.pval = overall.pval,
       Kr.df = Kr.df,
       pval.df = pval.df,
       pids.images = pids.images,
       images.missing.cells = images.missing.cells,
       n.images.association = n.images.association,
       arguments = arguments)
}
