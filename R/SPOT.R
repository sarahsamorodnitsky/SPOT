build_K_matrix <- function(data, radii, use.K, homogeneous, marked, cell.type, K.diff, nlarge, print.progress = TRUE) {

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
    if (print.progress) progress(i/(length(image.ids)/100))

    # Save the current image ID
    id.i <- image.ids[i]

    # Should the point process marked?
    if (!marked) {

      # Save the data
      data.xy <- data %>%
        dplyr::filter(id == id.i) %>%
        dplyr::select(x, y)

      # Convex hull for image
      w <- convexhull.xy(data.xy$x, data.xy$y)

      # Convert to a point process
      data.ppp <- as.ppp(data.xy, W=w)

      # Apply Ripley's K/Besag's L
      if (homogeneous) {

        # Ripley's K
        if (use.K) {
          Ki <- Ki.hom <- Kest(data.ppp, r = radii, nlarge = nlarge, correction = "isotropic")
        }

        # Besag's L
        if (!use.K) {
          Ki <- Ki.hom <- Lest(data.ppp, r = radii, nlarge = nlarge, correction = "isotropic")
        }
      }

      if (!homogeneous) {

        # Ripley's K
        if (use.K) {
          Ki <- Ki.inhom <- Kinhom(data.ppp, r = radii, nlarge = nlarge, correction = "isotropic")
        }

        # Besag's L
        if (!use.K) {
          Ki <- Li.inhom <- Linhom(data.ppp, r = radii, nlarge = nlarge, correction = "isotropic")
        }
      }
    }

    if (marked) {

      # Filter to this id
      data.types <- data %>% dplyr::filter(id == id.i) %>% dplyr::select(x, y, type)

      # Convert the types to a factor
      data.types$type <- as.factor(data.types$type)

      # Convert to a ppp
      w <- convexhull.xy(data.types$x, data.types$y)
      data.types.ppp <- as.ppp(data.types, W = w, marks = data.types$type)

      # Subset to specific cell type
      data.types.ppp.subset <- subset(data.types.ppp, marks %in% cell.type)

      # Are there 2 cell types?
      bivariate <- length(cell.type) == 2

      # Ripely's K/Besag's L for a specific cell type
      if (!bivariate) {

        # Ripley's K
        if (use.K) {

          if (homogeneous) {
            Ki <- Ki.type <- Kest(data.types.ppp.subset, r = radii, nlarge = nlarge, correction = "isotropic")
          }

          if (!homogeneous) {
            Ki <- Ki.type <-  Kinhom(data.types.ppp.subset, r = radii, nlarge = nlarge, correction = "isotropic")
          }
        }

        # Besag's L
        if (!use.K) {

          if (homogeneous) {
            Ki <- Ki.type <-  Lest(data.types.ppp.subset, r = radii, nlarge = nlarge, correction = "isotropic")
          }

          if (!homogeneous) {
            Ki <- Ki.type <-  Linhom(data.types.ppp.subset, r = radii, nlarge = nlarge, correction = "isotropic")
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
              Ki <- Ki.type <-  Kcross(data.types.ppp.subset, i = cell.type[1], j = cell.type[2],
                                       r = radii, nlarge = nlarge, correction = "isotropic")
            }

            if (!homogeneous) {
              Ki <- Ki.type <-  Kcross.inhom(data.types.ppp.subset, i = cell.type[1], j = cell.type[2],
                                             r = radii, nlarge = nlarge, correction = "isotropic")
            }
          }

          if (!use.K) {

            if (homogeneous) {
              Ki <- Ki.type <-  Lcross(data.types.ppp.subset, i = cell.type[1], j = cell.type[2],
                                       r = radii, nlarge = nlarge, correction = "isotropic")
            }

            if (!homogeneous) {
              Ki <- Ki.type <-  Lcross.inhom(data.types.ppp.subset, i = cell.type[1], j = cell.type[2],
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


spot <- function(data, radii, outcome, censor = NULL, model.type, use.K = TRUE,
                                 homogeneous = TRUE, adjustments = NULL, marked = FALSE,
                                 cell.type = NULL, K.diff = FALSE, pick.roi = "all", seed = 12345,
                                 nlarge = 10000, prop.nonzero = 0.2, print.progress = TRUE,
                                 use.max.radius.only = FALSE) {

  # Save the arguments
  arguments <- match.call()

  # -------------------------------------
  # Calculate the K(r)s or L(r)s
  # -------------------------------------

  # Summarize the images within each PID
  pids.images <- data %>%
    dplyr::select(PID, id) %>%
    distinct()

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
  pids.images <- pids.images %>% filter(id %in% Kr.df$id)

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
      Kr.df <- Kr.df %>% filter(selected.ROI == 1)

      # Filter to the number of images remaining
      pids.images <- pids.images %>% filter(id %in% Kr.df$id)

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

  # Would we like to use just one radius?
  if (use.max.radius.only) {
    if (paste0("radius.", max(radii)) %in% radii.to.save) {
      radii.to.save <- paste0("radius.", max(radii))
    }

    if (!((paste0("radius.", max(radii)) %in% radii.to.save))) {
      stop(paste0("Maximum radius was removed because fewer than ", prop.nonzero*100, "% of K(r) or L(r) were non-zero at this radius."))
    }
  }

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
                                          distinct(),
                                        Kr.subset,
                                        by = "PID")

      if (model.type == "survival") {
        # Save the model equation
        model.equation <- reformulate(c("ripley", adjustments), "Surv(get(outcome), get(censor))")

        # Run the model
        res <- coxph(model.equation, data = outcome.data, control = coxph.control(iter.max = 100))
      }

      if (model.type == "logistic") {
        # Save the model equation
        model.equation <- reformulate(c("ripley", adjustments), "get(outcome)")

        # Run the model
        res <- glm(model.equation, family = binomial(), data = outcome.data)
      }

      # Were random effects used?
      # used.marginal.model <- FALSE

      # # Association testing
      # if (model.type == "survival") {
      #
      #   # Fit a Cox PH model without random effects
      #   if (!any(duplicated(pids.images$PID))) {
      #
      #     # Save the model equation
      #     model.equation <- reformulate(c("ripley", adjustments), "Surv(get(outcome), get(censor))")
      #
      #     # Run the model
      #     coxph.control(iter.max = 100)
      #     res <- coxph(model.equation, data = outcome.data)
      #   }
      #
      #   # Fit a Cox PH model with random effects
      #   if (any(duplicated(pids.images$PID))) {
      #
      #     # Change value of used.random.effects
      #     used.marginal.model <- TRUE
      #
      #     # Save the model equation
      #     # model.equation <- reformulate(c("ripley", adjustments, "(1|PID)"), "Surv(get(outcome), get(censor))")
      #     model.equation <- reformulate(c("ripley", adjustments), "Surv(get(outcome), get(censor))")
      #
      #     # Run the model
      #     # coxme.control(iter.max = 100)
      #     # res <- coxme(model.equation, data = outcome.data)
      #     res <- coxph(model.equation, data = outcome.data, cluster = PID)
      #   }
      # }
      #
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
  overall.pval <- overall.pval.rounding <- NA

  # Combine results with and without rounding
  if (!is.logical(pval.df$pval)) {
    pval.df$pval.round <- pval.df$pval
    pval.df$pval.round[pval.df$pval >= 0.5] <- plyr::round_any(pval.df$pval[pval.df$pval >= 0.5], 0.1, f = floor)

    # Combine results (no rounding)
    overall.pval <- ACAT::ACAT(pval.df$pval)
    overall.pval.rounding <- ACAT::ACAT(pval.df$pval.round)
  }

  # Return
  list(overall.pval = overall.pval,
       overall.pval.rounding = overall.pval.rounding,
       Kr.df = Kr.df,
       pval.df = pval.df,
       pids.images = pids.images,
       images.missing.cells = images.missing.cells,
       n.images.association = n.images.association,
       radii = radii,
       adjustments = adjustments,
       homogeneous = homogeneous,
       marked = marked,
       cell.type = cell.type,
       seed = seed,
       arguments = arguments)
}
