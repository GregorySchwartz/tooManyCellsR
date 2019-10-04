#' Write a Matrix to a folder.
#'
#' This function will write a Matrix from the Matrix library to a temporary
#' directory containing matrix.mtx, genes.tsv, barcodes.tsv, and optionally a
#' labels.csv file.
#' @param mat The input Matrix with gene row names and cell barcode column
#'   names.
#' @param labels The input labels data frame with item (cell barcodes) and label
#'   (whatever labels you want to give them, such as tissue of origin, celltype,
#'   etc.) columns. Optional.
#' @return None
#' @export
#' @examples
#' input <- system.file("extdata", "mat.csv", package="TooManyCellsR")
#' df = read.csv(input, row.names = 1, header = TRUE)
#' mat = Matrix::Matrix(as.matrix(df), sparse = TRUE)
#' writeMatrixFiles(mat)

writeMatrixFiles = function(mat, labels = NULL) {

  # Output directory.
  tmp = tempdir()

  # Expression matrix.
  Matrix::writeMM(mat, paste0(tmp, "/matrix.mtx"))

  # Data frame of genes.
  utils::write.table(data.frame(x = rownames(mat), y = rownames(mat))
            , file = paste0(tmp, "/genes.tsv")
            , sep = "\t"
            , row.names = FALSE
            , col.names = FALSE
            , quote = FALSE
              )

  # Data frame of cell barcodes.
  utils::write.table(colnames(mat)
            , file = paste0(tmp, "/barcodes.tsv")
            , sep = "\t"
            , row.names = FALSE
            , col.names = FALSE
            , quote = FALSE
              )

  # Write labels file.
  if(!(is.null(labels))) {
    utils::write.table(labels
              , file = paste0(tmp, "/labels.csv")
              , sep = ","
              , row.names = FALSE
              , col.names = TRUE
              , quote = FALSE
                )
  }

  return(tmp)

}

#' Different error for importing data.
#'
#' This function will fail gracefully instead of stopping the program with an
#' error.
#' @param f The function to use.
#' @param file The input file to be read.
#' @return The imported data frame or NULL if an error occurred.
#' @export
#' @examples
#' input <- system.file("extdata", "mat.csv", package="TooManyCellsR")
#' fail = tryFunc(read.csv, "fail.csv")
#' fail
#' success = tryFunc(read.csv, input)
#' success
tryFunc = function(f, file) {
  return(tryCatch(f(file)
                , error = function(e) { print(paste0(file, " not found, ignoring import.")) }
                , warning = function(w) {}
                  )
         )
}

#' Import some 'too-many-cells make-tree' outputs into a data frame.
#'
#' This function will import some of the files resulting from a 'too-many-cells
#' make-tree' run into R as data frames. Does not include cluster list. Look at
#' the main tooManyCells function for the cluster list.
#' @param dir The output directory of a 'too-many-cells' run.
#' @return A list of each output. Reads the following files, see
#'   \url{https://gregoryschwartz.github.io/too-many-cells/} for more details:
#'   "dendrogram.svg", "clumpiness.pdf", "projection.pdf",
#'   "label_projection.pdf", "clumpiness.csv", "cluster_info.csv",
#'   "node_info.csv", and "cluster_diversity.csv".
#' @export
#' @examples
#' input <- system.file("extdata", "mat.csv", package="TooManyCellsR")
#' inputLabels <- system.file("extdata", "labels.csv", package="TooManyCellsR")
#' df = read.csv(input, row.names = 1, header = TRUE)
#' mat = Matrix::Matrix(as.matrix(df), sparse = TRUE)
#' labelsDf = read.csv(inputLabels, header = TRUE)
#' # Here we draw this small toy example with no filter or normalization, and
#' # decrease the size of the branches and increase the size of the leaf nodes.
#' # With non-toy real world single cell data, these options should not be
#' # necessary.
#' \dontrun{
#' tooManyCells( mat, labels = labelsDf
#'             , args = c( "make-tree"
#'                       , "--no-filter"
#'                       , "--normalization", "NoneNorm"
#'                       , "--draw-max-node-size", "40"
#'                       , "--draw-max-leaf-node-size", "70"
#'                       )
#'             )
#' res = importResults("out")
#' plot(res$treePlot, axes = FALSE)
#' }

importResults = function(dir = "out") {

  treePlotRes = tryFunc(imager::load.image, paste0(dir, "/dendrogram.svg"))
  clumpinessPlotRes = tryFunc(imager::load.image, paste0(dir, "/clumpiness.pdf"))
  projectionPlotRes = tryFunc(imager::load.image, paste0(dir, "/projection.pdf"))
  labelProjectionPlotRes = tryFunc(imager::load.image, paste0(dir, "/label_projection.pdf"))

  clumpinessRes = tryFunc(utils::read.csv, paste0(dir, "/clumpiness.csv"))
  clusterInfoRes = tryFunc(utils::read.csv, paste0(dir, "/cluster_info.csv"))
  nodeInfoRes = tryFunc(utils::read.csv, paste0(dir, "/node_info.csv"))
  clusterDiversityRes = tryFunc(utils::read.csv, paste0(dir, "/cluster_diversity.csv"))

  res = list( treePlot            = treePlotRes
            , clumpinessPlot      = clumpinessPlotRes
            , projectionPlot      = projectionPlotRes
            , labelProjectionPlot = labelProjectionPlotRes
            , clumpiness          = clumpinessRes
            , clusterInfo         = clusterInfoRes
            , nodeInfo            = nodeInfoRes
            , clusterDiversity    = clusterDiversityRes
            )

  return(res)

}

#' Execute 'too-many-cells'.
#'
#' This function will run 'too-many-cells' on a Matrix.
#'   Requires 'TooManyCells' to be installed (follow instructions at
#'   \url{https://gregoryschwartz.github.io/too-many-cells/} ).
#' @param mat The input Matrix with gene row names and
#'   cell barcode column names.
#' @param args The arguments to give to the command line program. See
#'   \url{https://gregoryschwartz.github.io/too-many-cells/} for more information.
#'   Defaults to "make-tree".
#' @param labels The input labels data frame with item (cell barcodes) and label
#'   (whatever labels you want to give them, such as tissue of origin, celltype,
#'   etc.) columns. Optional.
#' @param output The output folder for the 'too-many-cells' process. Defaults to
#'   "out".
#' @param prior The location of the tree that was already made (previous
#'   'too-many-cells' output) so quick visual or pruning changes can be made
#'   without remaking the tree (can potentially save hours).
#' @param docker If using 'too-many-cells' with docker, use this argument as the
#'   command to call. For instance, if version 0.2.1.0 was pulled from Docker
#'   Hub, set to "gregoryschwartz/too-many-cells:0.2.1.0".
#' @param mounts Additional directories to mount if needed for docker.
#'   The 'prior' argument will automatically mount if specified.
#' @return A list of each output, including the stdout. Reads the
#'   following files, see \url{https://gregoryschwartz.github.io/too-many-cells/} for
#'   more details: "dendrogram.svg", "clumpiness.pdf", "projection.pdf",
#'   "label_projection.pdf", "clumpiness.csv", "cluster_info.csv",
#'   "node_info.csv", and "cluster_diversity.csv".
#' @export
#' @examples
#' input <- system.file("extdata", "mat.csv", package="TooManyCellsR")
#' inputLabels <- system.file("extdata", "labels.csv", package="TooManyCellsR")
#' df = read.csv(input, row.names = 1, header = TRUE)
#' mat = Matrix::Matrix(as.matrix(df), sparse = TRUE)
#' labelsDf = read.csv(inputLabels, header = TRUE)
#' # Here we draw this small toy example with no filter or normalization, and
#' # decrease the size of the branches and increase the size of the leaf nodes.
#' # With non-toy real world single cell data, these options should not be
#' # necessary.
#' \dontrun{
#' res = tooManyCells( mat, labels = labelsDf
#'                   , args = c( "make-tree"
#'                             , "--no-filter"
#'                             , "--normalization", "NoneNorm"
#'                             , "--draw-max-node-size", "40"
#'                             , "--draw-max-leaf-node-size", "70"
#'                             )
#'                   )
#' plot(res$treePlot, axes = FALSE)
#' res$stdout
#' res$nodeInfo
#' }

tooManyCells = function( mat
                       , args = c("make-tree")
                       , labels = NULL
                       , output = "out"
                       , prior = NULL
                       , docker = NULL
                       , mounts = c()
                       ) {

  # Convert to absolute paths.
  dir.create(output, recursive = TRUE)
  output = normalizePath(output)

  # Write matrix.
  tmpdir = writeMatrixFiles(mat, labels)

  # Determine input matrix location.
  autoArgs = c("--matrix-path", tmpdir)

  # Determine the output location.
  autoArgs = c(autoArgs, "--output", output)

  # Determine if labels are needed.
  if(!is.null(labels)) {
    autoArgs = c(autoArgs, "--labels-file", paste0(tmpdir, "/labels.csv"))
  }

  # Determine if output already exists.
  if(!is.null(prior)) {
    prior = normalizePath(prior)
    autoArgs = c(autoArgs, "--prior", prior)
  }

  # Run 'too-many-cells'.
  stdout = if(is.null(docker)) {
             system2("too-many-cells", args = c(args, autoArgs), stdout = TRUE)
           } else {
             dockerArgs = c("run", "-i", "--rm"
                           , "-v" , paste0(tmpdir, ":", tmpdir)
                           , "-v", paste0(output, ":", output)
                           , if(is.null(prior)){ NULL } else { c("-v", paste0(dirname(prior), ":", dirname(prior))) }
                           , unlist(sapply(mounts, function(mount) c("-v", paste0(mount, ":", mount))))
                           , docker
                           )
             system2("docker", args = c(dockerArgs, args, autoArgs), stdout = TRUE)
           }

  # Get the output as a data frame.
  stdoutDf = utils::read.csv(text = stdout)

  # Plot tree if needed.
  if("make-tree" %in% args) {
    p = imager::load.image(paste0(output, "/dendrogram.svg"))
    graphics::plot(p, axes = FALSE)
  }

  # Clean up input files.
  suppressWarnings(file.remove(c( paste0(tmpdir, c("/matrix.mtx"))
                                , paste0(tmpdir, c("/barcodes.tsv"))
                                , paste0(tmpdir, c("/genes.tsv"))
                                , paste0(tmpdir, c("/labels.csv"))
                                )))

  # Import results.
  res = importResults(output)

  return(append(list(stdout = stdoutDf), res))

}
