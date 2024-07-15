#=======================================================================================================
#
# Load the internal collection
#
#=======================================================================================================

SCSBrainCellTypeCollection = function(organism = "human", useHomology = TRUE, 
                      addOldOrganismToSetNames = FALSE, namePattern = ".convertedFrom.%o",
                      addOldOrganismToSetDescriptions = FALSE,                          
                      descriptionPattern = " (Converted from %o.)")
{
  x = load(system.file("extdata/SCSBrainCellTypeCollection.rda", mustWork = TRUE, 
            package = "SCSBrainCellTypeCollection"));
  if (x!="SCSBrainCellTypeCollection")
     stop("Internal error: incorrect file content. Sorry!");

  # Convert the collection to the output organism. Since convertGeneSet does check whether new organism is
  # different from the current for each set, there's little overhead if 'organism' is already the same as
  # the organisms for which the gene sets were saved.
  if (is.null(organism)) return(get(x));
  if (is.na(organism)) return(get(x));

  suppressWarnings(convertCollectionToOrganism(get(x), organism = organism, useHomology = useHomology,
                   addOldOrganismToSetNames = addOldOrganismToSetNames,
                   namePattern = namePattern,
                   addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
                   descriptionPattern = descriptionPattern));
}

