#=======================================================================================================
#
# Load the internal collection
#
#=======================================================================================================

HuntingtonDiseaseWGCNACollection = function(organism = "human", useHomology = TRUE, 
                      addOldOrganismToSetNames = FALSE, namePattern = ".convertedFrom.%o",
                      addOldOrganismToSetDescriptions = FALSE,                          
                      descriptionPattern = " (Converted from %o.)")
{
  internalCollection.JAM = NULL;
  x = load(system.file("extdata/WGCNA.HD.PL.collection.rda", mustWork = TRUE, 
            package = "HuntingtonDiseaseWGCNACollection"));
  if (x!="WGCNA.HD.PL.collection")
     stop("Internal error: incorrect file content. Sorry!");

  # Convert the collection to the output organism. Since convertGeneSet does check whether new organism is
  # different from the current for each set, there's little overhead if 'organism' is already the same as
  # the organisms for which the gene sets were saved.
  if (is.null(organism)) return(get(x));
  if (is.na(organism)) return(get(x));

  convertCollectionToOrganism(get(x), organism = organism, useHomology = useHomology,
                   addOldOrganismToSetNames = addOldOrganismToSetNames,
                   namePattern = namePattern,
                   addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
                   descriptionPattern = descriptionPattern);
}

