
#=====================================================================================================
#
# Retrieve Illumina 450k annotation
#
#=====================================================================================================

getIllumina450kAnnotation = function(addEntrez = TRUE)
{
  env = new.env();
  x = load(system.file("extdata", "Illumina450kManifest.rda", package = "anRichment"),
           envir = env);
  if (x!="Illumina450kManifest.short") 
    stop("Internal error: inconsistency in Illumina annotation content. Sorry!")

  out = get(x, envir = env);      

  entrez = entrezFromMultipleSymbols(out$UCSC_RefGene_Name);
  cbind(out, Entrez = entrez);
}

