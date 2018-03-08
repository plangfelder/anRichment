#===========================================================================================
#
# userListEnrichment adapted to anRichment
#
#==========================================================================================

userListEnrichment = function(
    # Input gene classes (modules etc)
    geneR,
    labelR,

    # Input options
    omitCategories = "grey",
    inputIsEntrez = NULL,
    organism = "human",

    # Optional user-supplied gene lists
    fnIn = NULL, catNmIn = fnIn,

    # Which internal gene sets to use
    useBrainLists = FALSE,
    useBloodAtlases = FALSE,
    useStemCellLists = FALSE,
    useBrainRegionMarkers = FALSE,
    useImmunePathwayLists = FALSE,
    usePalazzoloWang = FALSE,

    # Output options
    outputGenes = FALSE,
    outputEntrez = NULL,
    outputCorrectedPvalues = TRUE,
    minGenesInCategory = 1,
    nameOut = "enrichment.csv"
)
{

  if (is.null(inputIsEntrez))
  {
    if (is.numeric(geneR)) 
    {
      inputIsEntrez = TRUE
    } else {
      tmp = suppressWarnings(as.numeric(geneR));
      if (sum(is.na(tmp)) > length(tmp)/2) inputIsEntrez = FALSE else inputIsEntrez = TRUE;
    }
  } 

  if (is.null(outputEntrez)) outputEntrez = inputIsEntrez;

  
  if (length(geneR) != length(labelR))
     stop("geneR and labelR must have same number of elements.")
  if (length(catNmIn) < length(fnIn)) {
     catNmIn = c(catNmIn, fnIn[(length(catNmIn) + 1):length(fnIn)])
     write("WARNING: not enough category names.  \n\t\t\t   Naming remaining categories with file names.",
            "")
    }
  if (is.null(fnIn) & (! (useBrainLists | useBloodAtlases | useStemCellLists | useBrainRegionMarkers
                              | useImmunePathwayLists | usePalazzoloWang)) )
        stop("Either enter user-defined lists or set one of the use_____ parameters to TRUE.")


   glIn = NULL
   if (length(fnIn)>0)
   {
      for (i in 1:length(fnIn))
      {
        ext = substr(fnIn[i], nchar(fnIn[i]) - 2, nchar(fnIn[i]))
        if (ext == "csv") {
            datIn = read.csv(fnIn[i])
            if (colnames(datIn)[2] == "Gene") {
                datIn = datIn[, 2:3]
            }
            else {
                datIn = datIn[, 1:2]
            }
        }
        else {
            datIn = scan(fnIn[i], what = "character", sep = "\n")
            datIn = cbind(datIn[2:length(datIn)], datIn[1])
        }
        colnames(datIn) = c("Gene", "Category")
        datIn[, 2] = paste(datIn[, 2], catNmIn[i], sep = "__")
        glIn = rbind(glIn, datIn)
      }
      if (inputIsEntrez)
      {
        entrez = as.numeric(as.character(glIn[, 1]));
      } else
        entrez = convert2entrez(organism = organism, symbol = glIn[, 1], ignoreCase = TRUE)
 
      userCollection = collectionFromGeneLists(entrez, setNames = glIn[, 2], organism = organism)
  } else
    userCollection = NULL;

  geneEntrez.all = if (inputIsEntrez) as.numeric(as.character(geneR)) else 
                     convert2entrez(organism = organism, symbol = geneR, ignoreCase = TRUE);

  fin = is.finite(geneEntrez.all)

  if (all(!fin)) stop("None of the gene symbols could be converted to Entrez IDs.");

  geneEntrez = geneEntrez.all[fin];
  labels = labelR[fin];

  
  internalColl = internalCollection(organism = organism);

  searchTerms = character(0);
  if (useBrainLists) searchTerms = c(searchTerms, "BrainLists");
  if (useBloodAtlases) searchTerms = c(searchTerms, "BloodAtlases");
  if (useStemCellLists) searchTerms = c(searchTerms, "StemCellLists");
  if (useBrainRegionMarkers) searchTerms = c(searchTerms, "BrainRegionMarkers");
  if (useImmunePathwayLists) searchTerms = c(searchTerms, "ImmunePathways");
  if (usePalazzoloWang) searchTerms = c(searchTerms, "PWLists");

  
  if (length(searchTerms)>0)
  {
    iColl.use = subsetCollection(internalColl, tags = searchTerms);
  } else
    iColl.use = NULL;

  if (is.null(iColl.use) && is.null(userCollection))
     stop("Either enter user-defined lists or set one of the use_____ parameters to TRUE.");

  # mergeCollections automatically drops NULL arguments.
  collection = mergeCollections(iColl.use, userCollection);

  geneListsIn.all = geneLists(collection);

  # This code is taken from Jeremy Miller's userListEnrichment

    uniqueEntrez = sort(unique(geneEntrez))
    nIn.all = length(uniqueEntrez);

  # Shorten geneListsIn to remove lists that don't overlap with uniqueEntrez

    overlapSizes = lapply(geneListsIn.all, function(l1, l2) length(intersect(l1, l2)), uniqueEntrez);
    keepInLists = overlapSizes > 0
    nIn.all = length(geneListsIn.all);
    geneListsIn = geneListsIn.all[keepInLists];
    catsIn = dataSetNames(collection)[keepInLists];

    typeIn = sapply(collection$dataSets[keepInLists], 
       function(ds) 
       {
          l = length(ds$groups);
          if (l>1) ds$groups[2] else if (l==1) ds$groups[1] else "Uncategorized";
       });

    catsR = sort(unique(labelR))
    omitCategories = c(omitCategories, "background")      
    catsR = catsR[!is.element(catsR, omitCategories)]
    lenAll = length(uniqueEntrez)
    nCols.pValues = 5;
    nIn = length(geneListsIn)
    nR = length(catsR);
    nComparisons = nR * nIn
    nComparisons.all = nR*nIn.all
    index = 1; 
    nOverlap = rep(0, nComparisons);
    pValues = rep(1, nComparisons);
    ovGenes = vector(mode = "list", length = nComparisons);
    isI = matrix(FALSE, lenAll, nIn);
    for (i in 1:nIn)
    {
      isI[, i] = is.element(uniqueEntrez, geneListsIn[[i]]);
    }
    for (r in 1:length(catsR)) 
    {
      isR  = is.element(uniqueEntrez,geneEntrez[(labelR == catsR[r])])
      for (i in 1:length(catsIn)) 
      { 
        isI.1 = isI[, i];
        lyn  = sum(isR&(!isI.1))
        lny  = sum(isI.1&(!isR))
        lyy  = sum(isR&isI.1)
        gyy  = uniqueEntrez[isR&isI.1]
        lnn  = lenAll - lyy - lyn - lny
        pv   = fisher.test(matrix(c(lnn,lny,lyn,lyy), 2, 2), alternative = "greater")$p.value
        nOverlap[index] = lyy;
        pValues[index] = pv;
        ovGenes[[index]] = gyy;
        index = index + 1
      }
    }

    # If needed, convert output entrez identifiers to gene symbols
    if (!outputEntrez)
    {
       # Add a missing value to each of the overlap lists to make sure all are at least length 1
       ovGenes.x = lapply(ovGenes, c, NA);
       ov.entrez = unlist(ovGenes.x);
       ov.symbols = convert2symbol(entrez = ov.entrez, organism = organism);
       index = unlist(mapply(function(i, l) rep(i, length(l)), 
                      1:length(ovGenes), ovGenes.x, SIMPLIFY = FALSE))
       ovGenes = tapply(ov.symbols, index, function(x) sort(x[!is.na(x)]));
    }

    results = list(pValues = data.frame(InputCategories = rep(catsR, rep(nIn, nR)),
                                        UserDefinedCategories = rep(catsIn, nR),
                                        Type = rep(typeIn, nR),
                                        NumOverlap = nOverlap,
                                        Pvalues = pValues,
                                        CorrectedPvalues = ifelse(pValues * nComparisons > 1, 1,
                                                                  pValues * nComparisons)),
                   ovGenes = ovGenes);
    namesOv = paste(results$pValues$InputCategories, "--", results$pValues$UserDefinedCategories);
    names(results$ovGenes) = namesOv
    if (outputCorrectedPvalues) {
        results$sigOverlaps = results$pValues[results$pValues$CorrectedPvalues < 0.05, c(1, 2, 3, 6)]
    } else {
        results$sigOverlaps = results$pValues[results$pValues$Pvalues < 0.05, c(1, 2, 3, 5)]
        write("Note that outputted p-values are not corrected for multiple comparisons.", 
            "")
    }
    results$sigOverlaps = results$sigOverlaps[order(results$sigOverlaps[, 4]), ];
    row.names(results$sigOverlaps) = NULL;
    
    rSig  = results$sigOverlaps
    nSig = nrow(rSig);
    if (nSig > 0)
    {
       rCats = paste(rSig$InputCategories,"--",rSig$UserDefinedCategories)
       rNums <- rep(0, nSig);
       rGenes <- rep("", nSig);
       for (i in 1:nSig)
       {
         rGn    = results$ovGenes[[which(names(results$ovGenes)==rCats[i])]]
         rNums[i]  = length(rGn);
         rGenes[i] = paste(rGn,collapse=", ");
       }
       rSig$NumGenes = rNums
       rSig$CategoryGenes = rGenes
       rSig = rSig[rSig$NumGenes>=minGenesInCategory,]
       if(!outputGenes) rSig = rSig[,1:4]
       results$sigOverlaps = rSig
       write.csv(results$sigOverlaps, file = nameOut, row.names = FALSE)
    }
 
    write(paste(length(namesOv), "comparisons were successfully performed."), 
        "")
    return(results)
}
  

  

