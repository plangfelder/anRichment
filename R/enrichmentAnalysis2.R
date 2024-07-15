#===================================================================================================
#
# mapCollectionIDs
#
#===================================================================================================

mapCollectionIDs = function(collection, new2old, keepOrder = TRUE)
{
  new = new2old[, 1];
  old = new2old[, 2];
  collection$dataSets = lapply(collection$dataSets, function(set) 
  {
    newIDs = new[old %in% set$data$Entrez];  # these may contain duplicates
    index1 = which(new %in% newIDs);  # this takes care of the duplicates
    index.old = match(old[index1], set$data$Entrez);
    if (keepOrder) 
    {
      order.new = order(index.old);
      newIDs = newIDs[order.new];
      index1 = index1[order.new];
      index.old = index.old[order.new];
    }
    set$data = data.frame(Entrez = newIDs, set$data[index.old, -1]);
    set;
  })
  collection;
}

#===================================================================================================
#
# Helper functions
#
#===================================================================================================


matchList = function(lst, ref, nomatch = NA_integer_, incomparables = NA, setAttributes = TRUE)
{
  lst2ref = match(unlist(lst), ref, nomatch = nomatch, incomparables = incomparables);
  out1 = tapply(lst2ref, rep(1:length(lst), sapply(lst, length)), as.integer);
  if (any(sapply(lst, length)==0))
  {
    out = vector(mode = "list", length = length(lst));
    out[sapply(lst, length) > 0] = out1;
  } else 
    out = out1;
  if (length(out)!=length(lst))
  {
    printFlush("matchList: mismatch of lengths of input and output. Dropping into browser.");
    browser();
  }
  if (setAttributes)
  {
    #tryCatch( {
    attributes(out) = attributes(lst);
    dim(out) = dim(lst);
    names(out) = names(lst); 
    #}, error = function(e) browser(e));
  }
  out;
}

.removeMissing = function(x) x[!is.na(x)]

#===================================================================================================
#
# enrichmentAnalysis
#
#===================================================================================================
# Main enrichment analysis function. 

enrichmentAnalysis.general = function(
               classLabels = NULL, 
               identifiers = NULL,
               classSuffixes = NULL,
               # Alternative input specification
               active = NULL,
               inactive = NULL,
               activeNames = names(active),
               # Optional map from active to inactive sets, may be useful to save memory
               active2inactive = NULL,

               # Reference collection: should be appropriate to the organism
               refCollection,

               # If data IDs and collections IDs aren't the same, a map from class/active to collection. 
               # The map can in principle be many-to-many.
               classID2setID = NULL,

               # optional sub-selection
               useGroups = NULL,

               # Specify what to use as the background list
               useBackground = c("intersection", "given", "reference", "custom"),
               customBackground = NULL,
               removeMissing = TRUE,

               # Included evidence levels
               useEvidence = "all",
               
               # Various calculation options
               removeDuplicatesInDifferentClasses = TRUE,
               ignoreLabels = NULL,
               alternative = "greater",  
                  ## usually only this makes sense for enrichment but sometimes one wants to know whether the observed
                  ## overlap differs from the expected one, not just whether it is greater. 

               # How much output to generate
               nBestDataSets = 10, 
               threshold = 0.05,
               thresholdType = c("Bonferroni", "FDR", "nominal"), 
              
               # Calculate multiple testing corrections? 
               getBonferroniCorrection = TRUE,
               getFDR = TRUE,

               # Return overlap features in the enrichment table?
               getOverlapIDs = TRUE,
               getOverlapSymbols = FALSE,
               ID2symbol = NULL,
               maxReportedOverlapGenes = 50,

               # Output options
               entrySeparator = "|",
               groupSeparator = entrySeparator,
               geneSeparator = entrySeparator,
               classColName = "class",

               # Return gene set details?
               getDataSetDetails = TRUE,

               # Diagnostic message options
               verbose = 1, indent = 0,

               # Undocumented arguments whose effect may change in the future
               ... )
{

   spaces = indentSpaces(indent);
   if (verbose > 0)
      printFlush(base::paste(spaces, "enrichmentAnalysis: preparing data.."));
   sAF = options("stringsAsFactors")
   options(stringsAsFactors = FALSE);
   on.exit(options(stringsAsFactors = sAF[[1]]), TRUE)

   spaces = indentSpaces(indent);

   useBackground = match.arg(useBackground);
   thresholdType = match.arg(thresholdType);

   if (!.isCollection(refCollection))
     stop("'refCollection' must be a valid collection.");

   organisms = lapply(refCollection$dataSets, getElement, "organism");
   same = sapply(organisms, function(o1, o2) sameOrganism(o1[1], o2[1]), organisms[[1]]);
   if (!all(same)) 
     stop("All data sets within 'refCollection' must correspond to the same organism.\n",
          "Function 'convertCollectionToOrganism' can be used to convert collection to an organism.");
   organism = organisms[[1]] [1];

   ignoreLabels = c(ignoreLabels, NA);

   if (length(useGroups)==0) 
   {
     useGroups = knownGroups(refCollection);
   } else
     refCollection = subsetCollection(refCollection, useGroups); 
   impliedGroups1 = impliedGroups(refCollection$groups);

   otherArgs = list(...);

   # Collect data set information and lists of IDs
   collectionIDLists = geneLists(refCollection, evidence = useEvidence, simplify = FALSE);
   dataSetIDs = sapply(refCollection$dataSets, getElement, "ID");
   dataSetNames = sapply(refCollection$dataSets, getElement, "name");
   shortSetNames = sapply(refCollection$dataSets, getElement, "shortName");
   dataSetGroups = sapply(refCollection$dataSets, function(ds) 
            base::paste(unlist(impliedGroups1[ds$groups]), collapse = groupSeparator) );

   allCollectionIDs = unique(unlist(collectionIDLists));
   if (!is.null(classID2setID))
   {
     if (is.null(dim(classID2setID)) )
        stop("'classID2setID' must be a 2-dimensional object with (at least) 2 columns.");
     class2setIDHasDuplicates = any(duplicated(classID2setID[, 1]));
     doMap = TRUE
   } else {
     doMap = FALSE;
     class2setIDHasDuplicates = FALSE;
   }

   if (!is.null(classID2setID))
   { 
     allCollectionIDs.mapped = classID2setID[classID2setID[, 2] %in% allCollectionIDs, 1] 
   } else
     allCollectionIDs.mapped = allCollectionIDs; 

   if (useBackground=="custom")
   {
     if (is.null(customBackground))
       stop("When 'useBackground' is 'custom', the custom background must be specified in 'customBackground'.")
     referenceBackground = customBackground;
   } else
     referenceBackground = allCollectionIDs.mapped;

   if (getOverlapSymbols && is.null(ID2symbol))
     stop("When 'getOverlapSymbols' is TRUE, 'ID2symbol translation table must be supplied.");

   if (getOverlapSymbols)
     backgroundSymbols = .translateUsingTable(referenceBackground, ID2symbol);

   # Check and prepare lists of active and inactive IDs

   if (!is.null(classLabels))
   {
     classLabels = as.matrix(classLabels);
     nActiveMetasets = ncol(classLabels);
     activeSizes.asGiven = lapply(1:nActiveMetasets, function(i)
     {
       out = table(classLabels[, i]);
       n1 = names(out);
       keep = which(!n1 %in% ignoreLabels)
       out = as.numeric(out)[keep];
       setNames(out, n1[keep]);
     })

     if (is.null(classSuffixes))
     {
       if (nActiveMetasets > 1) {
         classSuffixes = colnames(classLabels);
         if (is.null(classSuffixes))
            classSuffixes = spaste(" (label set ", prependZeros(1:nActiveMetasets), ")");
       } else 
         classSuffixes = "";
     }
     if (length(classSuffixes)!=nActiveMetasets)
       stop("When given, length of 'classSuffixes' must equal the number of label sets.");

     activeSizes.asGiven = unlist( .mymapply(function(lst, name) 
                                                setNames(lst, spaste(names(lst), name)), 
                                             activeSizes.asGiven, classSuffixes));
       
     if (removeDuplicatesInDifferentClasses)
        classLabels = as.matrix(apply(classLabels, 2, .removeDuplicatesWithDifferentLabels,
                                identifiers = identifiers, ignoredLabels = ignoreLabels, missingLabel = NA));
     active0 = apply(classLabels, 2, function(labels1)
     {
       out = tapply(identifiers, labels1, unique);  # Note to self: tapply removes missing 'labels1' automatically.
       out = out[ !names(out) %in% ignoreLabels ];
       out;
     });
     active = unlist( .mymapply(function(lst, name) setNames(lst, spaste(names(lst), name)), active0, classSuffixes), 
                      recursive = FALSE);
     activeNames = names(active);
     activeSizes.asGiven = activeSizes.asGiven[activeNames];

     if (ncol(classLabels)==1 && !any(is.na(classLabels))) 
     {
       inactive = list(identifiers);
       active2inactive = rep(1, length(active));
     } else {
       inactive = unlist(apply(classLabels, 2, function(labels1)
       {
         out1 = identifiers[!is.na(labels1)];
         list(out1);
       }), recursive = FALSE);
       active2inactive = unlist(lapply(1:ncol(classLabels), function(col) rep(col, length(active0[[col]]))));
     }
   } else {
     if (is.null(active))
        stop("Either 'classLabels' or 'active' must be supplied.");
     if (is.atomic(active)) active = list(active);
     activeSizes.asGiven = sapply(active, length);
   }
     
   nActive = length(active);
   if (is.null(activeNames))
     activeNames = spaste("active.", prependZeros(1:length(active), nchar(length(active))));

   active = lapply(active, .removeMissing);

   if (is.null(inactive))
   { 
     if (useBackground %in% c("intersection", "given")) 
     {
       stop(base::paste0("When 'active' list is specified, 'inactive' must be specified as well.\n",
                   "  Alternatively, you can set 'useBackground' to 'reference', 'allOrgGenes' or 'custom'\n",
                   "  BUT be aware that resulting enrichment statistics may be inflated.")); 
     } else
       inactive = referenceBackground;
   }
   privateInactive = !is.atomic(inactive)
   if (privateInactive)
   {
     if (is.null(active2inactive))
     {
       if (length(inactive)!=length(active)) 
         stop("When 'inactive' is a list, it must have the same length as 'active' or 'active2inactive' must ",
              "specify the mapping.");
       active2inactive = 1:length(active);
     } else {
       if (length(active2inactive)!=length(active))
         stop("When given, length of 'active2inactive' must equal length of 'active'.");
       if (any( replaceMissing(active2inactive < 1 | active2inactive > length(inactive), TRUE)))
         stop("When given all entries in 'active2inactive' must be between 1 and the number of components in 'inactive'.");
     }
   } else {
     inactive = list(inactive);
     active2inactive = rep(1, nActive);
   }

   nInactive = length(inactive);

   # Extend inactive to include all active genes within the active sets that refer to the same inactive. 
   ### Note: this makes a potentially somewhat restrictive assumption that the background for all active sets that refer
   ### to the same inactive is the same, namely the union of all the active sets and the referred to inactive.
   inactive.ext = lapply(1:nInactive, function(ina) 
       .removeMissing(unique(unlist( c(inactive[ina], active[active2inactive==ina])))));

   if (length(activeNames)!=length(active))
     stop("Length of 'activeNames' (", length(activeNames), ") must equal length of 'active' (", length(active), ").")

   #identIsInCollection = lapply(identifiers , `%in%`, referenceBackground);
    #### Note to self: identifiers may not be supplied and is currently not defined in all cases.
   #if (sum(unlist(identIsInCollection))==0)
   #  stop("There are no common identifiers in the supplied gene sets and the input 'identifiers'.\n",
   #       "   Please make sure the 'refCollection' and 'identifiers' are Entrez identifiers and\n",
   #       "   correspond to the same organism.");

   nDataSets = length(collectionIDLists)
   dataSetSizes = sapply(collectionIDLists, length);
   if (all(dataSetSizes==0))
     warning(immediate. = TRUE, 
             "Gene lists are empty (possibly because of restricting evidence).");

   useCompiledCode = FALSE;
   if ("useCompiledCode" %in% names(otherArgs)) useCompiledCode = otherArgs$useCompiledCode;

   # Calculate the necessary backgrounds and overlaps
   if (verbose > 0) printFlush(base::paste(spaces, " ..calculating overlaps.."));

   if (useCompiledCode)
   {
     # Convert all IDs to integer indices relative to maximalBackgound. The variables will have suffix .f (.factor).
     if (verbose > 1) printFlush(paste(spaces, "   (using compiled code)"));
     maximalBackground = .removeMissing(unique(c(unlist(inactive.ext), referenceBackground)))
     collectionIDLists.f = matchList(collectionIDLists, if (doMap) allCollectionIDs else maximalBackground)
     referenceBackground.f = match(referenceBackground, maximalBackground, incomparables = NA);
     active.f = matchList(active, maximalBackground);
     classBg.f = matchList(inactive.ext, maximalBackground);
     if (!is.null(classID2setID))
     { 
       classID2setID.ff = cbind( match(classID2setID[, 1], maximalBackground, incomparables = NA), 
                                 match(classID2setID[, 2], allCollectionIDs, incomparables = NA));
       stopifnot(all.equal(maximalBackground[classID2setID.ff[, 1]], allCollectionIDs[ classID2setID.ff[, 2]]))
     } else {
       classID2setID.ff = matrix(numeric(0), 0, 2);
     }
     iss = .intersectSizesForEnrichment(active.f, classBg.f, as.integer(active2inactive),
               collectionIDLists.f, as.logical(doMap),
               class2ref = classID2setID.ff, 
               class2setIDHasDuplicates = as.logical(class2setIDHasDuplicates),
               returnMappedSets = FALSE,
               bgTypeR = useBackground, mappedRefBg = as.integer(referenceBackground.f));
     setSizes = iss$effectiveSetSizes;
     rownames(setSizes) = dataSetIDs;
     classSizes = iss$effectiveClassSizes;
     names(classSizes) = activeNames;
     overlapSizes = iss$overlapSizes;
     effectiveBg = iss$effectiveBackgrounds;
     nBackgroundGenes = sapply(effectiveBg, length);
     active2bg = iss$class2effectiveBg + 1;
     nBg = length(effectiveBg);
     sets.mapped = NULL;
   } else {
     if (useBackground=="given") {
        effectiveBg = inactive.ext;
        active2bg = active2inactive;
     } else if (useBackground=="intersection") {
        effectiveBg = lapply(inactive.ext, intersect, referenceBackground);
        active2bg = active2inactive;
     } else if (useBackground %in% c("reference", "custom")) {
        effectiveBg = list(referenceBackground);
        active2bg = rep(1, nActive);
     } else 
        stop("Unhandled 'useBackground'.");
     nBg = length(effectiveBg);
     effectiveBg = lapply(effectiveBg, .removeMissing);
     classes.bg = .mymapply(match, active, effectiveBg[active2bg], incomparables = NA);
     classes.bg = lapply(classes.bg, .removeMissing)
     if (doMap)
     {
       if (class2setIDHasDuplicates) {
         collectionIDLists.flat = .indexedFlattenedList(collectionIDLists);
         mapIndex = classID2setID[, 2] %in% collectionIDLists.flat[[2]];
         sets.mapped = tapply(classID2setID[mapIndex, 1], collectionIDLists.flat[mapIndex], unique);
         sets.mapped = lapply(collectionIDLists, function(set) unique(classID2setID[ classID2setID[, 2] %in% set, 1]));
       } else 
         sets.mapped = lapply(collectionIDLists, function(set) classID2setID[ classID2setID[, 2] %in% set, 1]);
     } else 
       sets.mapped = collectionIDLists;
     setSizes = vector(mode = "list", length = nBg);
     overlapSizes = matrix(0, nDataSets, nActive);
     for (bg in 1:nBg) 
     {
         sets.bg = matchList(sets.mapped, effectiveBg[[bg]]);
         for (cl in which(active2bg==bg))
         {
           overlapSizes[, cl] = sapply(sets.bg, function(st) sum(st%in% classes.bg[[cl]]));
              ### FIXME: this could potentially be speeded up by using matchList 
              ### combined with sapply(x, function() sum(!is.na(x))
         }
         setSizes[[bg]] = sapply(sets.bg, function(x) sum(!is.na(x)));
     }
     setSizes = do.call(cbind, setSizes);
     classSizes = sapply(classes.bg, function(x) sum(!is.na(x))); 
     nBackgroundGenes = sapply(effectiveBg, length)
   }

   # Calculate enrichment p-values

   if (nBestDataSets > nDataSets) nBestDataSets = nDataSets;
   nComparisons = nActive * nDataSets;
   enrichmentIsValid = nBackgroundGenes > 0
   # these will hold the main results
   enrichment = matrix(1, nDataSets, nActive);
   if (any(nBackgroundGenes==0))
   {
     warning(immediate. = TRUE, 
        base::paste0("enrichmentAnalysis: number of effective background genes for some of the classes is zero.\n",
               "      This means that after restricting to the appropriate background type,\n",
               "      there are no genes left."));
     for (bg in which(nBackgroundGenes==0))
       enrichment[ , active2bg[[bg]]] = NA;
   }
  
   countsInDataSet = overlapSizes;
   enrichmentRatio.mat = matrix(0, nDataSets, nActive); ## a.k.a. effect size

   for (ac in 1:nActive)
   {
     bg1 = active2bg[ac];
     enrichmentRatio.mat[, ac] = countsInDataSet[, ac] * nBackgroundGenes[ bg1 ]/
                                   (setSizes[, bg1] * classSizes[ac]);
     enrichment[, ac] = .fisherPvalue( nCommon = overlapSizes[, ac],
                                      nInCat = setSizes[, bg1],
                                      nOutCat = nBackgroundGenes[bg1] - setSizes[, bg1],
                                      nInClass = rep(classSizes[ac], nDataSets), alternative = alternative);
   }

   # Calculate multiple testing corrections

   enrichment.Bonf = enrichment * nComparisons;
   enrichment.Bonf[ enrichment.Bonf > 1] = 1;
   #enrichment.FDR = try(qvalue(c(enrichment))$qvalues, silent = TRUE);
   enrichment.FDR = try(p.adjust(c(enrichment), method = "fdr"), silent = TRUE);
   if (inherits(enrichment.FDR, "try-error"))
   {
     warning(immediate. = TRUE,
         base::paste0("enrichmentAnalysis: FDR calculation failed.",
                if (thresholdType=="FDR") base::paste0("Will use nominal p-values\n", 
                    "                    to restrict returned terms.") else ""));
     if (thresholdType=="FDR") thresholdType = "nominal";
     enrichment.FDR = rep(NA, length(enrichment));
   } 
   dim(enrichment.FDR) = dim(enrichment);

   dimnames(enrichment) = dimnames(enrichment.Bonf) = dimnames(countsInDataSet) = 
       dimnames(enrichment.FDR) = dimnames(enrichmentRatio.mat) = list (dataSetIDs, activeNames);

   # Put together output tables and lists

   activeSizes = sapply(active, length);

   enrichment.threshold = switch(thresholdType, 
              nominal = enrichment,
              Bonferroni = enrichment.Bonf,
              FDR = enrichment.FDR);
   enrichmentTable = NULL;


   if (getDataSetDetails) dataSetDetails = list(); 
   if (any(enrichmentIsValid) && (!is.null(threshold) || nBestDataSets > 0))
   {
     if (verbose > 1) 
       printFlush(base::paste(spaces, "   ..putting together terms with highest enrichment significance.."));
     order = apply(enrichment, 2, order);
     dim(order) = dim(enrichment);

     for (cl in 1:nActive)
     {
        bg1 = active2bg[cl];
        if (getDataSetDetails) dataSetDetails[[cl]] = list();
        if (!is.null(threshold))
        {
           reportTerms = c(1:nDataSets)[enrichment.threshold[, cl] < threshold];
        } else 
           reportTerms = numeric(0);
        if (nBestDataSets > 0)
           reportTerms = unique(c(reportTerms, c(1:nDataSets)[order[1:nBestDataSets, cl]]));
        nRepTerms = length(reportTerms);
        if (nRepTerms > 0)
        {
           # Order the terms by increasing enrichment p-value
           reportTerms = reportTerms[ order(enrichment[ reportTerms, cl ]) ];

           if (getOverlapIDs | getOverlapSymbols | getDataSetDetails)
           {
              reportLists = collectionIDLists[reportTerms];
              if (doMap) {
                 if (is.null(sets.mapped))
                 {
                   reportLists = lapply(reportLists, function(set) unique(classID2setID[ classID2setID[, 2] %in% set, 1]));
                 } else
                   reportLists = sets.mapped[reportTerms]
              }
              overlapIDs0 = lapply(reportLists, intersect, intersect(active[[cl]], effectiveBg[[ bg1 ]]));
              if (getOverlapSymbols)
                 overlapSymbols0 = lapply(overlapIDs0, function(x)
                        backgroundSymbols[ match(x, referenceBackground, incomparables = NA)])
           }
           if (getOverlapIDs | getOverlapSymbols)
           {
              lengths = sapply(overlapIDs0, length);
              long = lengths>maxReportedOverlapGenes;
              overlapIDs = overlapIDs0;
              if (any(long)) for (i in which(long))
                overlapIDs[[i]] = paste0("(More than ", maxReportedOverlapGenes, " overlapping genes)");
              if (any(!long))
              {
                index = which(!long)
                if (getOverlapSymbols)
                {
                   if (getOverlapIDs) {
                     overlapIDs[index] = mapply(function(x, y) paste0(x, " (", y, ")"),
                                            overlapIDs[index], overlapSymbols0[index], SIMPLIFY = FALSE)
                   } else 
                     overlapIDs[index] = overlapSymbols0[index];
                }
                overlapIDs[index] = sapply(overlapIDs[index], function(s) base::paste(unique(s), collapse= geneSeparator));
             }
             overlapIDs = sapply(overlapIDs, identity);
           } else overlapIDs = NULL;
           expectedFrac = setSizes[reportTerms, bg1 ]/nBackgroundGenes[bg1]
           observedFrac = countsInDataSet[reportTerms, cl]/classSizes[cl];
           enrichmentRatio1 = observedFrac/expectedFrac;
           enrTab0 = list(class = rep(activeNames[cl], nRepTerms), 
                     rank = seq(length.out = nRepTerms),
                     dataSetID = dataSetIDs[reportTerms],
                     dataSetName = dataSetNames[reportTerms],
                     inGroups = dataSetGroups[reportTerms],
                     pValue = enrichment[reportTerms, cl],
                     Bonferroni = if (getBonferroniCorrection) enrichment.Bonf[reportTerms, cl] else NULL,
                     FDR = if (getFDR) enrichment.FDR[reportTerms, cl] else NULL,
                     nCommonGenes = countsInDataSet[reportTerms, cl],
                     fracOfEffectiveClassSize = observedFrac,
                     expectedFracOfEffectiveClassSize = expectedFrac,
                     enrichmentRatio = enrichmentRatio1,
                     classSize.asGiven = rep(activeSizes.asGiven[cl], nRepTerms),
                     validClassSize = rep(activeSizes[cl], nRepTerms),
                     effectiveClassSize = rep(classSizes[cl], nRepTerms), 
                     fracOfEffectiveSetSize = 
                          countsInDataSet[reportTerms, cl]/setSizes[reportTerms, bg1],
                     effectiveSetSize = setSizes[reportTerms, bg1],
                     shortDataSetName = shortSetNames[reportTerms],
                     overlapGenes = overlapIDs);
           l1 = sapply(enrTab0, length);
           #enrTab = as.data.frame(enrTab0[l1>0], row.names = enrTab0$rank);
           enrTab = try(as.data.frame(enrTab0[l1>0], row.names = enrTab0$rank));
           if (inherits(enrTab, "try-error")) browser();
           names(enrTab)[1] = classColName;
           rownames(enrTab) = NULL;
           if (is.null(enrichmentTable))
           {
              enrichmentTable = enrTab;
           } else 
              enrichmentTable = rbind(enrichmentTable, enrTab);

           if (getDataSetDetails)
           {
             for (rci in seq(length.out = nRepTerms))
             {
                gs = reportTerms[ rci ];
                dataSetDetails[[cl]][[rci]] = list(
                        dataSetID = enrTab$dataSetID[ rci ],
                        dataSetName = enrTab$dataSetName[rci],
                        dataSetDescription = refCollection$dataSets[[gs]]$description,
                        dataSetGroups = unlist(impliedGroups1[refCollection$dataSets[[gs]]$groups]),
                        enrichmentP = enrTab$pValue[rci],
                        commonGeneIDs = overlapIDs0[[rci]],
                        commonGeneSymbols = if (getOverlapSymbols) overlapSymbols0[[rci]] else NULL);
             }
             names(dataSetDetails[[cl]]) = enrTab$dataSetID;
             names(dataSetDetails)[cl] = activeNames[cl];
           }
        }
      }
   }
        
   list(enrichmentIsValid = enrichmentIsValid,
        enrichmentTable = enrichmentTable,
        pValues = enrichment,
        Bonferroni = if (getBonferroniCorrection) enrichment.Bonf else NULL,
        FDR = if (getFDR) enrichment.FDR else NULL,
        countsInDataSet = countsInDataSet, 
        enrichmentRatio = enrichmentRatio.mat,
        effectiveBackgroundSize = nBackgroundGenes,
        effectiveDataSetSizes = setSizes,
        effectiveClassSizes = classSizes,
        effectiveClass2bg = active2bg,
        dataSetDetails = if (getDataSetDetails) dataSetDetails else NULL); 
}

#===============================================================================================================
#
# Wrapper to emulate the old enrichmentAnalysis function
#
#===============================================================================================================

enrichmentAnalysis.Entrez = function(
               classLabels = NULL,
               identifiers = NULL,
               classSuffixes = NULL,
               # Alternative input specification
               active = NULL,
               inactive = NULL,
               activeNames = names(active),
               # Optional map from active to inactive sets, may be useful to save memory
               active2inactive = NULL,

               # Reference collection: should be appropriate to the organism
               refCollection,

               # optional sub-selection
               useGroups = NULL,

               # Specify what to use as the background list
               useBackground = c("intersection", "given", "reference", "allOrgGenes"),
               removeMissing = TRUE,

               # Included evidence levels
               useEvidence = "all",

               # Various calculation options
               removeDuplicatesInDifferentClasses = TRUE,
               ignoreLabels = NULL,
               alternative = "greater",

               # How much output to generate
               nBestDataSets = 10,
               threshold = 0.05,
               thresholdType = c("Bonferroni", "FDR", "nominal"),

               # Calculate multiple testing corrections? 
               getBonferroniCorrection = TRUE,
               getFDR = TRUE,

               # Return overlap features in the enrichment table?
               getOverlapEntrez = TRUE,
               getOverlapSymbols = FALSE,
               maxReportedOverlapGenes = 50,
               ID2symbol = NULL,

               # Output options
               entrySeparator = "|",
               groupSeparator = entrySeparator,
               geneSeparator = entrySeparator,
               classColName = "class",

               # Return gene set details?
               getDataSetDetails = TRUE,

               # Diagnostic message options
               verbose = 1, indent = 0,

               # Undocumented arguments whose effect may change in the future
               ... )
{

   organisms = lapply(refCollection$dataSets, getElement, "organism");
   same = sapply(organisms, function(o1, o2) sameOrganism(o1[1], o2[1]), organisms[[1]]);
   if (!all(same))
     stop("All data sets within 'refCollection' must correspond to the same organism.\n",
          "Function 'convertCollectionToOrganism' can be used to convert collection to an organism.");

   organism = organisms[[1]] [1];

  if (getOverlapSymbols && is.null(ID2symbol)) 
  {
    ID2symbol = geneAnnotationFromEntrez(allOrgGenes(organism = organism), organism = organism)[c("Entrez", "Symbol")];
  }


  enrichmentAnalysis.general(
               classLabels,
               identifiers,
               classSuffixes = classSuffixes,
               # Alternative input specification
               active = active,
               inactive = inactive,
               active2inactive = active2inactive,
               activeNames = activeNames,

               refCollection = refCollection,

               useGroups = useGroups,

               useBackground = useBackground,
               removeMissing = removeMissing,
               customBackground = if (useBackground=="allOrgGenes") allOrgGenes(organism) else NULL,

               useEvidence = useEvidence,

               removeDuplicatesInDifferentClasses = removeDuplicatesInDifferentClasses,
               ignoreLabels = ignoreLabels,

               nBestDataSets = nBestDataSets,
               threshold = threshold,
               thresholdType = thresholdType,

               getBonferroniCorrection = getBonferroniCorrection,
               getFDR = getFDR,

               getOverlapIDs = getOverlapEntrez,
               getOverlapSymbols = getOverlapSymbols,
               ID2symbol = ID2symbol,
               maxReportedOverlapGenes = maxReportedOverlapGenes,

               entrySeparator = entrySeparator,
               groupSeparator = groupSeparator,
               geneSeparator = geneSeparator,
               classColName = classColName,
               getDataSetDetails = getDataSetDetails,
               verbose = verbose, indent = indent, ... )
}


