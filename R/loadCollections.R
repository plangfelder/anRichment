#=======================================================================================================
#
# Utility functions
#
#=======================================================================================================
.mymapply = function(..., SIMPLIFY = FALSE) base::mapply(..., SIMPLIFY = SIMPLIFY)

.loadAsList = function(file)
{
  env = new.env();
  load(file = file, envir = env);
  as.list(env, all.names = TRUE);
};

.capitalize = function(s)
{
  if (substring(s, 1, 1)%in% letters) substring(s, 1, 1) = LETTERS[ match(substring(s, 1, 1), letters)];
  s;
}

#=======================================================================================================
#
# Loading of pre-made collections
#
#=======================================================================================================

.loadCollection = function(file, organism = "human", useHomology = TRUE,
                      addOldOrganismToSetNames = FALSE, namePattern = ".convertedFrom.%o",
                      addOldOrganismToSetDescriptions = FALSE,
                      descriptionPattern = " (Converted from %o.)")
{
  organism = organismLabels(organism)[1];  # This ensures organism is valid...
  collection = .loadAsList(system.file(file, mustWork = TRUE, package = "anRichment"))[[1]]

  # Convert the collection to the output organism. Since convertGeneSet does check whether new organism is
  # different from the current for each set, there's little overhead if 'organism' is already the same as
  # the organisms for which the gene sets were saved.
  if (length(organism)==0) return(collection);
  if (is.na(organism)) return(collection);

  convertCollectionToOrganism(collection, organism = organism, useHomology = useHomology,
                   addOldOrganismToSetNames = addOldOrganismToSetNames,
                   namePattern = namePattern,
                   addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
                   descriptionPattern = descriptionPattern);
}

#=======================================================================================================
#
# Loading of pre-made collections
#
#=======================================================================================================

internalCollection = function(organism = "human", useHomology = TRUE, 
                      addOldOrganismToSetNames = FALSE, namePattern = ".convertedFrom.%o",
                      addOldOrganismToSetDescriptions = FALSE,                          
                      descriptionPattern = " (Converted from %o.)",
                      trimmed = TRUE, ...)
{
  .loadCollection( 
     file = if (trimmed) "extdata/trimmedInternalCollection.rda" else 
                         "extdata/internalCollection.JAM.rda", 
     organism = organism, useHomology = useHomology,
     addOldOrganismToSetNames = addOldOrganismToSetNames,
     namePattern = namePattern,
     addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
     descriptionPattern = descriptionPattern);
}

HDSigDBCollection = function(organism = "human", useHomology = TRUE,
                      addOldOrganismToSetNames = FALSE, namePattern = ".convertedFrom.%o",
                      addOldOrganismToSetDescriptions = FALSE,
                      descriptionPattern = " (Converted from %o.)", ...)
{
  organism = organismLabels(organism)[1];
  if (organism %in% c("human", "mouse"))
  {
     fileOrg = organism
  } else 
     fileOrg = "human";

  .loadCollection(file = WGCNA::spaste("extdata/HDSigDB-", fileOrg, ".rda"),
        organism = organism, useHomology = useHomology,
        addOldOrganismToSetNames = addOldOrganismToSetNames,
        namePattern = namePattern,
        addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
        descriptionPattern = descriptionPattern);
}


HuntingtonsDiseaseGeneExpressionCollection = function(organism = "human", useHomology = TRUE,
                      addOldOrganismToSetNames = FALSE, namePattern = ".convertedFrom.%o",
                      addOldOrganismToSetDescriptions = FALSE,
                      descriptionPattern = " (Converted from %o.)")
{
  HDSigDBCollection(organism = organism, useHomology = useHomology, 
          addOldOrganismToSetNames = addOldOrganismToSetNames,
          namePattern = namePattern,
          addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
          descriptionPattern = descriptionPattern)
}



BioSystemsCollection = function(organism = "human", ...)
{
  org.ext = organismLabels(organism)[2];
  .loadAsList(system.file(paste0("extdata/BioSystemsCollection-", make.names(.capitalize(tolower(org.ext))), ".rda"),
           mustWork = TRUE, package = "anRichment"))[[1]];
}

MillerAIBSCollection = function(organism = "human", useHomology = TRUE,
                      addOldOrganismToSetNames = FALSE, namePattern = ".convertedFrom.%o",
                      addOldOrganismToSetDescriptions = FALSE,
                      descriptionPattern = " (Converted from %o.)", ...)
{
  .loadCollection(file = "extdata/MillerAIBSCollection.rda",
        organism = organism, useHomology = useHomology,
        addOldOrganismToSetNames = addOldOrganismToSetNames,
        namePattern = namePattern,
        addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
        descriptionPattern = descriptionPattern);
}

YangLiteratureCollection = function(organism = "human", useHomology = TRUE,
                      addOldOrganismToSetNames = FALSE, namePattern = ".convertedFrom.%o",
                      addOldOrganismToSetDescriptions = FALSE,
                      descriptionPattern = " (Converted from %o.)", ...)
{
  .loadCollection(file = "extdata/extendedCustomCollection.rda",
        organism = organism, useHomology = useHomology,
        addOldOrganismToSetNames = addOldOrganismToSetNames,
        namePattern = namePattern,
        addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
        descriptionPattern = descriptionPattern);
}

HDTargetDBCollection = function(organism = "human", useHomology = TRUE,
                      addOldOrganismToSetNames = FALSE, namePattern = ".convertedFrom.%o",
                      addOldOrganismToSetDescriptions = FALSE,
                      descriptionPattern = " (Converted from %o.)",
                      traceable = TRUE, ...)
{
  .loadCollection(file = spaste("extdata/HDTargetDB.collection.", if (traceable) "traceable." else "", "rda"),
        organism = organism, useHomology = useHomology,
        addOldOrganismToSetNames = addOldOrganismToSetNames,
        namePattern = namePattern,
        addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
        descriptionPattern = descriptionPattern);
}

PhenopediaCollection = function(organism = "human", useHomology = TRUE,
                      addOldOrganismToSetNames = FALSE, namePattern = ".convertedFrom.%o",
                      addOldOrganismToSetDescriptions = FALSE,
                      descriptionPattern = " (Converted from %o.)",
                      ...)
{
  .loadCollection(file = "extdata/PhenopediaCollection.rda",
        organism = organism, useHomology = useHomology,
        addOldOrganismToSetNames = addOldOrganismToSetNames,
        namePattern = namePattern,
        addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
        descriptionPattern = descriptionPattern);
}


allCollections = function(
  merge = FALSE,
  organism = "human", useHomology = TRUE,
  addOldOrganismToSetNames = FALSE, namePattern = ".convertedFrom.%o",
  addOldOrganismToSetDescriptions = FALSE,
  descriptionPattern = " (Converted from %o.)",

  trimInternal = TRUE,
  traceableHDTargetDB = TRUE,

  buildExternal = FALSE,
  genomicSpacings = 5e6,

  MSigDBxml)
{
  collections = list(
    YangLiterature = YangLiteratureCollection(organism = organism, useHomology= useHomology,
        addOldOrganismToSetNames = addOldOrganismToSetNames,
        namePattern = namePattern,
        addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
        descriptionPattern = descriptionPattern,
        trimmed = trimInternal),
    HDSigDB = HDSigDBCollection(organism = organism, useHomology= useHomology,
        addOldOrganismToSetNames = addOldOrganismToSetNames,
        namePattern = namePattern,
        addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
        descriptionPattern = descriptionPattern),
    HDTargetDB = HDTargetDBCollection(organism = organism, useHomology= useHomology,
        addOldOrganismToSetNames = addOldOrganismToSetNames,
        namePattern = namePattern,
        addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
        descriptionPattern = descriptionPattern,
        traceable = traceableHDTargetDB),
    internal = internalCollection(organism = organism, useHomology= useHomology,
        addOldOrganismToSetNames = addOldOrganismToSetNames,
        namePattern = namePattern,
        addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
        descriptionPattern = descriptionPattern, 
        trimmed = trimInternal),
    MillerAIBS = MillerAIBSCollection(organism = organism, useHomology= useHomology,
        addOldOrganismToSetNames = addOldOrganismToSetNames,
        namePattern = namePattern,
        addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
        descriptionPattern = descriptionPattern,
        trim = trimInternal),
    BioSystems = BioSystemsCollection(organism = organism),
    Phenopedia = PhenopediaCollection(organism = organism, useHomology= useHomology,
        addOldOrganismToSetNames = addOldOrganismToSetNames,
        namePattern = namePattern,
        addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
        descriptionPattern = descriptionPattern));

  if (buildExternal)
  {
    collections = c(collections, 
      list(GO = buildGOcollection(organism = organism),
           MSigDB = if (!is.null(MSigDBxml)) 
             MSigDBCollection(MSigDBxml, organism = organism, useHomology= useHomology,
              addOldOrganismToSetNames = addOldOrganismToSetNames,
              namePattern = namePattern,
              addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
              descriptionPattern = descriptionPattern) else NULL,
           genomicPosition = if (length(genomicSpacings) > 0)
             genomicPositionCollection(organism = organism, spacings = genomicSpacings,
                                     namePrefix = "Genomic position", 
                                     dataSource = "UCSC genome genome, via Bioconductor",
                                     overlapFactor = 2,
                                     membershipBy= "start", 
                                     useUnits = "Mb",
                                     unit = 1e6) else NULL));
  }
  collections = collections[sapply(collections, function(x) !is.null(x))];
  if (merge) do.call(mergeCollections, collections) else collections;
}
    
#=======================================================================================================
#
# Build a collection of GO terms
#
#=======================================================================================================

.conditionalNumeric = function(x)
{
  n = suppressWarnings(as.numeric(x));
  if (any( is.na(n) & !is.na(x))) x else n;
}

buildGOcollection = function(
  organism,
  termNames = NULL,
  includeOffspring = TRUE, 
  strict = TRUE,
  verbose = 2, indent = 0)
{
   knownOrganisms = organismLabels();
   thisOrgLabels = organismLabels(organism);
   orgInd = match(thisOrgLabels[2], knownOrganisms[,2])
   orgShorthand = knownOrganisms$shorthand[orgInd];
   
   spaces = indentSpaces(indent);

   orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6));
   if (includeOffspring) 
   {
      goColumn = "GOALL"
      evidenceColumn = "EVIDENCEALL";
   } else {
      goColumn = "GO";
      evidenceColumn = "EVIDENCE";
   }

   IDColumn = c(rep("ENTREZID", 4), "ORF", rep("ENTREZID", 6))[orgInd];

   missingPacks = NULL;
   packageName = base::paste("org.", orgShorthand, orgExtensions[orgInd], ".db", sep="");
   if (!require(packageName, character.only = TRUE))
     missingPacks = c(missingPacks, packageName);

   #if (!require(GO.db))
   #  missingPacks = c(missingPacks, "GO.db");

   if (!is.null(missingPacks))
     stop(base::paste("Could not load the requisite package(s)",
           base::paste(missingPacks, collapse = ", "), ". Please install the package(s)."))

   if (verbose > 0)
     printFlush(base::paste(spaces, "GOcollection: loading annotation data..."));

   orgDb = get(packageName);
   orgGOInfo = select(orgDb, keys = keys(orgDb, keytype = IDColumn), keytype = IDColumn, 
                      columns = c(goColumn));

   orgGOInfo = orgGOInfo[!is.na(orgGOInfo[[goColumn]]), ];

   #orgGOlists = eval(parse(text = base::paste("select(", packageName, ", columns = ",
   #                                        reverseMap[orgInd],", keys = keys(", packageName, "))", sep = "")));
   n = nrow(orgGOInfo);
   index = tapply(1:n, orgGOInfo[[goColumn]], identity);
   orgGOEntrez = lapply(index, function(i) orgGOInfo[[IDColumn]] [i]);
   orgGOEvidence = lapply(index, function(i) orgGOInfo[[evidenceColumn]] [i]);
   orgGOids = names(orgGOEntrez);
   nTerms = length(orgGOids);

   goKeys = keys(GO.db);
   goInfo = select(GO.db, columns = c("GOID", "ONTOLOGY", "TERM", "DEFINITION"), keys = goKeys)
   dbGoIDs = goInfo$GOID;
   dbGoOntologies = goInfo$ONTOLOGY;
   dbGoTerm =goInfo$TERM;
   dbGoDescriptions = goInfo$DEFINITION;

   ### This was old code, here for reference
   #goInfo = eval(parse(text = "AnnotationDbi::as.list(GO.db:::GOTERM)"));
   #if (length(goInfo) > 0)
   #{
   #   dbGoIDs = as.character(sapply(goInfo, GOID));
   #   dbGoOntologies = as.character(sapply(goInfo, Ontology));
   #   dbGoTerm = as.character(sapply(goInfo, Term));
   #   dbGoDescriptions = as.character(sapply(goInfo, Definition));
   #} else {
   #   dbGoIDs = "";
   #}

   GO.org2GO = match(orgGOids, dbGoIDs);

   if (any(is.na(GO.org2GO)))
   {
     # Drop those organism GO lists that do not map to the GO database
     keep = is.finite(GO.org2GO);
     orgGOids = orgGOids[keep];
     nTerms = sum(keep);
     GO.org2GO = match(orgGOids, dbGoIDs);
   }
     
   orgOntologies = dbGoOntologies[GO.org2GO];
   orgIDs = dbGoIDs[GO.org2GO];
   orgTerms = dbGoTerm[GO.org2GO];
   orgDescriptions = dbGoDescriptions[GO.org2GO];

   if (is.null(termNames)) termNames = orgTerms;

   term2orgTerm = match(tolower(termNames), tolower(orgTerms));
   if (any(is.na(term2orgTerm)))
   {
      missingTermNames = term2orgTerm[is.na(term2orgTerm)];
      term2orgID = match(missingTermNames, dbGoIDs);
      if (any(is.na(term2orgID)))
      {
         if (strict) 
         {
           stop("GOGenesInCategory: the following terms were found ",
                "in neither GO names nor GO term IDs:\n",
                base::paste(missingTermNames[is.na(term2orgID)], collapse = ", "), 
                "\nTo turn off strict checking, use 'strict = FALSE'.")
         } else 
           warning(base::paste0("Warning in GOGenesInCategory: the following terms were found ",
                           "in neither GO names nor GO term IDs:\n",
                           base::paste(missingTermNames[is.na(term2orgID)], collapse = ", ")),
                   immediate. = TRUE);
      }
     
      term2orgTerm[is.na(term2orgTerm)] = term2orgID;
   }


   # Note: we have to drop those termNames that have no match in GO databases since it does not make sense to
   # include empty sets in the collection.

   keepInputTerms = is.finite(term2orgTerm);
   termNames = termNames[keepInputTerms];
   nInputTerms = sum(keepInputTerms);

   term2orgTerm = term2orgTerm[keepInputTerms];

   dataSets = vector(mode="list", length = nInputTerms);
   names(dataSets) = termNames;

   termIDs = orgIDs[term2orgTerm];
   termNames.GO = orgTerms[term2orgTerm];
   termOntologies = orgOntologies[term2orgTerm];
   termDescriptions = orgDescriptions[term2orgTerm];

   knownEvidence = knownEvidenceCodes()$evidenceCode;

   if (verbose > 0)
   {
      cat(base::paste(spaces, " ..preparing term lists..."));
      if (nInputTerms > 10) pind = initProgInd();
   }

   step = floor(nInputTerms/100);
   for (t in 1:nInputTerms) 
   {
      c = term2orgTerm[t];
      te = as.character(orgGOEvidence[[c]]); # Term evidence codes
      unknownEvidence = is.na(match(te, knownEvidence));
      if (any(unknownEvidence))
      {
        printFlush(WGCNA::spaste("Have the following unknown evidence codes: \n   ",
                   base::paste(unique(te[unknownEvidence]), collapse = ", ")))
        te[unknownEvidence] = "other";
      }
      tc = 
      dataSets[[t]] = newGeneSet(geneEntrez = .conditionalNumeric(orgGOEntrez[[c]]),
                       geneEvidence = as.factor(te),
                       geneSource = "GO",
                       ID = termIDs[t],
                       name = termNames.GO[t],
                       description = termDescriptions[t],
                       source = "GO",
                       groups = c("GO", base::paste0("GO.", termOntologies[t])),
                       organism = thisOrgLabels,
                       internalClassification = c("GO", termOntologies[t], termIDs[t]),
                       lastModified = Sys.Date());
 
      if (nInputTerms > 10 && verbose > 0 && (t%%step==0)) pind = updateProgInd(t/nInputTerms, pind);
   }

   if (nInputTerms > 10 && verbose > 0)
   {
      pind = updateProgInd(1, pind);
      printFlush("");
   }

   groups = list( newGroup(name = "GO", description = "Gene Ontology", 
                            source = "http://www.geneontology.org"),
                       newGroup(name = "GO.BP", description = "Gene Ontology - Biological Proces",
                            source = "http://www.geneontology.org", parents = "GO"),
                       newGroup(name = "GO.CC", description = "Gene Ontology - Cellular Compartment",
                            source = "http://www.geneontology.org", parents = "GO"),
                       newGroup(name = "GO.MF", description = "Gene Ontology - Molecular Function",
                            source = "http://www.geneontology.org", parents = "GO"));

   newCollection(dataSets = dataSets, groups = groups);
}


#=========================================================================================================
#
# Build MSigDB collection from the MSigDB XML file
#
#=========================================================================================================


# Convert the molecular signature database MSigDB, in its XML or simple text form, to an anRichment
# collection. Will also try to add the opposite 

buildMSigDBCollection = function(file, MSDBVersion = "5.0", excludeCategories = c("C1", "C2", "C5"))
{
   t = xmlTreeParse(file);
   msdb1 = t$doc$children$MSIGDB
   n = length(msdb1);
   mSigCategories = matrix(
     c("H", "Hallmark gene sets",
       "C1", "positional gene sets",
       "C2", "curated gene sets",
       "C3", "motif gene sets",
       "C4", "computational gene sets",
       "C5", "GO gene sets",
       "C6", "oncogenic signatures",
       "C7", "immunologic signatures"),
      ncol = 2, byrow = TRUE);
   MSigGroups = character(0);
   if (FALSE)
   {
     fullDescriptions = character(0);
     for (s in 1:n)
     {
       msigNode = as.list(xmlToList(msdb1[[s]]));
       fullDescriptions[s] = msigNode$DESCRIPTION_FULL;
     }
   }

   # In the end the full descriptions are not that great.

   geneSets = list();
   index = 1;
   for (s in 1:n)
   {
     msigNode = as.list(xmlToList(msdb1[[s]]));
     if (msigNode$CATEGORY_CODE %in% excludeCategories) next;

     msigCat.ext = anRichmentMethods:::.translateUsingTable(msigNode$CATEGORY_CODE, mSigCategories);
     msigCat.rec = spaste("MSigDB ", msigNode$CATEGORY_CODE, ": ", msigCat.ext);
     msigSubCat.rec = if (msigNode$SUB_CATEGORY_CODE!="") 
                       paste(msigCat.rec, "-", msigNode$SUB_CATEGORY_CODE) else NULL;
     groups = c("MSigDB", msigCat.rec, msigSubCat.rec);

     MSigGroups = unique(c(groups, MSigGroups));

     geneSets[[index]] = newGeneSet(geneEntrez = strsplit(msigNode$MEMBERS_EZID, split = ",", fixed =TRUE)[[1]],
                geneEvidence = "other",
                geneSource = msigNode$CONTRIBUTOR,

                ID = spaste("MSigDB.", msigNode$SYSTEMATIC_NAME),
                name = spaste(msigNode$STANDARD_NAME, " (MSigDB)"),
                shortName = msigNode$STANDARD_NAME,
                description = msigNode$DESCRIPTION_BRIEF, 
                source = msigNode$CONTRIBUTOR,
                organism = msigNode$ORGANISM,
                groups = groups,          
                internalClassification = c(groups, spaste(msigNode$STANDARD_NAME, " (MSigDB)")),
                alternateNames = msigNode$HISTORICAL_NAMES,
                externalDB = spaste("Molecular Signatures Database (MSigDB) version ", MSDBVersion),
                externalAccession = msigNode$SYSTEMATIC_NAME,
                webLink = "http://software.broadinstitute.org/gsea/msigdb/index.jsp");
     index = index + 1;
   }

   groupList = .mymapply(newGroup, name = MSigGroups, description = MSigGroups,
                        MoreArgs = list(source = spaste("Molecular Signatures Database (MSigDB) version ", 
                                                        MSDBVersion)));

   MSigDBCollection = newCollection(dataSets = geneSets, groups = groupList);
   MSigDBCollection;
}

MSigDBCollection = function(file, MSDBVersion = "5.0", 
                      excludeCategories = c("C1", "C2", "C5"),
                      organism = "human", useHomology = TRUE,
                      addOldOrganismToSetNames = FALSE, namePattern = ".convertedFrom.%o",
                      addOldOrganismToSetDescriptions = FALSE,
                      descriptionPattern = " (Converted from %o.)",
                      verbose = 1, indent = 0, ...)
{
  spaces = indentSpaces(indent);
  if (verbose > 0) printFlush(spaste(spaces, "Building MSigDB collection from XML file.."));
  collection = buildMSigDBCollection(file = file, MSDBVersion = MSDBVersion, excludeCategories = excludeCategories);
  if (length(organism)==0) return(collection);
  if (is.na(organism)) return(collection);

  if (verbose > 0) printFlush(spaste(spaces, "Converting to organism `", organism, "'.."));
  suppressWarnings(convertCollectionToOrganism(collection, organism = organism, useHomology = useHomology,
                   addOldOrganismToSetNames = addOldOrganismToSetNames,
                   namePattern = namePattern,
                   addOldOrganismToSetDescriptions = addOldOrganismToSetDescriptions,
                   descriptionPattern = descriptionPattern));
}

#========================================================================================================
#
# create a collection of gene sets in which each gene set contains genes within a certain genomic region.
#
#========================================================================================================

genomicPositionCollection = function(
  organism, 
  spacings,
  # optional genome annotation table.
  genomeAnnotation = NULL,
  entrezColumn = NULL,
  chromosomeColumn = NULL,
  startBpColumn = NULL,
  endBpColumn = NULL,
  # Output options
  namePrefix = "Genomic position", 
  dataSource = "",
  overlapFactor = 2,
  membershipBy= c("start", "end", "overlap"),
  useUnits = "Mb", 
  unit = 1e6)
{
  organism = organismLabels(organism)[1]
  membershipBy = match.arg(membershipBy);
  if (!is.null(genomeAnnotation))
  {
    .checkCol = function(name, value, data)
    {
      if (length(value)!=1) stop("Length of ", name, " must be 1.");
      if (is.numeric(value))
      {
        if (value <1 || value > ncol(data)) 
           stop(name, " = ", value, " is not between 1 and 'ncol(data)' = ", ncol(data), ".");
      } else {
        if (!value %in% names(data))
          stop(formatLabels(spaste(
            name, " = ", value, " is not among the names of 'data' which are: ", paste(names(data), collapse = ", ")),
            maxCharPerLine = 100));
      }
    }
    .checkCol("entrezColumn", entrezColumn, genomeAnnotation);
    .checkCol("chromosomeColumn", chromosomeColumn, genomeAnnotation);
    .checkCol("startBpColumn", startBpColumn, genomeAnnotation);
    .checkCol("endBpColumn", endBpColumn, genomeAnnotation);
  } else {
    if (!organism %in% c("human", "mouse"))
      stop("The position collection is only available for human and mouse.");
    txNames = switch(organism, human = "TxDb.Hsapiens.UCSC.hg19.knownGene",
                               mouse = "TxDb.Mmusculus.UCSC.mm10.knownGene");
    txdb.lst = as.list(get(txNames));
    geneTxIds = txdb.lst$genes$tx_id;
    genomeAnnotation = data.frame(txdb.lst$transcripts[ match(geneTxIds, txdb.lst$transcripts$tx_id), ],
                     gene_id = txdb.lst$genes$gene_id);
    chromosomeColumn = "tx_chrom";
    endBpColumn = "tx_end";
    startBpColumn = "tx_start";
    entrezColumn = "gene_id";
  }

  chrLevels = sort(unique(as.character(genomeAnnotation[[chromosomeColumn]])));
  
  # process each chromosome separately
  dataSets = list();
  mainGroupName = namePrefix;
  groups = list(newGroup(mainGroupName, spaste(namePrefix, ": Genes by genomic location"),
                         dataSource));

  for (chr in chrLevels)
  {
    printFlush(spaste("Working on chromosome ", chr));
    chrGrpName = spaste(namePrefix, ", chromosome ", chr);
    chrGrp = newGroup(name = chrGrpName,
                      description = spaste(namePrefix, ": Genes on ", chr),
                      source = dataSource);

    keep = which(genomeAnnotation[[chromosomeColumn]]==chr);
    genomeAnnotation1 = genomeAnnotation[keep, ];
    maxPos = max(genomeAnnotation1[[endBpColumn]]);
    for (sp in spacings)
    {
      nSets = ceiling(overlapFactor * maxPos/sp);
      setStarts = c(0: (nSets-1)) * sp/overlapFactor;
      dataSets1 = lapply(1:nSets, function(i)
      {
        start1 = (i-1) * sp/overlapFactor+1;
        end1 = start1 + sp;
        if (membershipBy=="start") {
           keep.int = which(genomeAnnotation1[[startBpColumn]] >= start1 & genomeAnnotation1[[startBpColumn]] < end1);
        } else if (membershipBy=="end") {
           keep.int = which(genomeAnnotation1[[endBpColumn]] >= start1 & genomeAnnotation1[[endBpColumn]] < end1);
        } else
           keep.int = which(genomeAnnotation1[[endBpColumn]] >= start1 & genomeAnnotation1[[startBpColumn]] < end1)
  
        if (length(useUnits) > 0 && is.character(useUnits))
        {
           prettyStart = spaste(round(start1/unit, 1), " ", useUnits);
           prettyEnd = spaste(round(end1/unit, 1), " ", useUnits);
           prettySp = spaste(round(sp/unit, 1), " ", useUnits);
           setName = spaste(namePrefix, ", ", chr, ": ", prettySp, " at ", prettyStart);
           shortSetName = spaste(chr, ": ", prettySp, " starting at ", prettyStart);
        } else {
           prettyStart = formatC(start1, format = "d", big.mark = ",");
           prettyEnd = formatC(end1-1, format = "d", big.mark = ",")
           setName = spaste(namePrefix, ", ", chr, ": ", prettyStart, "--", prettyEnd); 
           shortSetName = spaste(chr, ": ", prettyStart, "--", prettyEnd); 

        }
        if (length(keep.int)>0)
        {
           newGeneSet(geneEntrez = genomeAnnotation[[entrezColumn]] [keep[keep.int]],
                      geneEvidence = "other", geneSource = dataSource,
                      ID = spaste(namePrefix, ".", chr, ".", i),
                      name = setName,
                      shortName = shortSetName,
                      description = spaste(namePrefix, ": Genes on chromosome ", shortSetName),
                      source = dataSource,
                      organism = organism,
                      internalClassification = c(mainGroupName, chrGrpName, setName),
                      groups = c(mainGroupName, chrGrpName),
                      externalDB = dataSource);
        } else NULL;
      });
      dataSets = c(dataSets, dataSets1);
    }
    groups = c(groups, list(chrGrp)); 
  }
  dataSets = dataSets[ sapply(dataSets, length) > 0];
  newCollection(dataSets = dataSets, groups = groups);
}

