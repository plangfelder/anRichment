
.mymapply = function(..., SIMPLIFY = FALSE) base::mapply(..., SIMPLIFY = SIMPLIFY)

.allClassLabels = list(dataSet = "PL-dataSet", geneSet = "PL-geneSet", 
                numericProperty = "PL-numericProperty", discreteProperty = "PL-discreteProperty",
                group = "PL-group", collection = "PL-collection");

.replaceMissing = function(x, replaceWith)
{
  if (missing(replaceWith))
  {
    if (is.logical(x)) {
      replaceWith = FALSE
    } else if (is.numeric(x)) {
      replaceWith = 0;
    } else if (is.character(x)) {
      replaceWith = ""
    } else stop("Need 'replaceWith'.");
  }
  x[is.na(x)] = replaceWith;
  x;
}

.equalOrMissing = function(x, y)
{
  out = .replaceMissing(x==y);
  out[is.na(x) & is.na(y)] = TRUE;
  out;
}

.listRep = function(data, n)
{
  out = list();
  if (n> 0) for (i in 1:n) out[[i]] = data;
  out;
}



#===================================================================================================
#
# Handling of evidence codes
#
#===================================================================================================

knownEvidenceCodes = function()
{
  codes.0 = c("EXP", "Inferred from Experiment", "Experimental", 
              "HTP", "Inferred from High Throughput Experiment", "High throughput experimental",
              "HDA", "Inferred from High Throughput Direct Assay", "High throughput experimental",
              "HMP", "Inferred from High Throughput Mutant Phenotype", "High throughput experimental",
              "HGI", "Inferred from High Throughput Genetic Interaction", "High throughput experimental",
              "HEP", "Inferred from High Throughput Expression Pattern", "High throughput experimental",
              "IDA", "Inferred from Direct Assay", "Experimental", 
              "IPI", "Inferred from Physical Interaction", "Experimental", 
              "IMP", "Inferred from Mutant Phenotype", "Experimental", 
              "IGI", "Inferred from Genetic Interaction", "Experimental", 
              "IEP", "Inferred from Expression Pattern", "Experimental", 
              "ISS", "Inferred from Sequence or Structural Similarity","Computational",  
              "ISO", "Inferred from Sequence Orthology","Computational",  
              "ISA", "Inferred from Sequence Alignment","Computational",  
              "ISM", "Inferred from Sequence Model","Computational",  
              "IGC", "Inferred from Genomic Context","Computational",  
              "IBA", "Inferred from Biological aspect of Ancestor","Computational",  
              "IBD", "Inferred from Biological aspect of Descendant","Computational",  
              "IKR", "Inferred from Key Residues","Computational",  
              "IMR", "Inferred from Missing Residues","Computational",  
              "IRD", "Inferred from Rapid Divergence","Computational",  
              "RCA", "Inferred from Reviewed Computational Analysis","Computational",  
              "TAS", "Traceable Author Statement", "Author statement", 
              "NAS", "Non-traceable Author Statement", "Author statement", 
              "IC", "Inferred by Curator", "Curator statement",
              "ND", "No biological Data available", "Curator statement",
              "IEA", "Inferred from Electronic Annotation", "Automatically assigned",
              "NR", "Not Recorded ", "Obsolete",
              "other", "other", "other");

  codes.mat = matrix(codes.0, ncol = 3, byrow = TRUE);
  colnames(codes.mat) = c("evidenceCode", "evidenceDescription", "evidenceType");
  codes.df = data.frame(codes.mat, numericCode = c(1:nrow(codes.mat)));
  codes.df;
}

# For a given set of evidence identifiers, return all unique matched numeric codes from all known evidence
# codes

matchEvidenceCode = function(evidence)
{
  if (length(evidence)==0)
    stop("'evidence' must be non-empty. Use evidence='all' to get genes with any evidence.");

  evidence = unique(evidence);

  ec = knownEvidenceCodes();
  nr = nrow(ec);
  codes.0 = tolower(unlist(ec[, c(1:3)]));
  if (length(evidence)==1 && evidence=="all")
  {
    row = 1:nrow(ec);
  } else {
    ev2codes = match(tolower(evidence), codes.0);
    if (any(is.na(ev2codes)))
      stop("Some entries in 'evidence' are not recognized:\n    ",
           base::paste(evidence[is.na(ev2codes)], collapse = ", "));
    ind = which(codes.0 %in% tolower(evidence))
    row = sort(unique(floor((ind-1) %% nr) + 1));
  }
  ec$numericCode[row];
}

 
#======================================================================================================
#
# Handling of organism labels
#
#======================================================================================================

organismLabels = function(organism = NULL)
{
  knownOrganisms = c("human", "Homo sapiens", "Hs",
                     "mouse", "Mus musculus", "Mm",
                     "rat", "Rattus norvegicus", "Rn",
                     "malaria", "Plasmodium falciparum", "Pf",
                     "yeast", "Saccharomyces cerevisiae", "Sc",
                     "fly", "Drosophila melanogaster", "Dm",
                     "bovine", "Bos taurus", "Bt",
                     "worm", "Caenorhabditis elegans", "Ce",
                     "canine", "Canis familiaris", "Cf",
                     "zebrafish", "Danio rerio", "Dr",
                     "chicken", "Gallus gallus", "Gg",
                     "mosquito", "Anopheles gambiae", "Ag",
                     "monkey", "Macaca mulatta", "Mmu",
                     "chimp", "Pan troglodytes", "Pt",
                     "pig", "Sus scrofa", "Ss")

  knownOrganisms.mat = matrix(knownOrganisms, ncol = 3, byrow = TRUE);
  colnames(knownOrganisms.mat) = c("commonName", "scientificName", "shorthand")
  knownOrganisms.df = as.data.frame(knownOrganisms.mat);
  if (is.null(organism))
    return(knownOrganisms.df)
  pos = charmatch(tolower(organism), tolower(knownOrganisms));
  if (is.na(pos))
    stop("Organism ", organism, " is not recognized. Recognized variants are:\n",
         base::paste( base::paste("   ", apply(knownOrganisms.df, 1, base::paste, collapse = ", "), collapse = "\n")))

  if (pos==0)
    stop("Organism ", organism, 
         " was multiply matched. Please supply more characters for unambiguous match.")

  row = floor((pos-1)/3) + 1;
  as.character(knownOrganisms.mat[row, ]);
}

organismShorthand = function(organism)
{
  if (is.null(organism)) stop("'organism' must be given.");
  organismLabels(organism)[3];
}

sameOrganism = function(org1, org2)
{
  sh1 = organismShorthand(org1);
  sh2 = organismShorthand(org2);

  sh1==sh2;
}

#======================================================================================================
#
# newGeneSet, newGeneProperty
#
#======================================================================================================

.checkDate = function(date, format="", origin = "1970-1-1")
{
  if (length(date)==0) return(NULL);
  if (is.character(date)) date = as.Date(date, format);
  if (is.numeric(date)) date = as.Date(date, origin = origin);
  as.Date(date);
}
  

.newDataSet = function(ID, name, shortName, description, source, organism, internalClassification, groups,
                       lastModified, format,
                       alternateNames,
                       externalDB, externalAccession,
                       webLink, 
                       weightIndex = NA)
{
  if (length(organism)!=3)
  {
    if (length(organism)!=1)
      stop("'organism' must be a single character string or a set of 3 standard labels.")
    organism = organismLabels(organism);
  }

  lastModified = .checkDate(lastModified, format = format);

  if (is.na(ID) || ID=="") stop("'ID' cannot be missing or empty.");
  if (is.na(name) || name=="") stop("'name' cannot be missing or empty.");
  
  out = list(ID = ID,
             name = name,
             shortName = shortName,
             description = description,
             source = source,
             organism = organism,
             internalClassification = internalClassification,
             groups = groups,
             lastModified = lastModified,
             alternateNames = alternateNames,
             externalDB = externalDB,
             externalAccession = externalAccession,
             webLink = webLink,
             weightIndex = weightIndex);
  class(out) = c(.allClassLabels$dataSet, class(out));
  out;
}

newGeneSet = function(geneEntrez, geneEvidence, geneSource, 
                      ID, name, shortName = name, 
                      description, source, organism, 
                      internalClassification, groups,
                      lastModified = Sys.Date(),
                      format = "%Y-%m-%d",
                      alternateNames = character(0),
                      externalDB = "",
                      externalAccession = "",
                      webLink = "")
{
  evidenceCodes = knownEvidenceCodes();

  n = length(geneEntrez);
      
  geneEvidence = as.character(geneEvidence);
  numEvidence = match(geneEvidence, evidenceCodes$evidenceCode);

  if (any(is.na(numEvidence))) 
    stop("Some entries in 'geneEvidence' were not recognized (or are missing).");

  if (length(geneEvidence)==1) geneEvidence = rep(geneEvidence, n);
  if (length(geneSource)==1) geneSource = rep(geneSource, n);

  if (length(geneEvidence)!=n) 
     stop("'geneEvidence' must be a scalar or a vector of length to that of 'geneEntrez'.");
  
  if (length(geneSource)!=n) 
     stop("'geneSource' must be a scalar or a vector of length to that of 'geneEntrez'.");

  base = .newDataSet(ID = ID, name = name, shortName = shortName, 
                     description = description,source = source, organism = organism, 
                     internalClassification = internalClassification, groups = groups, 
                     lastModified = lastModified, format = format,
                     alternateNames = alternateNames,
                     externalDB = externalDB,
                     externalAccession = externalAccession,
                     webLink = webLink);
  out = c( list(data = data.frame(Entrez = geneEntrez,
                                  evidence = factor(as.character(geneEvidence)),
                                  source = factor(as.character(geneSource)))),
           base,
           type = .allClassLabels$geneSet);
  class(out) = c(.allClassLabels$geneSet, class(base)); 
  out;
}


.completeName = function(name, completions, leadingSep = ".", trailingSep = "", insertAt = "end")
{
  if (length(name)==length(completions)) 
  {
    return(name)
  } else if (length(name)==1)
  {
    if (insertAt=="end")
    {
       return(WGCNA::spaste(name, leadingSep, completions, trailingSep))
    } else if (insertAt=="start")
    {
       return(WGCNA::spaste(leadingSep, completions, trailingSep, name))
    } else 
       stop("'insertAt' must be 'start' or 'end'");
  } else
    stop("'name' must either be a scalar or it must have the same length as 'completions'.");
}

.isDiscrete = function(x, maxDiscreteLevels)
{
  if (!is.numeric(x)) return(TRUE);
  length(unique(x)) <= maxDiscreteLevels;
}

.extendVector = function(x, currentID, newID, missingValue = NA)
{
  out = rep(missingValue, length(newID));
  if (length(x) > 0)
  {
    current2new = match(currentID, newID);
    if (any(is.na(current2new)))
       stop("Some elements in 'currentID' are not present in 'newID'.");
    out[match(currentID, newID)] = x;
  }
  out;
}

.extendGeneProperty = function(dataSet, currentID, newID, missingValue = NA)
{
  if (.isGeneProperty(dataSet))
  {
    dataSet$data = .extendVector(dataSet$data, currentID, newID, missingValue);
  }

  dataSet;
}

.checkOrExtend = function (x, n) 
{
  if (length(x)==1)
  {
    x = rep(x, n);
  } else if (length(x)!=n) stop("'x' must either have length 1 or length 'n'.");
  x;
}


# This function returns a list of "bare" gene property data sets. These do not constitute a collection
# since they miss identifiers, weights, group descriptions etc.

.newProperty.base = function(properties, weightIndex, discreteIndex, maxDiscreteLevels = 5,
                             IDs, names, shortNames, descriptions, source, organism,
                             internalClassificationBase, internalClassificationExtras,
                             groupBase, groupExtras, 
                             lastModified, format, 
                             alternateNames,
                             externalDB,
                             externalAccession,
                             webLink,
                             simplify = FALSE)
{
  properties = as.data.frame(properties);
  nSets = ncol(properties)
  geneDataSets = list();
  if (is.null(discreteIndex)) 
  { 
    discrete.log = sapply(properties, .isDiscrete, maxDiscreteLevels = maxDiscreteLevels);
  } else if (!is.logical(discreteIndex)) {
    discrete.log = rep(FALSE, nSets);
    discrete.log[discreteIndex] = TRUE;
  }
  for (set in 1:nSets)
  {
    base = .newDataSet(ID = IDs[set], 
                       name = names[set],
                       shortName = shortNames[set],
                       description = descriptions[set],
                       source = source,
                       organism = organism,
                       internalClassification = c(internalClassificationBase,
                                                  internalClassificationExtras[[set]]),
                       groups = c(groupBase, groupExtras[[set]]),
                       lastModified = lastModified[set],
                       format = format,
                       alternateNames = alternateNames[set],
                       externalDB = externalDB[set],
                       externalAccession = externalAccession[set],
                       webLink = webLink[set],
                       weightIndex = weightIndex);
    type = if (discrete.log[set]) .allClassLabels$discreteProperty else .allClassLabels$numericProperty;
    geneDataSets[[set]] = c(list( data = if (discrete.log[set]) as.factor(properties[, set]) 
                                          else properties[, set],
                                 type = type),
                           base);
    class(geneDataSets[[set]]) = c( type, class(base));
  }
  if (simplify && nSets==1) geneDataSets = geneDataSets[[1]];
  geneDataSets;
}

                                
# New gene property returns a collection of the properties since weights and identifiers are a part of a
# collection. Optionally the property can be added to an existing collection.


newGeneProperty = function(
        identifiers, properties, weights = NULL, discreteIndex = NULL, maxDiscreteLevels = 5,
        propertyNames = colnames(properties),
        IDs, names, shortNames, descriptions, source, organism, 
        internalClassificationBase, internalClassificationExtras,
        groupBase, groupExtras, 
        lastModified = Sys.Date(), 
        format = "%Y-%m-%d",
        alternateNames = list(character(0)),
        externalDB = "",
        externalAccession = "",
        webLink = "",
        newGroups = NULL, collection = NULL)
{

  stop("This function has not been properly tested and probably contains bugs. Sorry!");
  # FIXME: 
  # 1. It seems we always add a unit weight column which may not be necessary, especially when turning data
  # frames into a collection
  # 2. in the .newGeneProperty.base function, the weight index seems to be treated as a scalar.

  nGenes = length(identifiers);
  nProps = ncol(properties);

  IDs = .completeName(IDs, propertyNames);
  names = .completeName(names, propertyNames);
  descriptions = .completeName(descriptions, propertyNames, leadingSep = "(", trailingSep = ")",
                               insertAt = "start");

  lastModified = .checkOrExtend(lastModified, nProps);
  alternateNames = .checkOrExtend(alternateNames, nProps);
  externalDB = .checkOrExtend(externalDB, nProps);
  externalAccession = .checkOrExtend(externalAccession, nProps);
  webLink = .checkOrExtend(webLink, nProps);
  
  baseWeights = as.matrix(rep(1, nGenes));
  if (is.null(weights))
  {
    newWeights = baseWeights
  } else 
    newWeights = cbind(baseWeights, weights);

  if (is.null(collection))
  {
    collection = newCollection(identifiers = identifiers, 
                               weights = newWeights,
                               groups = newGroups)
    weightIndex = ncol(newWeights);
  } else {
    combinedIdentifiers = union(collection$identifiers, identifiers)
    
    collection$dataSets = lapply(collection$dataSets, .extendGeneProperty, collection$identifiers, 
                                 combinedIdentifiers);

    properties = apply(properties, 2, .extendVector, identifiers, combinedIdentifiers);
    collWeights = as.matrix(apply(collection$weights, 2, .extendVector, collection$identifiers, 
                                  combinedIdentifiers));
    
    if (!is.null(weights))
    {
      collWeights = cbind(collWeights, .extendVector(weights, identifiers, combinedIdentifiers))
      weightIndex = ncol(collWeights);
    } else
      weightIndex = 1;

    collection$weights = collWeights;
    collection$identifiers = combinedIdentifiers
    if (!is.null(newGroups))
      collection = do.call(addToCollection, c(list(collection = collection), newGroups));

    
  }
    
  newDataSets = .newProperty.base(properties,  weightIndex = weightIndex, 
                             discreteIndex = discreteIndex, maxDiscreteLevels = maxDiscreteLevels,
                             IDs = IDs, names = names, shortNames = shortNames, descriptions = descriptions, 
                             source = source, organism = organism,
                             internalClassificationBase = internalClassificationBase, 
                             internalClassificationExtras = internalClassificationExtras,
                             groupBase = groupBase, groupExtras = groupExtras, 
                             lastModified = lastModified, format = format, 
                             alternateNames = alternateNames,
                             externalDB = externalDB,
                             externalAccession = externalAccession,
                             webLink = webLink,
                             simplify = FALSE)

  collection$dataSets = c(collection$dataSets, newDataSets);

  collection;
}
  
#======================================================================================================
#
# newGroup
#
#======================================================================================================

newGroup = function(
  name, 
  description="", 
  source = "", 
  alternateNames = character(0), 
  parents = character(0))
{
  if (is.na(name) || name=="") stop("Name cannot be missing or an empty character string.");
  out = list(name = name, description = description, source = source,
             alternateNames = alternateNames, parents = parents);
  class(out) = c(.allClassLabels$group, class(out));
  out;
}

.standardGroupNameTable = function(groups)
{
  lst = lapply(groups, function(grp) 
       cbind(alternateName = grp$alternateNames, standardName = rep(grp$name, length(grp$alternateNames))));
  as.data.frame(do.call(rbind, lst));
}

#======================================================================================================
#
# mergeGroups
#
#======================================================================================================

mergeGroups = function(collection, groupNames)
{
  if (!.isCollection(collection)) stop("Given 'collection' is not a valid collection.");
  if (length(groupNames)<2) stop("At least two different groups must be specified in 'groupNames'.");
  in2groups = match(groupNames, collection$groups$name);
  if (any(is.na(in2groups))) 
    stop("The following given 'groupNames' were not found among groups in the given 'collection':\n   ",
         base::paste(groupNames[is.na(in2groups)], collapse = "\n   "));

  keepGroup = groupNames[1];
  dropGroups = groupNames[-1];

  collection$groups = collection$groups[-in2groups[-1], ];

  collection$dataSets = lapply(collection$dataSets, function(ds)
   {
      replace = ds$groups %in% dropGroups;
      ds$groups[replace] = keepGroup;
      ds$groups = unique(ds$groups);
      ds;
   });

  collection;
}

#======================================================================================================
#
# newCollection
#
#======================================================================================================

newCollection = function(dataSets = list(), 
                        groups = list(), 
                        identifiers = numeric(0), 
                        weights = matrix(nrow = 0, ncol = 1),
                        ...)
{
  out = list(dataSets = dataSets,
             groups = groups,
             identifiers = identifiers,
             weights = weights,
             evidenceCodes = knownEvidenceCodes());
  class(out) = c(.allClassLabels$collection, class(out));
  addToCollection(out, ...);
}

#======================================================================================================
#
# addToCollection
#
#======================================================================================================

.is = function(object, class, namesToCheck)
{
  return (inherits(object, class) && all(all(namesToCheck %in% names(object))));
}

.isDataSet = function(object)
{
  expectedNames = c("ID", "name", "description", "source", "organism", "internalClassification", 
                    "groups");
  .is(object, .allClassLabels$dataSet, expectedNames);
}


.isGeneSet = function(object)
{
  .is(object, .allClassLabels$geneSet, "data") && .isDataSet(object);
}

.isNumericProperty = function(object)
{
  .is(object, .allClassLabels$numericProperty, "data") && .isDataSet(object);
}

.isDiscreteProperty = function(object)
{
  .is(object, .allClassLabels$discreteProperty, "data") && .isDataSet(object);
}

.isGeneProperty = function(object)
{
  .isNumericProperty(object) || .isDiscreteProperty(object)
}

.isGroup = function(object)
{
  .is(object, .allClassLabels$group, names(newGroup("__")));
}

.isCollection = function(object)
{
  .is(object, .allClassLabels$collection, names(newCollection()));
}

addToCollection = function(collection, ...)
{
  dots = list(...)
  if (length(dots)==0) return(collection)

  if (is.null(names(dots))) names(dots) = rep("", length(dots));

  isGeneSet = sapply(dots, .isGeneSet);
  isGroup = sapply(dots, .isGroup);
  isCollection = sapply(dots, .isCollection);
  neither = !isGeneSet & !isGroup & !isCollection;
  if (any(neither))
    stop("Some given objects are neither valid gene sets nor valid groups nor valid collections. \n",
         "  Names and positions of the offending objects: \n   ",
         base::paste( base::paste0(names(dots)[neither], " (", which(neither), ")"), collapse = ", "));
  if (any(isGeneSet))
     collection$dataSets = c(collection$dataSets, dots[which(isGeneSet)]);
  if (any(isGroup))
     collection$groups = c(collection$groups, dots[which(isGroup)]);
  if (any(isCollection))
     collection = do.call(mergeCollections,c(list(collection), dots[isCollection]))
  collection;
}

#======================================================================================================
#
# .removeDuplicates
#
#======================================================================================================

.removeDuplicates = function(object, nameComponent = "name", stopOnDup = FALSE, warnDups = FALSE)
{
  names = lapply(object, getElement, nameComponent);
  duplicated = duplicated(names);
  if (any(duplicated))
  {
    if (stopOnDup)
    {
      stop("Some elements of 'object' are duplicated. Duplicated names:\n", 
           base::paste(names[duplicated], collapse = ", "))
    }
    if (warnDups)
    {
      warning(.immediate = TRUE, 
           ".removeDuplicates: Some elements of 'object' are duplicated. Duplicated names:\n",
           base::paste(names[duplicated], collapse = ", "))
    }
  }

  object[!duplicated];
}
  
                             

#======================================================================================================
#
# mergeCollections
#
#======================================================================================================

mergeCollections = function(...,
   stopOnDuplicates = FALSE)
{
  args = list(...);

  keep = sapply(args, length)>0;
  args = args[keep];

  if (length(args)==0) return(newCollection());

  isCollection = sapply(args, .isCollection);
  if (any(!isCollection))
    stop("Some given objects are not valid collections. \n",
         "  Names and positions of the offending objects: \n   ",
         base::paste( base::paste0(names(args)[!isCollection], " (", which(!isCollection), ")"), collapse = ", "));

  # Create a list of identifiers

  combinedIdentifiers = multiUnion( lapply(args, getElement, "identifiers"));

  # Extend weights for all collections

  combWeights.list = lapply(args, 
                       function(x) 
                       { 
                         apply(x$weights[, -1, drop = FALSE], 2, .extendVector, 
                               x$identifiers, combinedIdentifiers)
                       })

  combWeights = do.call(cbind, c(list(rep(1, length(combinedIdentifiers))), combWeights.list))

  # Extend gene properties

  combDataSets = do.call("c", lapply(args, function(x)
                       {
                         lapply(x$dataSets, .extendGeneProperty, x$identifiers, combinedIdentifiers)
                       }));

  # Build the new collection
 
  newCollection(dataSets =.removeDuplicates(combDataSets, stopOnDup = stopOnDuplicates),
                      groups = .removeDuplicates(do.call("c", lapply(args, getElement, "groups")),
                                 nameComponent = "name", warnDups = FALSE, stopOnDup = FALSE),
                      identifiers = combinedIdentifiers,
                      weights = combWeights);
}

#============================================================================================
#
# geneLists
#
#============================================================================================

# Isolate the gene lists from dataSets component of collection. Only keep genes that match the evidence
# and keep each gene only once. If no gene matches evidence, the corresponding component will be a numeric
# vector of length 0 (which is different from NULL).

# For data sets that are not gene sets, the returned gene list is numeric(0).

.geneLists = function(
   collection, evidence, tags, 
   matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                      "groupAlternateNames", "nameAndAlternates", "groupsAndAlternates"),
   searchType, invertSearch, exactMatch = TRUE, fixed = TRUE, ignore.case = TRUE,
   namesFrom = "ID", simplify = FALSE,
   firstDate = NULL, lastDate = NULL,
   dateFormat = "%Y-%m-%d", invertDateSearch = FALSE)
{
  if (length(tags) > 0)
  {
    index = .matchDataSets(collection, tags = tags, 
                           matchComponents = matchComponents, 
                           searchType = searchType, invert = invertSearch,
                           exactMatch = exactMatch, fixed = fixed, ignore.case = ignore.case,
                           firstDate = firstDate, lastDate = lastDate, dateFormat = dateFormat,
                           invert.date = invertDateSearch);
    if (length(index)==0)
    {
      warning(".geneLists: no data sets matched the given 'tags' and search type.");
      return(list());
    }
    if (!isTRUE(all.equal(index, 1:length(collection$dataSets))))  collection$dataSets = collection$dataSets[index];
  }

  dataType = sapply(collection$dataSets, getElement, "type");
  nDataSets = length(dataType);
  out = vector(mode = "list", length = nDataSets);
  geneSetIndex = which(dataType==.allClassLabels$geneSet);
  if (evidence=="all") 
  {
     geneLists = lapply(collection$dataSets[geneSetIndex], function(x) { unique(x$data$Entrez) });
  } else {
     keepNumCodes = matchEvidenceCode(evidence);
     geneLists = lapply(collection$dataSets[geneSetIndex], function(dataSet) 
              { unique(dataSet$data$Entrez [ dataSet$data$evidence %in% keepNumCodes] ) } );
  }
  out[geneSetIndex] = geneLists;
  for (ds in c(1:nDataSets)[-geneSetIndex]) out[[ds]] = numeric(0);
  if (length(namesFrom) > 0 && (is.character(namesFrom) && namesFrom!="")) 
  {
     names(out) = sapply(collection$dataSets, getElement, namesFrom);
  } else 
     names(out) = NULL;
  if (simplify && length(out)==1) out = out[[1]];
  out;
}

geneLists = function(collection, evidence = "all", 
                     tags = NULL, 
                      matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                      "groupAlternateNames", "nameAndAlternates", "groupsAndAlternates"),
                     searchType = c("any", "all"), invertSearch = FALSE,
                     exactMatch = TRUE, fixed = TRUE, ignore.case = TRUE,
                     firstDate = NULL, lastDate = NULL,
                     dateFormat = "%Y-%m-%d", invertDateSearch = FALSE,
                     namesFrom = c("ID", "name"),
                     simplify = TRUE)
{

  namesFrom = match.arg(namesFrom);
  .geneLists(collection, evidence = evidence, tags = tags, 
             matchComponents = matchComponents, searchType = searchType, 
             invertSearch = invertSearch, exactMatch = exactMatch, fixed = fixed, 
             ignore.case = ignore.case,
             firstDate = firstDate, lastDate = lastDate, dateFormat = dateFormat,
             invertDateSearch = invertDateSearch, 
             namesFrom = namesFrom, simplify = simplify)
}



#============================================================================================
#
# allGenes
#
#============================================================================================

# Gets all genes with specified evidence and puts them in a single vector. Useful for getting backgrounds.
  
allGeneSetGenes = function(
   collection, evidence = "all",
   tags = NULL, 
   matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                      "groupAlternateNames", "nameAndAlternates", "groupsAndAlternates"),
   searchType = c("any", "all"), invertSearch = FALSE,
   exactMatch = TRUE, fixed = TRUE, ignore.case = TRUE,
   firstDate = NULL, lastDate = NULL,
   dateFormat = "%Y-%m-%d", invertDateSearch = FALSE)
{
  if (!.isCollection(collection))
    stop("The argument 'collection' must be an anRichment collection.");

  geneLists.1 = .geneLists(collection, evidence = evidence, tags = tags,
                           matchComponents = matchComponents, searchType = searchType, 
                           invertSearch = invertSearch, exactMatch = exactMatch, fixed = fixed,
                           ignore.case = ignore.case,
                           firstDate = firstDate, lastDate = lastDate, dateFormat = dateFormat,
                           invertDateSearch = invertDateSearch);
  unique(unlist(geneLists.1));
}

allDataSetGenes = function(
  collection, evidence = "all", tags = NULL, 
  matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                      "groupAlternateNames", "nameAndAlternates", "groupsAndAlternates"),
  searchType = c("any", "all"), invertSearch = FALSE,
  exactMatch = TRUE, fixed = TRUE, ignore.case = TRUE,
  firstDate = NULL, lastDate = NULL,
  dateFormat = "%Y-%m-%d", invertDateSearch = FALSE)
{
  unique(c(allGeneSetGenes(collection, evidence, tags = tags, 
                           matchComponents = matchComponents, searchType = searchType,
                            invertSearch = invertSearch, exactMatch = exactMatch, fixed = fixed,
                            ignore.case = ignore.case,
                            firstDate = firstDate, lastDate = lastDate, dateFormat = dateFormat,
                            invertDateSearch = invertDateSearch)), 
                           collection$identifiers);
}

#============================================================================================
#
# dataSetNames and dataSetIDs
#
#============================================================================================

dataSetNames = function(collection, groups = NULL, dataSets = NULL)
{
  if (!is.null(groups) || !is.null(dataSets)) 
    collection = subsetCollection(collection, tags = c(groups, dataSets));

  out = sapply(collection$dataSets, getElement, "name");
  if (length(out)==0) out = character(0);

  as.character(out);
}

shortDataSetNames = function(collection, groups = NULL, dataSets = NULL)
{
  if (!is.null(groups) || !is.null(dataSets)) 
    collection = subsetCollection(collection, tags = c(groups, dataSets));

  out = sapply(collection$dataSets, getElement, "shortName");
  if (length(out)==0) out = character(0);

  as.character(out); 
}


dataSetIDs = function(collection, groups = NULL, dataSets = NULL)
{
  if (!is.null(groups) || !is.null(dataSets)) 
    collection = subsetCollection(collection, tags = c(groups, dataSets));

  out = sapply(collection$dataSets, getElement, "ID");
  if (length(out)==0) out = character(0);

  as.character(out);
}

#===========================================================================================
#
# Handling of groups
#
#===========================================================================================

.knownGroups = function(dataSets, sortBy = c("alphabet", "size"), reverseSort = FALSE)
{
  sortBy = match.arg(sortBy);
  groups.1 = unlist(lapply(dataSets, getElement, "groups"));
  if (length(groups.1)==0) return(character(0));
  sizes = table(groups.1)
  if (sortBy=="size")
  {
    order = order((as.numeric(reverseSort) - 0.5) * sizes);
    names(sizes)[order];
  } else
    sort(names(sizes), decreasing = reverseSort);
}

knownGroups = function(collection, sortBy = c("alphabet", "size"), reverseSort = FALSE )
{
  #sortBy = match.arg(sortBy);
  .knownGroups(collection$dataSets, sortBy, reverseSort);
}
 
nDataSets = function(collection)
{
  length(collection$dataSets);
}

.weightIndex = function(dataSets)
{
  out = sapply(dataSets, getElement, "weightIndex");
  if (length(out)==0) out = numeric(0);
  out;
}

#============================================================================================
#
# .matchDataSets
#
#============================================================================================

.indexedTags = function(collection, impliedGroups, 
                   matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                                     "groupAlternateNames",
                                     "nameAndAlternates", "groupsAndAlternates"))
{
  nGS = length(collection$dataSets)

  if ("nameAndAlternates" %in% matchComponents) 
    matchComponents = unique(c("name", "alternateNames", setdiff(matchComponents, "nameAndAlternates")));

  if ("groupsAndAlternates" %in% matchComponents) 
    matchComponents = unique(c("groups", "groupAlternateNames", setdiff(matchComponents, "groupsAndAlternates")));

  if ("groupAlternateNames" %in% matchComponents)
  {
     if (replaceMissing(match("groups", matchComponents), 10000) > match("groupAlternateNames", matchComponents))
       stop("When 'matchComponents' contains \"groupAlternateNames\", it must also contain \"groups\"\n",
            " __before__ \"groupAlternateNames\".");
  }

  setIndex = numeric(0);
  tags = character(0);
  for (mc in matchComponents)
  {
    if (mc=="groupAlternateNames")
    {
      alternateList = lapply(collection$groups, .getElement1, "alternateNames");
      names(alternateList) = sapply(collection$groups, getElement, "name");
      groupAlternates = lapply(groupsWithParents, function(grp) unlist(alternateList[grp]));
      tags1 = groupAlternates;
    } else {
      tags1 = lapply(collection$dataSets, getElement, mc);
      if (mc=="groups")
      {
        index = lapply(tags1, match, names(impliedGroups));
        groupsWithParents = .mymapply(function(g, ind1)
        {
          fin = is.finite(ind1);
          unique(c(g, unlist(impliedGroups[ ind1[fin]])))
        }, tags1, index);
        tags1 = groupsWithParents;
      }
    }
    index1 = c( rep(c(1:nGS), sapply(tags1, length)));
    tags = c(tags, unlist(tags1));
    setIndex = c(setIndex, index1);
  }
  data.frame(index = setIndex, tags = tags);
}

.getElement1 = function(lst, name)
{
  if (name %in% names(lst)) return(lst[[name]])
  character(0);
}

impliedGroups = function(groupList, queryGroups = NULL, includeSelf = TRUE,
                   get = c("parents", "children"))
{
  get = match.arg(get);
  groupNames = as.character(sapply(groupList, getElement, "name"));
  groupNames.ext = .indexedFlattenedList(lapply(groupList, getElement, "alternateNames"));
  groupNames.x = rbind(.indexedFlattenedList(groupNames), groupNames.ext);
  nGroups = length(groupNames);
  if (length(queryGroups)==0) queryGroups = groupNames;

  if (any(!queryGroups %in% groupNames.x$data))
    stop("Some of 'queryGroups' are not among given groups:\n",
         formatLabels( paste(queryGroups[!queryGroups %in% groupNames.x$data], collapse = ", "),
                       split = ", ", maxCharPerLine = 80), "\n\n");

  parents = lapply(groupList, .getElement1, "parents");
  parentIndex = lapply(parents, .translateUsingTable, groupNames.x[c(2, 1)]);
  .mymapply(function(i, name, parents) 
    if (any(is.na(i))) stop("Invalid parent specification for group ", name, ".\n  Parent group(s) ", 
                    paste(parents[is.na(i)], collapse = "; "), " does (do) not exist."),
    parentIndex, groupNames, parents);

  if (get=="parents") {
    impliedIndex1 = parentIndex;
  } else if (get=="children") {
    impliedIndex1 = .listRep(numeric(0), nGroups);
    for (g in 1:nGroups) for (pg in parentIndex[[g]])
      impliedIndex1[[ pg ]] = c(impliedIndex1[[pg]], g);
  }

  impliedIndex = lapply(match(queryGroups, groupNames), identity)
  lengths = sapply(impliedIndex, length);
  names(impliedIndex) = queryGroups;
  change = TRUE;
  while (change)
  {
    impliedIndex = lapply(impliedIndex, function(ii) 
    {
      add = unlist(impliedIndex1[ii]);
      if (any(add==ii[1]))
      {
        browser()
        stop("Loop in parent relationships: index ", ii[1], " appeared again.");
      }
      unique(c(ii, add))
    });
    newLengths = sapply(impliedIndex, length);
    if (all(lengths==newLengths)) {
      change = FALSE
    } else
      lengths = newLengths;
  }
  
  if (!includeSelf) impliedIndex = lapply(impliedIndex, `[`, -1);
  impliedGroups = lapply(impliedIndex, function(i) groupNames[i]);
  impliedGroups;
}

# return indices of data sets that match any or all of tags (according to searchType)

.matchDataSets = function(collection, tags, searchType = c("any", "all"), 
                          matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                                     "groupAlternateNames",
                                     "nameAndAlternates", "groupsAndAlternates"),
                          invert = FALSE,
                          exactMatch = TRUE, fixed = TRUE, ignore.case = TRUE,
                          firstDate = NULL, lastDate = NULL, 
                          dateFormat = "%Y-%m-%d", invert.date = FALSE,
                          impliedGroups = NULL)
{
  if (is.null(impliedGroups)) 
    impliedGroups = impliedGroups(collection$groups, queryGroups = NULL, includeSelf = TRUE,
                                  get = "parents");
  matches = 1:length(collection$dataSets);
  tags.use = tags;
  if (ignore.case) tags.use = tolower(tags);
  if (!is.null(tags))
  { 
    searchType = match.arg(searchType);
    knownTagIndex = .indexedTags(collection, matchComponents = matchComponents, 
                                 impliedGroups = impliedGroups);
    if (ignore.case) knownTagIndex$tags = tolower(knownTagIndex$tags);
    if (exactMatch)
    {
       if (searchType=="any")
       {
          matchedSets = unique(knownTagIndex$index[ knownTagIndex$tags %in% tags.use]);
       } else  
          matchedSets = lapply(tags.use, function(x, knownTagIndex)
                                 knownTagIndex$index[ x== knownTagIndex$tags], knownTagIndex);
    } else {
      matchedSets = lapply(tags.use, function(x, knownTagIndex, fixed, ignore.case) 
            knownTagIndex$index[ grep(x,  knownTagIndex$tags, fixed = fixed, ignore.case = if (fixed) FALSE else ignore.case) ], 
            knownTagIndex, fixed, ignore.case)
    }
    if (!exactMatch || searchType!="any") matchedSets = lapply(matchedSets, unique);
    if (searchType=="all")
    {
      matches = multiIntersect(matchedSets)
    } else {
      matches = multiUnion(matchedSets);
    }
    if (invert) matches = setdiff(1:length(collection$dataSets), matches);
  }

  if (!is.null(firstDate) | !is.null(lastDate))
  {
    firstDate = .checkDate(firstDate, format = dateFormat);
    lastDate = .checkDate(lastDate, format = dateFormat);
    setDates = lapply(collection$dataSets, getElement, "lastModified");
    dateMatches = which(sapply(setDates, 
      function(d, first, last) (is.null(first) || first<=d) && (is.null(last) || last>=d),
      firstDate, lastDate));
    if (invert.date) dateMatches = setdiff( 1:length(collection$dataSets), dateMatches);
    matches = intersect(matches, dateMatches);
  }

  matches;
}

#============================================================================================
#
# addDataSetsToGroup
#
#============================================================================================

addDataSetsToGroup = function(
  collection, dataSetTags, groups, 
  matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                      "groupAlternateNames", "nameAndAlternates", "groupsAndAlternates"),
  searchType = c("any", "all"),
  invertSearch = FALSE, exactMatch = TRUE, fixed = TRUE, ignore.case = TRUE,
  firstDate = NULL, lastDate = NULL, 
  dateFormat = "%Y-%m-%d", invertDateSearch = FALSE,
  verbose = 1)
{
  if (!.isCollection(collection)) 
    stop("'collection' must be a valid collection.");

  # If groups were supplied, convert them to names
  groups = sapply(groups, function(x) { if (.isGroup(x)) x$name else x });

  # Check that all supplied 'groups' are character
  groupIsChar = sapply(groups, function(x) { inherits(x, "character") } ); 
  if (!all(groupIsChar)) 
    stop("Supplied 'groups' must be either group objects (as returned by 'newGroup()')\n",
         "or character strings giving the names of the groups.");

  # Find data sets that match the given tags
  index = .matchDataSets(collection, dataSetTags,
                         matchComponents = matchComponents,
                         searchType = searchType, invert = invertSearch,
                         exactMatch = exactMatch, fixed = fixed, ignore.case = ignore.case,
                         firstDate = firstDate, lastDate = lastDate, dateFormat = dateFormat,
                         invert.date = invertDateSearch);

  nMatched = length(index);
  if (verbose > 0)
    printFlush(WGCNA::spaste("addDataSetsToGroup: given tags matched ", nMatched, " data sets."));

  if (nMatched==0) return(collection);

  # Add the groups to the matched data sets
  for (set in index)
    collection$dataSets[[set]]$groups = unique(c( collection$dataSets[[set]]$groups, groups));

  collection;
}



#============================================================================================
#
# subsetCollection
#
#============================================================================================


# Restrict collection to only those that match given tags

subsetCollection = function(
   collection, tags = NULL, 
   matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                      "groupAlternateNames", "nameAndAlternates", "groupsAndAlternates"),
   searchType = c("any", "all"), invertSearch = FALSE,
   exactMatch = TRUE, fixed = TRUE, ignore.case = TRUE,
   firstDate = NULL, lastDate = NULL,
   dateFormat = "%Y-%m-%d", invertDateSearch = FALSE)
{
  if (is.null(tags) && is.null(firstDate) && is.null(lastDate)) return(collection);

  implied = impliedGroups(collection$groups, includeSelf = TRUE);

  keepIndex = .matchDataSets(collection, tags, 
                             matchComponents = matchComponents,
                             searchType = searchType, invert = invertSearch,
                             exactMatch = exactMatch, fixed = fixed, ignore.case = ignore.case,
                             firstDate = firstDate, lastDate = lastDate, dateFormat = dateFormat,
                             invert.date = invertDateSearch, 
                             impliedGroups = implied)
  outGroups0 = .knownGroups(collection$dataSets[keepIndex]);
  grp2implied = match(outGroups0, names(implied));
  fin = is.finite(grp2implied);
  outGroups = unique(c(outGroups0, unlist(implied[fin]))); 
  availableGroups = sapply(collection$groups, getElement, "name");
  keepGroups = match(outGroups, availableGroups)
  if (any(is.na(keepGroups)))
  {
    warning(immediate. = TRUE,
       "subsetCollection: data sets contain references to groups that cannot be found\n",
       "  in the group information of the given 'collection'.");
    keepGroups = keepGroups[is.finite(keepGroups)];
  }
  outColl = newCollection(dataSets = collection$dataSets[keepIndex],
                  identifiers = collection$identifiers,
                  weights = collection$weights, groups = collection$groups[keepGroups]);

  # Remove weights that are not being used anymore
  usedWeights = c(1, sort(.unique.nm(.weightIndex(outColl$dataSets))))
  translationTable = cbind(usedWeights, 1:length(usedWeights));
  outColl$weights = collection$weights[, usedWeights, drop = FALSE];

  nDataSets = length(outColl$dataSets);
  if (nDataSets > 0)
  {
    weightIndices = .weightIndex(outColl$dataSets);
    changeIndex = which(is.finite(weightIndices));
    changeTo = .translate(weightIndices, translationTable);
    for (set in changeIndex)
      outColl$dataSets[[set]]$weightIndex = changeTo[set]
  }

  outColl;
}

#===================================================================================================
#
# Utility functions for checking similarity of phrases
#
#===================================================================================================

.splitToCharacters = function(s, len = 1)
{
  n = nchar(s);
  out = character(n-len+1);
  for (i in 1:n-len+1)
    out[i] = substring(s, i,i+len-1);
  out;
}

.table.allLevels = function(x, levels)
{
  sapply(levels, function(l, x) sum(l==x), x);
}

.freqSimilarity = function(f1, f2)
{
  w = pmean(f1, f2);
  1-weighted.mean( abs(f1-f2)/(f1+f2), w);
}

.phraseSimilarity1 = function(phrases, maxTuple = 3)
{

  n = length(phrases)
  out = matrix(NA, n, n);
  diag(out)= 1;
  dimnames(out) = list(phrases, phrases);
  if (n==1) return(out);

  # Drop case:
  phrases.l = tolower(phrases)

  # Remove all white spaces from the start, end, and reduce multi-white space to one
  phrases.l = sapply(phrases, function(s) sub("^ +", "", sub(" +$", "", gsub("  +", " ", s))));

  tuples = lapply(phrases.l, function(s)
    lapply(1:maxTuple, function(l, s) .splitToCharacters(s, l), s));

  levels = lapply(tuples, lapply, unique);

  frequencies = mapply(function(tuples1, levels1)
     mapply(.table.allLevels, tuples1, levels1, SIMPLIFY = FALSE),
     tuples, levels, SIMPLIFY = FALSE);

  tupleSim = rep(0, maxTuple)
  tupleWeight = rep(0, maxTuple)
  for (p1 in 1:(n-1)) for (p2 in (p1+1):n)
  {
    for (l in 1:maxTuple)
    {
      commonLevels = union(levels[[p1]] [[l]], levels[[p2]] [[l]]);
      commonFreq1 = rep(0, length(commonLevels));
      commonFreq1[match(levels[[p1]] [[l]], commonLevels)] = frequencies[[p1]] [[l]];
      commonFreq2 = rep(0, length(commonLevels));
      commonFreq2[match(levels[[p2]] [[l]], commonLevels)] = frequencies[[p2]] [[l]];
      tupleWeight[l] = length(commonLevels);
      tupleSim[l] = .freqSimilarity(commonFreq1, commonFreq2);
    }
    out[p1, p2] = out[p2, p1] = weighted.mean(tupleSim, tupleWeight);
  }
  out;
}

.phraseSimilarity = function(phrases1, phrases2=NULL, maxTuple = 3)
{
  if (is.null(phrases2)) return(.phraseSimilarity1(phrases1, maxTuple));
  n1 = length(phrases1)
  n2 = length(phrases2)
  out = matrix(NA, n1, n2);
  dimnames(out) = list(phrases1, phrases2);

  # Drop case:
  phrases1.l = tolower(phrases1)
  phrases2.l = tolower(phrases2)

  # Remove all white spaces from the start, end, and reduce multi-white space to one
  phrases1.l = sapply(phrases1.l, function(s) sub("^ +", "", sub(" +$", "", gsub("  +", " ", s))));
  phrases2.l = sapply(phrases2.l, function(s) sub("^ +", "", sub(" +$", "", gsub("  +", " ", s))));

  tuples1 = lapply(phrases1.l, function(s)
    lapply(1:maxTuple, function(l, s) .splitToCharacters(s, l), s));
  tuples2 = lapply(phrases2.l, function(s)
    lapply(1:maxTuple, function(l, s) .splitToCharacters(s, l), s));

  levels1 = lapply(tuples1, lapply, unique);
  levels2 = lapply(tuples2, lapply, unique);

  frequencies1 = mapply(function(tuples1, levels1)
     mapply(.table.allLevels, tuples1, levels1, SIMPLIFY = FALSE),
     tuples1, levels1, SIMPLIFY = FALSE);
  frequencies2 = mapply(function(tuples1, levels1)
     mapply(.table.allLevels, tuples1, levels1, SIMPLIFY = FALSE),
     tuples2, levels2, SIMPLIFY = FALSE);
  
  
  tupleSim = rep(0, maxTuple)
  tupleWeight = rep(0, maxTuple) 
  for (p1 in 1:n1) for (p2 in 1:n2)
  {
    for (l in 1:maxTuple)
    {
      commonLevels = union(levels1[[p1]] [[l]], levels2[[p2]] [[l]]);
      commonFreq1 = rep(0, length(commonLevels));
      commonFreq1[match(levels1[[p1]] [[l]], commonLevels)] = frequencies1[[p1]] [[l]];
      commonFreq2 = rep(0, length(commonLevels));
      commonFreq2[match(levels2[[p2]] [[l]], commonLevels)] = frequencies2[[p2]] [[l]];
      tupleWeight[l] = length(commonLevels);
      tupleSim[l] = .freqSimilarity(commonFreq1, commonFreq2);
    }
    out[p1, p2] = weighted.mean(tupleSim, tupleWeight);
  }
  out;
}


#===================================================================================================
#
# checkGroups
#
#===================================================================================================

# Check the consistency of group information in a collection. Specifically, check that all groups mentioned in
# data sets are also entered in groups and vice-versa. Use fuzzy matching to suggest possible matches for
# inconsistencies.

checkGroups = function(collection, similarityCut = 0.5, ignore.case = TRUE) #, considerParents = TRUE)
{
  children = impliedGroups(collection$groups, get = "children");
  if (ignore.case) children.l = lapply(children, tolower) else childen.l = children

  dataSetGroups = .knownGroups(collection$dataSets);
  if (ignore.case) dataSetGroups.l = tolower(dataSetGroups) else dataSetGroups.l = dataSetGroups;
  groupNames = names(children);
  if (ignore.case) groupNames.l = tolower(groupNames) else groupNames.l = groupNames;

  dataInGroups = dataSetGroups.l %in% groupNames.l;
  groupsInData = sapply(children.l, function(grps) any(grps %in% dataSetGroups.l));

  notInGroups = dataSetGroups[!dataInGroups];
  notInData = groupNames[!groupsInData];

  if (length(notInGroups)>0 || length(notInData) > 0)
  {
    similarity = .phraseSimilarity(dataSetGroups, groupNames)
    closeMatchIndex.dataSets = lapply(as.data.frame(t(similarity[!dataInGroups, , drop = FALSE])), function(s)
      {
         order = order(s, decreasing = TRUE);
         order[ s[order] >=similarityCut];
      });
    closeMatchNames.dataSets = sapply(closeMatchIndex.dataSets, function(i, names)
                 base::paste(names[i], collapse = "|"), groupNames);

    closeMatchIndex.groups = lapply(as.data.frame(similarity[, !groupsInData, drop = FALSE]), function(s)
      {
         order = order(s, decreasing = TRUE);
         order[ s[order] >=similarityCut];
      });
    closeMatchNames.groups = sapply(closeMatchIndex.groups, function(i, names)
                 base::paste(names[i], collapse = "|"), dataSetGroups);

    notInGroups.df = data.frame(dataSetGroup = notInGroups,
                                closeGroupMatches = closeMatchNames.dataSets);
    notInData.df = data.frame(group = notInData,
                                closeDataSetMatches = closeMatchNames.groups);

    rownames(notInData.df) = rownames(notInGroups.df) = NULL;
    out = list(allOK = FALSE, dataSetGroupsWithoutRecords = notInGroups.df,
               groupsWithoutMemberDataSets = notInData.df);
  } else out = list(allOK = TRUE, 
                    dataSetGroupsWithoutRecords = data.frame(dataSetGroup = character(0),
                                                             closeGroupMatches = character(0)),
                    groupsWithoutMemberDataSets = data.frame(group = character(0),
                                                             closeDataSetMatches = character(0)));

  out;
}  


# Drop those elements of "groups" that tag no data sets
dropUnreferencedGroups = function(collection, ignore.case = TRUE, verbose = 1, indent = 0)
{
  spaces = indentSpaces(indent);
  children = impliedGroups(collection$groups, get = "children");
  if (ignore.case) children.l = lapply(children, tolower) else children.l = children;
  dataSetGroups = .knownGroups(collection$dataSets);
  if (ignore.case) dataSetGroups.l = tolower(dataSetGroups) else dataSetGroups.l = dataSetGroups;
  groupNames = names(children);
  if (ignore.case) groupNames.l = tolower(groupNames) else groupNames.l = groupNames;
  groupsInData = sapply(children.l, function(grps) any(grps %in% dataSetGroups.l));
  dataInGroups = dataSetGroups.l %in% groupNames.l;

  if (any(!dataInGroups)) 
    warning("Some tags in data sets have no corresponding group. Run checkGroups for details.");

  if (verbose > 0 && sum(!groupsInData) > 0)
    printFlush(spaste(spaces, "Removing ", sum(!groupsInData), " groups:\n",
            formatLabels(paste( spaste("'", groupNames[!groupsInData], "'"), collapse = ", "),
                         maxCharPerLine = 80, split = "', ")));

  collection$groups = collection$groups[groupsInData];
  collection;
}






#===================================================================================================
#
# renameGroup
#
#===================================================================================================

renameGroup = function(collection, oldGroup, newGroup)
{
  collection$dataSets = lapply(collection$dataSets, function(ds1, old, new)
       {
         ds1$groups[ ds1$groups==old] = new;
         ds1;
       }, oldGroup, newGroup);
  collection$groups = lapply(collection$groups, function(g1, old, new)
       {
         if (g1$name==old) g1$name = new;
         g1$alternateNames[g1$alternateNames==old] = new;
         g1$parents[g1$parents==old] = new;
         g1;
       }, oldGroup, newGroup)
  collection;
}

renameMultipleGroups = function(collection, old2newTable)
{
  collection$dataSets = lapply(collection$dataSets, function(ds1)
       {
         ds1$groups = .translateUsingTable(ds1$groups, old2newTable, keepUntranslated = TRUE);
         ds1;
       });
  collection$groups = lapply(collection$groups, function(g1)
       {
         g1$name = .translateUsingTable(g1$name, old2newTable, keepUntranslated = TRUE)
         g1$alternateNames = .translateUsingTable(g1$alternateNames, old2newTable, keepUntranslated = TRUE)
         g1$parents = .translateUsingTable(g1$parents, old2newTable, keepUntranslated = TRUE)
         g1;
       })

  collection;
}

standardizeGroupReferences = function(collection)
{
  translationTable = do.call(rbind, lapply(collection$groups, function(.grp)
    data.frame(X = c(.grp$name, .grp$alternateNames), Y = .grp$name)));
  rownames(translationTable) = NULL
  collection$dataSets = lapply(collection$dataSets, function(.set)
  {
    .set$groups = .translateUsingTable(.set$groups, translationTable, keepUntranslated = TRUE);
    .set;
  })
  collection;
}
  
  
#===================================================================================================
#
# collection2dataFrames
#
#===================================================================================================

# Convert a collection to data frame(s). The main frame will have columns for all of the data set
# information. In addition, there will be a list containing gene sets, and a list containing the
# properties.

# Since each gene set contains, in addition to the actual genes, also evidence and source, it is in general
# difficult to return all this information within the data set information. 

# The idea is to return data that can be easily saved in plain text tables. 

# This function turns a list of several components into a data frame. Vectors are turned into a single
# character string by collapsing all entries separated by sep.

.list2dataFrame = function(lst, components, sep, organismComp)
{
  if (length(lst)==0) return(NULL);

  lst.2 =  lapply(components, function(com, lst)
                  {
                    if (com=="organism")
                    {
                       as.character(sapply( lapply(lst, getElement, com), `[[`, organismComp))
                    } else if (com=="lastModified") {
                       as.character(sapply( lapply(lst, getElement, com), as.character));
                    } else {
                       xx = try(as.character(sapply( lapply(lst, getElement, com), base::paste, collapse = sep)))
                       if (inherits(xx, "try-error")) browser();
                       xx;
                    }
                  }, lst)
  names(lst.2) = NULL;
  out = as.data.frame(lst.2);
  names(out) = components;
  rownames(out) = NULL;
  out;
}

.addIDsToGeneData = function(geneSet)
{
  nGenes = nrow(geneSet$data);
  data.frame(ID = rep(geneSet$ID, nGenes),
             GeneSetName = rep(geneSet$name, nGenes),
             as.data.frame(lapply(geneSet$data, as.character)));
}

geneSetInformation = function(
  collection, tags = NULL, 
  matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                      "groupAlternateNames", "nameAndAlternates", "groupsAndAlternates"),
  searchType = c("any", "all"), invertSearch = FALSE,
  exactMatch = TRUE, fixed = TRUE, ignore.case = TRUE,
  firstDate = NULL, lastDate = NULL,
  dateFormat = "%Y-%m-%d", invertDateSearch = FALSE,
  sep = "|", organismName = c("common", "scientific", "short"))
{
  organismName = match.arg(organismName);
  organismComp = match(organismName, c("common", "scientific", "short"));

  index1 = .matchDataSets(collection, tags, 
                          matchComponents = matchComponents,
                          searchType = searchType, invert = invertSearch,
                          exactMatch = exactMatch, fixed= fixed, ignore.case = ignore.case,
                          firstDate = firstDate, lastDate = lastDate, dateFormat = dateFormat,
                          invert.date = invertDateSearch);
  if (length(index1)==0) return(NULL);

  index = index1[which(sapply(collection$dataSets[index1], .isGeneSet))];

  if (length(index) > 0)
  {
    dataSets.keep = collection$dataSets[index];
    sizes = sapply(dataSets.keep, function(x) length(unique(x$data$Entrez)));

    dataSetComponents.gs = c("ID", "name", "shortName", "description", "source",
                             "organism", "internalClassification", "groups", "lastModified",
                             "alternateNames", "externalDB", "externalAccession",
                             "webLink");

    geneSetIDF = cbind(.list2dataFrame(dataSets.keep, components = dataSetComponents.gs,
                                       sep = sep, organismComp = organismComp), nGenes = sizes);
    rownames(geneSetIDF) = NULL;
  } else
    geneSetIDF = NULL

  geneSetIDF;
}

collection2dataFrames = function(collection,
                                 sep = "|",
                                 propertyNamesFrom = "ID",
                                 organismName = c("common", "scientific", "short"))
{
  if (!.isCollection(collection))
     stop("Given 'collection' does not appear to be a valid collection.");

  organismName = match.arg(organismName);
  organismComp = match(organismName, c("common", "scientific", "short"));

  isGeneSet = sapply(collection$dataSets, .isGeneSet);
  isGeneProperty = sapply(collection$dataSets, .isGeneProperty);

  dataSetComponents.gs = c("ID", "name", "shortName", "description", "source", 
                        "organism", "internalClassification", "groups", "lastModified",
                        "alternateNames", "externalDB", "externalAccession",
                        "webLink");
  dataSetComponents.gp = c(dataSetComponents.gs, "weightIndex");

  if (length(propertyNamesFrom)!=1 || !(propertyNamesFrom %in% dataSetComponents.gs))
    stop("'propertyNamesFrom' must be a single character string specifying a valid component\n",
         "  of a gene data set. Useful values are 'ID', 'name' and 'shortName'.");

  if (any(isGeneSet))
  {
    dataSetDF.geneSets = .list2dataFrame(collection$dataSets[isGeneSet], components = dataSetComponents.gs, 
                                       sep = sep, organismComp = organismComp);
    geneSet.DFlist = lapply(collection$dataSets[isGeneSet], .addIDsToGeneData);
    geneSet.DF = do.call(rbind, geneSet.DFlist);
  } else {
    dataSetDF.geneSets = NULL;
    geneSet.DF = NULL;
  }

  if (any(isGeneProperty))
  {
    dataSetDF.geneProperties = .list2dataFrame(collection$dataSets[isGeneProperty], 
                                               components = dataSetComponents.gp, sep = sep);
    geneProperty.DF = as.data.frame(
                         lapply(collection$dataSets[isGeneProperty], getElement, "data"));
    names(geneProperty.DF) = make.names(sapply(collection$dataSets[isGeneProperty], getElement,
                                        propertyNamesFrom))
    rownames(geneProperty.DF) = collection$identifiers;
  } else {
    dataSetDF.geneProperties = NULL;
    geneProperty.DF = NULL;
  }

  weights.DF = collection$weights;

  # Create a data frame of group information

  
  groupComponents = c("name", "description", "source", "alternateNames", "parents")
  group.DF = .list2dataFrame(collection$groups, groupComponents, sep = sep);

  list(geneSetInfo = dataSetDF.geneSets,
       dataPropertyInfo = dataSetDF.geneProperties,
       geneSetContent = geneSet.DF,
       genePropertyContent = geneProperty.DF,
       genePropertyWeights = weights.DF,
       evidenceCodes = collection$evidenceCodes,
       groupInfo = group.DF);
}


geneSetsToHumanReadable = function(collection, sep = ", ")
{
  content = character(nDataSets(collection));
  oldOrg = "";
  for (gs in 1:nDataSets(collection))
  {
    entrez = collection$dataSets[[gs]]$data$Entrez;
    org = collection$dataSets[[gs]]$organism[1];
    if (org!=oldOrg)
    {
      orgSymbols = orgObject(organism = org, objectName = "SYMBOL", convertToList = TRUE, uniqueOnly = TRUE,
                      warnNonUnique = FALSE)
      oldOrg = org;
    }
    symbols = orgSymbols[match(as.character(entrez), names(orgSymbols))];
    content[gs] = base::paste(symbols, collapse = sep);
  }
  data.frame(setName = as.character(dataSetNames(collection)),
      setDescription = as.character(lapply(collection$dataSets, getElement, "description")),
      source = as.character(lapply(collection$dataSets, getElement, "source")),
      `Gene symbols` = content);
}




#================================================================================================
#
# collectionFromDataFrames
#
#================================================================================================

.matchNames = function(names, df, ignoreCase, stopOnMissing = TRUE, stopOnNULL = FALSE)
{
  if (length(names)==0)
    if (stopOnNULL) stop("'names' is missing.") else return(NULL)
  if (is.numeric(names)) return(names);
  names.ic = if (ignoreCase) toupper(names) else names;

  df.names.ic = if (ignoreCase) toupper(names(df)) else names(df);

  names2df = match(names.ic, df.names.ic);
  if (stopOnMissing && any(is.na(names2df)))
    stop("The following columns were not found:\n",
         base::paste(names[is.na(names2df)], collapse = ", "));

  names2df;
}

.strsplit2 = function(what, split, removeLeadingAndTrailingSpaces, prune)
{
  if (removeLeadingAndTrailingSpaces)
 {
    pattern = WGCNA::spaste(" *", split, " *");
  } else
    pattern = split;

  s = strsplit(as.character(what), split = pattern, fixed = FALSE);
  if (prune) 
    s = lapply(s, function(x) x[x!=""]);

  s
}
 
#Test: .strsplit2(c("aaa| bbb | ccc ", "|aa"), split = "\\|", TRUE, TRUE)

.fillContent = function(df, ..., ignoreCase)
{
  args = list(...)
  nArgs = length(args)

  df.names.ic = if (ignoreCase) toupper(names(df)) else names(df);
  nr = nrow(df);
  for (a in 1:nArgs)
  {
    name = if (ignoreCase) toupper(args[[a]] [1]) else args[[a]] [1];
    if (!name %in% df.names.ic)
    {
      printFlush(WGCNA::spaste("Column ", args[[a]] [1], " not found in the given data frame.\n",
                        "  Will create the column with default values."));
      df1 = data.frame(rep(args[[a]] [2], nr));
      names(df1) = name;
      df = cbind(df, df1);
    }
  }

  df;
}

groupsFromDataFrame = function(
   groupDF = NULL,
   groupDF.nameCol = "Name",
   groupDF.descriptionCol = "Description",
   groupDF.sourceCol = "Source",
   groupDF.alternateNamesCol = "AlternateNames",
   groupDF.parentsCol = "Parents",

   # general options
   ignoreCase = TRUE,
   multiEntrySep = "\\|",
   ignoreSpacesAroundSep = TRUE,
   pruneEmptyEntries = TRUE)
{
    names1 = c(groupDF.nameCol, groupDF.descriptionCol, groupDF.sourceCol, groupDF.alternateNamesCol,
               groupDF.parentsCol);
    names2df = .matchNames(names1, groupDF, ignoreCase);

    groups = mapply(newGroup, name = groupDF[, names2df[1]],
                    description = groupDF[, names2df[2]],
                    source = .strsplit2(groupDF[, names2df[3]], split = multiEntrySep,
                                  removeLeadingAndTrailingSpaces = ignoreSpacesAroundSep,
                                  prune = pruneEmptyEntries),
                    alternateNames = .strsplit2(groupDF[, names2df[4]], split = multiEntrySep,
                                  removeLeadingAndTrailingSpaces = ignoreSpacesAroundSep,
                                  prune = pruneEmptyEntries),
                    parents = .strsplit2(groupDF[, names2df[5]], split = multiEntrySep,
                                  removeLeadingAndTrailingSpaces = ignoreSpacesAroundSep,
                                  prune = pruneEmptyEntries),
                    SIMPLIFY = FALSE);
}


collectionFromDataFrames = function(
   # Gene set information
   geneSetInfoDF = NULL,
   geneSetInfoDF.IDcol = "ID",
   geneSetInfoDF.nameCol = "Name",
   geneSetInfoDF.shortNameCol = "ShortName",
   geneSetInfoDF.descriptionCol = "Description",
   geneSetInfoDF.organismCol = "Organism",
   geneSetInfoDF.sourceCol = "Source",
   geneSetInfoDF.groupsCol = "Groups",
   geneSetInfoDF.internalClassificationCol = "InternalClassification",
   geneSetInfoDF.lastModifiedCol = "LastModified",
   lastModifiedFormat = "%Y-%m-%d",
   geneSetInfoDF.alternateNamesCol = "AlternateNames",
   geneSetInfoDF.externalDBCol = "ExternalDB",
   geneSetInfoDF.externalAccessionCol = "ExternalAccession",
   geneSetInfoDF.webLinkCol = "webLink",

   # Gene set information
   geneSetContentDF = NULL,
   geneSetContentDF.nameCol = "GeneSetName",
   geneSetContentDF.entrezCol = "Entrez",
   geneSetContentDF.evidenceCol = "Evidence",
   geneSetContentDF.sourceCol = "Source",

   # property data set information
   propertyInfoDF = NULL,
   propertyInfoDF.IDcol = geneSetInfoDF.IDcol,
   propertyInfoDF.nameCol = geneSetInfoDF.nameCol,
   propertyInfoDF.shortNameCol = geneSetInfoDF.shortNameCol,
   propertyInfoDF.descriptionCol = geneSetInfoDF.descriptionCol,
   propertyInfoDF.organismCol = geneSetInfoDF.organismCol,
   propertyInfoDF.sourceCol = geneSetInfoDF.sourceCol,
   propertyInfoDF.groupsCol = geneSetInfoDF.groupsCol,
   propertyInfoDF.internalClassificationCol = geneSetInfoDF.internalClassificationCol,
   propertyInfoDF.lastModifiedCol = geneSetInfoDF.lastModifiedCol,
   propertyInfoDF.alternateNamesCol = geneSetInfoDF.alternateNamesCol,
   propertyInfoDF.externalDBCol = geneSetInfoDF.externalDBCol,
   propertyInfoDF.externalAccessionCol = geneSetInfoDF.externalAccessionCol,
   propertyInfoDF.webLinkCol = geneSetInfoDF.webLinkCol,
   propertyInfoDF.weightIndexCol = "weightIndex",

   # Gene property information
   propertyContentDF = NULL,
   propertyContentDF.entrezCol = NULL,
   propertyContentDF.colnamesMatch = "ID",

   # Weight information
   weightDF = NULL,

   # information about groups
   groupDF = NULL,
   groupDF.nameCol = "Name",
   groupDF.descriptionCol = "Description",
   groupDF.sourceCol = "Source",
   groupDF.alternateNamesCol = "AlternateNames",
   groupDF.parentsCol = "Parents",

   # general options
   ignoreCase = TRUE,
   multiEntrySep = "\\|",
   ignoreSpacesAroundSep = TRUE,
   pruneEmptyEntries = TRUE
)
{

  if (is.null(geneSetInfoDF) && is.null(propertyInfoDF))
    stop("At least on of 'geneSetInfoDF' or 'propertyInfoDF' must be given and non-NULL.");

  # dataSet information for gene sets

  if (!is.null(geneSetInfoDF))
  {
    if (is.null(geneSetContentDF))
       stop("When 'geneSetInfoDF' is non-NULL, valid 'geneSetContentDF' must be given as well.");

    geneSetInfoDF = .fillContent( geneSetInfoDF,
                                  c(geneSetInfoDF.lastModifiedCol, as.character(Sys.Date())),
                                  c(geneSetInfoDF.alternateNamesCol, ""),
                                  c(geneSetInfoDF.externalDBCol, ""),
                                  c(geneSetInfoDF.externalAccessionCol, ""),
                                  c(geneSetInfoDF.webLinkCol, ""),
                                  ignoreCase = ignoreCase);

    names1 = c(geneSetInfoDF.IDcol, 
               geneSetInfoDF.nameCol, 
               geneSetInfoDF.shortNameCol, 
               geneSetInfoDF.descriptionCol, 
               geneSetInfoDF.organismCol, 
               geneSetInfoDF.sourceCol,
               geneSetInfoDF.groupsCol,  
               geneSetInfoDF.internalClassificationCol,
               geneSetInfoDF.lastModifiedCol, 
               geneSetInfoDF.alternateNamesCol,
               geneSetInfoDF.externalDBCol, 
               geneSetInfoDF.externalAccessionCol, 
               geneSetInfoDF.webLinkCol);
    
    names2df = .matchNames(names1, geneSetInfoDF, ignoreCase);
    multiCols = match(c(geneSetInfoDF.groupsCol, geneSetInfoDF.internalClassificationCol,
                        geneSetInfoDF.alternateNamesCol), names1);
    geneDataSet.lst = as.list(geneSetInfoDF[, names2df]);
    geneDataSet.lst[multiCols] = lapply(geneDataSet.lst[multiCols], .strsplit2, split = multiEntrySep,
                                  removeLeadingAndTrailingSpaces = ignoreSpacesAroundSep,
                                  prune = pruneEmptyEntries);

    # Check that geneDataSet names are unique - this is a limitation stemming from the fact that gene sets
    # are matched to geneDataSets by name, not by ID.

    geneDataSet.presentNames = geneSetInfoDF[, names2df[2]];
    if (any(duplicated(geneDataSet.presentNames)))
    {
      stop("Gene set names in 'geneSetInfoDF' must be unique. The following names are duplicated:\n",
           base::paste(geneDataSet.presentNames[ duplicated(geneDataSet.presentNames)], collapse = ", "));
    }

    # Preprocess geneSetContentDF.
    names2 = c(geneSetContentDF.nameCol, geneSetContentDF.entrezCol, 
               geneSetContentDF.evidenceCol, geneSetContentDF.sourceCol);
    geneSet.names2df = .matchNames(names2, geneSetContentDF, ignoreCase);

    geneSet.presentNames = sort(unique(as.character(geneSetContentDF[, geneSet.names2df[1]])));

    # For now require that the names of gene sets in geneSetInfoDF and geneSetContentDF agree perfectly

    if (!isTRUE(all.equal(geneSet.presentNames, sort(geneDataSet.presentNames))))
    {
      common = intersect(geneSet.presentNames, geneDataSet.presentNames);
      stop("Gene set names in 'geneSetInfoDF' and 'geneSetContentDF' do not match.\n",
           "Gene set names present in 'geneSetInfoDF' and missing in 'geneSetContentDF':\n   ",
           base::paste(geneDataSet.presentNames[ !geneDataSet.presentNames%in%common], collapse = "\n   "),
           "\nGene set names present in 'geneSetContentDF' and missing in 'geneSetInfoDF':\n   ",
           base::paste(geneSet.presentNames[ !geneSet.presentNames%in%common], collapse = "\n   "));
    }

    geneSetContentDF.lst = tapply(1:nrow(geneSetContentDF), geneSetContentDF[, geneSet.names2df[1]],
                           function(index, df) df[index, ], geneSetContentDF[, geneSet.names2df[-1]], 
                           simplify = FALSE)

    order = match(geneDataSet.presentNames, geneSet.presentNames);
    if (any(is.na(order)))
      stop("Internal error: some entries of 'order' are missing.");

    geneSetContentDF.lst.ord = geneSetContentDF.lst[order];
    geneSetContentDF.lst.entrez = lapply(geneSetContentDF.lst.ord, `[[`, 1);
    geneSetContentDF.lst.evidence = lapply(geneSetContentDF.lst.ord, `[[`, 2 );
    geneSetContentDF.lst.source = lapply(geneSetContentDF.lst.ord, `[[`, 3 );

    # Create the gene sets
    geneSets = mapply(newGeneSet, geneEntrez = geneSetContentDF.lst.entrez,
                      geneEvidence = geneSetContentDF.lst.evidence,
                      geneSource = geneSetContentDF.lst.source,
                      ID = geneDataSet.lst[[1]],
                      name = geneDataSet.lst[[2]],
                      shortName = geneDataSet.lst[[3]],
                      description = geneDataSet.lst[[4]],
                      organism = geneDataSet.lst[[5]],
                      source = geneDataSet.lst[[6]],
                      groups = geneDataSet.lst[[7]],
                      internalClassification = geneDataSet.lst[[8]],
                      lastModified = geneDataSet.lst[[9]],
                      alternateNames = geneDataSet.lst[[10]],
                      externalDB = geneDataSet.lst[[11]],
                      externalAccession = geneDataSet.lst[[12]],
                      webLink = geneDataSet.lst[[13]],
                    MoreArgs = list( format = lastModifiedFormat),
                    SIMPLIFY = FALSE);
  } else
    geneSets = list();

  # Process groups
  if (!is.null(groupDF))
  {
    groups = groupsFromDataFrame(
       groupDF = groupDF,
       groupDF.nameCol = groupDF.nameCol,
       groupDF.descriptionCol = groupDF.descriptionCol,
       groupDF.sourceCol = groupDF.sourceCol,
       groupDF.alternateNamesCol = groupDF.alternateNamesCol,
       groupDF.parentsCol = groupDF.parentsCol,
       ignoreCase = ignoreCase,
       multiEntrySep = multiEntrySep,
       ignoreSpacesAroundSep = ignoreSpacesAroundSep,
       pruneEmptyEntries = pruneEmptyEntries);
  } else 
    groups = list();

  # Create the collection

  coll = newCollection(dataSets = geneSets, groups = groups)

  # dataSet information for gene properties
  if (!is.null(propertyInfoDF))
  {
    if (is.null(propertyContentDF))
       stop("When 'propertyInfoDF' is non-NULL, valid 'propertyContentDF' must be given as well.");


    propertyInfoDF = .fillContent( propertyInfoDF,
                                  c(propertyInfoDF.lastModifiedCol, as.character(Sys.Date())),
                                  c(propertyInfoDF.alternateNamesCol, ""),
                                  c(propertyInfoDF.externalDBCol, ""),
                                  c(propertyInfoDF.externalAccessionCol, ""),
                                  c(propertyInfoDF.webLinkCol, ""),
                                  ignoreCase = ignoreCase);


    names1 = c(propertyInfoDF.IDcol,
               propertyInfoDF.nameCol,
               propertyInfoDF.shortNameCol,
               propertyInfoDF.descriptionCol,
               propertyInfoDF.organismCol,
               propertyInfoDF.sourceCol,
               propertyInfoDF.groupsCol,
               propertyInfoDF.internalClassificationCol,
               propertyInfoDF.lastModifiedCol,
               propertyInfoDF.alternateNamesCol,
               propertyInfoDF.externalDBCol,
               propertyInfoDF.externalAccessionCol,
               propertyInfoDF.webLinkCol,
               propertyInfoDF.weightIndexCol);

    names2df = .matchNames(names1, propertyInfoDF, ignoreCase);
    multiCols = match(c(propertyInfoDF.groupsCol, propertyInfoDF.internalClassificationCol,
                        propertyInfoDF.alternateNamesCol), names1);
    propertyInfo.lst = as.list(propertyInfoDF[, names2df]);
    propertyInfo.lst[multiCols] = lapply(propertyInfo.lst[multiCols], .strsplit2, split = multiEntrySep,
                                  removeLeadingAndTrailingSpaces = ignoreSpacesAroundSep,
                                  prune = pruneEmptyEntries);

    propertyDataSet.presentIDs = propertyInfoDF[, names2df[1]];

    # Preprocess propertyContentDF.
    geneProperty.names2df = .matchNames(propertyContentDF.entrezCol, propertyContentDF, ignoreCase);
    identifiers = propertyContentDF[, geneProperty.names2df];

    propertyContentDF.r = propertyContentDF[, -geneProperty.names2df]
    geneProperty.presentIDs = sort(unique(names(propertyContentDF.r)));

    # For now require that the names of gene sets in propertyInfoDF and propertyContentDF agree perfectly

    if (!isTRUE(all.equal(geneProperty.presentIDs, sort(propertyDataSet.presentIDs))))
    {
      common = intersect(geneProperty.presentIDs, propertyDataSet.presentIDs);
      stop("Gene set IDs in 'propertyInfoDF' and 'propertyContentDF' do not match.\n",
           "Gene set IDs present in 'propertyInfoDF' and missing in 'propertyContentDF':\n   ",
           base::paste(propertyDataSet.presentIDs[ !propertyDataSet.presentIDs%in%common], collapse = "\n   "),
           "\nGene set IDs present in 'propertyContentDF' and missing in 'propertyInfoDF':\n   ",
           base::paste(geneProperty.presentIDs[ !geneProperty.presentIDs%in%common], collapse = "\n   "));
    }

    order = match(geneProperty.presentIDs, propertyDataSet.presentIDs);
    if (any(is.na(order)))
      stop("Internal error: some entries of 'order' are missing.");

    # Create the gene properties
    nProp = nrow(propertyInfoDF);
    coll = newGeneProperty(
              identifiers = identifiers,
              properties = propertyContentDF.r[ , order],
              weights = weightDF, 
              propertyNames = rep("", nProp),
              IDs = propertyInfo.lst[[1]],
              names = propertyInfo.lst[[2]],
              shortNames = propertyInfo.lst[[3]],
              descriptions = propertyInfo.lst[[4]],
              source = propertyInfo.lst[[6]],
              organism = propertyInfo.lst[[5]],
              groupBase = "",
              groupExtras = propertyInfo.lst[[7]],
              internalClassificationBase = "",
              internalClassificationExtras = propertyInfo.lst[[8]],

              lastModified = propertyInfo.lst[[9]],
              alternateNames = propertyInfo.lst[[10]],
              externalDB = propertyInfo.lst[[11]],
              externalAccession = propertyInfo.lst[[12]],
              webLink = propertyInfo.lst[[13]],

              collection = coll);
  }

  coll;
  
}

#===========================================================================================
#
# simplistic collection from a gene list
#
#==========================================================================================
# Assumes that column 1 are gene Entrez IDs, column 2 are gene set names.

collectionFromGeneLists = function(identifiers, setNames = NULL, organism, userName = "user",
                                   sortSets = TRUE)
{
  if (!is.null(setNames))
  {
     if (length(identifiers)!=length(setNames))
        stop("If given, 'setNames' must have the same length as 'identifiers'.")
  } else {
    # Assume identifiers is a list with one component per set whose names give the set names
    if (is.null(names(identifiers)))
      stop("If 'setNames' is not given, 'identifiers' must have valid names that will be used as set names.");
    setNames = unlist(mapply(function(name, set) rep(name, length(set)), 
                             names(identifiers), identifiers, 
                             SIMPLIFY = FALSE));
    identifiers = unlist(identifiers);
  }
  identifiers.all = identifiers;
  sets.all = setNames;
  keep = !is.na(identifiers.all);

  identifiers = identifiers.all[keep];
  sets = sets.all[keep];

  setNames = unique(sets);
  if (sortSets) setNames = sort(setNames);
  nSets = length(setNames);
  geneSets = list();

  IDs = WGCNA::spaste("User.", prependZeros(1:nSets));
  for (set in 1:nSets)
  {
    geneSets[[set]] = newGeneSet(geneEntrez = identifiers[sets==setNames[set]],
                                 geneEvidence = "other",
                                 geneSource = userName,
                                 ID = IDs[set], name = setNames[set], 
                                 description = WGCNA::spaste(setNames[set], " (supplied by ", userName, ")"),
                                 source = userName, organism = organism,
                                 internalClassification = c(userName, setNames[set]),
                                 groups = userName, lastModified = Sys.Date());
  }
  userGroups = list(newGroup(userName, WGCNA::spaste("Supplied by ", userName), userName));
  newCollection(geneSets, userGroups);
}


#================================================================================================
#
# Collection from classes (e.g., WGCNA modules)
#
#================================================================================================

collectionFromClasses = function(
  identifiers, 
  labels,
  organism,
  shortAnalysisName,
  
  # Data frame with the following components: 
  classDescription, 
  classCol = "class",
  nameCol,
  shortNameCol,
  descriptionCol,
  ignoreClasses = 0,
  evidence = "other",
  source = shortAnalysisName,
  IDBase, IDBaseSep = ".",
  groups,
  internalClassification,
  lastModified = Sys.Date(),
  format = "%Y-%m-%d")
{
  classes = sort(unique(labels));
  classes = setdiff(classes, ignoreClasses);
  
  nClasses = length(classes);
  width = min(3, nchar(nClasses));
  geneSets = lapply(classes, function(cl)
  {
    row = match(cl, classDescription[[classCol]]);
    if (is.na(row)) stop("Cannot find class ", cl, " among entries in column ", classCol);
    
    newGeneSet(
        geneEntrez = identifiers[labels==cl],
        geneEvidence = evidence, 
        geneSource = source,
        organism = organism,
        ID = WGCNA::spaste(IDBase, IDBaseSep, prependZeros(match(cl, classes), width)),
        name = WGCNA::spaste(classDescription[[nameCol]] [row], " [", shortAnalysisName, "]"),
        description = WGCNA::spaste(classDescription[[descriptionCol]] [row]),
        source = source,
        groups = sapply(groups, getElement, "name"),
        internalClassification = internalClassification,
        lastModified = lastModified,
        format = format);
  })

  newCollection(dataSets = geneSets, groups = groups);
}

# Adapt enrichmentLabels for use as classDescription in collectionFromClasses

enrichmentLabelsToClassDescription = function(eLabels,
   shortAnalysisName,
   fullAnalysisName,
   shortNameFirst = TRUE)
{

  stripLeadingPart = function(s, sep, emptyReplacement = "")
  {
    if (length(s)>1) return(sapply(s, stripLeadingPart, sep));

    colonPos = regexpr(sep, s);
    if (colonPos==-1) s1 = "" else s1 = substring(s, colonPos+nchar(sep));
    s1 = sub("^  *", "", sub("  *$", "", s1));
    if (s1=="") s1 = emptyReplacement;
    s1;
  }

  if (shortNameFirst) 
  {
     eLabels$enrichmentLabel = WGCNA::spaste(shortAnalysisName, " ", eLabels$enrichmentLabel)
  } else
     eLabels$enrichmentLabel = WGCNA::spaste(eLabels$enrichmentLabel, " [", shortAnalysisName, "]")

  eLabels$description = WGCNA::spaste("Module ", eLabels$class, " determined in ", fullAnalysisName,
                         ". Highest enriched terms (CT=Cell Type, BR=Brain Region, GO=Gene Ontology): ",
                         stripLeadingPart(eLabels$enrichmentLabel, "):", "(none)"));
  eLabels;
}

#================================================================================================
#
# Utility functions for enrichment analysis
#
#================================================================================================

# This function removes duplicate identifiers whose labels are not all the same, ignoring 'ignoredLabels'.
# This function assumes that labels are a single vector, and sets the labels of removed identifiers to NA.

.removeDuplicatesWithDifferentLabels = function(labels, identifiers, ignoredLabels, 
                                      missingLabel)
{
  #notIgnored = !(labels %in% c(ignoredLabels, missingLabel));
  #labels.ni = labels[notIgnored];
  #ids.ni = identifiers[notIgnored];

  duplicatedIDs = setdiff(unique(identifiers[duplicated(identifiers)]), c(ignoredLabels, missingLabel))
  if (length(duplicatedIDs)==0) return(labels);
  duplicated = identifiers %in% duplicatedIDs;

  labels.dup = labels[ duplicated ];    
  ids.dup = identifiers[ duplicated ];
  different = !sapply(tapply(labels.dup, ids.dup, function(x) all(x==x[1], na.rm = TRUE)), identity);   
  if (any(different))
  {
    removeIDs = names(different)[different];
    labels[ identifiers %in% removeIDs ] = missingLabel;
  }
  labels;
}


# Calculates overlap p-values

.fisherPvalue = function(nCommon, nInCat, nOutCat, nInClass, alternative = c("greater", "two.sided", "less"))
{
  alternative = match.arg(alternative) 
  nComparisons = length(nCommon);
  if (alternative %in% c("greater", "two.sided")) {
    result.g = rep(1, nComparisons)
    calculate = which(nCommon > 0);
    if (length(calculate) > 0) 
      result.g[calculate] = phyper(nCommon[calculate] - 1, nInCat[calculate], nOutCat[calculate], nInClass[calculate],  
            lower.tail = FALSE)
  } 
  if (alternative %in% c("less", "two.sided")) {
      result.l = phyper(nCommon, nInCat, nOutCat, nInClass, lower.tail = TRUE)
  } 
  if (alternative=="greater") {
    result = result.g;
  } else if (alternative=="less") {
    result = result.l;
  } else if (alternative=="two.sided") {
    result = 2*pmin(result.l, result.g);
    result[result > 1] = 1;
  } else stop(paste0("Unrecognized 'alternative' ", alternative));
  result;
}


#.oddsRatio = function(nCommon, nInCat, nOutCat, nInClass)
#{
  # Calculating odds ratio:
  # nCommon 		nInClass - nCommon
  # nInCat - nCommon	nOutCat - nInClass + nCommon

  # sample estimate OR = ad/bc


#===================================================================================================
#
# enrichmentAnalysis
#
#===================================================================================================
   
# Main enrichment analysis function. 

# Note: removal of duplicates is mandatory since I am assuming that each entrez code is present only once.

enrichmentAnalysis = function(
               classLabels, 
               identifiers,
               # Alternative input specification
               active = NULL,
               inactive = NULL,
               activeNames = names(active),

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

               # safe-to-use missing label (more convenient than NA)
               safeToUseMissingLabel = -394879358L,

               # How much output to generate
               nBestDataSets = 10, 
               threshold = 0.05,
               thresholdType = c("Bonferroni", "FDR", "nominal"), 
              
               # Calculate multiple testing corrections? 
               getBonferroniCorrection = TRUE,
               getFDR = TRUE,

               # Return overlap genes in the enrichment table?
               getOverlapEntrez = TRUE,
               getOverlapSymbols = FALSE,
               maxReportedOverlapGenes = 50,

               # Output options
               entrySeparator = "|",
               groupSeparator = entrySeparator,
               geneSeparator = entrySeparator,
               classColName = "class",
               combineEnrichmentTables = !is.null(active),
               combinedMultipleTestingCorrection = TRUE,

               # Return gene set details?
               getDataSetDetails = TRUE,

               # Diagnostic message options
               verbose = 1, indent = 0 )
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

   if (length(useGroups)==0) 
   {
     useGroups = knownGroups(refCollection);
   } else
     refCollection = subsetCollection(refCollection, useGroups); 

   impliedGroups1 = impliedGroups(refCollection$groups);

   if (useBackground=="allOrgGenes")
   {
     referenceBackground = allOrgGenes(organism)
   } else
     referenceBackground = allGeneSetGenes(refCollection, evidence = useEvidence);

   if (getOverlapSymbols)
     backgroundSymbols = convert2symbol(entrez = referenceBackground, organism = organism)

   if (!is.null(active))
   {
     if (is.null(inactive))
     { 
       if (useBackground %in% c("intersection", "given")) 
       {
         stop(base::paste0("When 'active' list is specified, 'inactive' must be specified as well.\n",
                     "  Alternatively, you can set 'useBackground' to 'reference' or 'allOrgGenes'\n",
                     "  BUT be aware that resulting enrichment statistics may be inflated.")); 
       } else
         inactive = referenceBackground;
     }
     privateInactive = !is.atomic(inactive)
     if (privateInactive)
     {
       if (length(inactive)!=length(active)) 
         stop("When 'inactive' is a list, it must have the same length as 'active'.");
     }
     identifiers = sort(unique(c(unlist(active), unlist(inactive))))
     if (is.atomic(active)) active = list(active);
     if (is.null(activeNames)) 
       activeNames = spaste("active.", prependZeros(1:length(active), nchar(length(active))));
     classLabels = mapply(function(a1, n1) c("inactive", n1) [ as.numeric(identifiers %in% a1) + 1],
                    active, activeNames, SIMPLIFY = FALSE);
     if (privateInactive)
       classLabels = mapply(function(lab1, active1, inactive1) 
       {
         lab1[!identifiers %in% c(active1, inactive1)] = safeToUseMissingLabel;
         lab1;
       }, classLabels, active, inactive, SIMPLIFY = FALSE);

     classLabels = do.call(cbind, classLabels);
     ignoreLabels = "inactive";
   }

   if (length(safeToUseMissingLabel)!=1)
     stop("'safeToUseMissingLabel' must by a single identifier (preferably an integer number).");

   if (is.numeric(classLabels)) 
   {
     if (!is.finite(safeToUseMissingLabel)) 
       stop("'safeToUseMissingLabel' must be finite.");
   } else {
     if (is.na(safeToUseMissingLabel))
       stop("'safeToUseMissingLabel' cannot be NA.");
   }

   classLabels = as.matrix(classLabels);
   classLabels[is.na(classLabels)] = safeToUseMissingLabel

   # Remove any identifiers that are missing

   keep = !is.na(identifiers);
   classLabels = classLabels[keep, , drop = FALSE];
   identifiers = identifiers[keep];

   nSets = ncol(classLabels);
   nGivenRaw = nrow(classLabels);

   if (removeDuplicatesInDifferentClasses)
      classLabels = as.matrix(apply(classLabels, 2, .removeDuplicatesWithDifferentLabels,
                              identifiers = identifiers,  ignoredLabels = ignoreLabels,
                                    missingLabel = safeToUseMissingLabel));

   entrezLists = geneLists(refCollection, evidence = useEvidence, simplify = FALSE);
   dataSetIDs = sapply(refCollection$dataSets, getElement, "ID");
   dataSetNames = sapply(refCollection$dataSets, getElement, "name");
   shortSetNames = sapply(refCollection$dataSets, getElement, "shortName");
   dataSetGroups = sapply(refCollection$dataSets, function(ds) 
            base::paste(unlist(impliedGroups1[ds$groups]), collapse = groupSeparator) );

   nDataSets = length(entrezLists)
   dataSetSizes = sapply(entrezLists, length);

   if (all(dataSetSizes==0))
     warning(immediate. = TRUE, 
             "Gene lists are empty (possibly because of restricting evidence).");

   if (nBestDataSets > nDataSets) nBestDataSets = nDataSets;

   common = intersect(referenceBackground, identifiers);

   if (length(common)==0)
     stop("There are no common identifiers in the supplied gene sets and the input 'identifiers'.\n",
          "   Please make sure the 'refCollection' and 'identifiers' are Entrez identifiers and\n",
          "   correspond to the same organism.");

   identIsInCollection = identifiers %in% referenceBackground;

   ### stuff defined here is undefined if any(classLabels==safeToUseMissingLabel) 
   ### 
   if (!any(classLabels==safeToUseMissingLabel))
   {
     validIdentifiers.0 = unique(identifiers);
     identIsInCollection.unique.0 = validIdentifiers.0 %in% referenceBackground;
     if (useBackground %in% c("intersection", "given")) 
     {
        if (verbose > 2)
          printFlush(base::paste(spaces, "   ..restricting collection to given identifiers.."));
        effEntrezLists.0 = lapply(entrezLists, intersect, validIdentifiers.0);
        effDataSetSizes.0 = sapply(effEntrezLists.0, length);
     } else {
        effEntrezLists.0 = entrezLists;
        effDataSetSizes.0 = dataSetSizes;
     }
     anyInvalid = FALSE;
   } else {
     anyInvalid = TRUE;
   }

   setResults = list();

   ignoreLabels = c(ignoreLabels, safeToUseMissingLabel);
   for (set in 1:nSets)
   {
      if (verbose > 0)
        printFlush(base::paste(spaces, " ..working on label set", set, ".."));
      labelLevels = levels(factor(classLabels[, set]));

      keep = !(labelLevels %in% ignoreLabels);
      if (sum(keep)==0)
          stop("No class labels were kept after removing labels that are supposed to be ignored.");
      labelLevels = labelLevels[keep]
      nLabelLevels = length(labelLevels);

      nComparisons = nLabelLevels * nDataSets;

      modCodes = list();

      if (useBackground %in% c("intersection", "reference", "allOrgGenes"))
      {
        for (ll in 1:nLabelLevels)
           modCodes[[ll]] = unique(identifiers[classLabels[, set]==labelLevels[ll] & identIsInCollection]);
      } else if (useBackground == "given") {
        for (ll in 1:nLabelLevels)
           modCodes[[ll]] = unique(identifiers[classLabels[, set]==labelLevels[ll]])
      }
      modSizes.eff = sapply(modCodes, length)

      # May need to restrict the refCollection gene sets to genes in common with identifiers

      validLabels = (classLabels[, set]!=safeToUseMissingLabel)
      if (anyInvalid)
      {
        validIdentifiers = unique(identifiers[validLabels]);
        identIsInCollection.unique = validIdentifiers %in% referenceBackground;
        if (useBackground %in% c("intersection", "given")) 
        {
           if (verbose > 2)
             printFlush(base::paste(spaces, "   ..restricting collection to given identifiers.."));
           effEntrezLists = lapply(entrezLists, intersect, validIdentifiers);
           effDataSetSizes = sapply(effEntrezLists, length);
        } else {
           effEntrezLists = entrezLists;
           effDataSetSizes = dataSetSizes;
        }
      } else {
        validIdentifiers = validIdentifiers.0;
        identIsInCollection.unique = identIsInCollection.unique.0;
        effEntrezLists = effEntrezLists.0;
        effDataSetSizes = effDataSetSizes.0;
      }

      # Calculate the total number of background genes
      nBackgroundGenes = switch(useBackground,
            intersection = sum(identIsInCollection.unique),
            given = length(validIdentifiers),
            reference = length(referenceBackground),
            allOrgGenes = length(referenceBackground) );

      enrichmentIsValid = nBackgroundGenes > 0
      if (nBackgroundGenes==0)
      {
        warning(immediate. = TRUE, 
           base::paste0("enrichmentAnalysis: number of effective background genes for set ", set, " is zero.\n",
                  "      This means that after restricting to the appropriate background type,\n",
                  "      there are no genes left."));
      }

      # these will hold the main results
      enrichment = matrix(if (enrichmentIsValid) 1 else NA, nDataSets, nLabelLevels);
      countsInDataSet = matrix(0, nDataSets, nLabelLevels);
      enrichmentRatio.mat = matrix(0, nDataSets, nLabelLevels); ## a.k.a. effect size

      keepSets = which(effDataSetSizes > 0);
      effNDataSets = length(keepSets);
      if (effNDataSets > 0)
      {
        entrezLists.2 = effEntrezLists[keepSets];
        dataSetSizes.2 = effDataSetSizes[keepSets];

        countsInDataSet.1 = matrix(0, effNDataSets, nLabelLevels);

        if (verbose > 2)
           printFlush(base::paste(spaces, "   ..calculating overlaps.."));

        for (ll in 1:nLabelLevels)
          countsInDataSet.1[, ll] = sapply(entrezLists.2, function(x) sum(x%in% modCodes[[ll]]));
             ### Note: entrezLists.2 must not contain duplicates for this to be valid.

        countsInDataSet[keepSets, ] = countsInDataSet.1;

        # Flatten countsInDataSet to a vector
        dim(countsInDataSet.1) = NULL;
        dataSetSizes.ext = rep(dataSetSizes.2, nLabelLevels);
        modSizes.ext = rep(modSizes.eff, each = effNDataSets);

        enrichmentRatio1 = countsInDataSet.1 * nBackgroundGenes / (dataSetSizes.ext * modSizes.ext)
        enrichmentRatio.mat[keepSets, ] = enrichmentRatio1;

        if (verbose > 1)
           printFlush(base::paste(spaces, "   ..calculating enrichments (this may take a while).."));

        enrichment[keepSets, ] = .fisherPvalue( nCommon = countsInDataSet.1,
                                      nInCat = dataSetSizes.ext,
                                      nOutCat = nBackgroundGenes - dataSetSizes.ext,
                                      nInClass = modSizes.ext);
      }
    

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
          dimnames(enrichment.FDR) = dimnames(enrichmentRatio.mat) = list (dataSetIDs, labelLevels);

      modSizes = table(classLabels[ !(classLabels[, set] %in% ignoreLabels), set]);

      enrichment.threshold = switch(thresholdType, 
              nominal = enrichment,
              Bonferroni = enrichment.Bonf,
              FDR = enrichment.FDR);
      enrichmentTable = NULL;

      if (getDataSetDetails) dataSetDetails = list(); 
      if (enrichmentIsValid && (!is.null(threshold) || nBestDataSets > 0))
      {
        if (verbose > 1) 
          printFlush(base::paste(spaces, "   ..putting together terms with highest enrichment significance.."));
        order = apply(enrichment, 2, order);
        dim(order) = dim(enrichment);

        for (ll in 1:nLabelLevels)
        {
           if (getDataSetDetails) dataSetDetails[[ll]] = list();
           if (!is.null(threshold))
           {
              reportTerms = c(1:nDataSets)[enrichment.threshold[, ll] < threshold];
           } else 
              reportTerms = numeric(0);
           if (nBestDataSets > 0)
              reportTerms = unique(c(reportTerms, c(1:nDataSets)[order[1:nBestDataSets, ll]]));
           nRepTerms = length(reportTerms);

           if (nRepTerms > 0)
           {
              # Order the terms by increasing enrichment p-value
              reportTerms = reportTerms[ order(enrichment[ reportTerms, ll ]) ];

              if (getOverlapEntrez | getOverlapSymbols)
              {
                 overlapEntrez = lapply(effEntrezLists[reportTerms], intersect, modCodes[[ll]]);
                 lengths = sapply(overlapEntrez, length);
                 long = lengths>maxReportedOverlapGenes;
                 if (any(long)) for (i in which(long))
                   overlapEntrez[[i]] = WGCNA::spaste("(More than ", maxReportedOverlapGenes, " overlapping genes)");
                 if (any(!long))
                 {
                   index = which(!long)
                   if (getOverlapSymbols)
                   {
                      overlapSymbols = lapply(overlapEntrez[index], function(x)
                             backgroundSymbols[ match(x, referenceBackground)])
                      if (getOverlapEntrez) {
                        overlapEntrez[index] = mapply(function(x, y) WGCNA::spaste(x, " (", y, ")"),
                                               overlapEntrez[index], overlapSymbols, SIMPLIFY = FALSE)
                      } else 
                        overlapEntrez[index] = overlapSymbols;
                   }
                   overlapEntrez[index] = sapply(overlapEntrez[index], base::paste, collapse= geneSeparator);
                }
                overlapEntrez = sapply(overlapEntrez, identity);
              } else overlapEntrez = NULL;
              expectedFrac = effDataSetSizes[reportTerms]/nBackgroundGenes
              observedFrac = countsInDataSet[reportTerms, ll]/modSizes.eff[ll];
              enrichmentRatio1 = observedFrac/expectedFrac;
              enrTab0 = list(class = rep(labelLevels[ll], nRepTerms), 
                        rank = seq(length.out = nRepTerms),
                        dataSetID = dataSetIDs[reportTerms],
                        dataSetName = dataSetNames[reportTerms],
                        inGroups = dataSetGroups[reportTerms],
                        pValue = enrichment[reportTerms, ll],
                        Bonferroni = if (getBonferroniCorrection) enrichment.Bonf[reportTerms, ll] else NULL,
                        FDR = if (getFDR) enrichment.FDR[reportTerms, ll] else NULL,
                        nCommonGenes = countsInDataSet[reportTerms, ll],
                        fracOfEffectiveClassSize = observedFrac,
                        expectedFracOfEffectiveClassSize = expectedFrac,
                        enrichmentRatio = enrichmentRatio1,
                        classSize = rep(modSizes[ll], nRepTerms),
                        effectiveClassSize = rep(modSizes.eff[ll], nRepTerms), 
                        fracOfEffectiveSetSize = 
                             countsInDataSet[reportTerms, ll]/effDataSetSizes[reportTerms],
                        effectiveSetSize = effDataSetSizes[reportTerms],
                        shortDataSetName = shortSetNames[reportTerms],
                        overlapGenes = overlapEntrez);
              l1 = sapply(enrTab0, length);
              enrTab = as.data.frame(enrTab0[l1>0], row.names = enrTab0$rank);
              #enrTab = try(as.data.frame(enrTab0[l1>0], row.names = enrTab0$rank));
              #if (inherits(enrTab, "try-error")) browser();
              names(enrTab)[1] = classColName;
              rownames(enrTab) = NULL;
              if (is.null(enrichmentTable))
              {
                 enrichmentTable = enrTab;
              } else 
                 enrichmentTable = rbind(enrichmentTable, enrTab);

              if (getDataSetDetails)
              {
                identifiers.valid = identifiers;
                identifiers.valid[ classLabels[, set] %in% ignoreLabels] = NA;
                for (rci in seq(length.out = nRepTerms))
                {
                   gs = reportTerms[ rci ];
                   geneCodes = intersect(modCodes[[ll]], effEntrezLists[[gs]])
                   dataSetDetails[[ll]][[rci]] = list(
                           dataSetID = enrTab$dataSetID[ rci ],
                           dataSetName = enrTab$dataSetName[rci],
                           dataSetDescription = refCollection$dataSets[[gs]]$description,
                           dataSetGroups = unlist(impliedGroups1[refCollection$dataSets[[gs]]$groups]),
                           enrichmentP = enrTab$pValue[rci],
                           commonGeneEntrez = geneCodes,
                           commonGenePositions = match(geneCodes, identifiers.valid));
                }
                names(dataSetDetails[[ll]]) = make.names(enrTab$dataSetID);
                names(dataSetDetails)[ll] = make.names(labelLevels[ll]);
              }
           }
         }
      }
           
      setResults[[set]] = list(enrichmentIsValid = enrichmentIsValid,
                               enrichmentTable = enrichmentTable,
                               pValues = enrichment,
                               Bonferroni = if (getBonferroniCorrection) enrichment.Bonf else NULL,
                               FDR = if (getFDR) enrichment.FDR else NULL,
                               countsInDataSet = countsInDataSet, 
                               enrichmentRatio = enrichmentRatio.mat,
                               effectiveBackgroundSize = nBackgroundGenes,
                               effectiveDataSetSizes = effDataSetSizes,
                               effectiveClassSizes = modSizes.eff,
                               dataSetDetails = if (getDataSetDetails) dataSetDetails else NULL);
                               
   }

     

   
   if (nSets==1)
   {
      c(setResults[[1]],
        list( identifierIsInCollection = identIsInCollection,
              effectiveClassLabels = classLabels));
           
   } else {
      if (combineEnrichmentTables && nSets > 1)
      {
        keep = sapply(setResults, getElement, "enrichmentIsValid");
        enrTable1 = do.call(rbind, lapply(setResults[keep], getElement, "enrichmentTable"));
        countsInDataSet = do.call(cbind, lapply(setResults[keep], getElement, "countsInDataSet"))
        enrichmentRatio = do.call(cbind, lapply(setResults[keep], getElement, "enrichmentRatio"))
        p1 = do.call(cbind, lapply(setResults[keep], getElement, "pValues"));
        # Columns are classes, rows are gene sets
        n = prod(dim(p1));
        pb1 = n*p1;
        pb1[pb1>1] = 1;
        Bonf1 = if (getBonferroniCorrection) pb1 else NULL;
        if (getFDR)
        {
          FDR1 = p.adjust(c(p1), method = "fdr");
          dim(FDR1) = dim(p1);
          dimnames(FDR1) = dimnames(p1);
        } else
          FDR1 = NULL;
        if (combinedMultipleTestingCorrection)
        {
          enrTable1$Bonferroni = enrTable1$pValue * n;
          enrTable1$Bonferroni[enrTable1$Bonferroni > 1] = 1;
          if (getFDR) {
            tableFDR = numeric(nrow(enrTable1));
            for (col in 1:ncol(p1))
            {
              tabRows = which(enrTable1$class==colnames(p1)[col]);
              if (any(!enrTable1$dataSetID[ tabRows] %in% rownames(FDR1)))
              {
                stop("Internal error when filling common table FDR.");
              }
              tableFDR[tabRows] = FDR1[as.character(enrTable1$dataSetID[ tabRows]), col];
            }
            enrTable1$FDR = tableFDR
          }
        }
        out = list(enrichmentTable = enrTable1,
                   pValues = p1,
                   Bonferroni = Bonf1,
                   FDR = FDR1,
                   countsInDataSet = countsInDataSet);
      } else out = list();
      c(out, list(setResults = setResults,
            identifierIsInCollection = identIsInCollection, effectiveClassLabels = classLabels));
   }
}

#======================================================================================================
#
# orgObject
#
#======================================================================================================

# This is a convenience function to retrieve a specified object from an appropriate org.Xx.eg.db package
orgExtension = function(organism)
{
  short = organismShorthand(organism);
  if (short=="Sc") 
  {
    ext = "sgd"
  } else
    ext = "eg"
  ext;
}

orgObject = function(organism, objectName, convertToList = TRUE, uniqueOnly = FALSE,
                     warnNonUnique = TRUE)
{
  short = organismShorthand(organism);
  ext = base::paste0(".", orgExtension(organism))
  fullName = base::paste0("org.", short, ext, objectName);
  lib = base::paste0("org.", short, ext, ".db");
  require(lib, character.only = TRUE);
  obj = get(fullName);
  if (convertToList)
  {
    mapped = mappedkeys(obj);
    obj = AnnotationDbi::as.list(obj[mapped]);
    if (uniqueOnly)
    {
      lens = sapply(obj, length);
      keep = lens==1;
      if (warnNonUnique & any(!keep)) warning(immediate. = TRUE,
         base::paste0("orgObject: some elements of ", objectName, " for organism ", organism, 
                " are either empty or non-unique."));
      obj = unlist(obj[keep]);
    }
  }
  obj;
}

processAnnotationObject = function(object, uniqueOnly = TRUE,
                               warnNonUnique = FALSE, reference = NULL)
{
  m = mappedkeys(object);
  xx = as.list(object[m]);

  if (!is.null(reference))
  {
    xx = c(xx, UnknownID = list(NA));
    ref2xx = match(reference, names(xx));
    ref2xx[is.na(ref2xx)] = length(xx);

    xx = xx[ref2xx]
  }

  lengths = sapply(xx, length);
  if (any(lengths!=1) && uniqueOnly)
  {
    if (is.null(reference))
    {
      xx = xx[lengths==1];
    } else
      xx[lengths!=1] = NA;
    if (warnNonUnique)
      warning("getAnnotationObject: Non-unique entries have been dropped.");
  }

  if (uniqueOnly) xx = unlist(xx);

  xx;
}

# Convenience function to create a gene annotation data frame from Entrez IDs.

geneAnnotationFromEntrez = function(entrez, organism, includePosition = FALSE)
{
  short = organismShorthand(organism);
  ext = base::paste0(".", orgExtension(organism))
  lib = base::paste0("org.", short, ext, ".db");
  if (!require(lib, character.only = TRUE))
    stop("The required package ", lib, " could not be loaded. Please install it and try again.");

  symbolObj = get(base::paste0("org.", short, ext, "SYMBOL"));
  nameObj = get(base::paste0("org.", short, ext, "GENENAME"));
  out = data.frame(Entrez = entrez,
             Symbol = processAnnotationObject(symbolObj, reference = entrez),
             Name = processAnnotationObject(nameObj, reference = entrez));
  if (includePosition)
  {
    chrObj = get(base::paste0("org.", short, ext, "CHR"))
    posObj = get(base::paste0("org.", short, ext, "CHRLOC"))
    out = cbind(out, Chr = processAnnotationObject(chrObj, reference = entrez),
                Loc = processAnnotationObject(posObj, reference = entrez));
  }
  out;
}




#======================================================================================================
#
# entrez2symbol
#
#======================================================================================================

# Here entrez can be a (nested) list of vectors
convert2symbol = function(entrez, organism)
{
  if (length(entrez)==0) return(character(0));
  e2s = orgObject(organism, "SYMBOL", convertToList = TRUE, uniqueOnly = TRUE);
  if (is.atomic(entrez))
  {
    entrez = list(entrez);
    atomic = TRUE;
  } else
    atomic = FALSE;
  out = lapply(entrez, .translateUsingTable, data.frame(Entrez = as.numeric(names(e2s)), Symbol = as.character(e2s)));
  if (atomic) out = out[[1]];
  out;
}

if (FALSE)
{
  #test

  xx = orgObject("mouse", "ALIAS2EG", convertToList = TRUE, uniqueOnly = FALSE, warnNonUnique = FALSE)

  convert2entrez(organism = "mouse", symbol = 
         list(c("vlcad", "Acadvl", "AI196007", "Bcd-1", "Bcd1"), list(c("Abpa27", "Sal-1", "Abpa"), 
         c("AI413825", "Abc", "Abc2", "D2H0S1474", "Abca2"))));
}

convert2entrez = function(organism, symbol = NULL, refSeq = NULL, ignoreCase = FALSE,
                          dropNonConverted = FALSE, useAlias = TRUE,
                          noMatch = NA, multipleMatches = NA,
                          annotation = NULL,
                          symbolCol = "symbol",
                          refSeqCol = "refSeq",
                          aliasCol = "synonyms",
                          aliasSep = "; ",
                          aliasSplit.fixed = TRUE,
                          entrezCol = "Entrez",
                          ...)
{
  .noMatch1 = -43237387L;
  .multiMatch1 = -84763262L;
  #if (useAlias)
  #  printFlush("FYI: symbol matching includes aliases. To get old behaviour, use 'useAlias = FALSE'.");

  if (!is.null(symbol))
  {
    identifiers = symbol;
    objName = "SYMBOL2EG";
    sourceCol = symbolCol;
    sourceCol.input = "symbolCol";
  } else if (!is.null(refSeq)) {
    identifiers = refSeq;
    objName = "REFSEQ2EG";
    sourceCol = refSeqCol;
    sourceCol.input = "refSeqCol";
  } else stop("Either `symbol' or 'refSeq' must be given.");

  if (length(identifiers)==0) return(numeric(0));

  if (!is.null(annotation))
  {
    if (!sourceCol %in% names(annotation)) stop("'", sourceCol.input, "' must be among 'names(annotation)'.");
    if (!aliasCol %in% names(annotation)) stop("'aliasCol' must be among 'names(annotation)'.");
    if (!entrezCol %in% names(annotation)) stop("'entrezCol' must be among 'names(annotation)'.");
  }

  args = list(...);
  if ("s2e" %in% names(args))
  {
    s2e = args$s2e
  } else if (is.null(annotation))
    {
      s2e = orgObject(organism, objName, convertToList = TRUE, uniqueOnly = FALSE, warnNonUnique = FALSE);
    } else {
      keep = !is.na(annotation[[sourceCol]])
      s2e = setNames(annotation[[entrezCol]][keep], annotation[[sourceCol]] [keep]);
    }

  if ("s2e2" %in% names(args))
  {
    s2e2 = args$s2e2;
  } else {
    if (!is.null(symbol) && useAlias)
    {
      if (is.null(annotation))
      { 
        s2e2 = orgObject(organism, "ALIAS2EG", convertToList = TRUE, uniqueOnly = FALSE, warnNonUnique = FALSE)
      } else {
        split = strsplit(annotation[[aliasCol]], split = aliasSep, fixed = aliasSplit.fixed)
        aliasIndex = .indexedFlattenedList(split);
        aliasIndex = aliasIndex[ !is.na(aliasIndex$data), ];
        s2e2 = tapply(annotation[[entrezCol]] [ aliasIndex$index], aliasIndex$data, identity);
      }
      ref2 = names(s2e2);
      if (ignoreCase) ref2 = toupper(ref2);
    } else s2e2 = NULL;
  }

  if (!is.atomic(identifiers)) { 
    if (!is.null(symbol))
    {
      out = lapply(identifiers, function(id) convert2entrez(organism = organism,
               symbol = id, refSeq = NULL,
               ignoreCase = ignoreCase, dropNonConverted = dropNonConverted,
               useAlias = useAlias, noMatch = noMatch, multipleMatches = multipleMatches,
               annotation = annotation, symbolCol = symbolCol, refSeqCol = refSeqCol,
               aliasCol = aliasCol, aliasSep = aliasSep, aliasSplit.fixed = aliasSplit.fixed, 
               entrezCol = entrezCol,
               s2e = s2e, s2e2 = s2e2))
    } else 
      out = lapply(identifiers, function(id) convert2entrez(organism = organism,
               symbol = NULL, refSeq = id,
               ignoreCase = ignoreCase, dropNonConverted = dropNonConverted,
               useAlias = useAlias, noMatch = noMatch, multipleMatches = multipleMatches,
               annotation = annotation, symbolCol = symbolCol, refSeqCol = refSeqCol,
               aliasCol = aliasCol, aliasSep = aliasSep, aliasSplit.fixed = aliasSplit.fixed,
               entrezCol = entrezCol,
               s2e = s2e))
    return(out);
  } 

  ref = names(s2e);
  if (ignoreCase) 
  {
    identifiers = toupper(identifiers);
    ref = toupper(ref);
  }
    
  if (!is.null(s2e2))
  {
    ref2 = names(s2e2);
    if (ignoreCase) ref2 = toupper(ref2);
  }

  ident2ref = match(identifiers, ref);
  lens = sapply(s2e, length);
  unique = .replaceMissing(lens[ident2ref]==1);
  multiple = .replaceMissing(lens[ident2ref]>1);
  out = rep(.noMatch1, length(identifiers));
  out[multiple] = .multiMatch1;
  out[unique] = unlist(s2e[ident2ref[unique]]); 

  nonConverted = !identifiers %in% ref;
  if (!is.null(symbol) && any(nonConverted) && useAlias)
  {
    id2ref2 = match(identifiers[nonConverted], ref2);
    lens2 = sapply(s2e2, length);
    unique = .replaceMissing(lens2[id2ref2]==1);
    multiple = .replaceMissing(lens2[id2ref2]>1);
    out[nonConverted][unique] = unlist(s2e2[id2ref2[unique]])
    out[nonConverted][multiple] = .multiMatch1;
  }
  if (dropNonConverted) {
    out = out[out!=.noMatch1]; 
  } else 
    out[out==.noMatch1] = noMatch;

  out[.replaceMissing(out==.multiMatch1)] = multipleMatches;
  out
}

entrezFromMultipleSymbols = function(symbol, split = ";", organism = "human", ignoreCase = FALSE)
{
  symbol.split = lapply(strsplit(symbol, split = split), unique)
  symbols.indexed = .indexedFlattenedList(symbol.split);

  entrez.indexed = convert2entrez(organism = organism, symbol = symbols.indexed$data,
                                    useAlias = FALSE, noMatch = NA, multipleMatches = -1,
                                    ignoreCase = ignoreCase);
  entrez.indexed.alias = convert2entrez(organism = organism, symbol = symbols.indexed$data,
                                    useAlias = TRUE, noMatch = NA, multipleMatches = -1,
                                    ignoreCase = ignoreCase);

  .findEntrez = function(entrez1)
  {
    # Remove all missing 
    finEntrez = !is.na(entrez1);
    entrez1 = entrez1[finEntrez];
    symbols1 = symbols.indexed[finEntrez, ]
    var = tapply(entrez1, symbols1$index, var, na.rm = TRUE)
    var[is.na(var)] = 0;
    keep = as.numeric(names(var)[var==0]);
    entrez1[entrez1==-1] = NA;
    entrez.all = rep(NA, length(symbol));
    keep2index = match(keep, symbols1$index)
    entrez.all[keep] = entrez1[keep2index];
    entrez.all;
  }

  entrez = .findEntrez(entrez.indexed);
  entrez.alias = .findEntrez(entrez.indexed.alias);
  entrez[is.na(entrez)] = entrez.alias[is.na(entrez)];
  entrez;
}


#======================================================================================================
#
# Default backgrounds
#
#======================================================================================================

allOrgGenes = function(organism)
{
  obj = orgObject(organism = organism, objectName = "CHR", convertToList = TRUE, uniqueOnly = FALSE)
  invisible(as.numeric(names(obj)));
}

#======================================================================================================
#
# Translation of Entrez IDs between organisms
#
#======================================================================================================

.translateUsingTable = function(x, translationTable, keepUntranslated = FALSE)
{
  if (is.atomic(x))
  {
    out = translationTable[ match(x, translationTable[, 1]), 2]
    if (keepUntranslated) out[is.na(out)] = x[is.na(out)];
    out;
  } else lapply(x, .translateUsingTable, translationTable, keepUntranslated = keepUntranslated);
}


entrezHomologs = function(orgFrom, orgTo, useHomology = TRUE, version = NULL)
{
  if (is.null(version))
  {
    homologyFile = "Homology-allOrganisms-2019.rda"
  } else if (version < 2019.03 & version >= 2014.04) 
  {
    homologyFile = "Homology-allOrganisms.rda"
  } else if (version < 2014.04) {
    homologyFile = "Homology-Hs-Mm.rda"
  } else homologyFile = "Homology-allOrganisms-2019.rda"
  shortFrom = organismShorthand(orgFrom);
  shortTo = organismShorthand(orgTo);
  homology = NULL;
  x = load(system.file("extdata", homologyFile, package = "anRichment"));

  if (! ("homology"%in% x)) stop("Internal error: file ", homologyFile,
                                 " does not contain the required variable 'homology'.");
  columnOrg = gsub("Entrez.", "", names(homology));
  colFrom = match(shortFrom, columnOrg);
  colTo = match(shortTo, columnOrg);
  if (is.na(colFrom) || is.na(colTo))
  {
    useHomology = FALSE
    warning(WGCNA::spaste("Homology information is not available for translating between ", orgFrom, "\n",
                  "   and ", orgTo, ". Will match genes by name. "));
  }

  if (useHomology)
  {
    map = homology[, c(colFrom, colTo)]
    keep = rowSums(is.na(map))==0;
    map = map[keep, ];
  } else {
    # Convert entrez identifiers to gene symbols
    egSymbol.from = toupper(orgObject(shortFrom, "SYMBOL", convertToList = TRUE, uniqueOnly = TRUE,
                            warnNonUnique = FALSE));
    egSymbol.to = toupper(orgObject(shortTo, "SYMBOL", convertToList = TRUE, uniqueOnly = TRUE,
                            warnNonUnique = FALSE));
    common = intersect(egSymbol.from, egSymbol.to);
    common2from = match(common, egSymbol.from);
    common2to = match(common, egSymbol.to);
    map = cbind(as.numeric(names(egSymbol.from)[common2from]), 
                as.numeric(names(egSymbol.to)[common2to]));
  }
  rownames(map) = NULL;
  colnames(map) = WGCNA::spaste("Entrez.", c(shortFrom, shortTo));
  map;
}


# If entrezMap is supplied, it is assumed that first column corresponds to orgFrom and second to orgTo.

mapEntrez = function(entrez, orgFrom, orgTo, entrezMap = NULL, useHomology = TRUE, version = NULL)
{
  if (is.null(entrezMap))
    entrezMap = entrezHomologs(orgFrom, orgTo, useHomology = useHomology, version = version); 

  out = .translateUsingTable(entrez, entrezMap);
  if (!any(is.finite(out)))
  {
     #browser()
     warning("mapEntrez: None of the genes in the existing gene lists could be mapped.")
  }

  out;
}

#======================================================================================================
#
# Convert a gene sets and collections to another organism
#
#======================================================================================================

# Convert a gene set

.convertDataSetToOrganism = function(dataSet, organism, addOldOrganismToName = FALSE, 
                                    namePattern = ".convertedFrom.%o",
                                    addOldOrganismToDescription = FALSE, 
                                    descriptionPattern = " (Converted from %o.)")
{
  if (sameOrganism(dataSet$organism[1], organism)) return(dataSet);

  if (addOldOrganismToName)
    dataSet$name = WGCNA::spaste(dataSet$name, gsub("%o", dataSet$organism[1], namePattern, fixed = TRUE));

  if (addOldOrganismToDescription)
    dataSet$description = WGCNA::spaste(dataSet$description, 
                               gsub("%o", dataSet$organism[1], descriptionPattern, fixed = TRUE));

  dataSet$organism = organismLabels(organism);
  dataSet;
}

convertGeneSetToOrganism = function(geneSet, organism, 
                                    entrezMap = NULL,
                                    useHomology = TRUE,
                                    addOldOrganismToName = FALSE, 
                                    namePattern = ".convertedFrom.%o",
                                    addOldOrganismToDescription = FALSE, 
                                    descriptionPattern = " (Converted from %o.)")
{

  if (sameOrganism(geneSet$organism[1], organism)) return(geneSet);

  entrez = geneSet$data$Entrez;
  newEntrez = mapEntrez(entrez, orgFrom = geneSet$organism[1], orgTo = organism,
                        entrezMap = entrezMap, useHomology = useHomology);
  keep = !is.na(newEntrez);
  geneSet$data$Entrez = newEntrez;
  geneSet$data = geneSet$data[keep, , drop = FALSE];

  .convertDataSetToOrganism(geneSet, organism, 
                           addOldOrganismToName = addOldOrganismToName,
                           namePattern = namePattern,
                           addOldOrganismToDescription = addOldOrganismToDescription,
                           descriptionPattern = descriptionPattern);
                           
}
  

# Convert gene identifiers in gene sets and gene properties from their original organism (whatever that may
# be) to the new given organism.

convertCollectionToOrganism = function(collection, organism, useHomology = TRUE, 
          addOldOrganismToSetNames = FALSE, namePattern = ".convertedFrom.%o",
          addOldOrganismToSetDescriptions = FALSE, 
          descriptionPattern = " (Converted from %o.)",
          dropEmptySets = TRUE)
{
  isGeneSet = sapply(collection$dataSets, .isGeneSet);
  isProperty = sapply(collection$dataSets, .isGeneProperty);

  if (any(!isGeneSet & !isProperty))
    stop("Internal error: some data sets in the collection are neither gene sets not gene properties.");

  # First gene sets. Each gene set can in principle correspond to its own organism.

  orgShort = organismShorthand(organism);

  if (any(isGeneSet))
  {
    # Get the shorthand of set organisms
    setOrgs = sapply(collection$dataSets[isGeneSet], function(x) x$organism[3]); 
    orgLevels = unique(setOrgs)
    # Only convert sets whose organism is different from 'organism'
    orgLevels = orgLevels[orgLevels!=orgShort];
    if (length(orgLevels) > 0)
    {
      entrezMaps = lapply(orgLevels, entrezHomologs, orgTo = organism, useHomology = useHomology);
      for (o in 1:length(orgLevels))
      {
         sets1 = which(isGeneSet)[setOrgs==orgLevels[o] ];
         collection$dataSets[sets1] = lapply(collection$dataSets[sets1], convertGeneSetToOrganism,
                                         organism = organism, entrezMap = entrezMaps[[o]], 
                                         addOldOrganismToName = addOldOrganismToSetNames, 
                                         namePattern = namePattern,
                                         addOldOrganismToDescription = addOldOrganismToSetDescriptions,
                                         descriptionPattern = descriptionPattern
                                         );
      }
    }
  }
  if (dropEmptySets)
  {
    lens = sapply(collection$dataSets[isGeneSet], function(ds) nrow(ds$data));
    keep = rep(TRUE, nDataSets(collection));
    keep[isGeneSet] [lens==0] = FALSE;
    collection$dataSets = collection$dataSets[keep];
  }
    
  if (any(isProperty))
  {
    # Gene properties all have the same single set of identifiers associated with them, so they must all
    # correspond to a single organism. Take it from the first gene property.
    oldOrg = collection$dataSets[[which(isProperty)[1]]]$organism[1];
    newIdentifiers = mapEntrez(collection$identifiers, orgFrom = oldOrg, orgTo = organism, 
                                  useHomology = useHomology);
    keep = !is.na(newIdentifiers);
    if (!any(keep)) 
       stop("None of the identifiers could be mapped. Please check that the original organism is correct.");
    collection$identifiers = collection$newIdentifiers[keep];
    collection$weights = collection$weights[keep, , drop = FALSE];
    collection$dataSets[isProperty] = lapply(collection$dataSets[isProperty], 
              function(gp, index, org)
              {
                gp$data = gp$data[index];
                #gp$organism = org;
                .convertDataSetToOrganism(gp, organism = org, 
                                       addOldOrganismToName = addOldOrganismToSetNames,
                                       namePattern = namePattern,
                                       addOldOrganismToDescription = addOldOrganismToSetDescriptions,
                                       descriptionPattern = descriptionPattern);
              }, keep, organism);
  }

  collection;
}

#===========================================================================================================
#
# Utility functions
#
#===========================================================================================================

.unique.nm = function(x)
{
  unique(x[!is.na(x)])
}

.translate = function (data, dictionary) 
{
  translated = dictionary[match(data, dictionary[, 1]), 2]
  attributes(translated) = attributes(data)
  translated
}


#=====================================================================================================
#
# enrichmentLabels
#
#=====================================================================================================

.matchCol = function(col, colNames)
{
  name = as.character(match.call(expand.dots = FALSE)$col)
  if (!is.numeric(col))
    col = match(col, colNames);
  if (is.na(col)) stop("Undefined or mis-specified '", name, "'.")
  if (col<1 | col>length(colNames)) 
    stop("'", name, "' is out of range: must be between 1 and the number of columns in input.");
  col
}


enrichmentLabels = function(
  # Enrichment table and the relevant columns
  enrichmentTable, 
  classCol = "class", 
  pValueCol = "Bonferroni", 
  setNameCol = "dataSetName", 
  groupCol = "inGroups", 
  fractionCol = "fracOfEffectiveClassSize", 
  classSizeCol = "classSize", 

  classLabels = NULL,

  # Restrictions on terms to be included
  pValueThreshold = 0.05, 
  minSize = 0.1, 

  # Which groups to include as separate entries
  focusOnGroups = c("all", "GO.BP", "GO.MF", "GO.CC", "Brain region markers", "Cell type markers"), 
  excludeGroups = NULL,
  groupShortNames = NULL, 
  numericClassLabels = TRUE, 

  entrySeparator = "|",

  # Optional: collection or at least group information (component of collection) if implied groups are to be used

  collection = NULL,
  groupList = NULL,

  # What to include in the one-line summary
  includeClassSize = TRUE, 
  outputClassPrefix = if (includeClassSize) "M." else "", 
  orderBySignificance = TRUE,
  includeFraction = TRUE, 
  fractionDigits = 2)
{

  colNames = colnames(enrichmentTable);
  classCol = .matchCol(classCol, colNames);
  pValueCol = .matchCol(pValueCol, colNames);
  setNameCol = .matchCol(setNameCol, colNames);
  groupCol = .matchCol(groupCol, colNames); 
  fractionCol = .matchCol(fractionCol, colNames); 
  classSizeCol = .matchCol(classSizeCol, colNames); 

  nFocus = length(focusOnGroups)

  if (nFocus ==0) stop("A group (or groups) to focus on must be given in 'focusOnGroups'.")

  nRow = nrow(enrichmentTable);

  if (nRow==0) return(NULL);

  groupSplit = strsplit(enrichmentTable[, groupCol], split = entrySeparator, 
                        fixed = TRUE);

  if (!is.null(collection)) groupList = collection$groups;
  if (!is.null(groupList))
  {
    uniqueGroups = unique(unlist(groupSplit));
    existing = unique(c(sapply(groupList, getElement, "name"), unlist(lapply(groupList, getElement, "alternateNames"))));
    uniqueGroups = intersect(uniqueGroups, existing);
    implied = impliedGroups(groupList, queryGroups = uniqueGroups, get = "parents");
    groupSplit = lapply(groupSplit, function(gs)
    {
      gs2 = intersect(gs, names(implied));
      implied1 = implied[gs2];
      unique(c(gs, unlist(implied1)));
    });
  }
  lengths = sapply(groupSplit, length);
  groupIndex = data.frame(group = unlist(groupSplit), row = rep(c(1:nRow), lengths));
  dropRows = unique(groupIndex$row[ groupIndex$group %in% excludeGroups]);
  groupIndex = groupIndex[ !(groupIndex$row %in% dropRows), ]

  allClasses = sort(unique(enrichmentTable[, classCol]));
  if (is.null(classLabels))
  {
    classSizes = enrichmentTable[ match(allClasses, enrichmentTable[, classCol]), classSizeCol];
  } else {
    tab = table(as.vector(classLabels));
    classSizes = as.numeric(tab[match(allClasses, names(tab))])
  }

  nClasses = length(allClasses);
  bestLabels = matrix("", nClasses, nFocus);
  bestP = matrix(NA, nClasses, nFocus);
  fraction = matrix(NA, nClasses, nFocus);
  colnames(bestLabels) = colnames(bestP) = colnames(fraction) = focusOnGroups;
  for (f in 1:nFocus)
  {
    if (focusOnGroups[f]=="all") 
    {
      rowIndex = groupIndex$row;
    } else {
      rowIndex = groupIndex$row [ groupIndex$group == focusOnGroups[f] ];
    }
      
    table = enrichmentTable[rowIndex, ];
    p = as.numeric(as.character(table[, pValueCol]));

    # If minSize was specified, prefer sets that are larger than threshold and pass the p-value
    # threshold.
    if (is.null(minSize)) 
    {
      preference = rep(1, nrow(table))
    } else {
      preference = -as.numeric((table[, fractionCol] > minSize) & ( p<=pValueThreshold) )
    }

    order = order(table[, classCol], preference, p )
    enr.ord = table[order, ];
    class = enr.ord[, classCol];

    n = length(class);
    bestIndex = c(1, c(2:n)[ class[-1]!=class[-n] ]);
    classes1 = class[bestIndex];
    this2all = match(classes1, allClasses);
    bestP[ this2all, f]  = as.numeric(as.character(enr.ord[bestIndex, pValueCol]))
    bestLabels[ this2all, f] = enr.ord[bestIndex, setNameCol];
    fraction[this2all, f] = enr.ord[bestIndex, fractionCol];
  }
  bestInfo = .interleave(lapply(list(bestLabels, signif(bestP, 2), round(fraction, 2)), as.data.frame),
                  nameBase = rep("", 3));
  colnames(bestInfo) = WGCNA::spaste( rep(c("highestEnrichedSet.", "pValue.", "fraction."), nFocus),
                               rep(focusOnGroups, rep(3, nFocus)));

  #numCols = c(grep("fraction", colnames(bestInfo)), grep("pValue", colnames(bestInfo)));
  #bestInfo[, numCols] = apply(bestInfo[, numCols, drop = FALSE], 2, 
  #                            function(x) {signif(as.numeric(as.character(x)), 2)})

  bestP[is.na(bestP)] = 2;
  fraction[is.na(fraction)] = -1;

  #bestLabels.keep = WGCNA::spaste(bestLabels, "; ");
  bestLabels.keep = bestLabels;
  if (includeFraction)
    bestLabels.keep = WGCNA::spaste(bestLabels, " (", round(fraction, fractionDigits), ")");

  dim(bestLabels.keep) = dim(bestLabels);
  bestLabels.keep[ bestP > pValueThreshold | fraction < minSize ] = "";

  separators = ifelse( rowSums(bestLabels.keep!="") > 0, ": ", "");
  enrlabels0 = WGCNA::spaste(separators, sapply(1:nClasses, function(cl)
                      {
                        bl = bestLabels.keep[cl, ];
                        p1 = bestP[cl, ];
                        dup = duplicated(bl) | bl=="";
                        if (!is.null(groupShortNames))
                          bl = WGCNA::spaste(bl, " (", groupShortNames, ")")
                        if (orderBySignificance)
                        {
                           order1 = order(p1[!dup]);
                        } else 
                           order1 = 1:sum(!dup);
                        base::paste(bl[!dup][order1], collapse = "; ")
                      } ));
  sizeString = if (includeClassSize) WGCNA::spaste(" (", classSizes, ")") else rep("", nClasses);
  enrlabels = sub("; $", "", WGCNA::spaste(outputClassPrefix, allClasses, sizeString, enrlabels0));
  enrlabels = gsub(") (", "; ", enrlabels, fixed = TRUE);

  out = data.frame(enrichmentLabel = enrlabels, 
                   class = allClasses, classSize = classSizes, bestInfo, bestLabel = bestLabels);
  if (numericClassLabels)
    out = out[order(as.numeric(as.character(allClasses))), ]

  out; 
}

.interleave = function(matrices, nameBase = names(matrices), sep = ".", baseFirst = TRUE,
                      doInterleave = TRUE, check.names = TRUE)
{
  nMats = length(matrices)
  nCols = ncol(matrices[[1]]);

  dims = lapply(matrices, dim);

  if (any(nameBase!=""))
  {
    if (baseFirst)
    {
       for (m in 1:nMats) colnames(matrices[[m]]) = spaste(nameBase[m], sep, colnames(matrices[[m]]));
    } else {
       for (m in 1:nMats) colnames(matrices[[m]]) = spaste(colnames(matrices[[m]]), sep, nameBase[m]);
    }
  }
  matrices = lapply(matrices, function(mat) if (is.list(mat)) as.data.frame(mat) else mat)

  if (doInterleave)
  {
    out = as.data.frame(lapply(1:nCols,
                               function(index, matrices)
                                  as.data.frame(lapply(matrices,
                                            function(x, i) x[, i, drop = FALSE], index), check.names = FALSE),
                               matrices), check.names = check.names);
  } else
    out = as.data.frame(do.call(cbind, .removeListNames(matrices)));

  if (!is.null(rownames(matrices[[1]]))) rownames(out) = make.unique(rownames(matrices[[1]]))
  out;
}



#====================================================================================================
#
# Gene set clustering functions
#
#====================================================================================================

dataSetMatrix = function(collection, background = NULL,
                          includeProperties = TRUE,
                          includeGeneSets = TRUE,
                          dropZeroOverlapSets = FALSE,
                          namesFrom = c("ID", "name"))
{
  if (!.isCollection(collection))
    stop("Input 'collection' must be a valid collection object.");

  namesFrom = match.arg(namesFrom);
  setNames = as.character(sapply(collection$dataSets, getElement, namesFrom));
  if (is.null(background)) background = allDataSetGenes(collection);
  nGenes = length(background);
  nSets = nDataSets(collection);

  isNumProp = sapply(collection$dataSets, inherits, "PL-numericProperty");
  isDiscProp = sapply(collection$dataSets, inherits, "PL-discreteProperty");
  isGeneSet = sapply(collection$dataSets, inherits, "PL-geneSet");
  includeSets = (isGeneSet & includeGeneSets) | ( (isNumProp | isDiscProp) & includeProperties);
  includeIndex = which(includeSets);
  nSets = length(includeIndex);

  out = matrix(0, nGenes, nSets);

  id2background = match(collection$identifiers, background);
  fin.identifiers = is.finite(id2background);
  id2background.fin = id2background[fin.identifiers];
  any.id2background.fin = any(id2background.fin);
  keep = rep(FALSE, nSets);
  for (set in 1:nSets)
  {
    s = includeIndex[set];
    if (isGeneSet[s])
    {
       set2background = match(collection$dataSets[[s]]$data$Entrez, background);
       fin = is.finite(set2background);
       if (any(fin)) 
       {
          out[ set2background[fin], set] = 1;
          keep[set] = TRUE;
       }
    } else if (any.id2background.fin) {
       keep[set] = TRUE;
       out[-id2background.fin, set] = NA;
       out[ id2background.fin, set] = collection$dataSets[[s]]$data[fin.identifiers];
    }
  }
  colnames(out) = setNames;
  rownames(out) = background;
  if (dropZeroOverlapSets) out = out[, keep];     
  out;
}

.dataSetSimilarity = function(dsMat, 
                              type = c("product", "cosine", "cov", "cor"),
                              productNormalization = c("min", "mean", "max"),
                              TOMType = c("none", "min", "mean"),
                              verbose = 0, indent = 0
                              )
{
  type = match.arg(type);
  productNormalization = match.arg(productNormalization)
  TOMType = match.arg(TOMType);

  spaces = indentSpaces(indent);

  if (TOMType!="none" && type=="cov")
    stop("TOM cannot be used in conjunction with covariance (cov) 'type'.");

  if (type=="cor" || type=="cosine")
  { 
    sim = cor(dsMat, use = 'p', cosine = type=="cosine");
  } else if (type=="cov") {
    sim = cov(dsMat, use = 'p');
  } else {
    fin.dsMat = !is.na(dsMat);
    dsMat[!fin.dsMat] = 0;
    counts1 = t(dsMat^2 * fin.dsMat) %*% fin.dsMat;
    productNormFnc = match.fun(WGCNA::spaste("p", productNormalization));
    counts = do.call(productNormFnc, list(counts1, t(counts1))); 
    products = t(dsMat) %*% dsMat;
    sim = products/counts;
    sim[!is.finite(sim)] = 0;
  }

  if (TOMType!="none")
    sim = TOMsimilarity(sim, TOMType = "unsigned", TOMDenom = TOMType, verbose = verbose, indent = indent)

  sim;
}

dataSetSimilarity = function(
    collection,

    # background against which to calculate obverlaps. defaults to all genes in the collection.
    background = NULL,

    # Subsetting optiona
    tags = NULL, 
    matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                      "groupAlternateNames", "nameAndAlternates", "groupsAndAlternates"),

    searchType = c("any", "all"), invertSearch = FALSE,
    exactMatch = TRUE, fixed = TRUE, ignore.case = TRUE,

    # Special subsetting options
    includeProperties = TRUE,
    includeGeneSets = TRUE,
    dropZeroOverlapSets = FALSE,

    namesFrom = c("ID", "name"),

    # Alternative to all of the above: pre-calculated data set matrix
    dataSetMat = NULL,

    # Distance options
    type = c("product", "cosine", "cov", "cor"),
    productNormalization = c("min", "mean", "max"),
    TOMType = c("none", "min", "mean"),
    verbose = 1,
    indent = 0
  )
{

  spaces = indentSpaces(indent);

  if (is.null(dataSetMat))
  {
    if (!.isCollection(collection)) 
      stop("Input 'collection' must be a valid collection object.");
  
  
    collection = subsetCollection(collection, tags = tags, 
                                  matchComponents = matchComponents, searchType = searchType,
                                  invertSearch = invertSearch, exactMatch = exactMatch, fixed = fixed,
                                  ignore.case = ignore.case);

    if (is.null(background)) background = allDataSetGenes(collection);
    if (verbose > 0) 
      printFlush(WGCNA::spaste(spaces, "dataSetSimilarity: creating membership matrix.."));

    dataSetMat = dataSetMatrix(collection, background = background,
                            includeProperties = includeProperties,
                            includeGeneSets = includeGeneSets,
                            dropZeroOverlapSets = dropZeroOverlapSets,
                            namesFrom = namesFrom);
  }
  
  if (verbose > 0) 
    printFlush(WGCNA::spaste(spaces, "dataSetSimilarity: calculating similarity matrix.."));

  .dataSetSimilarity(dataSetMat, type = type, productNormalization = productNormalization,
                     TOMType = TOMType);
}


#======================================================================================================
#
# Supporting functions for matching subcollections to main collection
#
#======================================================================================================

subCollectionIndex = function(subCollection, combinedCollection, combinedIDNames = NULL,
       stopOnUnmatched = TRUE)
{
  if (is.null(combinedIDNames))
     combinedIDNames = spaste(dataSetIDs(combinedCollection), dataSetNames(combinedCollection));

  subIDNames = spaste(dataSetIDs(subCollection), dataSetNames(subCollection));

  out = match(subIDNames, combinedIDNames);
  if (any(is.na(out)))
    (if (stopOnUnmatched) stop else warning)(spaste(
       "Some ID-name combinations in 'subCollection' could not be matched\n",
       "to ID-name combinations in 'combinedCollection'."));
  out;
}

subCollectionRowsInEnrichmentTable = function(enrichmentTable, subCollection)
{  
  if (!is.null(enrichmentTable$dataSetID))
  {
    subIDNames = spaste(dataSetIDs(subCollection), dataSetNames(subCollection))
    etIDNames = spaste(enrichmentTable$dataSetID, enrichmentTable$dataSetName);
  } else {
    subIDNames = dataSetNames(subCollection)
    etIDNames = enrichmentTable$dataSetName;
  }
  which(etIDNames %in% subIDNames);
}

splitEnrichmentResultsBySubcollections = function(enr, combinedCollection, collections,
           dropColumns = character(0))
{
   if (length(enr)==0) return(NULL);
   tables1 = enr[c("countsInDataSet", "pValues")];
   tablesByCollection = lapply(collections, function(coll)
   {
      if (length(enr$enrichmentTable) > 0)
      {
         rows = subCollectionRowsInEnrichmentTable(enr$enrichmentTable, coll);
         etab = enr$enrichmentTable[rows, ];
         # Rank the enrichment rank
         rankCol = which(names(etab)=="rank");
         if (length(rankCol)!=1) {
           printFlush("splitEnrichmentResultsBySubcollections: Internal error. Dropping into browser.")
           browser();
         }
         classLevels = unique(etab$class);
         rank2 = etab$rank;
         for (cl in classLevels)
           rank2[etab$class==cl] = rank(etab$rank[etab$class==cl]);
         etab2 = cbind( etab[, 1:(rankCol-1), drop = FALSE],
                        rank = rank2,
                        rankAcrossAllCollections = etab$rank,
                        etab[, -c(1:rankCol)]);
         etab2 = etab2[, !colnames(etab2) %in% dropColumns];
      } else
         etab2 = NULL;
      index = subCollectionIndex(coll, combinedCollection);
      if (!is.null(tables1[[1]])) 
      {
        index2 = match(dataSetIDs(coll), rownames(tables1[[1]]))
        if(!all(.replaceMissing(index==index2)))
        {
           printFlush("splitEnrichmentResultsBySubcollections: found inconsistency.\n Possible",
                      "cause: duplicated IDs of collection gene or data sets.");
           browser();
        }
      }
      out = try(c( list(enrichmentTable = etab2),
         lapply(tables1, function(t) if (is.null(t)) NULL else t[index, ])));
      if (inherits(out, "try-error"))
      {
         printFlush("splitEnrichmentResultsBySubcollections: found error, dropping into browser.");
         browser();
      }
      out;
   });
   tablesByCollection;
}

splitEnrichmentTableBySubcollections = function(enrTab, combinedCollection, collections,
           dropColumns = character(0))
{
   if (length(enrTab)==0) return(NULL);
   if (is.null(nrow(enrTab)) || nrow(enrTab)==0) return(enrTab);
   tablesByCollection = lapply(collections, function(coll)
   {
      rows = subCollectionRowsInEnrichmentTable(enrTab, coll);
      etab = enrTab[rows, ];
      # Rank the enrichment rank
      rankCol = which(names(etab)=="rank");
      classLevels = unique(etab$class);
      rank2 = etab$rank;
      for (cl in classLevels) 
        rank2[etab$class==cl] = rank(etab$rank[etab$class==cl]);
      etab2 = cbind( etab[, 1:(rankCol-1), drop = FALSE],
                     rank = rank2,
                     rankAcrossAllCollections = etab$rank,
                     etab[, -c(1:rankCol)]);
      etab2[, !colnames(etab2) %in% dropColumns];
   });
   tablesByCollection;
}


#=========================================================================================================
#
# Functions specific to enrichment analysis for data where multiple probes (variables) map to one gene,
# especially for methylation data
#
#=========================================================================================================

.removeListNames = function(lst) 
{
  names(lst) = NULL;
  lst;
}

.indexedFlattenedList = function(lst)
{  
  n = length(lst);
  lengths = sapply(lst, length);
  index = do.call("c", .mymapply(rep, 1:n, lengths));
  data.frame(index = index, data = unlist(.removeListNames(lst)));
}

.multiGrepl = function(patterns, x, ...)
{
  mat = do.call(cbind, lapply(patterns, function(p) as.numeric(grepl(p, x, ...))));
  rowSums(mat)>0;
}

.selectRows = function(data,
   keepValues = NULL,        
   keepValueColumns = NULL,  
   keepPatterns = NULL,      
   keepPatternColumns = NULL,
   ignoreCase = TRUE)
{

  if (ignoreCase) trafo = toupper else trafo = identity;

  keepRows = rep(TRUE, nrow(data));
  if (!is.null(keepValues))
  {
    if (is.atomic(keepValues))
    {
       if (length(keepValueColumns) !=1) 
         stop("'keepValues' should be a list with one component per entry in 'keepValueColumns'.");
       keepValues = list(keepValues);
    }
    if (length(keepValues)!=length(keepValueColumns))
       stop("Length of the 'keepValues' list must equal length of 'keepValueColumns'.")
    for (rc in 1:length(keepValueColumns))
    {
      col = match(keepValueColumns[rc], colnames(data))
      if (is.na(col)) 
        stop("Column ", keepValueColumns[rc], " not found among the columns of 'data':\n", 
             base::paste(colnames(data), collapse = ", "));
      keepRows = .replaceMissing(keepRows & (trafo(data[, col]) %in% trafo(keepValues[[rc]])));
      if (sum(keepRows, na.rm = TRUE)==0)
        stop("No rows were left after restriction applied to ",  keepValueColumns[rc]);
    }
  }

  if (!is.null(keepPatterns))
  {
    #browser();
    if (is.atomic(keepPatterns))
    { 
       if (length(keepPatternColumns) !=1)
         stop("'keepPatterns' should be a list with one component per entry in 'keepPatternColumns'.");
       keepPatterns = list(keepPatterns);
    }
    if (length(keepPatterns)!=length(keepPatternColumns))
       stop("Length of the 'keepPatterns' list must equal length of 'keepPatternColumns'.")
    for (rc in 1:length(keepPatternColumns))
    {
      col = match(keepPatternColumns[rc], colnames(data))
      if (is.na(col))
        stop("Column ", keepPatternColumns[rc], " not found among the columns of 'data':\n",
             base::paste(colnames(data), collapse = ", "));
      keepRows = .replaceMissing(keepRows & .multiGrepl(keepPatterns[[rc]], data[, col], 
                                                      ignore.case = ignoreCase));
      if (sum(keepRows, na.rm = TRUE)==0)
        stop("No rows were left after restriction applied to ",  keepPatternColumns[rc]);
    }
    printFlush(WGCNA::spaste("Number of retained probes: ", sum(keepRows)));
  }

  keepRows;
}

# Pick one representative probe for each gene from data.

representativeProbes = function(
   mtmat, 

   # Mapping of variables to genes: either Entrez...
   entrez = NULL,

   #... or an annotation table with other necessary information
   annot = NULL,
   annotIDCol = "Name",
   entrezCol = NULL,
   symbolCol = "UCSCRefGeneName",
   symbolSep = ";",
   organism = "human",

   # Optional restrictions on CpG probes
   keepValues = NULL,
   keepValueColumns = NULL,
   keepPatterns = NULL, 
   keepPatternColumns = NULL,
   ignoreCase = TRUE,

   # Options governing selection of representatives
   selectionMethod = "maxVariance",
   selectionStatisticFnc = NULL,
   useGroupHubs = TRUE,
   connectivityPower = 6, 
   minProportionPresent = 1,
   statisticFncArguments = list(),
   adjacencyArguments = list(),

   # Consensus options 
   consensusQuantile = 0,
   calibration = "full quantile",

   getRepresentativeData = TRUE,
   dataToReduce = NULL,
   verbose = 2, indent = 0)
{
  if (isMultiData(mtmat, strict = FALSE))
  {
     if (!isMultiData(mtmat, strict = TRUE))
        stop("If 'mtmat' is a multiData structure, it must have the same columns in each set.");
     multi = TRUE;
  } else {
     multi = FALSE;
     mtmat = multiData(mtmat);
  }
  nVars = checkSets(mtmat)$nGenes;

  if (!is.null(annot))
  {
    nameCol = .matchNames(annotIDCol, annot, ignoreCase = FALSE, stopOnMissing = TRUE);
    annot = annot[ match(mtd.colnames(mtmat), annot[, nameCol]), ]

    selectRows = .selectRows(annot, keepValues = keepValues, keepValueColumns = keepValueColumns,
                                keepPatterns = keepPatterns, keepPatternColumns,
                                ignoreCase = ignoreCase);
  } else
    selectRows = rep(TRUE, nVars);

  if (is.null(entrez)) 
  {
     if (is.null(annot))
       stop("If 'entrez' is not given, 'annot' must be given.");

     eCol = .matchNames(entrezCol, annot, ignoreCase = FALSE, stopOnNULL = FALSE, stopOnMissing = TRUE);
     if (is.null(eCol))
     {
        entrez = entrezFromMultipleSymbols(annot[ match(mtd.colnames(mtmat), annot[, nameCol]),
                                                .matchNames(symbolCol, annot, ignoreCase = FALSE,
                                                            stopOnNULL = TRUE, stopOnMissing = TRUE)],
                                        split = symbolSep, organism = organism);
     } else
       entrez = annot[ match(mtd.colnames(mtmat), annot[, nameCol]), eCol];
  }

  if (length(entrez)!=nVars)
    stop("If given, the length of 'entrez' must equal the number of variables (columns) in 'mtmat'.");

  keepProbes = is.finite(entrez) & selectRows;

  mtmat.keep = mtd.subset(mtmat, , keepProbes);
      
  repData = consensusRepresentatives(mtmat.keep, 
                          group = entrez[keepProbes],
                          colID = mtd.colnames(mtmat.keep),
                          consensusQuantile = consensusQuantile,
                          calibration = calibration,
                          method = selectionMethod,
                          selectionStatisticFnc = selectionStatisticFnc,
                          statisticFncArguments = statisticFncArguments,
                          adjacencyArguments = adjacencyArguments,
                          useGroupHubs = useGroupHubs,
                          connectivityPower = connectivityPower,
                          minProportionPresent = minProportionPresent,
                          getRepresentativeData = getRepresentativeData,
                          verbose = verbose, indent = indent);

  if (!is.null(dataToReduce))
  {
    repIndex = match(repData$representatives[, 2], mtd.colnames(mtmat));
    if (is.null(dim(dataToReduce))) {
       reducedData = dataToReduce[repIndex];
       names(reducedData) = repData$representatives[, 1];
    } else {
       reducedData = dataToReduce[repIndex, ];
       colnames(reducedData) = repData$representatives[, 1];
    }
  } else 
    reducedData = NULL; 

  list(representatives = repData$representatives,
       varSelected = repData$varSelected,
       representativeData = repData$representativeData,
       reducedData = reducedData);
}

#==============================================================================================================
#
# Accessor function for standardized groups
#
#==============================================================================================================

standardizedGroups.df = function()
{
  x = load(system.file("extdata", "standardizedGroups-general.rda", package = "anRichment"))
  get(x);
}

standardizedGroups = function(groups.df = standardizedGroups.df())
{
  groupsFromDataFrame(groups.df) 
}
