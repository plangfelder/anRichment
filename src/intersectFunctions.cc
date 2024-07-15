/*
 *
 * Functions for fast calculation of intersects and intersect sizes of lists of vectors
 *
 * No fancy algorithms here, speedup is from collecting many calculations into one C++ call 
 * and from having pre-sorted the vectors.
 *
 * Code by Peter Langfelder, licensed under GPL-2
 *
 */
  

#include <Rcpp.h>
#include <vector>

#include "matchFunctions.cc"

using namespace Rcpp;
using namespace std;


// Note: this function is not exported because it returns C indexing (0-based).
List match_int(const List & ls1, const IntegerVector & v2, const IntegerVector & incomparables)
{
  size_t nSets = ls1.size();
  List out(nSets);
  HashData_Int hashTable(v2, NA_INTEGER, incomparables);

  for (size_t set = 0; set<nSets; set++)
  {
    IntegerVector inter2 = hashTable.HashLookup(ls1[set]);
    out[set] = inter2;
  }
  return out;
}

// [[Rcpp::export(name = ".match_int_C")]]
List match_int_C(const List & ls1, const IntegerVector & v2, const IntegerVector & incomparables)
{
  size_t nSets = ls1.size();
  List out(nSets);
  HashData_Int hashTable(v2, NA_INTEGER, incomparables);

  for (size_t set = 0; set<nSets; set++)
  {
    IntegerVector inter2 = hashTable.HashLookup(ls1[set]);
    for (size_t i=0; i<inter2.size(); i++) if (inter2[i]!=NA_INTEGER) inter2[i]++;
    out[set] = inter2;
  }
  return out;
}

// [[Rcpp::export(name = ".intersect_int")]]
List intersect_int(const List & ls1, const IntegerVector & v2, const IntegerVector & incomparables)
{
  size_t nSets = ls1.size();
  List out(nSets);
  HashData_Int hashTable(v2, NA_INTEGER, incomparables);

  for (size_t set = 0; set<nSets; set++)
  {
    IntegerVector inter2 = hashTable.intersect(ls1[set]);
    out[set] = inter2;
  }
  return out;
}

List intersect_int_list(const List & lst, const List & ref, const IntegerVector & lst2ref, 
                        const IntegerVector & incomparables)
{
  size_t nSets = lst.size(), nRef = ref.size();
  
  // Check inputs for consistency
  if (lst.size() != lst2ref.size())
    stop("intersect_int_list: arguments lst and lst2ref have inconsistent lengths.");

  // Prepare hash tables for each component of ref
  vector <HashData_Int> hashTable(nRef);

  for (size_t r=0; r < nRef; r++)
  {
     IntegerVector refV = ref[r];
     hashTable[r].init(refV, NA_INTEGER, incomparables);
  }

  // Do the matching
  List out(nSets);
  for (size_t set = 0; set<nSets; set++)
  {
    size_t refIndex = lst2ref[set];
    if (refIndex >= nRef) 
       stop("intersect_int_list: Invalid entry in lst2ref for set index " + to_string(set) + ": " + to_string(refIndex)
                           + "\n   note: maximum valid value is " + to_string(nRef));
    IntegerVector inter2 = hashTable[ lst2ref[set]].intersect(lst[set]);
    out[set] = inter2;
  }
  return out;
}

// [[Rcpp::export(name = ".intersectSize_int")]]
IntegerVector intersectSize_int(const List & ls1, const IntegerVector & v2, const IntegerVector & incomparables)
{
  size_t nSets = ls1.size();
  IntegerVector out(nSets);
  HashData_Int hashTable(v2, NA_INTEGER, incomparables);

  for (size_t set = 0; set<nSets; set++)
  {
    out[set] = hashTable.intersectSize(ls1[set]);
  }
  return out;
}

// Helper function: count the number of non-missing entries

size_t nNonmissing(IntegerVector v)
{
  size_t count = 0;
  for (size_t i=0; i<v.length(); i++) if (v[i]!=NA_INTEGER) count++;
  return count;
}

  

/*==================================================================================================
 *
 * The big enrichment utility function to calculate background, set and intersect sizes with respect 
 * to a given background. 
 *
 * All classes and sets need to be given as integers with respect to background which is assumed to 
 * be 1:backgroundSize.
 *
 * The identifiers in classes and in reference sets need not be the same; they will be mapped using 
 * a user-supplied mapping. The only requirement is that one class identifier maps to at most one 
 * reference identifier, i.e., the map is a function from class IDs to reference IDs (not necessarily 
 * 1-1).
 *
 * Inputs are: 
 *  . a list with class members for each class
 *  . a list with backgrounds
 *  . an integer vector with mapping from class to bakground
 *  . a list with reference set members
 *  . map between class and reference IDs: a two-column integer matrix
 *  . reference background: either mapped union of all ref sets or a custom set. Important: 
 *    the reference background has to refer to class members, not set members since the custom 
 *    background in the calling R function contains class IDs, not set IDs.
 *
 * =================================================================================================*/

// [[Rcpp::export(name=".intersectSizesForEnrichment")]]
List intersectSizesForEnrichment(
   const List & classes, 
   const List & classBg,
   const IntegerVector & class2bg,

   const List & refSets, 
   const LogicalVector & doMap,
   const IntegerMatrix & class2ref,
   const LogicalVector & class2setIDHasDuplicates,

   const CharacterVector & bgTypeR,
   const IntegerVector & mappedRefBg,

   const LogicalVector & returnMappedSets)
{

  size_t nClasses = classes.length();
  size_t nBg = classBg.length();
  size_t nSets = refSets.length();
  size_t nMap = class2ref.nrow();

  int useMap = doMap[0];

  const string knownBgTypes[4] = {"custom", "given", "intersection", "reference" };
  enum CBgTypes {custom, given, intersection, reference};
  const size_t nBgTypes = 4;

  int bgType = -1;
  for (size_t i=0; i<nBgTypes; i++) if (bgTypeR[0]==knownBgTypes[i]) bgType = i;

  if (bgType==-1)
    stop("Invalid bgTypeR given.");
  
  if (useMap && (nMap==0))
    stop("Empty (zero-row) ID translation matrix. Use doMap = FALSE to suppress mapping.");

  if (class2bg.length() < nClasses)
    stop("Length of class2bg must equal length of classes");

  // Check validity of class2bg. Note: it is assumed to use R indexing which is 1-based. Will convert it to 0-based later.
  for (size_t i=0; i<nClasses; i++) 
  {
    if (class2bg[i] < 1 || class2bg[i] >nBg) 
      stop("Invalid entry in class2bg: " + to_string(class2bg[i] ) + " is not between 1 and " + to_string(nBg) +
           " (number of backgrounds).");
  }

  // Create a vector of incomparables. Here hard-coded to NA.
  IntegerVector incomparables(1, NA_INTEGER);

  // Map the reference sets using the supplied map, if requested.
  IntegerVector class2ref_class = class2ref(_, 0),
                class2ref_ref = class2ref(_, 1);

  List mappedRefSets(nSets);
  if (useMap)
  {
    List set2map = match_int(refSets, class2ref_ref, incomparables);
    // Map the reference set to new IDs
    for (size_t set=0; set < nSets; set++)
    {
      IntegerVector map1 = set2map[set];
      IntegerVector mappedSet(nNonmissing(map1));
      size_t indx = 0;
      for (size_t i = 0; i < map1.length(); i++) if (map1[i]!=NA_INTEGER)
      {
        mappedSet[indx] = class2ref_class[ map1[i] ];
        indx++;
      }
      if (class2setIDHasDuplicates)
      {
        HashData_Int mappedSetHash(mappedSet, NA_INTEGER, incomparables, 1);
        size_t nUnique = mappedSetHash.nUnique(), indx = 0;
        IntegerVector uniqueMappedSet(nUnique);
        for (size_t i = 0; i < mappedSet.size(); i++) if (mappedSetHash.isUnique(i)) 
        {
          uniqueMappedSet[indx] = mappedSet[i];
          indx++;
        }
        mappedRefSets[set] = uniqueMappedSet;
      } else
        mappedRefSets[set] = mappedSet;
    }
  } else 
    for (size_t set=0; set < nSets; set++)
       mappedRefSets[set] = refSets[set];

  List effectiveBg;
  IntegerVector class2effectiveBg(nClasses);

  switch(bgType)
  {
     case given:
        effectiveBg = classBg;
        class2effectiveBg = clone(class2bg);
        break;
     case custom:
     case reference:
        effectiveBg.push_back(mappedRefBg);
        fill_n(class2effectiveBg.begin(), nClasses, 1); // Use 1-based indexing
        break;
     case intersection:
        effectiveBg = intersect_int(classBg, mappedRefBg, incomparables);
        class2effectiveBg = clone(class2bg);
        break;
     default: 
        stop("intersectSizesForEnrichment_C: Unhandled bgType " + to_string(bgType));
  } 
    
  // Convert R's 1-based indexing to C 0-based indexing.
  for (size_t i=0; i<nClasses; i++)
  {
    class2effectiveBg[i] = class2effectiveBg[i] - 1;
  }

  size_t nCombinedBg = effectiveBg.size();

  // Restrict classes to the respective classBg
  List classesInBg = intersect_int_list(classes, effectiveBg, class2effectiveBg, incomparables);
  // Calculate effective class sizes 
  IntegerVector effectiveClassSizes(nClasses);
  for (size_t cl=0; cl < nClasses; cl++) 
  {
    IntegerVector tmp = classesInBg[cl];
    effectiveClassSizes[cl] = tmp.length();
  }

  // Calculate effective set sizes and overlaps with respect to each background.
  IntegerMatrix effectiveSetSizes(nSets, nCombinedBg);
  for (size_t bg=0; bg < nCombinedBg; bg++)
  {
     IntegerVector bgv1 = effectiveBg[bg];
     IntegerVector ess1 = intersectSize_int(mappedRefSets, bgv1, incomparables);
     for (size_t set = 0; set < nSets; set++) 
       effectiveSetSizes(set, bg) = ess1[set];
  }
  
  IntegerMatrix overlapSizes(nSets, nClasses);
  for (size_t cl = 0; cl < nClasses; cl++)
  {
    IntegerVector class1 = classesInBg[cl];
    IntegerVector os1 = intersectSize_int(mappedRefSets, class1, incomparables);
    for (size_t set=0; set<nSets; set++) overlapSizes(set, cl) = os1[set];
  }

  List out;
  out["effectiveClassSizes"] = effectiveClassSizes,
  out["effectiveSetSizes"] = effectiveSetSizes,
  out["overlapSizes"] = overlapSizes,
  out["effectiveBackgrounds"] = effectiveBg,
  out["class2effectiveBg"] = class2effectiveBg;
  if (returnMappedSets[0]) out["mappedRefSets"] = mappedRefSets;

  return out;
}
    
