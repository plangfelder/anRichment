
# anRichment
# Annotation and Enrichment Functions 

Peter Langfelder, Jeremy A. Miller, Steve Horvath 

With contributions from Jim Wang, Mike Palazzolo, X. William Yang, CHDI, Rancho
Biosciences and others

### Contact
Peter (dot) Langfelder (at) gmail (dot) com


### R package anRichment 

Package anRichment is an add-on package for the [R statististical language and 
environment](https://www.r-project.org). The aim of this package is to provide tools for
enrichment calculations, access to multiple reference gene sets such as GO, KEGG, Reactome as well as
more specialized collections of gene sets collected from published literature, tools for creating
custom collections of reference gene sets and others.

### Caution: experimental package!

This package should be considered experimental. The "API" (such as it is), i.e., the function names,
arguments and returned values, can and probably will change. Furthermore, multiple functions (for example, building GO
collection and various identifier mapping functions) rely on annotation packages from
[Bioconductor](https://www.bioconductor.org) that are periodically updated, making it likely that running the same
calculations using different versions of these packages will produce (hopefully only slightly) different results.

### R tutorials

A [set of tutorials](Tutorials/index.html)
that illustrate various aspects of this annotation and enrichement package are available.

### Download and installation

Although the package can be installed under older R versions,
we strongly recommend using the newest R
version because the external annotation packages (GO, organism-specific genome annotation) are continuously updated.

__Simple installation__

The package can be installed using for example the function `github_install` from R package devtools. If devtools is not
installed in your R, install it first using 

    install.packages(devtools)

Then load the devtools package and use `github_install`:

    library("devtools");
    install_github("anRichment", username = "plangfelder");

A potential alternative is using the function `githubinstall` from the package also named githubinstall. The installation
and use is similar:

    # If needed, install the package first
    install.packages("githubinstall")
    library("githubinstall");
    githubinstall("anRichment");


__Manual installation__

One can also download the anRichment package directly from github and install it manually. This will require that all
dependencies are installed beforehand. For example, one may execute the following commands in R:


    needs = c("AnnotationDbi", "GO.db", "org.Hs.eg.db", "XML", "WGCNA",
              "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Mmusculus.UCSC.mm10.knownGene");
    installDeps = unlist(lapply(needs, function(pkg)
      if (require(pkg, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)) character(0) else pkg));
    if (length(installDeps)> 0)
    {
     cat("Installing dependencies...\n");
     if (!require("BiocManager")) {
       install.packages("BiocManager")
       library("BiocManager");
     }
     BiocManager::install(installDeps, suppressUpdates = TRUE, suppressAutoUpdate = TRUE);
     cat(".. done installing dependencies.\n");
    };


If this code generates package installation errors, investigate availability of the prerequisites for your version of R.
Contact Peter Langfelder if you encounter any other errors or bugs in the code.


### Problems installing or using the package

Most installation errors are caused by missing dependencies, potentially because they are unavailable for a particular
(usually old) R version. Please try to install the newest (or at least not very old, say no more than 2-3 years old)
version of R and all needed packages, then try installing anRichment again.

The package code itself may and probably does contain bugs causing errors or perhaps incorrect calculation results. If you
come across one, please first try the same code with the newest version of anRichment. If that fails to fix the problem, 
email Peter Langfelder. 

### Acknowledgments

This package originated as an attempt to merge Jeremy Miller's function `userListEnrichment` and the gene
lists it provides, and Peter Langfelder's `GOenrichmentAnalysis`. Both these functions were provided as part
of the WGCNA R package.

The package is currently maintained by Peter Langfelder and contains gene sets compiled by multiple
contributors. At the risk of putting out an incomplete list, these include Jeremy A. Miller, Michael
Palazzolo, Jim Wang, X. William Yang and his group (in particular, Jeff Cantle and Nan Wang),
Dan Salomon and his group, Jeff Aaronson and the team at Rancho Biosciences. PL
gratefully ackowledges the support of CHDI Foundation.


