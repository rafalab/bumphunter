\name{bumphunternews}
\title{Bumphunter News}
\encoding{UTF-8}

\section{Version 1.9}{
  \itemize{
    \item Fixed bug in the bootstrap method for bumphunter. The
  coefficient for which inference was calculated was always 2 instead of
  the coefficient coef passed as an argument.
  }
}

\section{Version 1.3}{
  \itemize{
    \item Fixed NAMESPACE issues.
    \item nearestgene() now loads TT.rda into the local environment
    \item known_transcripts(), and consequently matchGenes(), was
    yielding mostly NA for the Refseq / 'annotation' column (iff
    there were multiple refseq ids for the given Entrez/Gene id)
    \item Fixed data/TT.rda, as the new value of known_transcripts().
    \item Cleaned up NEWS.
    \item Added GitHub URL to DESCRIPTION.
    \item Cleaned up duplicates in Suggests:.
  }
}

\section{Version 1.1}{
  \itemize{
    \item Added NEWS file.
    \item Fixed a bug related to >= for numerics.
    \item Added smoothing using a gaussian kernal as implemented in the
    locfit package through the function locfitByCluster.
    \item Added closeSockets for cleanup for doParallel on Windows.
    \item More bugfixes for windows; now using foreachCleanup().
    \item Added a 'bumps' class and print method.
    \item annotateNearest / regionMatch now give NA annotations for
    queries with no nearest subject (perhaps because the seqname
    is missing from the subject). Previously this was taken to be
    mistaken input and a hard error raised.
    \item Speedup of fwer computations using foreach.
    \item Added boundedClusterMaker.
    \item bug fix to internal function .getModT (which are not used in
    the main bumphunter functions).  Now the t-statistics returned a correct.
  }
}

