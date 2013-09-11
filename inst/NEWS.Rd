\name{bumphunternews}
\title{Bumphunter News}
\encoding{UTF-8}

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

