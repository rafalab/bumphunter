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
  }
}

