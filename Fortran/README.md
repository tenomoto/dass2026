# README for Fortran

## Makefile

### Mac

No editing is required.

### Windows

Edit Makefile to point BLAS/LAPACK as follows.

    LIBS = -L"C:/rtools45/x86_64-w64-mingw32.rcrt/lib" -lblas -llapack
    # Mac
    #LIBS = -framework accelerate

### Other platforms

Install BLAS/LAPACK or an alternative (e.g. openblas) if not available.
Point these libraries in `LIBS`.

## Visualization

An R scirpt is provided to visualize the binary output.

    % ./var
    % Rscript plot.R var
    % ./enkf
    % Rscript plot.R enkf


