for prg in calcF CalcLaplacian copier time_derivative; do
    echo Compile ${prg}
    mpicxx -o ${prg} `/opt/adios2/bin/adios2-config --cxx-flags` ${prg}.cpp `/opt/adios2/bin/adios2-config --cxx-libs`
done

