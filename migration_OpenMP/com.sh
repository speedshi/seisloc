FC=ifort

if [ "$FC" = "ifort" ]; then
   FCFLAG="-warn"
else
   FCFLAG="-ffree-line-length-none"
fi


cd ../funclib/;
${FC} -c paramod.F90;
${FC} -c *.F90;
cd ../migration_OpenMP/;

${FC} ${FCFLAG} ../funclib/*.o migrationloc.F90 -o migrationloc -I../funclib/ -fopenmp -m64 
