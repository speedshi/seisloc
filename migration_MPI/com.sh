FC=ifort
FCMPI=mpif90

if [ "$FC" = "ifort" ]; then
   FCFLAG="-warn"
else
   FCFLAG="-ffree-line-length-none"
fi


cd ../funclib/;
${FC} -c paramod.F90;
${FC} -c *.F90;
cd ../migration_MPI/;

${FCMPI} ${FCFLAG} ../funclib/*.o migrationloc.F90 -o migrationloc -I../funclib/