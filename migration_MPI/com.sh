cd ../funclib/;
ifort -c paramod.F90;
ifort -c *.F90;
cd ../migration_MPI/;

mpif90 ../funclib/*.o migrationloc.F90 -o migrationloc -I../funclib/