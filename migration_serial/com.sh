cd ../funclib/;
ifort -c paramod.F90;
ifort -c *.F90;
cd ../migration_serial/;

ifort ../funclib/*.o migrationloc.F90 -o migrationloc -I../funclib/