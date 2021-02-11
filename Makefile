#--------------------------------------------------------------------------------------------------#
# UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
# FACULTAD DE INGENIERÍA
# DIVISIÓN DE INGENIERÍA EN CIENCIAS DE LA TIERRA
# DEPARTAMENTO DE GEOFÍSICA
#
# Asignatura: Inversión de Datos Geofísicos
#   Profesor: Dr. Mauricio Nava Flores
#
# Archivo Makefile prueba
#
# Fecha de elaboración: 03/DIC/2020
# Autores:
#           Hernández Terán Oscar 
#           Jimenez Morales Nicole Amayrani
#           Martinez Reyes Xavier
#           Reyes Romero Alejandro
#--------------------------------------------------------------------------------------------------#
# Variables:
CF = gfortran
OC = -Wall -O3 -Wno-uninitialized
DE = /home/oscar/Dropbox/9Semestre/Inversion/ProyectoFinal/Códigos/Funciones

DExec = DirExe
DComp = DirComp
DResl = DirResl

VPATH = $(DE)

# Acciones:
all: dir comp
	$(CF) $(OC) $(DComp)/NumKind.o $(DComp)/bin2dec.o $(DComp)/dec2bin.o $(DComp)/mean.o $(DComp)/init_random_seed.o $(DComp)/AG_fxy.o \
	-o $(DExec)/exe

dir:
	- mkdir -p $(DExec) $(DComp) $(DResl)

comp: NumKind.o mean.o bin2dec.o dec2bin.o init_random_seed.o AG_fxy.f08
	$(CF) -I$(DComp) $(OC) -c AG_fxy.f08 -o $(DComp)/AG_fxy.o

NumKind.o: NumKind.f08
	$(CF) -J$(DComp) $(OC) -c $(DE)/NumKind.f08 -o $(DComp)/NumKind.o

bin2dec.o: bin2dec.f08
	$(CF) -J$(DComp) $(OC) -c $(DE)/bin2dec.f08 -o $(DComp)/bin2dec.o

dec2bin.o: dec2bin.f08
	$(CF) -J$(DComp) $(OC) -c $(DE)/dec2bin.f08 -o $(DComp)/dec2bin.o

mean.o: mean.f08
	$(CF) -J$(DComp) $(OC) -c $(DE)/mean.f08 -o $(DComp)/mean.o

init_random_seed.o: init_random_seed.f08
	$(CF) -J$(DComp) $(OC) -c $(DE)/init_random_seed.f08 -o $(DComp)/init_random_seed.o

exe:
	$(DExec)/exe

# Limpieza de directorio:
.PHONY: clean
clean:
	- rm -r $(DExec) $(DComp) $(DResl) exe
	clear
#--------------------------------------------------------------------------------------------------#
