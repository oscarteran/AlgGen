# -*- coding: utf-8 -*-
"""
Spyder Editor

Este es un programa en Python que sirve para la visualizacion de los datos
Sirve para ver:
    -Curvas de convergencia
    -Superficies aproximadas
    -Datos Geofísicos: Para esto hay que definir si los aarchivos traen encabezado o no
    
Si, le dedique más esfuerzo que el que le debí haber dedicado.
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import math
import pandas as pd

#Funciones no propias:
def truncate(number, decimals=0):
    """
    Returns a value truncated to a specific number of decimal places.
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer.")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more.")
    elif decimals == 0:
        return math.trunc(number)

    factor = 10.0 ** decimals
    return math.trunc(number * factor) / factor

#Funciones propias:
def Esferica(x,y):
    return x ** 2+y**2

def Rastrigin(x,y):
    return x ** 2 - 10 * np.cos(2 * np.pi * x) ** 2 + y ** 2 - 10 * np.cos(2 * np.pi * y) + 20

def mallaCuadrada(Min,Max):
    x = np.linspace(Min,Max,1000)
    y = np.linspace(Min,Max,1000)
    x,y= np.meshgrid(x,y)
    return x,y

def GrafContorno(x,y,z): #falta poner los datos 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set(title = 'Anomalía magnética', xlabel = 'x[Km]', ylabel = 'y[km]')
    cs = ax.contourf(x,y,z)
    c =fig.colorbar(cs)
    plt.show()
    
#Lectura de datos nDimensionales:
def LeerOpt(filename,case):  
    lines = len(open(filename).readlines())
    if case == 1:
        lines = lines - 1
    with open(filename) as f:
        line = f.readline()
    col = (len(line.split()))
    f.close()
    
    M = np.zeros((lines,col))
    cnt=0
    with open(filename,'r') as f:
        if case == 1:
            next(f)
        lin = f.readline()
        for i in range(lines):
            lin_sep = lin.split()
            #Por alguna razón para leer los datos de sísmica, me sale con col-1, para los demas es co
            #for j in range(col) 
            for j in range(col):
                M[i,j]=(float(lin_sep[j]))
            lin = f.readline()
            cnt= cnt + 1
        f.close()
    return M

def DefinirRuta(NombreArchivo):
    if os.path.exists(NombreArchivo):
        ruta =os.getcwd()
        filename= ruta +'/'+ NombreArchivo
    else:
        sys.exit('No existe el archivo')
    return filename

def grafSuperficie(x,y,z,x1,y1,z1):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    cs = ax.contourf(x,y,z, cmap='inferno')
    ax.scatter(x1,y1,c='g',marker ='.')
    c =fig.colorbar(cs)
    plt.suptitle('Esférica Q=50 Pc =1.0 Pm =0.1')
    plt.savefig('EsfQ50Pc1.0.png')
    plt.show()
    
def grafCurvConv(x,y):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set(title = 'Curva de Convergencia', xlabel = '# de generaciones', ylabel = 'Convergencia')
    ax.scatter(x,y,marker='.')
    plt.suptitle('Esférica Q=50 Pc =1.0 Pm =0.1')
    plt.grid()
    plt.savefig('ConvEsfQ50Pc1.0.png')
    plt.show()
    
def grafDatos1d(x,y,title,xlabel,ylabel):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if title == 'SEV':
        ax.set(title = title, xlabel = xlabel, ylabel = ylabel, xscale = 'log')
    else:
        ax.set(title = title, xlabel = xlabel, ylabel = ylabel, xscale = 'linear')
    ax.scatter(x, y, marker='.')
    plt.grid()
    plt.show()

def grafComp1d(x,y,x1,y1,title,xlabel,ylabel):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if title == 'SEV':
        ax.set(title = title, xlabel = xlabel, ylabel = ylabel, xscale = 'log')
    else:
        ax.set(title = title, xlabel = xlabel, ylabel = ylabel, xscale = 'linear')
    ax.plot(x, y, c = 'r',label= 'Datos Reales')
    ax.plot(x1, y1, c = 'b' ,label= 'Datos Invertidos')
    plt.legend()
    plt.grid()
    plt.show()
	

def MallaDatos(x,y,z):
    for i in range(len(x)):
        dx= truncate(np.abs(x[i+1]-x[i]),5)
        if (dx != 0):
            break
    for i in range(len(y)):
        dy= truncate(np.abs(y[i+1]-y[i]),5)
        if (dy != 0):
            break
    
    Nx = int(((np.max(x) - np.min(x)) / (dx)) + 1)
    Ny = int(((np.max(y) - np.min(y)) / (dy)) + 1)
    
    x.shape=(Nx,Nx)
    y.shape=(Ny,Ny)
    z.shape=(Nx,Ny)
    return x,y,z    



#Programa principal
#Se ingresa el archivo
f1 = DefinirRuta('Sismg_1.dat') #Datos observados
f = DefinirRuta('AnomOpt.dat') #Invertida

#Lectura de datos
M = LeerOpt(f1,2)
M1 = LeerOpt(f,2)
 

#M1 = LeerOpt(f1,2)#Se pone 1 si se quiere saltar el primer renglon

#Para graficar datos 1d:
#grafDatos1d(M[:,0], M[:,1], 'SEV', 'ab/2[m]', 'Rho Aparente')
#Para comparar resultados:
grafComp1d(M[:,0],M[:,1],M1[:,0],M1[:,1],'Sismograma','Tiempo[ms]','Distancia[m]')

#Para graficar datos 2d:
#x, y, z = MallaDatos(M[:,0], M[:,1], M[:,2])
#GrafContorno(x, y, z)
#Para curva de convergencia
#grafCurvConv(M[:,0],M[:,1])

#Análisis de rendimiento:
#[x,y]=mallaCuadrada(-5.12, 5.12) #Rastrigin
#[x,y]=mallaCuadrada(-5, 5) #Esférica
#Esferica
#z=Esferica(x, y)
#Rastrigin
#z=Rastrigin(x, y)

#grafSuperficie(x, y, z, M1[:,3], M1[:,4], M1[:,2])
