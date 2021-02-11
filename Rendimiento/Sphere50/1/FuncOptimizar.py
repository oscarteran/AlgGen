# -*- coding: utf-8 -*-
"""
Spyder Editor

Este es un programa en Python que sirve para la visualizacion de los datos
Si, le dedique más esfuerzo que el que le debí haber dedicado.
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

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
    ax.scatter(x,y,marker='.')
    plt.show()
    
#def grafDatos2d(x,y,z,title,xlabel,ylabel):

#Se ingresa el archivo
f = DefinirRuta('CurvaConv24.dat') #curva de convergencia
f1 = DefinirRuta('Resultados.dat')

#Lectura de datos
M = LeerOpt(f,2)
M1 = LeerOpt(f1,2)#Se pone 1 si se quiere saltar el primer renglon

#Para graficar datos 1d:
#grafDatos1d(M[:,0], M[:,1], 'Perfil', 't[s]', 'x[m]')
#Para graficar datos 2d:
    
#Para curva de convergencia
grafCurvConv(M[:,0],M[:,1])

#Análisis de rendimiento:
#[x,y]=mallaCuadrada(-5.12, 5.12) #Rastrigin
[x,y]=mallaCuadrada(-5, 5) #Esférica
#Esferica
z=Esferica(x, y)
#Rastrigin
#z=Rastrigin(x, y)

grafSuperficie(x, y, z, M1[:,3], M1[:,4], M1[:,2])