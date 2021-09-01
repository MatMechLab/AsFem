#!/usr/bin/env python3
'''
This script can read a 2D image and convert it into geo file
@Author: Yang Bai
@Date  : 2021.08.29 
@Usage : Image2Mesh.py mypic.jpg -range 1.0 2.0
'''
import matplotlib.pyplot as plt
from scipy.interpolate import interpolate
import numpy as np
import sys

minval=250;maxval=256

Dx=0.1;Dy=0.1  # pixel size <----> mesh size

imgname=sys.argv[1]
print('************ start to read image ...')


iInd=0
for args in sys.argv:
    iInd+=1
    if '-range' in args:
        break
if len(sys.argv)>=iInd+2:
    if '-range' in sys.argv[iInd-1]:
        minval=float(sys.argv[iInd+1-1])
        maxval=float(sys.argv[iInd+2-1])

iInd=0
for args in sys.argv:
    iInd+=1
    if '-dx' in args:
        break
if len(sys.argv)>=iInd+1:
    if '-dx' in sys.argv[iInd-1]:
        Dx=float(sys.argv[iInd+1-1])
        Dy=Dx

img=plt.imread(imgname)
print('*** image file is:%s'%(imgname))
print('*** threashhold value is: [%f,%f]'%(minval,maxval))
print('*** mesh size=%f'%(Dx))

def rgb2gray(rgb):
    return np.dot(rgb[...,:3], [0.299, 0.587, 0.144])

grayimg=rgb2gray(img)
print('*** image resolution=',np.shape(grayimg))

Nx=np.shape(grayimg)[0];Ny=np.shape(img)[1]

gray1d=grayimg.ravel()
print('*** nx=%g,ny=%g'%(Nx,Ny))
print('*** image size=(%f,%f)'%((Ny+1)*Dx,(Nx+1)*Dy))

Theta=-180.0/180.0*np.pi

X=np.zeros(Nx*Ny)
Y=np.zeros(Nx*Ny)
Value=np.zeros(Nx*Ny)
px=[];py=[];pv=[];tol=2.0e-2
for j in range(Ny):
    for i in range(Nx):
        x=i*Dx;y=j*Dy
        X[j*Nx+i]=y #x*np.cos(Theta)-y*np.sin(Theta)
        Y[j*Nx+i]=x #x*np.sin(Theta)+y*np.cos(Theta)
        Value[j*Nx+i]=grayimg[i,j]
        if grayimg[i,j]>=minval and grayimg[i,j]<=maxval:
            px.append(y)
            py.append(x)
            pv.append(grayimg[i,j])


plt.figure(1)
plt.imshow(grayimg,cmap = plt.get_cmap('gray'))

# plt.figure(2)
# plt.scatter(X,Y,Value)

plt.figure(3)
plt.scatter(px,py,pv)

px.append(0.0);py.append(0.0)
px.append((Ny+1)*Dx);py.append(0.0)
px.append((Ny+1)*Dx);py.append((Nx+1)*Dy)
px.append(0.0);py.append((Nx+1)*Dy)


filename=imgname[0:-4]+'.geo'
inp=open(filename,'w')
inp.write('SetFactory("OpenCASCADE");\n')
inp.write('dx=0.1;\n')
for i in range(len(px)):
    str='Point(%d)={%14.6e,%14.6e,0.0,dx};\n'%(i+1,px[i],py[i])
    inp.write(str)
inp.write('\n\n')

inp.close()

print('write result to %s'%(filename))


plt.show() 
