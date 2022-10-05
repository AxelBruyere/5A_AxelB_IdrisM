import numpy as np
import math
import cv2
from matplotlib import pyplot as plt

# Distance focale c270HD = 4.0mm


image_loin = cv2.imread('capture_mire_1.png')
image_proche = cv2.imread('capture_mire_0.png')

N = 7*7

tailleI = image_loin.shape
centre_Ix = tailleI[1]/2
centre_Iy = tailleI[0]/2
coord_px1 = cv2.findChessboardCorners(image_proche, (7, 7))
coord_px2 = cv2.findChessboardCorners(image_loin, (7, 7))
coord_px = []
coord_mm = []
coord_mx = 0
coord_my = 0

#Initialisation des matrices coord_px et coord_mm

for coord in coord_px1[1]:
    coord_px.append((coord[0][0], coord[0][1]))
for coord in coord_px2[1]:
    coord_px.append((coord[0][0], coord[0][1]))

print("Coord_px :", coord_px)

for i in range(N):
    coord_mm.append((coord_mx, coord_my, 0))
    coord_mx += 20
    if coord_mx == 140:
        coord_mx = 0
        coord_my += 20

coord_mx = 0
coord_my = 0

for i in range(N):
    coord_mm.append((coord_mx, coord_my, 120))
    coord_mx += 20
    if coord_mx == 140:
        coord_mx = 0
        coord_my += 20
U1 = []
U2 = []
A = []
for i in range(N*2):
    U1.append(coord_px[i][0]-centre_Ix)
    U2.append(coord_px[i][1]-centre_Iy)


for i in range(N*2):
    a1 = U2[i] * coord_mm[i][0]
    a2 = U2[i] * coord_mm[i][1]
    a3 = U2[i] * coord_mm[i][2]
    a4 = U2[i]
    a5 = -U1[i] * coord_mm[i][0]
    a6 = -U1[i] * coord_mm[i][1]
    a7 = -U1[i] * coord_mm[i][2]
    A.append((a1, a2, a3, a4, a5, a6, a7))
A_inv = np.linalg.pinv(A)

L = A_inv.dot(U1)

o2c = 1 / math.sqrt(L[4]*L[4]+L[5]*L[5]+L[6]*L[6])

Beta = o2c * math.sqrt(L[0]*L[0]+L[1]*L[1]+L[2]*L[2])
o2c = -o2c

o1c = L[3]*o2c/Beta
r11 = L[0]*o2c/Beta
r12 = L[1]*o2c/Beta
r13 = L[2]*o2c/Beta
r21 = L[4]*o2c
r22 = L[5]*o2c
r23 = L[6]*o2c

r3 = np.cross((r11, r12, r13), (r21, r22, r23))
r31 = r3[0]
r32 = r3[1]
r33 = r3[2]

gamma = - math.atan(r12/r11)
phi = - math.atan(r23/r33)
omega = math.atan(r13/(-r23 * math.sin(phi) + r33 * math.cos(phi)))

B = []
R = []
for i in range(N*2):
    B.append((U2[i], - (r21*coord_mm[i][0]+r22*coord_mm[i][1]+r23*coord_mm[i][2]+o2c)))
    R.append(-U2[i] * (r31*coord_mm[i][0]+r32*coord_mm[i][1]+r33*coord_mm[i][2]))

B_inv = np.linalg.pinv(B)
M = B_inv.dot(R)
f = 4

s2 = f / (M[1]/Beta)
s1 = f / M[1]
S = []
X = []
Mext = []
print("Centre en X ", centre_Ix)
print("Center en Y", centre_Iy)
print("*************************************")
print("r11", r11)
print("r12", r12)
print("r13", r13)
print("r21", r21)
print("r22", r22)
print("r23", r23)
print("*************************************")
print("r31", r31)
print("r32", r32)
print("r33", r33)
print("*************************************")
print("Beta", Beta)
print("oc1", o1c)
print("oc2", o2c)
print("oc3", M[0])
print("*************************************")
print("omega", omega)
print("phi", phi)
print("gamma", gamma)
print("*************************************")
print("f1", f/s1)
print("f2", f/s2)
print("*************************************")
print("s1", s1)
print("s2", s2)
print("*************************************")


S.append((f/s1, 0, centre_Ix, 0))
S.append((0, f/s2, centre_Iy, 0))
S.append((0, 0, 1, 0))

Mext.append((r11, r12, r13, o1c))
Mext.append((r21, r22, r23, o2c))
Mext.append((r31, r31, r33, M[0]))
Mext.append((0, 0, 0, 1))
Ualpha = np.asarray(S).dot(np.asarray(Mext))
red = [0,0,255]

for x in range(7):
    for y in range(7):
        X = (x*20, y*20, 0, 1)
        U_new = Ualpha.dot(np.transpose(np.asarray(X)))
        print(U_new[2])

        pixelx = int(U_new[0]/U_new[2])
        pixely = int(U_new[1]/U_new[2])

  
        

        
        image_modif = image_proche
        image_modif[pixely, pixelx] = red
        image_modif[pixely+1, pixelx] = red
        image_modif[pixely, pixelx+1] = red
        image_modif[pixely, pixelx-1] = red
        image_modif[pixely-1, pixelx] = red
#cv2.imshow('image', image_modif)
cv2.imwrite("result_pierre.png", image_modif)
print(X)











