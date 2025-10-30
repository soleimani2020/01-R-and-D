import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import random
from itertools import islice
import pandas as pd
import statistics
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from itertools import chain


# This function gets three points -on graphen - , generates two in plane vectors , cross product the vectors on the plane , returns 
# the vector prependicular to the plane 

def Graphene_Normal_Vector(p0,p1,p2):
    
    points = [p0,p1,p2]
    p0, p1, p2 = points
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    # These two vectors are in the plane
    ux, uy, uz = u = [x1-x0, y1-y0, z1-z0] #first vector
    vx, vy, vz = v = [x2-x0, y2-y0, z2-z0] #sec vector
    # the cross product is a vector normal to the plane
    u_cross_v = [uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx] #cross product
    point  = np.array(p1)
    normal = np.array(u_cross_v) # The a,b,c coefficients are obtained from a vector normal to the plane
    d = -point.dot(normal) #  dot product of the normal vector with any one of the point position vectors
    # ax+by+cz-d=0
    #print('plane equation:\n{:1.4f}x + {:1.4f}y + {:1.4f}z + {:1.4f} = 0'.format(normal[0], normal[1], normal[2], d))
    xx, yy = np.meshgrid(range(10), range(10))
    z = (-normal[0] * xx - normal[1] * yy - d) * 1. / normal[2]
    #plt3d = plt.figure().gca(projection='3d')
    # Specify x0, y0 and z0 coordinates of the arrow locations
    # Specify normal[0], normal[1] and normal[2] direction components of the normal vector 
    #plt3d.quiver(x0, y0, z0, normal[0], normal[1], normal[2], color="m")
    #plt3d.plot_surface(xx, yy, z)
    #plt3d.set_xlabel("X", color='red', size=18)
    #plt3d.set_ylabel("Y", color='green', size=18)
    #plt3d.set_zlabel("Z", color='b', size=18)
    #plt.show()
    return list([0,0,1])



def Octanol_vector(P,Q):
    POINTS = [P,Q]
    P, Q  = POINTS
    px, py, pz = P
    qx, qy, qz = Q
    distance = [px -qx, py - qy ,pz-qz ]
    norm = math.sqrt(distance[0] ** 2 + distance[1] ** 2+distance[2] ** 2)
    direction = [distance[0] / norm, distance[1] / norm,distance[2] / norm]
    return norm ,  direction , distance


# this function gets two vectors and calculates the angle between the two vectors 

def Angle(vector_1,vector_2):
    unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
    unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    angle_rad = np.arccos(dot_product)
    deg = np.rad2deg(angle_rad)
    return deg 







def One_Snapshot1(filename,Down,Up):
    with open(filename) as f:
        next_n_lines = list(islice(f,2))
        next_n_lines = list(islice(f, 42)) 
        np_lines = np.array(next_n_lines)    
        np_lines_lists = np.char.split(np_lines) 
        np_lines_array = np.array(np_lines_lists.tolist())
        raw_data = pd.DataFrame(np_lines_array)
        raw_data_oct1 = raw_data.astype({0:'str', 1:'str', 2:'str',3:'float64', 4:'float64', 5:'float64'})
        X=raw_data_oct1[3]
        Y=raw_data_oct1[4]
        Z=raw_data_oct1[5]
        TRJ_Oxygen=[]
        TRJ_CARBON=[]
        ANGLES=[]
        
        for i in range(0,len(raw_data_oct1[0])):
            if  Down < raw_data_oct1[5][i] < Up :
                if (raw_data_oct1[1][i] == "BHQd" ) :
                    Trj_O=[raw_data_oct1[3][i],raw_data_oct1[4][i],raw_data_oct1[5][i]]
                    Trj_C=[raw_data_oct1[3][i+41],raw_data_oct1[4][i+41],raw_data_oct1[5][i+41]]
                    My_Octanol_vector=Octanol_vector(Trj_C,Trj_O)[2]
                    My_Graphene_Normal_Vector=Graphene_Normal_Vector([0.740000,2.060000 ,67.860001 ],[3.220000,2.010000 ,68.160004 ],[0.840000, 4.950000,67.639999])
                    My_Angle=Angle(My_Graphene_Normal_Vector,My_Octanol_vector)
                    ANGLES.append(My_Angle)
        My_Angles=ANGLES
        return (My_Angles)
    




def Angle_calc_one_snapshot(filename,Down,Up):
    list1=One_Snapshot1(filename,Down,Up)

    return list1 , list1





def Angle_calc_all_snapshots(Down,Up):
    All_S=[]
    for it in range(0,3,1):   #(6000,30000,200)
        FILE = "conf"+str(it)+".gro"
        All_S.append(Angle_calc_one_snapshot(FILE,Down,Up)[0])
    ENSEMBLE_S=All_S
    AVE_ENSEMBLE_S=np.mean(ENSEMBLE_S)
    
    return AVE_ENSEMBLE_S , ENSEMBLE_S


A = Angle_calc_all_snapshots(Down=-100, Up=100)[1]

A = list(chain.from_iterable(A))

print(A)

mean = statistics.mean(A)
sd = statistics.stdev(A)

fig = plt.figure()

# Use a histogram to visualize the distribution
plt.hist(A, bins=20, density=True, alpha=0.6, color='g')

# Plot the normal distribution curve
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mean, sd)
plt.plot(x, p, 'k', linewidth=2)

plt.axvline(mean, color='grey', linestyle='--', linewidth=2, label="Mean: {:.2f}, Std: {:.2f}".format(mean, sd))
plt.legend(loc='upper right')
plt.savefig("AAA.png")

print("Program terminated successfully.")





