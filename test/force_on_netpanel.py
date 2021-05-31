import os
import sys
import numpy as np
import matplotlib.pyplot as plt


print(os.getcwd())
sys.path.append(os.path.join(os.getcwd()))  

import hydroModules as hdm
import enviromentModules as env


u=[1.0,0,0]

node_position=[]
for i in range(11):
    for j in range(11):
        node_position.append([0.0,round(0.1*i,2),round(-0.1*j,2)])
        
node_position=np.array(node_position)

# print(len(node_position))
net_panel=[]
for i in range(10):
    for j in range(10):
        net_panel.append([int(1+i+11*j),    int(1+i+1+11*j),
                          int(1+i+11*(j+1)),int(1+i+1+11*(j+1))])
# print(net_panel)
# print(len(net_panel))
# print(node_position[net_panel[1]])


def show_panel(nodes,element):
    ax = plt.figure().add_subplot(projection='3d')

    ax.scatter(nodes[:,0], nodes[:,1], nodes[:,2], color='r')
    for each in element:
        list_of_node=[each[0],each[1],each[3],each[2],each[0]]
        ax.plot(nodes[list_of_node][:,0],nodes[list_of_node][:,1],nodes[list_of_node][:,2],color='b')
       
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()


def coe(modle,inflow_angle,sn,re):
    drag_coefficient, lift_coefficient = 0, 0
    if modle == 'S1':  # aarsnes 1990
        drag_coefficient = 0.04 + (-0.04 + sn - 1.24 * pow(sn, 2) + 13.7 * pow(sn, 3)) * np.cos(inflow_angle)
        lift_coefficient = (0.57 * sn - 3.54 * pow(sn, 2) + 10.1 * pow(sn, 3)) * np.sin(2 * inflow_angle)
            
    if modle == 'S2':  # Loland 1991
        drag_coefficient = 0.04 + (-0.04 + 0.33 * sn + 6.54 * pow(sn, 2) - 4.88 * pow(sn, 3)) * np.cos(inflow_angle)
        lift_coefficient = (-0.05 * sn + 2.3 * pow(sn, 2) - 1.76 * pow(sn, 3)) * np.sin(2 * inflow_angle)
            
    if modle == 'S3':  # Kristiansen 2012
        a1 = 0.9
        a2 = 0.1
        b1 = 1.0
        b2 = 0.1
        reynolds_number = re
        cd_cylinder = -78.46675 \
                     + 254.73873 * np.log10(reynolds_number) \
                     - 327.88640 * pow(np.log10(reynolds_number), 2) \
                     + 223.64577 * pow(np.log10(reynolds_number), 3) \
                     - 87.922340 * pow(np.log10(reynolds_number), 4) \
                     + 20.007690 * pow(np.log10(reynolds_number), 5) \
                     - 2.4489400 * pow(np.log10(reynolds_number), 6) \
                     + 0.1247900 * pow(np.log10(reynolds_number), 7)
        cd_zero = cd_cylinder * (sn * (2 - sn)) / (2.0 * pow((1 - sn), 2))
        cn_45 = cd_cylinder * sn / (2.0 * pow((1 - sn), 2))

        cl_45=np.pi*cn_45/(8.0+cn_45)
        cl_zero=(0.5 * cd_zero -cl_45)/pow(2,0.5)
        drag_coefficient = cd_zero *(a1 * np.cos(inflow_angle) + a2 * np.cos(3 * inflow_angle))
        lift_coefficient = cl_zero *(b1 * np.sin(2 * inflow_angle) + b2 * np.sin(4 * inflow_angle))
                                     
    if modle == 'S4':  # Fridman 1973
        reynolds_number =re
        reynolds_star = reynolds_number / (2 * sn)
        coe_tangent = 0.1 * pow(reynolds_number, 0.14) * sn
        coe_normal = 3 * pow(reynolds_star, -0.07) * sn
        drag_coefficient = coe_normal * np.cos(inflow_angle) * pow(np.cos(inflow_angle), 2) + coe_tangent * np.sin(
            inflow_angle) * pow(
            np.sin(inflow_angle), 2)
        lift_coefficient = coe_normal * np.sin(inflow_angle) * pow(np.cos(inflow_angle), 2) + coe_tangent * np.cos(
            inflow_angle) * pow(
            np.sin(inflow_angle), 2)

    if modle == 'S5':  # Lee 2005 # polynomial fitting
        drag_coefficient = 0.556 * pow(inflow_angle, 7) - 1.435 * pow(inflow_angle, 6) - 2.403 * pow(
            inflow_angle, 5) + 11.75 * pow(inflow_angle, 4) - 13.48 * pow(inflow_angle, 3) + 5.079 * pow(
            inflow_angle, 2) - 0.9431 * pow(inflow_angle, 1) + 1.155
        lift_coefficient = -10.22 * pow(inflow_angle, 9) + 69.22 * pow(inflow_angle, 8) - 187.9 * pow(
            inflow_angle, 7) + 257.3 * pow(inflow_angle, 6) - 181.6 * pow(inflow_angle, 5) + 59.14 * pow(
            inflow_angle, 4) - 7.97 * pow(inflow_angle, 3) + 2.103 * pow(
            inflow_angle, 2) + 0.2325 * pow(inflow_angle, 1) + 0.0129

    return [drag_coefficient,lift_coefficient]
    



hdm.two_dimensional.row_water=1000.0
net=hdm.two_dimensional.netting('S3',net_panel,0.2,0.005)


fb=net.cal_buoy_force(node_position)
fh=net.force_on_element(node_position,u,)
f_node=net.distribute_force(node_position)

# f=net.cal_buoy_force(node_position)
# print('buoycy force is ')
# print(fb)
# print('total force')
# print(fb.sum(axis=0))


# print('hydrodynamic force is ')
# print(fh)
# print('total force')
# print(fh.sum(axis=0))

theta=np.linspace(0.0,np.pi/2.0-0.01,20)
f_=[]


cd_cl=[]


for th in theta:
    print(th)
    cd_cl.append(coe('S3',th,0.2,1000.0))
    # then rotate the nodes along Y axis
    tra_ma = np.array([[np.cos(th), 0,        -np.sin(th)],
                       [0,          1,        0],
                       [np.sin(th),          0,        np.cos(th)]])


    posi=np.dot(node_position,tra_ma)
    fh=net.force_on_element(posi,u,)
    f_.append(fh.sum(axis=0))
    # show_panel(posi,net_panel)
f_=np.array(f_)    
cd_cl=np.array(cd_cl)
fig,ax=plt.subplots(1,3)
ax[0].set_title('FX')
ax[0].plot(theta,f_[:,0])

ax[1].set_title('FZ')
ax[1].plot(theta,-f_[:,2])

ax[2].set_title('FZ/FX')
ax[2].plot(theta,-f_[:,2]/f_[:,0],label='FZ')
ax[2].plot(theta,cd_cl[:,1]/cd_cl[:,0])
# plt.legend()
plt.show()


# fig,ax=plt.subplots(1,3)
# for a in angles:
#     cd_cl.append(coe('S1',a,0.2,1000.0))
    
# cd_cl=np.array(cd_cl)
# ax[0].set_title('cd')
# ax[0].plot(angles,cd_cl[:,0])


# ax[1].set_title('cl')
# ax[1].plot(angles,cd_cl[:,1])


# ax[2].set_title('cl/cd')
# ax[2].plot(angles,cd_cl[:,1]/cd_cl[:,0])

# plt.show()