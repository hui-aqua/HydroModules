import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs


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
    
cd,cl=coe('S2',0,0.2,1000.0)
print(cd,cl)
angles=np.linspace(0.0,np.pi/2.0-0.01,20)
cd_cl=[]



fig,ax=plt.subplots(1,3)
for a in angles:
    cd_cl.append(coe('S1',a,0.2,1000.0))
    
cd_cl=np.array(cd_cl)
ax[0].set_title('cd')
ax[0].plot(angles,cd_cl[:,0])


ax[1].set_title('cl')
ax[1].plot(angles,cd_cl[:,1])


ax[2].set_title('cl/cd')
ax[2].plot(angles,cd_cl[:,1]/cd_cl[:,0])

plt.show()