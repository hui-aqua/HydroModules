"""
License
This file is part of UiS-Aqua.

UiS-Aqua is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

UiS-Aqua is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with UiS-Aqua. If not, see <https://www.gnu.org/licenses/>.


-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Plot figure(s)
please email: hui.cheng@uis.no \n
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.weight'] = 'regular'
plt.rcParams["mathtext.default"] = "it"
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams['image.cmap'] = 'summer'


g1 = 2
g2 = 1
gs = gridspec.GridSpec(g1, g2)  # Create 1x2 sub plots
def plotwaves(sea,time_slice,floder):
    print(sea)
    time_max = 3600  # [s]
    dt = 1
    time_frame = np.arange(0, time_max, dt)
    
    elevations = sea.get_elevations_with_time(np.array([0, 0, 0]), time_frame)
    position2 = np.zeros((200, 3))
    for i in range(200):
        position2[i] = [100-i, 0, 0]
    
    fig=plt.figure(figsize=(6.7, 6.5))
    t=time_slice
    ax = plt.subplot(gs[0, 0])
    plt.title("Elevation with time at position x=0, y=0")
    plt.plot(time_frame, elevations, color='b', linewidth=0.5)
    plt.xlabel("time (s)")
    plt.ylabel("surface elevation (m)")
    plt.xlim(0, 3600)
    plt.ylim(-5, 5)

    ax = plt.subplot(gs[1, 0])
    plt.title("Elevation at " + str(t) + " s")

    plt.plot([i for i in range(100)], sea.get_elevation_at_nodes(position2, t), color='b')
    plt.plot([-1000, 100], [0, 0], color='r', linewidth=0.5)
    plt.xlabel("X (m)")
    plt.ylabel("elevation (m)")
    plt.xlim(-100, 100)
    plt.ylim(-60, 6)

    plt.tight_layout()
    fig.savefig(floder, dpi=600)
    plt.close(fig)
