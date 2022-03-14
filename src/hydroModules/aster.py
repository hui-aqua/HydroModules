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
Any questions about this code,
please email: hui.cheng@uis.no \n
"""

import numpy as np

# Functions used by code_aster
def get_position_aster(table_aster):
    """
    Module public function.\n
    :param table_aster: A table from Code_Aster by ``POST_RELEVE_T`` command with NOM_CHAM=('DEPL')
    :return:  [np.array].shape=(N,3) Unit: [m]. A numpy array of all the nodes positions.
    """
    content = table_aster.EXTR_TABLE()
    original_x = content.values()['COOR_X']
    original_y = content.values()['COOR_Y']
    original_z = content.values()['COOR_Z']
    delta_x = content.values()['DX']
    delta_y = content.values()['DY']
    delta_z = content.values()['DZ']
    position = np.array([original_x, original_y, original_z]) + np.array([delta_x, delta_y, delta_z])
    return np.transpose(position)

def get_velocity_aster(table_aster):  # to get the velocity
    """
    Module public function.\n
    :param table_aster: A table from Code_Aster by ``POST_RELEVE_T`` command with NOM_CHAM=('VITE')
    :return:  [np.array].shape=(N,3) Unit: [m/s]. A numpy array of all the nodes velocities.
    """
    content = table_aster.EXTR_TABLE()
    velocity_x = content.values()['DX']
    velocity_y = content.values()['DY']
    velocity_z = content.values()['DZ']
    velocity = np.array([velocity_x, velocity_y, velocity_z])
    return np.transpose(velocity)

def get_acceleration_aster(table_aster):  # to get the acceration
    """
    Module public function.\n
    :param table_aster: A table from Code_Aster by ``POST_RELEVE_T`` command with NOM_CHAM=('ACCE')
    :return:  [np.array].shape=(N,3) Unit: [m/s2]. A numpy array of all the nodes velocities.
    """
    return get_velocity_aster(table_aster)


def get_displace_vector(table_aster):
    """
    Module public function.\n
    :param table_aster: A table from Code_Aster by ``POST_RELEVE_T`` command with NOM_CHAM=('DEPL')
    :return:  [np.array].shape=(N,3) Unit: [m]. A numpy array of all the nodes positions.
    """
    content = table_aster.EXTR_TABLE()
    delta_x = content.values()['DX']
    delta_y = content.values()['DY']
    delta_z = content.values()['DZ']
    return np.mean(np.np.transpose(np.array([delta_x, delta_y, delta_z])),axis=0)
