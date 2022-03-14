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
wave direction is x+
"""
import numpy as np
from numpy import pi


class Airywave:
    """
    Using Airy wave theory      \n
    Ref. DNV GL-RP205 Ver. 2008:P45
    Linear wave theory.
    """

    def __init__(self, wave_height=1.0, wave_period=10.0, water_depth=60.0, direction=0.0, phase=0.0):
        """
        :param wave_height: [float] | Unit: [m]. wave height.
        :param wave_period: [float] | Unit: [s]. wave period.
        :param water_depth: [float] | Unit: [m]. wave depth.
        :param direction: [float] | Unit: [degree]. direction of propagation, measured from the positive x-axis.
        :param phase: [float] | Unit: [degree]. phase.
        """
        self.gravity = 9.81
        self.wave_Height = wave_height
        self.wave_Period = wave_period
        self.water_Depth = water_depth
        self.wave_beta = pi * direction / 180.0
        self.wave_phase = pi * phase / 180.0

        # 1 Calculation
        # wave length
        alpha = [1, 0.666, 0.445, -0.105, 0.272]
        omega_ba = 4.0 * pow(pi, 2) * water_depth / \
            self.gravity / pow(wave_period, 2)
        f_omega = 0.0
        for index, item in enumerate(alpha):
            f_omega += item * pow(omega_ba, index)
        self.wave_Length = self.wave_Period * pow(self.gravity * self.water_Depth, 0.5) * pow(
            f_omega / (1 + omega_ba * f_omega), 0.5)
        # wave number
        self.wave_k = 2 * pi / self.wave_Length
        # angular frequency
        self.omega = pow(self.gravity * self.wave_k *
                         np.tanh(self.wave_k * self.water_Depth), 0.5)
        # phase velocity
        self.wave_phase_velocity = pow(
            self.gravity / self.wave_k * np.tanh(self.wave_k * self.water_Depth), 0.5)

        # 2 for easy calculation
        self.pi_h_t = pi * wave_height / wave_period
        self.pi_h_t_2 = 2 * wave_height * pow(pi / wave_period, 2)

    def __str__(self):
        """ Print the information of object. """
        s0 = 'The environment is airy wave condition and the specific parameters are:\n'
        s1 = 'water Depth= ' + str(self.water_Depth) + ' m\n'
        s2 = 'wave Period= ' + str(self.wave_Period) + ' s\n'
        s2_1 = 'wave number= ' + str(self.wave_k) + ' 1/m\n'
        s3 = 'wave Length= ' + str(self.wave_Length) + ' m\n'
        s4 = 'wave Height= ' + str(self.wave_Height) + ' m\n'
        s5 = 'wave phase velocity= ' + str(self.wave_phase_velocity) + ' m/s\n'
        s6 = 'wave direction:= ' + str(self.wave_beta) + ' degree\n'
        return s0 + s1 + s2 + s2_1 + s3 + s4 + s5 + s6

    def calc_theta(self, position, global_time):
        """
        A private function. \n
        :param position: [np.array].shape=(n,3) or [np.array].shape=(n,2) coordinates Unit: [m]. \n
        The coordinate of the point which you want to know the wave surface elevation. can be [x,y] or [x,y,z]
        :param global_time: [float] Unit: [s].
        :return: [float] Unit: [m]. The sea surface level in Z direction. At the targeted position.
        """
        if len(position.shape) == 1:
            # only one point
            return self.wave_k * (position[0] * np.cos(self.wave_beta) + position[1] * np.sin(
                self.wave_beta)) - self.omega * global_time + self.wave_phase
        elif len(position.shape) == 2:
            # a list of point
            return self.wave_k * (position[:, 0] * np.cos(self.wave_beta) + position[:, 1] * np.sin(
                self.wave_beta)) - self.omega * global_time + self.wave_phase

    def get_elevation(self, position, global_time):
        """
        A private function. \n
        :param position: [np.array].shape=(n,3) coordinates Unit: [m]. The position of the point which you want to know the wave surface elevation.
        :param global_time: [float] Unit: [s].
        :return: scale or [np.array].shape=(n,) Unit: [m]. The sea surface level in Z direction.
        """
        return self.wave_Height / 2 * np.cos(self.calc_theta(position, global_time))

    def get_velocity_with_time(self, position, global_time, irregularwaves=False):
        """
        :param position: [np.array].shape=(1,3) | Unit: [m]. The position of the point which you want to know the wave velocity
        :param global_time: [float] Unit: [s]. The time which you want to know the wave velocity
        :return:  [np.array].shape=(1,3) Unit: [m/s]. A numpy array of the velocity at the targeted point.
        """
        theta = self.calc_theta(position, global_time)
        eta = self.get_elevation(position, global_time)
        position_z = position[2]
        velocity = np.zeros((len(global_time), 3))
        velocity[:, 0] = np.cos(self.wave_beta) * self.pi_h_t * np.cosh(
            self.wave_k * (position_z + self.water_Depth)) * np.cos(theta) / np.sinh(self.wave_k * self.water_Depth)
        velocity[:, 1] = np.sin(self.wave_beta) * self.pi_h_t * np.cosh(
            self.wave_k * (position_z + self.water_Depth)) * np.cos(theta) / np.sinh(self.wave_k * self.water_Depth)
        velocity[:, 2] = self.pi_h_t * np.sinh(self.wave_k * (position_z + self.water_Depth)) * np.sin(theta) / np.sinh(
            self.wave_k * self.water_Depth)
        if not irregularwaves:
            velocity[position[2] > eta] = 0.0
        return velocity

    def get_acceleration_with_time(self, position, global_time, irregularwaves=False):
        """
        :param position: [np.array].shape=(1,3) | Unit: [m]. The position of the point which you want to know the wave acceleration.
        :param global_time: [float] Unit: [s]. The time which you want to know the wave velocity
        :return: [np.array].shape=(1,3) Unit: [m/s]. A numpy array of the acceleration at the targeted point.
        """
        theta = self.calc_theta(position, global_time)
        eta = self.get_elevation(position, global_time)
        position_z = position[2]
        acceleration = np.zeros((len(global_time), 3))
        acceleration[:, 0] = np.cos(self.wave_beta) * self.pi_h_t_2 * np.cosh(
            self.wave_k * (position_z + self.water_Depth)) * np.sin(theta) / np.sinh(self.wave_k * self.water_Depth)
        acceleration[:, 1] = np.sin(self.wave_beta) * self.pi_h_t_2 * np.cosh(
            self.wave_k * (position_z + self.water_Depth)) * np.sin(theta) / np.sinh(self.wave_k * self.water_Depth)
        acceleration[:, 2] = -self.pi_h_t_2 * np.sinh(self.wave_k * (position_z + self.water_Depth)) * np.cos(
            theta) / np.sinh(
            self.wave_k * self.water_Depth)
        if not irregularwaves:
            acceleration[position[2] > eta] = 0.0
        return acceleration

    def get_elevations_with_time(self, position, time_list):
        """
        Public function.\n
        :param position: [np.array].shape=(n,3) Unit: [m]. The position of one node
        :param time_list: [np.array].shape=(n,1) | Uint: [s]. The time sequence for getting the elevations
        :return: Get a list of elevations at one position with a time sequence \n
        """
        return self.get_elevation(position, time_list)

    def get_elevation_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point: [np.array].shape=(n,3) Unit: [m]. A list of node positions
        :param global_time: time [s] \n
        :return: Get a list of elevation at a list of point \n
        """
        return self.get_elevation(list_of_point, global_time)

    def get_velocity_at_nodes(self, list_of_point, global_time, irregularwaves=False):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        theta = self.calc_theta(list_of_point, global_time)
        eta = self.get_elevation(list_of_point, global_time)
        positions_z = list_of_point[:, 2]
        velocities = np.zeros((len(list_of_point), 3))
        velocities[:, 0] = np.cos(self.wave_beta) * self.pi_h_t * np.cosh(
            self.wave_k * (positions_z + self.water_Depth)) * np.cos(theta) / np.sinh(self.wave_k * self.water_Depth)
        velocities[:, 1] = np.sin(self.wave_beta) * self.pi_h_t * np.cosh(
            self.wave_k * (positions_z + self.water_Depth)) * np.cos(theta) / np.sinh(self.wave_k * self.water_Depth)
        velocities[:, 2] = self.pi_h_t * np.sinh(self.wave_k * (positions_z + self.water_Depth)) * np.sin(theta) / np.sinh(
            self.wave_k * self.water_Depth)
        if not irregularwaves:
            velocities[list_of_point[:, 2] > eta] = 0.0
        return velocities

    def get_acceleration_at_nodes(self, list_of_point, global_time, irregularwaves=False):
        """
        Public function.\n
        :param list_of_point: [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: time [s] \n
        :return: Get a list of acceleration at a list of point \n
        """
        theta = self.calc_theta(list_of_point, global_time)
        eta = self.get_elevation(list_of_point, global_time)
        positions_z = list_of_point[:, 2]
        accelerations = np.zeros((len(list_of_point), 3))
        accelerations[:, 0] = np.cos(self.wave_beta) * self.pi_h_t_2 * np.cosh(
            self.wave_k * (positions_z + self.water_Depth)) * np.sin(theta) / np.sinh(self.wave_k * self.water_Depth)
        accelerations[:, 1] = np.sin(self.wave_beta) * self.pi_h_t_2 * np.cosh(
            self.wave_k * (positions_z + self.water_Depth)) * np.sin(theta) / np.sinh(self.wave_k * self.water_Depth)
        accelerations[:, 2] = -self.pi_h_t_2 * np.sinh(self.wave_k * (positions_z + self.water_Depth)) * np.cos(theta) / np.sinh(
            self.wave_k * self.water_Depth)
        if not irregularwaves:
            accelerations[list_of_point[:, 2] > eta] = 0.0
        return accelerations


if __name__ == "__main__":
    pass
