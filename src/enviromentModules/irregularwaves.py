import numpy as np
from . import Airywave as wave

# import Airywave as wave
# import wave_spectrum as wsp


class summation:
    """
    Irregular random waves, representing a real sea state,
    can be modelled as a summation of sinusoidal wave components. 
    The simplest random wave model is the linear longcrested wave model given by
    DNV-RP-C205, Section 3.3.2.1

    """

    def __init__(self, waveSpectrum, water_depth, wave_direction):
        """
        Parameters
        ------------
        waveSpectrum: A n*2 array list of wave spectrum, the first colum is w, the second is the S(w)
        water_depth: water depth of the sea, assume flat sea floor. A position number | float | Unit [m]
        wave_direction: direction of wave propagation. | float | Unit [degree]
        """

        self.water_depth = water_depth
        self.list_of_waves = []
        d_fre = abs(waveSpectrum[1, 0]-waveSpectrum[0, 0])
        print(d_fre)
        for each in waveSpectrum:
            xi = np.sqrt(2 * d_fre * each[1])
            wave_period = 2 * np.pi / each[0]
            self.list_of_waves.append(
                wave.Airywave(xi * 2, wave_period, water_depth, wave_direction, np.random.uniform(0, 360)))

    def __str__(self):
        """ Print the information of the present object. """
        s0 = 'The environment is irregular waves condition and the specific parameters are:\n'
        s1 = 'significant wave height = ' + str(self.hs) + ' m\n'
        s2 = 'peak period= ' + str(self.tp) + ' s\n'
        s3 = 'Number of wave components is ' + str(len(self.list_of_waves))
        return s0 + s1 + s2 + s3

    # elevation

    def get_elevations_with_time(self, position, time_list):
        """
        Public function.\n
        :param position: [np.array].shape=(n,3) | Unit: [m]. The position of one node
        :param time_list: [np.array].shape=(n,1) | Uint: [s]. The time sequence for getting the elevations \n
        :return: Get a list of elevations at one position with a time sequence \n
        """
        wave_elevations = np.zeros((len(self.list_of_waves), len(time_list)))
        for index, each_wave in enumerate(self.list_of_waves):
            wave_elevations[index] = each_wave.get_elevation(
                position, time_list)
        return np.sum(wave_elevations, axis=0)

    def get_elevation_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point: [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: time [s] \n
        :return: Get a list of elevation at a list of point \n
        """
        wave_elevations = np.zeros(
            (len(self.list_of_waves), len(list_of_point)))
        for index, each_wave in enumerate(self.list_of_waves):
            wave_elevations[index] = each_wave.get_elevation_at_nodes(
                list_of_point, global_time)
        return np.sum(wave_elevations, axis=0)

    # velocity
    def get_velocity_with_time(self, position, time_list):
        """
        Public function.\n
        :param position: [np.array].shape=(n,3) | Unit: [m]. The position of one node
        :param time_list: [np.array].shape=(n,1) | Uint: [s]. The time sequence for getting the elevations \n
        :return: Get a list of elevations at one position with a time sequence \n
        """
        waves_velocities = np.zeros(
            (len(self.list_of_waves), len(time_list), 3))
        for index, each_wave in enumerate(self.list_of_waves):
            waves_velocities[index] = each_wave.get_velocity_with_time(
                position, time_list, irregularwaves=True)
        return np.sum(waves_velocities, axis=0)

    def get_velocity_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        node_velocity_set = np.zeros(
            (len(self.list_of_waves), len(list_of_point), 3))
        list_elevation = self.get_elevation_at_nodes(
            list_of_point, global_time)
        points = list_of_point.copy()
        # wheeler stretching method
        points[:, 2] = (list_of_point[:, 2] - list_elevation) * \
            self.water_depth / (self.water_depth + list_elevation)
        # linear stretching method
        # points[:, 2][points[:, 2] >= 0] = 0
        for index, each_wave in enumerate(self.list_of_waves):
            node_velocity_set[index] = each_wave.get_velocity_at_nodes(
                points, global_time, irregularwaves=True)
        velocities = np.sum(node_velocity_set, axis=0)
        # ensure the velocity is zero above the wave elevation
        velocities[list_of_point[:, 2] > list_elevation] = 0
        return velocities

    def get_acceleration_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of node positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        node_acceleration_set = np.zeros(
            (len(self.list_of_waves), len(list_of_point), 3))
        list_elevation = self.get_elevation_at_nodes(
            list_of_point, global_time)
        points = list_of_point.copy()
        # wheeler stretching method
        points[:, 2] = (list_of_point[:, 2] - list_elevation) * \
            self.water_depth / (self.water_depth + list_elevation)
        # linear stretching method
        # points[:, 2][points[:, 2] >= 0] = 0
        for index, each_wave in enumerate(self.list_of_waves):
            node_acceleration_set[index] = each_wave.get_acceleration_at_nodes(
                points, global_time, irregularwaves=True)
        accelerations = np.sum(node_acceleration_set, axis=0)
        accelerations[list_of_point[:, 2] > list_elevation] = 0
        return accelerations


if __name__ == "__main__":
    sea_state = irregular_sea(5, 8, 3, 60, 0)
    print(sea_state)
