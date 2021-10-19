import numpy as np
from . import Airywave as wave
from . import wave_spectrum as wsp


# import Airywave as wave
# import wave_spectrum as wsp


class irregular_sea:
    """
    The default wave spectra is JONSWAP
    ...
    ----------
    *NOTE: The maximum applied simulation time should be less than 3h.
    """

    # TODO change the input to wave spectra
    def __init__(self, significant_wave_height, peak_period, gamma, water_depth, wave_direction):
        """
        Parameters
        ------------
        significant_wave_height:significant wave height | float | Unit [m]
        peak_period:peak period | float | Unit [s]
        gamma:gamma | float | Unit [-]
        water_depth: water depth of the sea, assume flat sea floor. A position number | float | Unit [m]
        wave_direction: direction of wave propagation. | float | Unit [degree]
        """
        self.hs = significant_wave_height
        self.tp = peak_period
        self.water_depth = water_depth
        self.list_of_waves = []

        time_max = 3600 * 3  # 3h, We assume the simulations will not exceed 3h.
        fre_max = 3  # we assume the highest eigen frequency of studied structure is below 3 Hz.
        d_fre = 2 * np.pi / time_max  # get the resolution for frequency.
        fre_range = np.arange(d_fre, fre_max, d_fre)
        design_wave_spectra = wsp.jonswap_spectra(fre_range, significant_wave_height, peak_period, gamma)
        list_xi = np.sqrt(2 * d_fre * design_wave_spectra)
        for index, item in enumerate(list_xi):
            wave_period = 2 * np.pi / fre_range[index]
            self.list_of_waves.append(
                wave.Airywave(item * 2, wave_period, water_depth, wave_direction, np.random.uniform(0, 360)))

    def __str__(self):
        """ Print the information of the present object. """
        s0 = 'The environment is irregular waves condition and the specific parameters are:\n'
        s1 = 'significant wave height = ' + str(self.hs) + ' m\n'
        s2 = 'peak period= ' + str(self.tp) + ' s\n'
        s3 = 'Number of wave components is ' + str(len(self.list_of_waves))
        return s0 + s1 + s2 + s3

    ## elevation

    def get_elevations_with_time(self, position, time_list):
        """
        Public function.\n
        :param position: [np.array].shape=(n,3) | Unit: [m]. The position of one node
        :param time_list: [np.array].shape=(n,1) | Uint: [s]. The time sequence for getting the elevations \n
        :return: Get a list of elevations at one position with a time sequence \n
        """
        wave_elevations = np.zeros((len(self.list_of_waves), len(time_list)))
        for index, each_wave in enumerate(self.list_of_waves):
            wave_elevations[index] = each_wave.get_elevation(position, time_list)
        return np.sum(wave_elevations, axis=0)

    def get_elevation_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point: [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: time [s] \n
        :return: Get a list of elevation at a list of point \n
        """
        wave_elevations = np.zeros((len(self.list_of_waves), len(list_of_point)))
        for index, each_wave in enumerate(self.list_of_waves):
            wave_elevations[index] = each_wave.get_elevation_at_nodes(list_of_point, global_time)
        return np.sum(wave_elevations, axis=0)

    ## velocity
    def get_velocity_with_time(self, position, time_list):
        """
        Public function.\n
        :param position: [np.array].shape=(n,3) | Unit: [m]. The position of one node
        :param time_list: [np.array].shape=(n,1) | Uint: [s]. The time sequence for getting the elevations \n
        :return: Get a list of elevations at one position with a time sequence \n
        """
        waves_velocities = np.zeros((len(self.list_of_waves), len(time_list), 3))
        for index, each_wave in enumerate(self.list_of_waves):
            waves_velocities[index] = each_wave.get_velocity_with_time(position, time_list, irregularwaves=True)
        return np.sum(waves_velocities, axis=0)

    def get_velocity_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        node_velocity_set = np.zeros((len(self.list_of_waves), len(list_of_point), 3))
        list_elevation = self.get_elevation_at_nodes(list_of_point, global_time)
        points = list_of_point.copy()
        # wheeler stretching method
        points[:, 2] = (list_of_point[:, 2] - list_elevation) * self.water_depth / (self.water_depth + list_elevation)
        # linear stretching method
        # points[:, 2][points[:, 2] >= 0] = 0
        for index, each_wave in enumerate(self.list_of_waves):
            node_velocity_set[index] = each_wave.get_velocity_at_nodes(points, global_time, irregularwaves=True)
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
        node_acceleration_set = np.zeros((len(self.list_of_waves), len(list_of_point), 3))
        list_elevation = self.get_elevation_at_nodes(list_of_point, global_time)
        points = list_of_point.copy()
        # wheeler stretching method
        points[:, 2] = (list_of_point[:, 2] - list_elevation) * self.water_depth / (self.water_depth + list_elevation)
        # linear stretching method
        # points[:, 2][points[:, 2] >= 0] = 0
        for index, each_wave in enumerate(self.list_of_waves):
            node_acceleration_set[index] = each_wave.get_acceleration_at_nodes(points, global_time, irregularwaves=True)
        accelerations = np.sum(node_acceleration_set, axis=0)
        accelerations[list_of_point[:, 2] > list_elevation] = 0
        return accelerations


if __name__ == "__main__":
    sea_state = irregular_sea(5, 8, 3, 60, 0)
    print(sea_state)
