"""
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


class net2netWake:
    """
    A module that can be use to deal with net to net wake effect.\n
    **note1:** Can only apply to a single fish cage.
    **note2:** Only work for current velocity.
    **note3** Only work for the net element, not for cable/ pipes
    """

    def __init__(self, model_index, initial_node_position, hydro_element, direction, net_solidity=0):
        """
        :param model_index: [string] Unit: [-]. To indicate the model function, e.g.: 'factor-0.9', 'loland-0.5', 'hui-1'.
        :param initial_node_position: [np.array].shape=(N,3) Unit: [m]. The initial coordinates of N nodes in cartesian coordinate system.
        :param hydro_element: [list] Unit: [-]. A python list to indicate how the lines are connected.
        :param current_velocity: [np.array].shape=(1,3) or a [list] of three values Unit: [m/s]. The current velocity [ux,uy,uz] in cartesian coordinate system.
        :param origin: [np.array].shape=(1,3) or a [list] of three values Unit: [m]. The origin [x,y,z] for detecting the elements in the wake region. For a fish cage, the origin is usually sit in the floating collar.
        :param dw0: [float] Unit: [m]. The diameter of the physical net twines.
        :param net_solidity: [float] Unit: [-]. The solidity of netting.
        """

        self.positions = np.array(initial_node_position)

        if np.linalg.norm(direction) != 0:
            self.direction = np.array(direction)/np.linalg.norm(direction)
        else:
            self.direction = np.array([1.0, 0.0, 0.0])

        self.elements = hydro_element

        self.origin = np.mean(self.positions, axis=0)

        self.sn = net_solidity
        self.wake_type = str(model_index).split("-")[0]
        self.wake_value = str(model_index).split("-")[1]
        self.wake_element_indexes = self.get_element_in_wake()
        self.get_element_in_wake()
        self.reduction_factors = np.reshape([1.0]*len(self.elements), (-1, 1))
        self.update_reduction(initial_node_position)

        print("\n net2net wake effect is initialized.\n")

    def __str__(self):
        s0 = "The selected wake model is " + str(self.wake_type) + "\n"
        s1 = "The index of the nodes in the wake region is " + \
            str(self.wake_element_indexes) + "\n"
        S = s0 + s1
        return S

    def is_element_in_wake(self, one_element):
        """
        A private function to tell if an element is in wake region or not.\n
        :param one_element: [List] of three values Unit [-]. A list to contain the index of nodes in one element. The list should contain at least two values.
        :return: True if the element in the wave region.
        """
        element_center = np.array([0.0, 0.0, 0.0])  # must be a float
        for node in one_element:
            element_center += np.array(
                self.positions[int(node)]) / len(one_element)
        vector_point_to_origin = np.array(element_center - self.origin)
        if np.dot(vector_point_to_origin, self.direction) < 0:
            return True
        else:
            return False

    def get_element_in_wake(self):
        """
        A private function to go go through all the elements and find the elements in the wake region.\n
        :return: [List] Unit [-]. A list of "indexes of the elements" in the wake region.
        """
        elements_in_wake = []
        for index, element in enumerate(self.elements):
            if self.is_element_in_wake(element):
                elements_in_wake.append(index)
        return elements_in_wake

    def cal_factor1(self):
        """
        A private function calculate the flow reduction factor according to Loland (1991). r=1-0.14*CD (CD is the input value in the ```wakeModel```) \n
        :return: Flow reduction factor: [float] Unit [-].
        """
        cd_0 = self.sn-1.24*pow(self.sn, 2)+13.7*pow(self.sn, 3)
        return (1 - 0.46 * cd_0)

    def cal_factor2(self, element_index, node_position):
        """
        A private function calculate the flow reduction factor according to Hui et al., (2020).\n
        :param alf: [float] Unit [rad]. Definition: the angle between normal vector of a net panel and the flow direction.
        :return: Flow reduction factor: [float] Unit [-].
        """
        p1 = node_position[self.elements[element_index][0]]
        p2 = node_position[self.elements[element_index][1]]
        p3 = node_position[self.elements[element_index][2]]
        a1 = p1 - p2
        a2 = p1 - p3
        unit_normal_vector = np.cross(
            a1, a2) / np.linalg.norm(np.cross(a1, a2))
        cosine_alpha = np.dot(unit_normal_vector,
                              self.direction) / np.linalg.norm(self.direction)

        reduction_factor = (cosine_alpha + 0.05 - 0.38 *
                            self.sn) / (cosine_alpha + 0.05)
        return max(0, reduction_factor)

    def update_reduction(self, position):
        for i in self.wake_element_indexes:
            if self.wake_type in ['factor']:
                self.reduction_factors[i] = min(1.0, float(self.wake_value))
            elif self.wake_type in ['loland']:
                self.reduction_factors[i] = float(self.cal_factor1())
            elif self.wake_type in ['hui']:
                self.reduction_factors[i] = float(
                    self.cal_factor2(i, position))
            else:
                print("the selected wake type " +
                      str(self.wake_type) + " is not supported.")
                exit()
        print('updated current velocity reduction')

    def cal_volume_tetrahedron(self, point1, point2, point3):
        origin = self.origin
        origin[2] = 0.0

        vector_a = point1-origin
        vector_b = point2-origin
        vector_c = point3-origin

        return abs(np.dot(vector_a, (np.cross(vector_b, vector_c))))/6.0

    def cal_cage_volume(self, position):
        """[summary]

        Args:
            position ([type]): [description]
        """

        total_volume = 0.0
        for element in self.elements:
            # loop based on the hydrodynamic elements
            if len(element) == 3:
                total_volume += self.cal_volume_tetrahedron(position[element[0]],
                                                            position[element[1]],
                                                            position[element[2]])
            elif len(element) == 4:
                for i in range(4):
                    node_3 = [int(k) for k in element]
                    # delete the i node to make the square to a triangle
                    node_3.pop(i)
                    total_volume += 0.5*self.cal_volume_tetrahedron(position[node_3[0]],
                                                                    position[node_3[1]],
                                                                    position[node_3[2]])
            else:
                print(str(element)+'is not supported by the netting element')
                exit()

        return total_volume
    
    
    def cal_area_triangle(self,point1, point2, point3):
        v_ab=point1-point2
        v_ac=point1-point3
        return 0.5*np.linalg.norm(np.cross(v_ab,v_ac))
        
    def cal_cage_net_area(self,position):
        """[summary]

        Args:
            position ([type]): [description]
        """
        total_net_area = 0.0
        for element in self.elements:
            # loop based on the hydrodynamic elements
            if len(element) == 3:
                total_net_area += self.cal_area_triangle(position[element[0]],
                                                         position[element[1]],
                                                         position[element[2]])
            elif len(element) == 4:
                for i in range(4):
                    node_3 = [int(k) for k in element]
                    # delete the i node to make the square to a triangle
                    node_3.pop(i)
                    total_net_area += 0.5*self.cal_area_triangle(position[node_3[0]],
                                                                 position[node_3[1]],
                                                                 position[node_3[2]])
            else:
                print(str(element)+'is not supported by the netting element')
                exit()

        return total_net_area


if __name__ == "__main__":
    pass
