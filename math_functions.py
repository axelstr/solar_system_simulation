
class math_functions():
    def __init__(self):
        pass

    def norm(self,R):
        """Returns the norm of the 3D input array R."""
        norm_squared = 0
        for i in range(3):
            norm_squared += R[i]**2
        return norm_squared**(1/2)

    def vector_add(self, array_of_arrays):
        """Adds arrays indexwise to new array. All arrays in input must consist
        of floats or integers and be of equal lengths."""
        # assert correct input
        array_lengths = [len(array) for array in array_of_arrays]
        out_array_length = array_lengths[0]
        for length in array_lengths:
            assert length == out_array_length, 'Array lengths must be equal:\n'+str(array_of_arrays)

        # add vectors
        out = [0 for _ in range(out_array_length)]
        for array in array_of_arrays:
            for i in range(out_array_length):
                out[i] += array[i]
        return out

    def vector_scale(self, array, scalar):
        """Returns array with all numbers in array multiplied with scalar."""
        return [x*scalar for x in array]

    def vector_append(self, array_of_arrays):
        """Append arrays indexwise to new array. All arrays in input must consist
        of floats or integers and be of equal lengths."""
        # assert correct input
        array_lengths = [len(array) for array in array_of_arrays]
        out_array_length = array_lengths[0]
        for length in array_lengths:
            assert length == out_array_length, 'Array lengths must be equal:\n'+str(array_of_arrays)

        # add vectors
        out = [[] for _ in range(out_array_length)]
        for array in array_of_arrays:
            for i in range(out_array_length):
                out[i] = [out[i] ]
        return out
