# %%
import numpy

class Struct():
    '''
    A way to create a structure in python
    '''
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)  

    def __repr__(self):
        str = 'Struct fields:\n'
        for key, value in self.__dict__.items():
            if isinstance(value, numpy.ndarray):
                new_var = value.shape
                value = ['shape: ', new_var]
            str += '- {0}: {1}\n'.format(key, value)
        return str

if __name__ == '__main__':
    st = Struct()
    st.a = 67
    st.b = [1,2,3]
    st.c = 'Hello world!'
    st.D = numpy.random.randn(20,34,5)

    print(st)
# %%
