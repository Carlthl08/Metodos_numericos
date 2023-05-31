import streamlit as st
import pandas as pd
import numpy as np
import sympy as sy
from matplotlib import pyplot as plt
from sympy.plotting.plot import MatplotlibBackend, Plot
from sympy.plotting import plot3d,plot3d_parametric_line
import plotly as ply
import plotly.express as ex
import plotly.graph_objects as gro
from plotly.subplots import make_subplots





def get_sympy_subplots(plot:Plot):
    """
    It takes a plot object and returns a matplotlib figure object

    :param plot: The plot object to be rendered
    :type plot: Plot
    :return: A matplotlib figure object.
    """
    backend = MatplotlibBackend(plot)

    backend.process_series()
    backend.fig.tight_layout()
    return backend.plt


def forward_difference(f, x, h):
    """
    Aproxima la derivada de una función en un punto utilizando el método de diferencias hacia adelante.

    Args:
        f: La función para la cual se desea aproximar la derivada.
        x: El punto en el cual se quiere aproximar la derivada.
        h: El tamaño del incremento en x.

    Returns:
        La aproximación de la derivada en el punto dado.
    """
    return (f.subs(x, x + h) - f.subs(x, x)) / h



st.title(':blue[Derivación Numérica]')

r'''
**Derivación Numérica**

La derivación numérica es una técnica utilizada para aproximar la derivada de una función en un punto dado utilizando métodos computacionales. En lugar de calcular la derivada analíticamente, se utiliza una aproximación basada en diferencias finitas.

Existen varios métodos de derivación numérica, siendo los más comunes el método de diferencias hacia adelante, el método de diferencias hacia atrás y el método de diferencias centradas.

**Método de Diferencias Hacia Adelante**

El método de diferencias hacia adelante estima la derivada de una función $f(x)$ en un punto $x_0$ utilizando la siguiente fórmula:

$$
f'(x_0) \approx \frac{{f(x_0 + h) - f(x_0)}}{{h}}
$$

Donde $h$ es un pequeño incremento en $x$.

**Método de Diferencias Hacia Atrás**

El método de diferencias hacia atrás estima la derivada de una función \(f(x)\) en un punto \(x_0\) utilizando la siguiente fórmula:

$$
f'(x_0) \approx \frac{{f(x_0) - f(x_0 - h)}}{{h}}
$$

**Método de Diferencias Centradas**

El método de diferencias centradas estima la derivada de una función \(f(x)\) en un punto \(x_0\) utilizando la siguiente fórmula:

$$
f'(x_0) \approx \frac{{f(x_0 + h) - f(x_0 - h)}}{{2h}}
$$

Este método proporciona una mayor precisión en comparación con los métodos de diferencias hacia adelante y hacia atrás.

Es importante tener en cuenta que la elección del tamaño del incremento \(h\) es crucial para obtener una buena aproximación de la derivada. Si \(h\) es demasiado pequeño, se pueden introducir errores de redondeo y de truncamiento. Si \(h\) es demasiado grande, la aproximación puede ser imprecisa.

La derivación numérica es ampliamente utilizada en situaciones donde no es posible obtener la derivada analítica de una función o cuando se trabaja con datos discretos.


'''
st.subheader(':blue[Ejemplo]')

r'''

Supongamos que queremos aproximar la derivada de la función $f(x) = x^2$ en el punto \(x = 2\) utilizando el método de diferencias hacia adelante. Tomaremos un incremento \(h = 0.1\).

Aplicamos la fórmula del método de diferencias hacia adelante:

$$
f'(x_0) \approx \frac{{f(x_0 + h) - f(x_0)}}{{h}}
$$

Sustituimos los valores:

$$
f'(2) \approx \frac{{f(2 + 0.1) - f(2)}}{{0.1}}
$$

Calculamos $f(2 + 0.1)$ y $f(2)$:

$$
f(2 + 0.1) = (2 + 0.1)^2 = 4.41
$$
$$
f(2) = 2^2 = 4
$$

Sustituimos los valores en la fórmula:

$$
f'(2) \approx \frac{{4.41 - 4}}{{0.1}} = 4.1
$$

Por lo tanto, utilizando el método de diferencias hacia adelante con un incremento \(h = 0.1\), aproximamos la derivada de la función \(f(x) = x^2\) en el punto \(x = 2\) como \(f'(2) \approx 4.1\).


'''
st.subheader('Método')
xxs = st.text_input('Ingrese la función $f(x)$: ',value='(x - 1)**2')



fx = sy.parse_expr(xxs,transformations='all')
intstrr = ''

xx = st.number_input('Ingrese el punto a aproximar:',value=1)

hh = st.number_input('Ingrese el incremento h:',value=0.1)


method = forward_difference(fx,xx,hh)

st.latex(r'\frac{\partial f}{\partial x}'+sy.latex(fx)+' = '+str(method) )


if len(fx.free_symbols)<= 2:
    if len(fx.free_symbols) == 1:
        func = sy.lambdify(list(fx.free_symbols),fx)
        plo = gro.Figure()
        plo.add_trace(gro.Scatter(x=np.linspace(-10,10,1000),y=func(np.linspace(-10,10,1000))))
        st.plotly_chart(plo)
        p =sy.plot(fx,show=False)
        pl = get_sympy_subplots(p)



    if  len(fx.free_symbols) == 2:
        func = sy.lambdify(list(fx.free_symbols),fx)
        plo = gro.Figure()
        ran = np.linspace(-10,10,100)
        su = [[func(ran[xs],ran[ys]) for xs in range (0,len(ran)) ] for ys in range(0,len(ran))]
        plo.add_trace(gro.Surface(z=su))
        st.plotly_chart(plo)
        p =plot3d(fx,show=False)
        pl = get_sympy_subplots(p)




#COLOCA TU METODO AQUI y PASA LA  FUNCION ALOJADA EN fx
