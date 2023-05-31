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



def continium_minimun_quads_aprox(fx,interval,symb,degree):
    """
    Given a function, an interval, a symbol and a degree, it returns the polynomial that best approximates the function in
    the given interval

    :param fx: The function to be approximated
    :param interval: the interval in which the function is defined
    :param symb: the symbol that will be used to represent the variable in the function
    :param degree: The degree of the polynomial
    :return: The function that is the best aproximation of the given function in the given interval.
    """

    m = []


    for i in range(0,degree+1):
        aux = []
        for j in range(0,degree+1):
            aux.append(sy.integrate((symb**i)*(symb**j),(symb,interval[0],interval[1])))
        m.append(aux)

    #pprint(Matrix(m))


    b = []

    for i in range(0,degree+1):
        b.append(sy.integrate((symb**i)*fx,(symb,interval[0],interval[1])))

    #pprint(Matrix(b))

    sol = sy.Matrix(m).inv() * sy.Matrix(b)

    expr = 0

    for i in range(0,degree+1):
        expr = expr + (sol[i]*symb**i)

    #pprint(expr)


    p = sy.plot(fx,(symb,interval[0],interval[1]),show=False)
    p.append(sy.plot(expr,(symb,interval[0],interval[1]),show=False)[0])

    #p.show()


    return sy.expand(expr),get_sympy_subplots(p)





t = r'''
La aproximación continua de mínimos cuadrados es un método utilizado para encontrar una función continua que se ajuste de manera óptima a un conjunto de datos. A diferencia de la aproximación discreta, en la cual se busca ajustar una función a puntos de datos específicos, la aproximación continua se enfoca en ajustar una función a una distribución continua de datos.

El proceso básico de aproximación continua de mínimos cuadrados implica los siguientes pasos:

1. Definir el problema: Tener un conjunto de datos continuos que representan una distribución de valores en un rango determinado. Se busca encontrar una función continua que se ajuste a esta distribución de datos.

2. Seleccionar el tipo de función: Decidir qué tipo de función se utilizará para aproximarse a los datos continuos. Al igual que en la aproximación discreta, se pueden utilizar funciones polinómicas, exponenciales, logarítmicas, trigonométricas u otras funciones, dependiendo del problema y los datos disponibles.

3. Formular el problema de mínimos cuadrados: Definir una función de error que representa la diferencia entre los valores predichos por la función continua y los valores reales de los datos. La función de error suele ser la integral del cuadrado de las diferencias (residuos) entre los valores predichos y los valores reales a lo largo de todo el rango de datos.

4. Minimizar el error: Encontrar los coeficientes o parámetros de la función continua que minimizan la función de error. Esto se puede lograr mediante técnicas matemáticas como la optimización numérica, donde se busca el mínimo de la función de error utilizando métodos como el cálculo de variaciones, el método de los mínimos cuadrados ponderados o técnicas de ajuste de curvas.

5. Evaluar el ajuste: Una vez que se encuentran los coeficientes de la función continua que minimizan el error, se puede evaluar qué tan bien se ajusta la función a los datos continuos. Esto se puede hacer calculando métricas de ajuste como el coeficiente de determinación (R^2) o visualizando los datos junto con la función ajustada en un gráfico.

'''



e = r'''


Supongamos que tenemos los siguientes puntos de datos:

```
x = [1, 2, 3, 4, 5]
y = [2, 3, 5, 6, 8]
```

Queremos encontrar una función continua que se ajuste de manera óptima a estos datos.

**Definición del problema:** Tenemos un conjunto de datos continuos representados por los puntos $(x, y)$. Queremos encontrar una función continua que se ajuste a esta distribución de datos.

**Selección del tipo de función:** Utilizaremos una función polinómica de grado 2, es decir, una función de la forma $f(x) = ax^2 + bx + c$.

**Formulación del problema de mínimos cuadrados:** La función de error será la integral del cuadrado de las diferencias entre los valores predichos por la función polinómica y los valores reales de los datos a lo largo de todo el rango de datos.

**Minimización del error:** Encontraremos los coeficientes $a$, $b$ y $c$ que minimizan la función de error mediante técnicas de optimización numérica.

**Evaluación del ajuste:** Calcularemos métricas de ajuste, como el coeficiente de determinación ($R^2$), para evaluar qué tan bien se ajusta la función polinómica a los datos continuos.

'''


st.title(':blue[Aproximación continua de minimos cuadrados]')

st.write(t)

st.subheader(':blue[Ejemplo]')
st.write(e)

st.subheader('Método')
xxs = st.text_input('Ingrese la función $f(x)$: ',value='cos(pi*x)')



fx = sy.parse_expr(xxs,transformations='all')
intstrr = ''


fxxs = st.text_input('Ingrese el intervalo $[a,b]$: ',value='[-1,1]')


for t in fxxs:

    if t != '{' and t != '}' and t != '[' and t != ']' and t != '(' and t != ')' and t != ' ':
        intstrr = intstrr + t

interval = list(map(float,intstrr.split(',')))



degree = st.slider('Grado del polinomio de aproximación: ',1,10,value=2)

method = continium_minimun_quads_aprox(fx,interval,list(fx.free_symbols)[0],int(degree))

st.write('El polinomio esta dado por:')
st.latex('P_{'+str(degree)+'}(x)='+sy.latex(method[0]))




plo = gro.Figure()
func = sy.lambdify(list(fx.free_symbols)[0],fx)
aproxfunc = sy.lambdify(list(fx.free_symbols)[0],method[0])
plo.add_trace(gro.Scatter(x = np.linspace(interval[0],interval[1],1000),y=func(np.linspace(interval[0],interval[1],1000)),name='Función', marker_color='rgba(152, 0, 0, .8)'))
plo.add_trace(gro.Scatter(x=np.linspace(interval[0],interval[1],1000),y=aproxfunc(np.linspace(interval[0],interval[1],1000)),name='Aproximación',fill='tonexty'))
st.plotly_chart(plo)


