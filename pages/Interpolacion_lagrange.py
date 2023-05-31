import streamlit as st
import pandas as pd
import numpy as np
import sympy as sy
from matplotlib import pyplot as plt
from sympy.plotting.plot import MatplotlibBackend, Plot
from sympy.plotting import plot3d,plot3d_parametric_line
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

def li(v, i):
    """
    The function takes a list of numbers and an index, and returns the Lagrange interpolating polynomial for the list of
    numbers with the index'th number removed

    :param v: the list of x values
    :param i: the index of the x value you want to interpolate
    :return: the Lagrange interpolating polynomial for the given data points.
    """
    x = sy.symbols('x')

    s = 1
    st = ''
    for k in range(0,len(v)):
        if k != i:
            st = st + '((' + str(x) + '-'+ str(v[k])+')/('+str(v[i])+'-'+str(v[k])+'))'
            s = s*((x-v[k])/(v[i]-v[k]))

    return s

def Lagrange(v,fx):
    """
    It takes in a list of x values and a list of y values, and returns the Lagrange polynomial that interpolates those
    points

    :param v: list of x values
    :param fx: The function you want to interpolate
    :return: the Lagrange polynomial.
    """
    #print(v)
    #print(fx)
    lis = []
    for i in range(0,len(v)):
        lis.append(li(v,i))

    sums = 0

    for k in range(0,len(v)):
        sums = sums+(fx[k]*lis[k])

    #print(sums)

    sy.simplify(sums)

    sy.pprint(sums)

    p1 = sy.plot(sums,show=False)
    p2 = get_sympy_subplots(p1)
    p2.plot(v,fx,"o")
    #p2.show()
    return sy.expand(sums), p2,lis

st.title(':blue[Interpolación de Lagrange]')

r'''

**Interpolación de Lagrange**

La interpolación de Lagrange es un método utilizado para aproximar una función desconocida a través de un conjunto de puntos de datos conocidos. Este método construye un polinomio que pasa exactamente por todos los puntos de datos, lo que permite estimar el valor de la función en cualquier otro punto dentro del rango de datos.

La idea principal detrás de la interpolación de Lagrange es utilizar una combinación lineal de polinomios de Lagrange para construir el polinomio interpolante. Cada polinomio de Lagrange está asociado con un punto de datos y está diseñado para ser igual a 1 en ese punto y igual a 0 en todos los demás puntos.

La fórmula general del polinomio interpolante de Lagrange es:

$$
P(x) = \sum_{i=0}^{n} y_i \cdot L_i(x)
$$

Donde:
- \(P(x)\) es el polinomio interpolante.
- $n$ es el número de puntos de datos.
- $x$ es el punto en el que se desea estimar el valor de la función.
- $y_i$ son los valores de la función en los puntos de datos.
- $L_i(x)$ son los polinomios de Lagrange asociados con los puntos de datos, definidos por:

$
L_i(x) = \prod_{j=0, j \neq i}^{n} \frac{x - x_j}{x_i - x_j}
$

La fórmula del polinomio interpolante de Lagrange se obtiene sumando los productos de los valores de la función en los puntos de datos $y_i$ y los polinomios de Lagrange $L_i(x)$ evaluados en el punto \(x\).

La interpolación de Lagrange es ampliamente utilizada en diversas áreas, como análisis numérico y gráficos por computadora, y proporciona una aproximación suave de la función original a través de los puntos de datos.

'''

st.subheader(':blue[Ejemplo]')

r'''

Supongamos que tenemos los siguientes puntos de datos:

```
x = [1, 2, 4]
y = [2, 3, 1]
```

Queremos encontrar el polinomio interpolante de Lagrange que pase por estos puntos y utilizarlo para estimar el valor de la función en otro punto, por ejemplo, en \(x = 3\).

**Paso 1:** Definir los puntos de datos.

```
x = [1, 2, 4]
y = [2, 3, 1]
```

**Paso 2:** Calcular los polinomios de Lagrange.

Para \(x = 1\):
```
L0(x) = (x - 2)(x - 4) / (1 - 2)(1 - 4) = -x^2 + 6x - 8
```

Para \(x = 2\):
```
L1(x) = (x - 1)(x - 4) / (2 - 1)(2 - 4) = 3x^2 - 14x + 12
```

Para \(x = 4\):
```
L2(x) = (x - 1)(x - 2) / (4 - 1)(4 - 2) = -2x^2 + 9x - 6
```

**Paso 3:** Construir el polinomio interpolante.

```
P(x) = 2 * L0(x) + 3 * L1(x) + 1 * L2(x)
```

Simplificando, obtenemos:
```
P(x) = -x^2 + 6x - 8 + 3x^2 - 14x + 12 - 2x^2 + 9x - 6
P(x) = x^2 + x - 2
```

Por lo tanto, el polinomio interpolante de Lagrange es \(P(x) = x^2 + x - 2\).

**Paso 4:** Estimar el valor de la función en \(x = 3\).

Sustituyendo \(x = 3\) en el polinomio interpolante:
```
P(3) = 3^2 + 3 - 2 = 12
```

Por lo tanto, el valor estimado de la función en \(x = 3\) utilizando el polinomio interpolante de Lagrange es \(P(3) = 12\).

'''



st.subheader('Método')

filess = None
if filess != None:
    fi = pd.read_csv(filess)
    st.write('Los datos a interpolar son: ')
    st.write(fi)
    x = list(fi['x'])
    fx = list(fi['y'])
else:
    xxs = st.text_input('Ingrese los valores de $x_k$: ',value='{1,2,3,4}')

    xsstr = ''


    for i in xxs:

        if i != '{' and i != '}' and i != '[' and i != ']' and i != '(' and i != ')' and i != ' ':
            xsstr = xsstr + i

    fxxs = st.text_input('Ingrese los valores de $f(x_k)$: ',value='{-1,3,4,5}')

    x = list(map(float,xsstr.split(',')))
    intstrr = ''




    for t in fxxs:

        if t != '{' and t != '}' and t != '[' and t != ']' and t != '(' and t != ')' and t != ' ':
            intstrr = intstrr + t

    fx = list(map(float,intstrr.split(',')))


#st.write(x)
#st.write(fx)
#data = [x,fx]
#st.write(data)


method = Lagrange(x,fx)

st.write('_Los polinomios fundamentales de Lagrange estan dados por:_')
lli = r'''l_i(x) = \begin{cases}'''
for t in range(0,len(method[2])):
    lli = lli +'l_'+str(t)+r'='+sy.latex(sy.expand(method[2][t]))+r'\\'
lli = lli + r'\end{cases}'
st.latex(lli)
st.write('_El polinomio de Interpolacion está dado por:_')
st.latex(r'p_n(x) = \sum_{i=0}^{n} l_i(x)f(x_i)')
st.latex('p_n(x) =' + sy.latex(method[0]))

func = sy.lambdify(sy.symbols('x'),method[0])
funcdata = pd.DataFrame(dict(x=np.linspace(-10,10,1000),y=func(np.linspace(-10,10,1000))))

plo = gro.Figure()

plo.add_trace(gro.Scatter(x=np.linspace(-10,10,1000),y=func(np.linspace(-10,10,1000)),name='Interpolación'))
plo.add_trace(gro.Scatter(x=x,y=fx, marker_color='rgba(152, 0, 0, .8)',name='Datos'))
#plo.add_hline(y=0)
#plo.add_vline(x=0)
plo.update_layout(title='Grafica de la Interpolación')
st.plotly_chart(plo)


