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


def discrete_minimun_quads_aprox(xs,y,functionss,symbs):
    """
    Given a set of points $(x_i,y_i)$ and a set of functions (x)$, it finds the linear combination of the functions that
    best fits the points

    :param xs: list of x values
    :param y: the y values of the data points
    :param functionss: a list of functions that will be used to approximate the data
    :param symbs: the symbol that you want to use for the function
    :return: The expression of the function that best fits the data.
    """

    m = []


    for i in range(0,len(xs)):
        aux = []

        for j in range(0,len(functionss)):

            aux.append(functionss[j])
        m.append(aux)


    #pprint(Matrix(m))

    mev = []
    for i in range(0,len(m)):
        aux = []

        for j in range(0,len(m[0])):
            if len(m[i][j].free_symbols) > 0:
                aux.append(m[i][j].subs(symbs,xs[i]))
            else:
                aux.append(m[i][j])
        mev.append(aux)

    #pprint(Matrix(mev))

    mevT = sy.Matrix(mev).transpose()
    #pprint(mevT)

    a = mevT*sy.Matrix(mev)

    #pprint(a)

    b = mevT*sy.Matrix(y)

    #pprint(b)

    ainv = a.inv()

    xsol = ainv*b

    #pprint(xsol)


    expr = xsol[0]+xsol[1]*symbs


    p = sy.plot(expr,show=False)
    p2 = get_sympy_subplots(p)

    p2.plot(xs,y,"o")
    #p2.show()
    return sy.expand(expr),p2

t = r'''
La aproximación discreta de mínimos cuadrados es un método utilizado para encontrar una función que se ajuste de manera óptima a un conjunto de puntos de datos. En lugar de intentar encontrar una función exacta que pase por todos los puntos, el enfoque de mínimos cuadrados busca minimizar la suma de los cuadrados de las diferencias entre los valores predichos por la función y los valores reales de los datos.

El proceso básico de aproximación discreta de mínimos cuadrados implica los siguientes pasos:

1. Definir el problema: Tener un conjunto de datos discretos que consiste en pares de coordenadas (x, y) o datos (x_i, y_i) donde se desea encontrar una función que se ajuste a estos puntos.

2. Seleccionar el tipo de función: Decidir qué tipo de función se utilizará para aproximarse a los datos. Esto puede incluir funciones polinómicas, exponenciales, logarítmicas, trigonométricas u otras funciones dependiendo del problema y los datos disponibles.

3. Formular el problema de mínimos cuadrados: Definir una función de error que representa la diferencia entre los valores predichos por la función y los valores reales de los datos. La función de error suele ser la suma de los cuadrados de las diferencias (residuos) entre los valores predichos y los valores reales.

4. Minimizar el error: Encontrar los coeficientes o parámetros de la función que minimizan la función de error. Esto se puede lograr mediante técnicas matemáticas como la optimización numérica, donde se busca el mínimo de la función de error mediante métodos como el método de Newton, el gradiente descendente o la descomposición QR.

5. Evaluar el ajuste: Una vez que se encuentran los coeficientes de la función que minimizan el error, se puede evaluar qué tan bien se ajusta la función a los datos. Esto se puede hacer calculando métricas de ajuste como el coeficiente de determinación (R^2) o visualizando los datos junto con la función ajustada en un gráfico.

'''


e = r'''
Supongamos que tenemos los siguientes puntos de datos:

```
x = [1, 2, 3, 4, 5]
y = [2, 3, 5, 6, 8]
```

Queremos encontrar una función que se ajuste de manera óptima a estos puntos de datos.

**Definición del problema:** Tenemos un conjunto de puntos de datos representados por los pares $(x, y)$. Queremos encontrar una función que se ajuste a estos puntos.

**Selección del tipo de función:** Podemos utilizar una función polinómica de grado 1, es decir, una función de la forma $f(x) = ax + b$.

**Formulación del problema de mínimos cuadrados:** La función de error será la suma de los cuadrados de las diferencias entre los valores predichos por la función y los valores reales de los datos.

**Minimización del error:** Encontraremos los coeficientes $a$ y $b$ que minimizan la función de error mediante técnicas de optimización numérica.

**Evaluación del ajuste:** Podemos calcular métricas de ajuste, como el coeficiente de determinación ($R^2$), para evaluar qué tan bien se ajusta la función a los puntos de datos.

'''
st.title(':blue[Aproximación discreta de minimos cuadrados]')

st.write(t)

st.subheader(':blue[Ejemplo]')
st.write(e)


xxs = st.text_input('Ingrese los valores de $x_n$: ',value='{-1,1,3,4}')

xsstr = ''


for i in xxs:

    if i != '{' and i != '}' and i != '[' and i != ']' and i != '(' and i != ')' and i != ' ':
        xsstr = xsstr + i

fxxs = st.text_input('Ingrese los valores de $f(x_n)$: ',value='{6,1,11,3}')

x = list(map(float,xsstr.split(',')))
intstrr = ''




for t in fxxs:

    if t != '{' and t != '}' and t != '[' and t != ']' and t != '(' and t != ')' and t != ' ':
        intstrr = intstrr + t

fx = list(map(float,intstrr.split(',')))



funx = st.text_input('Ingrese las funciones $f_k(x)$ a considerar:',value='{1,x**2}')
funcstr = ''

for t in funx:

    if t != '{' and t != '}' and t != '[' and t != ']' and t != '(' and t != ')' and t != ' ':
        funcstr = funcstr +t

#st.write(funcstr)
funcs = []
for i in funcstr.split(','):
    funcs.append(sy.parse_expr(i,transformations='all'))

sym = list(funcs[0].free_symbols)

l = 0
while l < len(funcs):
    if len(funcs[l].free_symbols) != 0:
        sym = list(funcs[l].free_symbols)
        break
    l += 1

#st.write(str(sym))
method = discrete_minimun_quads_aprox(x,fx,funcs,sym[0])

st.write('La combinacion lineal que mejor se ajusta a los datos es:')
st.latex('f(x)='+sy.latex(method[0]))


func = sy.lambdify(sym[0],method[0])

plo = gro.Figure()
plo.add_trace(gro.Scatter(x=x,y=fx,name='Datos'))
plo.add_trace(gro.Scatter(x=np.linspace(min(x)-10,max(x)+10,1000),y=func(np.linspace(min(x)-10,max(x)+10,1000)),name='Aproximación'))

st.plotly_chart(plo)
