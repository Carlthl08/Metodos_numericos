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


def trapezoidal_rule(f, a, b, n):
    """
    Aproxima el valor de una integral definida utilizando el método del trapecio.

    Args:
        f: La función a integrar.
        a: El límite inferior de integración.
        b: El límite superior de integración.
        n: El número de subintervalos.

    Returns:
        El valor aproximado de la integral definida.
    """
    h = (b - a) / n  # Tamaño de cada subintervalo
    x = sy.symbols('x')
    points = [a + i * h for i in range(n + 1)]  # Puntos de división entre los subintervalos
    integral_sum = f.subs(x, a) + f.subs(x, b)  # Suma de los valores en los extremos de la integral

    for i in range(1, n):
        integral_sum += 2 * f.subs(x, points[i])  # Suma de los valores en los puntos interiores de los subintervalos

    integral_approx = (h / 2) * integral_sum  # Aproximación de la integral

    return integral_approx


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




st.title(':blue[Integración Numérica]')
r'''
**Integración Numérica**

La integración numérica es una técnica utilizada para aproximar el valor de una integral definida de una función utilizando métodos computacionales. Estos métodos se basan en dividir el intervalo de integración en subintervalos más pequeños y utilizar fórmulas de aproximación para estimar la contribución de cada subintervalo a la integral total.

Existen varios métodos de integración numérica, siendo los más comunes el método del trapecio, la regla del punto medio y la regla de Simpson.

**Método del Trapecio**

El método del trapecio aproxima la integral de una función $f(x)$ en un intervalo $[a, b]$ dividiendo el intervalo en $n$ subintervalos de igual longitud y aproximando cada subintervalo por un trapecio. La fórmula general para el método del trapecio es:

$$
\int_{a}^{b} f(x) \, dx \approx \frac{h}{2} \left( f(a) + 2\sum_{i=1}^{n-1} f(x_i) + f(b) \right)
$$

Donde $h$ es la longitud de cada subintervalo y $x_i$ son los puntos de división entre los subintervalos.

**Regla del Punto Medio**

La regla del punto medio aproxima la integral de una función $f(x)$ en un intervalo $[a, b]$ dividiendo el intervalo en $n$ subintervalos de igual longitud y evaluando la función en el punto medio de cada subintervalo. La fórmula general para la regla del punto medio es:

$$
\int_{a}^{b} f(x) \, dx \approx h \sum_{i=1}^{n} f\left(\frac{{x_{i-1} + x_i}}{2}\right)
$$

Donde $h$ es la longitud de cada subintervalo y $x_{i-1}$ y $x_i$ son los puntos de división entre los subintervalos.

**Regla de Simpson**

La regla de Simpson aproxima la integral de una función $f(x)$ en un intervalo $[a,b]$ dividiendo el intervalo en $n$ subintervalos de igual longitud y utilizando polinomios de segundo grado para aproximar cada subintervalo. La fórmula general para la regla de Simpson es:

$$
\int_{a}^{b} f(x) \, dx \approx \frac{h}{3} \left( f(a) + 4\sum_{i=1}^{\frac{n}{2}} f(x_{2i-1}) + 2\sum_{i=1}^{\frac{n}{2}-1} f(x_{2i}) + f(b) \right)
$$

Donde $h$ es la longitud de cada subintervalo y $x_{2i-1}$ y $x_{2i}$ son los puntos de división entre los subintervalos.

Estos métodos de integración numérica proporcionan una aproximación del valor de una integral definida, y la precisión de la aproximación depende del número de subintervalos utilizados y la regularidad de la función integrada.
'''

r'''

**Ejemplo: Método del Trapecio**

Aproximemos el valor de la integral definida

$$
\int_{0}^{1} e^{-x^2} \, dx
$$

utilizando el método del trapecio. Dividiremos el intervalo \([0, 1]\) en 4 subintervalos de igual longitud.

Aplicamos la fórmula del método del trapecio:

$$
\int_{a}^{b} f(x) \, dx \approx \frac{h}{2} \left( f(a) + 2\sum_{i=1}^{n-1} f(x_i) + f(b) \right)
$$

Sustituimos los valores:

$$
\int_{0}^{1} e^{-x^2} \, dx \approx \frac{1}{8} \left( e^0 + 2\left( e^{-\left(\frac{1}{4}\right)^2} + e^{-\left(\frac{2}{4}\right)^2} + e^{-\left(\frac{3}{4}\right)^2} \right) + e^{-1^2} \right)
$$

Calculamos los valores de \(e^{-\left(\frac{i}{4}\right)^2}\) para \(i = 1, 2, 3\):

$$
e^{-\left(\frac{1}{4}\right)^2} \approx 0.939413
$$
$$
e^{-\left(\frac{2}{4}\right)^2} \approx 0.882497
$$
$$
e^{-\left(\frac{3}{4}\right)^2} \approx 0.771051
$$

Sustituimos los valores en la fórmula:

$$
\int_{0}^{1} e^{-x^2} \, dx \approx \frac{1}{8} \left( 1 + 2(0.939413 + 0.882497 + 0.771051) + e^{-1^2} \right)
$$

Calculamos $e^{-1^2}$:

$$
e^{-1^2} \approx 0.367879
$$

Sustituimos el valor en la fórmula:

$$
\int_{0}^{1} e^{-x^2} \, dx \approx \frac{1}{8} \left( 1 + 2(0.939413 + 0.882497 + 0.771051) + 0.367879 \right)
$$

Realizamos los cálculos:

$$
\int_{0}^{1} e^{-x^2} \, dx \approx 0.746824
$$

Por lo tanto, utilizando el método del trapecio con 4 subintervalos, aproximamos el valor de la integral definida $\int_{0}^{1} e^{-x^2} \, dx$ como \(0.746824\).

'''
st.subheader('Método')
xxs = st.text_input('Ingrese la función $f(x)$: ',value='exp(-x)')



fx = sy.parse_expr(xxs,transformations='all')
intstrr = ''





fxxs = st.text_input('Ingrese el intervalo $[a,b]$: ',value='[0,1]')


for t in fxxs:

    if t != '{' and t != '}' and t != '[' and t != ']' and t != '(' and t != ')' and t != ' ':
        intstrr = intstrr + t

interval = list(map(float,intstrr.split(',')))

symbs = list(fx.free_symbols)

dxx = ''
integ = sy.integrate(fx, symbs)
for i in symbs:
    dxx = dxx +'d'+ str(i)

nn = st.number_input('Numero de subintervalos: ',value=1 )

st.latex(r'\int '+sy.latex(fx)+dxx+' = '+sy.latex(integ)+' + C')

method = trapezoidal_rule(fx,interval[0],interval[1],nn)

st.latex(r'\int_{'+str(interval[0])+'}^{'+str(interval[1])+'}'+sy.latex(fx)+dxx+' = '+str(method))

if len(fx.free_symbols)<= 2:
    if len(fx.free_symbols) == 1:
        func = sy.lambdify(list(fx.free_symbols),fx)
        plo = gro.Figure()
        plo.add_trace(gro.Scatter(x=np.linspace(interval[0],interval[1],1000),y=func(np.linspace(interval[0],interval[1],1000)),fill='tozeroy'))
        plo.update_layout(title='Integral de f'+str(tuple(fx.free_symbols))+' en el intervalo '+str(interval))
        st.plotly_chart(plo)


    if  len(fx.free_symbols) == 2:
        func = sy.lambdify(list(fx.free_symbols),fx)
        plo = gro.Figure()
        ran = np.linspace(-10,10,100)
        su = [[func(ran[xs],ran[ys]) for xs in range (0,len(ran)) ] for ys in range(0,len(ran))]
        plo.add_trace(gro.Surface(z=su))
        st.plotly_chart(plo)
        p =plot3d(fx,show=False)
        pl = get_sympy_subplots(p)

        st.pyplot(pl)


#COLOCA TU METODO AQUI y PASA LA  FUNCION ALOJADA EN fx
