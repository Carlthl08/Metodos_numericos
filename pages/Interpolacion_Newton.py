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

def diff_div(v,fx,order):
    """
    > The function takes in a list of values, a list of function values, and an order, and returns a list of divided
    differences

    :param v: the list of x values
    :param fx: the function you want to differentiate
    :param order: the order of the derivative you want to take
    :return: the difference quotient of the function f(x)
    """

    m = []

    for i in range(0,len(fx)):
        #print(fx[i])
        if i + 1 < len(fx) and i +order < len(v):
            #print(v[i+1],v[i],"/",fx[i+order]," ",fx[i])
            m.append((fx[i+1]-fx[i])/(v[i+order]-v[i]))
    return m

def divided_diff(fx,v):
    """
    The function takes in a list of x values and a list of f(x) values, and returns a list of lists of divided differences

    :param fx: the function to be interpolated
    :param v: list of x values
    :return: The divided difference table is being returned.
    """
    x = v
    nfx = fx
    m = []
    for i in range(0,len(v)-1):
        nx = diff_div(v,nfx,i+1)
        #print(nx)
        m.append(nx)
        nfx = nx

    #print(m)
    return m

def table_divided_dif(x,fx,diff):
    table = []
    for i in range(0,len(fx)):
        table.append([x[i],fx[i]])
        table.append([' ',' '])

    for i in range(1,len(table)):
        for j in range(0,len(diff)):
            for k in range(0,len(diff[j])):
                table[i].append(diff[j][k])


    return table


def print_divdiff(v):
    """
    The function `print_divdiff` prints a divided difference table given a list of values.

    :param v: The input parameter "v" is expected to be a list of lists, where the first list contains the x values and the
    second list contains the corresponding f(x) values. The remaining lists contain the divided differences of the function
    """
    table = []
    s = ''
    for i in range(0,len(v[0])):
        table.append([str(v[0][i] ),str(v[1][i]) ])
        table.append([])
    #print(table)
    for k in range(2,len(v)):
        aux = k-1
        for j in range(0,len(v[k])):
            for t in range(0,k):
                table[aux+j].append(' ')
            table[aux+j].append( str(v[k][j]))
            aux = aux + 1
    m = 0
    for i in table:
        if len(i) > m:
            m = len(i)

    for t in table:
        while len(t) != m:
            t.append('')

    return table
def Newton_interpolation(fx,v):
    """
    It takes in a list of x values and a list of f(x) values, and returns a polynomial that interpolates the points

    :param fx: a list of the function values
    :param v: list of x values
    :return: The function is being returned.
    """
    diff = divided_diff(fx,v)
    x = sy.symbols('x')

    expr = v[0]
    expr_reg = v[0]

    for i in range(0,len(diff)):
        s = diff[i][0]
        p = 1
        for k in range(0,len(v)):

            p = p*(x-v[k])
            #print(p, "p",k)
            if k == i:
                break
        s = s * p
        expr = expr + s

    for i in range(0,len(diff)):
        s = diff[i][-1]
        p = 1
        for k in range(len(v)-1,0,-1):

            p = p*(x-v[k])
            #print(p, "p",k)
            if k == i:
                break
        s = s * p
        expr_reg = expr_reg + s

    #pprint(expr)

    p = sy.plot(expr,(x,-10,10),show=False)
    #p.append(sy.plot(expr_reg,(x,-10,10),show=False)[0])
    p2 = get_sympy_subplots(p)
    p2.plot(v,fx,"o")



    #p2.show()
    difdiv = [v,fx]
    for i in diff:
        difdiv.append(i)

    return sy.expand(expr),p2, print_divdiff(difdiv),sy.expand(expr_reg)





st.title(':blue[Interpolación de Newton]')

r'''
**Interpolación de Newton**

La interpolación de Newton es un método utilizado para aproximar una función desconocida a través de un conjunto de puntos de datos conocidos. Este método utiliza polinomios de Newton para construir el polinomio interpolante, lo que permite estimar el valor de la función en cualquier otro punto dentro del rango de datos.

La idea principal detrás de la interpolación de Newton es utilizar diferencias divididas para construir el polinomio interpolante. Las diferencias divididas son las diferencias entre los valores de la función en los puntos de datos y se utilizan para calcular los coeficientes del polinomio interpolante.

El polinomio interpolante de Newton se define de la siguiente manera:

$$
P(x) = f[x_0] + f[x_0,x_1](x-x_0) + f[x_0,x_1,x_2](x-x_0)(x-x_1) + \ldots + f[x_0,x_1,\ldots,x_n](x-x_0)(x-x_1)\ldots(x-x_{n-1})
$$

Donde:
- $P(x)$ es el polinomio interpolante.
- $f[x_0]$, $f[x_0,x_1]$, $f[x_0,x_1,x_2]$, etc., son las diferencias divididas de la función.
- $x_0, x_1, x_2, \ldots, x_n$ son los puntos de datos conocidos.

Las diferencias divididas se calculan de la siguiente manera:

$$
f[x_i] = y_i
$$
$$
f[x_i,x_{i+1}] = \frac{f[x_{i+1}] - f[x_i]}{x_{i+1} - x_i}
$$
$$
f[x_i,x_{i+1},x_{i+2}] = \frac{f[x_{i+1},x_{i+2}] - f[x_i,x_{i+1}]}{x_{i+2} - x_i}
$$
y así sucesivamente.

'''

st.subheader('Ejemplo')
r'''

Supongamos que tenemos los siguientes puntos de datos:

```
x = [1, 2, 4]
y = [2, 3, 1]
```

Queremos encontrar el polinomio interpolante de Newton que pase por estos puntos y utilizarlo para estimar el valor de la función en otro punto, por ejemplo, en \(x = 3\).

**Paso 1:** Definir los puntos de datos.

```
x = [1, 2, 4]
y = [2, 3, 1]
```

**Paso 2:** Calcular las diferencias divididas.

Primero, calculamos las diferencias divididas de primer orden:
```
f[x0, x1] = (y1 - y0) / (x1 - x0) = (3 - 2) / (2 - 1) = 1
f[x1, x2] = (y2 - y1) / (x2 - x1) = (1 - 3) / (4 - 2) = -1
```

Luego, calculamos las diferencias divididas de segundo orden:
```
f[x0, x1, x2] = (f[x1, x2] - f[x0, x1]) / (x2 - x0) = (-1 - 1) / (4 - 1) = -0.67
```

**Paso 3:** Construir el polinomio interpolante.

Utilizamos las diferencias divididas para construir el polinomio interpolante de Newton:

```
P(x) = y0 + f[x0, x1](x - x0) + f[x0, x1, x2](x - x0)(x - x1)
     = 2 + 1(x - 1) + (-0.67)(x - 1)(x - 2)
```

Simplificando, obtenemos:
```
P(x) = 2 + (x - 1) - 0.67(x - 1)(x - 2)
     = 2 + x - 1 - 0.67(x^2 - 3x + 2)
     = -0.67x^2 + 3.34x - 3.34
```

Por lo tanto, el polinomio interpolante de Newton es \(P(x) = -0.67x^2 + 3.34x - 3.34\).

**Paso 4:** Estimar el valor de la función en \(x = 3\).

Sustituyendo \(x = 3\) en el polinomio interpolante:
```
P(3) = -0.67(3)^2 + 3.34(3) - 3.34
     = -6.03 + 10.02 - 3.34
     = 0.65
```

Por lo tanto, el valor estimado de la función en \(x = 3\) utilizando el polinomio interpolante de Newton es \(P(3) = 0.65\).
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

    fxxs = st.text_input('Ingrese los valores de $f(x_k)$: ',value='{1,2.14557,3.141592,4}')

    x = list(map(float,xsstr.split(',')))
    intstrr = ''




    for t in fxxs:

        if t != '{' and t != '}' and t != '[' and t != ']' and t != '(' and t != ')' and t != ' ':
            intstrr = intstrr + t

    fx = list(map(float,intstrr.split(',')))


#st.write(x)
#st.write(fx)



method = Newton_interpolation(fx,x)
try:
    st.write('Diferencias Divididas')

    st.table(method[2])
except:
    st.error('This is an error', icon="🚨")
st.write('_El polinomio de Interpolacion está dado por:_')
st.write('Progresivo: ')
st.latex('p_n(x) = ' +sy.latex(method[0]))
st.write('Regresivo:')
st.latex('p_n(x) = ' +sy.latex(method[-1]))



func = sy.lambdify(sy.symbols('x'),method[0])
funcdata = pd.DataFrame(dict(x=np.linspace(-10,10,1000),y=func(np.linspace(-10,10,1000))))
func2 = sy.lambdify(sy.symbols('x'),method[-1])
plo = gro.Figure()

plo.add_trace(gro.Scatter(x=np.linspace(-10,10,1000),y=func(np.linspace(-10,10,1000)),name='Interpolación'))
#plo.add_trace(gro.Scatter(x=np.linspace(-10,10,1000),y=func2(np.linspace(-10,10,1000)),name='Interpolación2'))
plo.add_trace(gro.Scatter(x=x,y=fx, marker_color='rgba(152, 0, 0, .8)',name='Datos'))

#plo.add_hline(y=0)
#plo.add_vline(x=0)
plo.update_layout(title='Grafica de la Interpolación')
st.plotly_chart(plo)

