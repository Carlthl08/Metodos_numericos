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


def spline_natural(fx,v):
    """
    It takes a list of x values and a list of y values, and returns a list of sympy expressions that represent the cubic
    spline interpolation of the data

    :param fx: list of f(x) values
    :param v: list of x values
    """

    inter = []
    fxinter = []

    hi =[]
    for i in range(0,len(v)-1):
        inter.append((v[i],v[i+1]))
        fxinter.append((fx[i],fx[i+1]))

    #print(inter)
    for i in range(0,len(inter)):
        hi.append(inter[i][1]-inter[i][0])

    m = np.zeros(len(v)**2).reshape(len(fx),len(fx))
    #print(hi)
    #print(m)
    for i in range(0,len(v)):
        for j in range(0,len(v)):
            if (i == j and i == 0 and j == 0) or (j == i and i == len(v)-1 and j == len(v)-1):
                m[i][j] = 1
                continue
            else:
                if (i == j):
                    m[i][j] = 2*(hi[i-1]+hi[i])
                    m[i][j-1] = hi[i-1]
                    m[i][j+1] = hi[i]

    b = np.zeros(len(v))

    for i in range(1,len(v)-1):
        b[i] = ((1/hi[i])*(fx[i+1]-fx[i]))-((1/hi[i-1])*(fx[i]-fx[i-1]))

    #print(m)
    #pprint(Matrix(b.transpose()))

    c = (sy.Matrix(m).inv())*sy.Matrix(b.transpose())
    #pprint(c)
    b = []

    for i in range(0,len(hi)):
        b.append(((fx[i+1]-fx[i])/hi[i])-((((2*c[i])+c[i+1])*hi[i])/3))

    #pprint(Matrix(b))

    d = []

    for i in range(0,len(hi)):
        d.append((c[i+1]-c[i])/(3*hi[i]))

    #pprint(Matrix(d))


    x = sy.symbols('x')
    spl = []
    for i in range(0,len(inter)):
        spl.append(fx[i]+ (b[i]*(x-v[i]))+(c[i]*((x-v[i])**2)) + (d[i]*((x-v[i])**3)))

    #pprint(Matrix(spl))



    p = sy.plot(spl[0], (x,inter[0][0],inter[0][1]),show=False)

    for i in range(1, len(spl)):
        paux = sy.plot(spl[i],(x,inter[i][0],inter[i][1]),show=False)
        p.append(paux[0])


    p2 = get_sympy_subplots(p)
    p2.plot(v,fx,"o")
    #p2.show()

    return spl, p2

def spline_sujeto(fx,v,fpx0,fpx1 ):
    """
    It takes a list of x values, a list of y values, and the first and second derivatives of the first and last points, and
    returns a plot of the cubic spline interpolation

    :param fx: the function values
    :param v: the x values of the points
    :param fpx0: the first derivative of the function at the first point
    :param fpx1: the derivative of the function at the last point
    """

    inter = []
    fxinter = []

    hi =[]
    for i in range(0,len(v)-1):
        inter.append((v[i],v[i+1]))
        fxinter.append((fx[i],fx[i+1]))

    #print(inter)
    for i in range(0,len(inter)):
        hi.append(inter[i][1]-inter[i][0])

    m = np.zeros(len(v)**2).reshape(len(fx),len(fx))
    #print(hi)
    #print(m)
    for i in range(0,len(v)):
        for j in range(0,len(v)):
            if (i == j and i == 0 and j == 0) :
                m[i][j] = 2*hi[i]
                m[i][j+1] = hi[i]
                continue
            elif (j == i and i == len(v)-1 and j == len(v)-1):
                m[i][j] = 2*hi[-1]
                m[i][j-1] = hi[-1]
                continue
            else:
                if (i == j):
                    m[i][j] = 2*(hi[i-1]+hi[i])
                    m[i][j-1] = hi[i-1]
                    m[i][j+1] = hi[i]

    b = np.zeros(len(v))
    b[0] = ((3/hi[0])*(fx[1]-fx[0]))- (3*fpx0)
    b[-1] = (3*fpx1)-((3/hi[-1])*(fx[-1]-fx[len(fx)-2]))

    for i in range(1,len(v)-1):
        b[i] = ((3/hi[i])*(fx[i+1]-fx[i]))-((3/hi[i-1])*(fx[i]-fx[i-1]))

    #print(m)
    #pprint(Matrix(b.transpose()))

    c = (sy.Matrix(m).inv())*sy.Matrix(b.transpose())
    #pprint(c)
    b = []

    for i in range(0,len(hi)):
        b.append(((fx[i+1]-fx[i])/hi[i])-((((2*c[i])+c[i+1])*hi[i])/3))

    #pprint(Matrix(b))

    d = []

    for i in range(0,len(hi)):
        d.append((c[i+1]-c[i])/(3*hi[i]))

    #pprint(Matrix(d))


    x = sy.symbols('x')
    spl = []
    for i in range(0,len(inter)):
        spl.append(fx[i]+ (b[i]*(x-v[i]))+(c[i]*((x-v[i])**2)) + (d[i]*((x-v[i])**3)))

    #pprint(Matrix(spl))



    p = sy.plot(spl[0], (x,inter[0][0],inter[0][1]),show=False)

    for i in range(1, len(spl)):
        paux = sy.plot(spl[i],(x,inter[i][0],inter[i][1]),show=False)
        p.append(paux[0])


    p2 = get_sympy_subplots(p)
    p2.plot(v,fx,"o")
    return spl,p2


st.title(':blue[Interpolación por Splines Cubicos]')

r'''

**Interpolación por Splines Cúbicos**

La interpolación por splines cúbicos es un método utilizado para aproximar una función desconocida a través de un conjunto de puntos de datos conocidos. Este método divide el rango de datos en segmentos más pequeños y ajusta polinomios cúbicos en cada segmento, garantizando que la función interpolante sea suave y tenga derivadas continuas.

La idea principal detrás de la interpolación por splines cúbicos es encontrar un conjunto de polinomios cúbicos que se ajusten a los puntos de datos y satisfagan ciertas condiciones de suavidad. Estas condiciones incluyen la continuidad de la función, la continuidad de la primera derivada y la continuidad de la segunda derivada en los puntos de intersección de los segmentos.

El polinomio spline cúbico en cada segmento se define de la siguiente manera:

$$
S_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3
$$

Donde:
- $S_i(x)$ es el polinomio spline cúbico en el segmento \(i\).
- $x_i$ es el punto de inicio del segmento \(i\).
- $a_i$, $b_i$, $c_i$ y $d_i$ son los coeficientes del polinomio.

Para determinar los coeficientes de los polinomios spline cúbicos, se deben resolver un sistema de ecuaciones lineales utilizando las condiciones de suavidad. Estas condiciones aseguran que la función interpolante sea continua y suave en todos los puntos de datos.

La interpolación por splines cúbicos es ampliamente utilizada debido a su capacidad para proporcionar una aproximación suave de la función original, manteniendo la continuidad y suavidad de las derivadas. También permite interpolar datos no equidistantes.

'''

st.subheader(':blue[Ejemplo]')

r'''
Supongamos que tenemos los siguientes puntos de datos:

```
x = [1, 2, 3, 4]
y = [3, 1, 2, 4]
```

Queremos encontrar los polinomios spline cúbicos que se ajusten a estos puntos y utilizarlos para estimar el valor de la función en otro punto, por ejemplo, en \(x = 2.5\).

**Paso 1:** Definir los puntos de datos.

```
x = [1, 2, 3, 4]
y = [3, 1, 2, 4]
```

**Paso 2:** Calcular los coeficientes de los polinomios spline cúbicos.

Para cada segmento \(i\) entre los puntos \(x_i\) y \(x_{i+1}\), se deben determinar los coeficientes \(a_i\), \(b_i\), \(c_i\) y \(d_i\) del polinomio spline cúbico \(S_i(x)\).

**Paso 3:** Construir el polinomio spline cúbico.

Una vez que se han determinado los coeficientes de los polinomios spline cúbicos, podemos construir el polinomio spline cúbico global utilizando los polinomios correspondientes a cada segmento.

En este caso, los polinomios spline cúbicos resultantes son:

Para el segmento entre \(x_0 = 1\) y \(x_1 = 2\):
```
S0(x) = 3 + (-2)(x - 1) + 0(x - 1)^2 + (1)(x - 1)^3
      = -2x + 4
```

Para el segmento entre \(x_1 = 2\) y \(x_2 = 3\):
```
S1(x) = 1 + (-2)(x - 2) + 3(x - 2)^2 + (-2)(x - 2)^3
      = 7x^3 - 27x^2 + 34x - 12
```

Para el segmento entre \(x_2 = 3\) y \(x_3 = 4\):
```
S2(x) = 2 + (1)(x - 3) + (-2)(x - 3)^2 + (1)(x - 3)^3
      = 5x^3 - 39x^2 + 99x - 64
```

Por lo tanto, el polinomio spline cúbico global es:

```
S(x) = -2x + 4 en el intervalo [1, 2]
     = 7x^3 - 27x^2 + 34x - 12 en el intervalo [2, 3]
     = 5x^3 - 39x^2 + 99x - 64 en el intervalo [3, 4]
```

**Paso 4:** Estimar el valor de la función en \(x = 2.5\).

Dado que \(x = 2.5\) se encuentra en el intervalo [2, 3], utilizamos el polinomio \(S1(x)\) para estimar el valor de la función:

```
S1(2.5) = 7(2.5)^3 - 27(2.5)^

2 + 34(2.5) - 12
        = 9.625
```

Por lo tanto, el valor estimado de la función en \(x = 2.5\) utilizando la interpolación por splines cúbicos es \(S(2.5) = 9.625\).
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

splinetype = st.selectbox('Ingrese el tipo de Spline a utilizar:',('Natural','Sujeto'),index=1)
dfxstr = ''
if splinetype =='Sujeto':
    dfxx = st.text_input('Ingrese el valor de las derivadas de el primer y ultimo termino:',value='{10,-10}')
    for t in dfxx:

        if t != '{' and t != '}' and t != '[' and t != ']' and t != '(' and t != ')' and t != ' ':
            dfxstr = dfxstr + t

    dfx = list(map(float,dfxstr.split(',')))

    #st.write(dfx)

if splinetype == 'Natural':
    method = spline_natural(fx,x)

if splinetype == 'Sujeto':
    method = spline_sujeto(fx,x,dfx[0],dfx[1])

st.write('Los Splines estan dados por:')

l = r'''f(x) = \begin{cases}'''


for i in method[0]:
    l = l + sy.latex(sy.expand(i)) + r'\\'


l = l + r'''\end{cases}'''

st.latex(l)

plo = gro.Figure()

for i in range(0, len(method[0])):
    spl = sy.lambdify(sy.symbols('x'),method[0][i])
    plo.add_trace(gro.Scatter(x=np.linspace(x[i],x[i+1],1000),y=spl(np.linspace(x[i],x[i+1],1000)),name='Spline '+str(i)))

plo.add_trace(gro.Scatter(x=x,y=fx,mode='markers',name='Puntos'))
plo.update_layout(title='Grafica de la Interpolación')
st.plotly_chart(plo)

