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


def parse_inputsys(inp):
    eq = []
    sfix = ''

    for i in inp:
        if i != '{' and i != '}' and i != '[' and i != ']' and i != '(' and i != ')' and i != ' ':
            sfix = sfix + i

    for t in sfix.split(','):
        eq.append(sy.parse_expr(t,transformations='all'))

    return eq

def get_maxsymbs(sys):
    maxs = 0
    s = []
    for i in sys:
        if max(maxs,len(i.free_symbols)) != maxs:
            s = list(i.free_symbols)

    return s

def jacobian(ff,symb):
    """
    It takes a vector of functions and a vector of symbols and returns the Jacobian matrix of the functions with respect to
    the symbols
    :param ff: the function
    :param symb: the symbols that are used in the function
    :return: A matrix of the partial derivatives of the function with respect to the variables.
    """
    m = []

    for i in range(0,len(ff)):
        aux  = []
        for j in range(0,len(symb)):
            aux.append(sy.diff(ff[i],symb[j]))
        m.append(aux)

    return np.array(m)

def eval_matrix(matrix , v,symb):
    """
    It takes a matrix, a list of symbols and a list of values, and returns the matrix with the symbols substituted by the
    values

    :param matrix: the matrix of the system of equations
    :param v: the vector of values for the variables
    :param symb: the symbols that will be used in the matrix
    :return: the matrix with the values of the variables substituted by the values of the vector v.
    """
    e = 0
    mm = []
    for i in range(0,len(matrix)):
        aux = []
        ev = []
        for k in range(0,len(symb)):
            ev.append((symb[k],v[k]))
        for j in range(len(matrix[i])):
            aux.append(matrix[i][j].subs(ev).evalf())
        mm.append(aux)
    return np.array(mm)

def evalVector(ff, x0,symb):
    """
    > Given a list of symbolic expressions, a list of values for the symbols, and a list of the symbols, evaluate the
    symbolic expressions at the given values

    :param ff: the vector of functions
    :param x0: initial guess
    :param symb: the symbols that are used in the symbolic expression
    :return: the value of the function at the point x0.
    """
    v = []
    for i in range(0,len(ff)):
        ev = []

        for k in range(0,len(x0)):
            ev.append((symb[k],x0[k]))

        v.append(ff[i].subs(ev).evalf())
    return np.array(v)

def NewtonMethod( ff, x0,symb ):
    """
    The function takes in a vector of functions, a vector of initial guesses, and a vector of symbols. It then calculates
    the Jacobian matrix, the Jacobian matrix evaluated at the initial guess, the inverse of the Jacobian matrix evaluated at
    the initial guess, the vector of functions evaluated at the initial guess, and then the Newton step.

    The function returns the Newton step.

    :param ff: the function we want to find the root of
    :param x0: initial guess
    :param symb: the symbols used in the function
    :return: The return value is the x_np1 value.
    """
    j = jacobian(ff,symb)
    #print("Jacobian Matrix")
    #pprint(Matrix(j))
    jev = sy.Matrix( eval_matrix(j,x0,symb))
    #print("J(",x0,")")
    #pprint(jev)

    jinv = jev.inv()
    #print("F(",x0,")")
    ffev = sy.Matrix(evalVector(np.transpose(ff),x0,symb))
    #print("J^-1(",x0,")*","F(",x0,")")
    mm = sy.Matrix(jinv)*ffev
    #pprint(mm)
    x_np1 = sy.Matrix(np.transpose(np.array(x0)))
    #pprint(x_np1-mm)
    return list(x_np1-mm)

def norm_inf(x_0,x_1):
    """
    > The function `norm_inf` takes two vectors `x_0` and `x_1` and returns the maximum absolute difference between the two
    vectors

    :param x_0: the initial guess
    :param x_1: the vector of the current iteration
    :return: The maximum difference between the two vectors.
    """
    a = [abs(x_1[i]-x_0[i]) for i in range(len(x_0))]
    return max(a)

def newton_method(ff,x_0,symbs,error,maxiter):
    """
    Given a function (x,y)$, a starting point $, and a list of symbols,
    the function will return the next point $ in the Newton's method sequence

    :param ff: the function to be minimized
    :param x_0: initial guess
    :param symbs: the symbols that we're using in the function
    :return: the final value of x_0, the list of x values, and the list of y values.
    """
    #pprint(Matrix(x_0))
    xs = []
    ys = []
    xns = [x_0]
    erros = []
    iterr = 0
    while True and iterr < maxiter:

        x_1 = NewtonMethod(ff,x_0,symbs)
        #print(x_1)
        ninf = norm_inf(x_0,x_1)
        erros.append(ninf)
        #print(ninf)

        x_0 = list(x_1)
        xns.append(tuple(x_0))
        xs.append(x_0[0])
        ys.append(x_0[1])
        if ninf < error:
            #print("Iteraciones: ",iterr)
            break
        iterr = iterr+1

    #print(x_0)
    return xns,xs,ys,erros


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


st.title(':blue[Metodo de Newton]')

r'''
**Método de Newton para Sistemas de Ecuaciones No Lineales**

El método de Newton es una técnica iterativa utilizada para encontrar soluciones de sistemas de ecuaciones no lineales. Este método se basa en la aproximación de una función no lineal por medio de su serie de Taylor truncada hasta el primer término, y luego se busca el punto donde esta aproximación se anula.

Supongamos que tenemos un sistema de ecuaciones no lineales:

$$
\begin{align*}
f_1(x_1, x_2, \ldots, x_n) &= 0 \\
f_2(x_1, x_2, \ldots, x_n) &= 0 \\
\ldots \\
f_n(x_1, x_2, \ldots, x_n) &= 0 \\
\end{align*}
$$

El método de Newton busca una solución aproximada para este sistema iterativamente. Comenzamos con una estimación inicial del punto de solución \((x_1^0, x_2^0, \ldots, x_n^0)\) y utilizamos la siguiente fórmula para mejorar la aproximación:

$$
\mathbf{X}^{k+1} = \mathbf{X}^k - J^{-1}(\mathbf{X}^k) \cdot \mathbf{F}(\mathbf{X}^k)
$$

Donde:
- $\mathbf{X}^k$ es la estimación actual del punto de solución en la iteración \(k\).
- $\mathbf{F} \mathbf{X}^k)$ es el vector de funciones $f_1, f_2, \ldots, f_n$ evaluado en $\mathbf{X}^k$.
- $J(\mathbf{X}^k)$ es la matriz Jacobiana evaluada en $\mathbf{X}^k$.
- $J^{-1}(\mathbf{X}^k)$ es la inversa de la matriz Jacobiana evaluada en $\mathbf{X}^k$.

El proceso se repite hasta que se alcanza la convergencia, es decir, cuando la diferencia entre dos iteraciones sucesivas es lo suficientemente pequeña.

Es importante tener en cuenta que el método de Newton puede no converger o converger a un mínimo local. Por lo tanto, es recomendable seleccionar una buena estimación inicial y verificar la convergencia mediante criterios establecidos.

El método de Newton para sistemas de ecuaciones no lineales es ampliamente utilizado debido a su eficiencia y rapidez de convergencia en muchos casos. Sin embargo, también tiene limitaciones y puede requerir ajustes adicionales en situaciones particulares.
'''



st.subheader(':blue[Ejemplo]')


r'''
Supongamos que queremos encontrar una solución para el siguiente sistema de ecuaciones no lineales:

$$
\begin{align*}
x^2 + y^2 &= 25 \\
x + y &= 7
\end{align*}
$$

Podemos aplicar el método de Newton para encontrar una aproximación de la solución. Comenzamos con una estimación inicial del punto de solución \((x^0, y^0)\) y utilizamos las siguientes ecuaciones para mejorar la aproximación:

$$
\begin{align*}
x^{k+1} &= x^k - \frac{{f_1(x^k, y^k)}}{{\frac{{\partial f_1}}{{\partial x}}(x^k, y^k)}} \\
y^{k+1} &= y^k - \frac{{f_2(x^k, y^k)}}{{\frac{{\partial f_2}}{{\partial y}}(x^k, y^k)}}
\end{align*}
$$

Donde:
- $x^k, y^k$ son las estimaciones actuales del punto de solución en la iteración $k$.
- $f_1(x, y) = x^2 + y^2 - 25$ es la primera ecuación del sistema.
- $f_2(x, y) = x + y - 7$ es la segunda ecuación del sistema.
- $\frac{{\partial f_1}}{{\partial x}}, \frac{{\partial f_1}}{{\partial y}}, \frac{{\partial f_2}}{{\partial x}}, \frac{{\partial f_2}}{{\partial y}}$ son las derivadas parciales de las funciones respecto a las variables correspondientes.

Podemos seleccionar una estimación inicial $(x^0, y^0) = (2, 5)$ y aplicar el método de Newton iterativamente hasta alcanzar la convergencia.

**Paso 1:** Estimación inicial
$(x^0, y^0) = (2, 5)$

**Paso 2:** Iteraciones
Aplicamos las ecuaciones iterativamente hasta alcanzar la convergencia:

- Iteración 1:
  $$
  \begin{align*}
  x^1 &= x^0 - \frac{{f_1(x^0, y^0)}}{{\frac{{\partial f_1}}{{\partial x}}(x^0, y^0)}} = 2 - \frac{{(2^2 + 5^2 - 25)}}{{(2x^0)}} = \frac{{19}}{{4}} \approx 4.75 \\
  y^1 &= y^0 - \frac{{f_2(x^0, y^0)}}{{\frac{{\partial f_2}}{{\partial y}}(x^0, y^0)}} = 5 - \frac{{(2 + 5 - 7)}}{{1}} = 5
  \end{align*}
  $$

- Iteración 2:
  $$

  x^2 = x^1 - \frac{{f_1(x^1, y^1)}}{{\frac{{\partial f_1}}{{\partial x}}(x^1, y^1)}} = \frac{{19}}{{4}} -\frac{{(\left(\frac{{19}}{{4}}\right)^2 + 5^2 - 25)}}{{\left(\frac{{19}}{{2}}\right)}} \approx 4.794 \\
  y^2 = y^1 - \frac{{f_2(x^1, y^1)}}{{\frac{{\partial f_2}}{{\partial y}}(x^1, y^1)}} = 5 - \frac{{(\frac{{19}}{{4}} + 5 - 7)}}{{1}} = 5

  $$

Continuamos con las iteraciones hasta que la diferencia entre dos iteraciones sucesivas sea lo suficientemente pequeña.

**Paso 3:** Convergencia
Verificamos la convergencia mediante criterios establecidos, como la tolerancia o el número máximo de iteraciones. Si la convergencia no se alcanza, podemos ajustar la estimación inicial y repetir el proceso.


En este ejemplo, después de varias iteraciones, obtenemos una estimación del punto de solución $(x, y) \approx (4.795,2.205)$ como una aproximación de la solución del sistema de ecuaciones no lineales.
'''





st.subheader('Método')
sys = st.text_input('Ingrese el sistema de ecuaciones ',value=r'{x^2+3*y*x+2,y**3*x**2-2*y**3-5}')
try:
    system = parse_inputsys(sys)
except :
    st.error('Error al introducir el sistema de ecuaciones', icon="🚨")

st.write('_El sistema de ecuaciones es:_')
psys ='F'+str(tuple(get_maxsymbs(system)))+ r'''
    =\begin{cases}

'''
for i in system:
    psys = psys + sy.latex(i)
    psys = psys + r' = 0\\'

psys = psys + r'\end{cases}'

st.latex(psys)

fx = sy.lambdify(list(get_maxsymbs(system)),system[0])
fx2 = sy.lambdify(list(get_maxsymbs(system)),system[1])




try:

    x,y = sy.symbols('x,y')
    st.write('_Grafica 3D del sistema_')

    plo1 = gro.Figure()

    ran = np.linspace(-10,10,100)
    su1 = [[fx(ran[xs],ran[ys]) for xs in range (0,len(ran)) ] for ys in range(0,len(ran))]
    su2 = [[fx2(ran[xs],ran[ys]) for xs in range (0,len(ran)) ] for ys in range(0,len(ran))]

    plo1.add_trace(gro.Surface(z=su1,name='Ecuación 1',opacity=.7))
    plo1.add_trace(gro.Surface(z=su2,name='Ecuación 2',opacity=.7))

    st.plotly_chart(plo1)


except Exception as excep:
    st.error(excep)
    st.error('Error al graficar', icon="🚨")


initaprx = st.text_input('Ingrese una aproximacion inicial $x_0$: ',value=[-1,1])

intaprox = []
intstr = ''




for i in initaprx:

    if i != '{' and i != '}' and i != '[' and i != ']' and i != ' ':
        intstr = intstr + i

try:
    st.write('La aproximacion inicial es: ')
    intaprox = list(map(int, intstr.split(',')))
    st.latex(sy.latex(sy.Matrix(list(intaprox))))
except:
    st.error('Error al introducir la aproximación inicial', icon="🚨")


err = st.text_input('Ingrese el error de tolerancia: ',value='0.00001')
try:
    st.write('El error de tolerancia es:', float(err))
except:
    st.error('Error al introducir el error de tolerancia', icon="🚨")


maxiter = st.slider('Maximo de Iteraciones',10,1000,10)


st.write('_Matrix Jacobiana_:')
symbs = (get_maxsymbs(system))


st.latex(sy.latex(sy.Matrix(jacobian(system,symbs))))

method = newton_method(system,list(intaprox),symbs,float(err),maxiter)

tabb = []

for i in range(0,len(method[1])):
    aux = list(method[0][i])
    aux.append(method[3][i])
    tabb.append(aux)

cols = list(map(str, list(symbs)))
cols.append('Error')
#st.write(tabb)
table = pd.DataFrame(tabb,columns=cols)

st.write(table)
try:


    #st.write(method[2])

    xs =[float(i) for i in method[1]]
    xy =[float(i) for i in method[2]]
    evalfx = [float(fx(i[0],i[1])) for i in method[0]]

    plo = gro.Figure()
    plo.add_trace(gro.Scatter3d(x=xs,y=xy,z=evalfx,name='Aproximaciones'))

    ranx = np.linspace(int(method[1][-1]-1),int(method[1][-1]+1),100)
    rany = np.linspace(int(method[2][-1]-1),int(method[2][-1]+1),100)
    su1 = [[fx(ranx[xs],rany[ys]) for xs in range (0,len(ranx)) ] for ys in range(0,len(rany))]
    su2 = [[fx2(ranx[xs],rany[ys]) for xs in range (0,len(ranx)) ] for ys in range(0,len(rany))]


    plo.add_trace(gro.Surface(z=su1,name='Ecuación 1',opacity=.8))
    plo.add_trace(gro.Surface(z=su2,name='Ecuación 2',opacity=.7))

    plo.update_layout(title='Grafica del Sistema')
    st.plotly_chart(plo)




except Exception as e:
    st.error(str(e))




