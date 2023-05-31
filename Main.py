import streamlit as st
import pandas as pd
import numpy as np
import sympy as sy
from matplotlib import pyplot as plt
from sympy.plotting.plot import MatplotlibBackend, Plot
from sympy.plotting import plot3d,plot3d_parametric_line
import plotly as ply
import base64

file1_ = open("./Im1.png", "rb")


#Portada
contents1 = file1_.read()
data_url1 = base64.b64encode(contents1).decode("utf-8")
file1_.close()

st.markdown(
    f'<img src="data:image/gif;base64,{data_url1}" alt="Escudo_Unam">',
    unsafe_allow_html=True,
)


st.title('''UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO\n

            ''',)

st.header('FACULTAD DE ESTUDIOS SUPERIORES ACATLÁN')

st.header('MÉTODOS NUMERICOS II')

st.subheader('''
PROFESOR\n
_Julio César Galindo López_
''')

st.subheader('''''')

st.subheader('''
        Presenta \n
            Edgar Ramírez Bautista
''')



