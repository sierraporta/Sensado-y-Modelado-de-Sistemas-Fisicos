# Sensado-y-Modelado-de-Sistemas-Fisicos

## Utils:
1. [Utils001_CurveFit.ipynb](https://github.com/sierraporta/Sensado-y-Modelado-de-Sistemas-Fisicos/blob/main/Utils001_CurveFit.ipynb) es un código útil para preparar ajuste de curvas experimentales.

## Dedicados
1. El primer código [Detecting_Points.ipynb](https://github.com/sierraporta/Sensado-y-Modelado-de-Sistemas-Fisicos/blob/main/Detecting_Points.ipynb) trata de un pequeño código básico para estraer características de una imagen, especialmente rata de extraer punto de una imagen usando blur y centroides. Esto es especialmnete útil para nuestro proyecto de cálculo de gravedad usando imágenes extraidas de un video de reproduccion de una pequeña esfera cayendo libremente.
2. Módulo de ingesta de video [Video decomposing.ipynb](https://github.com/sierraporta/Sensado-y-Modelado-de-Sistemas-Fisicos/blob/main/Video_Analysis.ipynb): fotogramas, tiempos y métricas básicas. Este bloque de **utilidad** prepara cualquier video para análisis cuantitativo. Es la **primera etapa común** para los proyectos de:
- **Caída libre** (estimación de $g$ a partir de $y(t)=y_0+v_0 t+\tfrac12 g t^2$),
- **Pelota rebotando** (energía y coeficiente de restitución),
- **Péndulo** (período $T$, amortiguamiento),
- **Viscosímetro de esfera** (velocidad terminal $v_t$ y $\mu$),
- y en general, **cualquier experimento** que parta de un video y necesite un *pipeline* reproducible.
3. **Detección de picos en un péndulo simple**. Este script genera datos simulados del movimiento de un péndulo simple, reconstruye el ángulo $θ(t)$ a partir de las posiciones $(x,y)$ y aplica un filtro de suavizado (Savitzky–Golay) para reducir el efecto del ruido. Posteriormente utiliza scipy.signal.find_peaks para identificar los máximos locales de la oscilación y graficarlos sobre la señal.
La visualización permite:
- Ver la dinámica oscilatoria del péndulo.
- Identificar los picos detectados automáticamente (en rojo).
- Usar la separación temporal entre picos para estimar el período experimental del péndulo.
4. **Resolviendo Ecuaciones Diferenciales**. Hay un Notebook: [intro_linear_ode_cooling.ipynb](https://github.com/sierraporta/Sensado-y-Modelado-de-Sistemas-Fisicos/blob/main/intro_linear_ode_cooling.ipynb). 
5. **Sistema Masa Resorte**. Si ya tienes $x(t_i)$ de la cámara (posiciones muestreadas en el tiempo), puedes estimar parámetros físicos ajustando directamente el modelo de ecuaciones diferenciales del oscilador. Trabaja siempre respecto al equilibrio, definiendo $y(t)=x(t)-x_{\mathrm{eq}}$ para que la gravedad desaparezca. El modelo libre y lineal con amortiguamiento viscoso es
\[
m\,\ddot y + b\,\dot y + k\,y = 0
\;\;\Longleftrightarrow\;\;
\ddot y + 2\zeta\omega_0\,\dot y + \omega_0^2\,y = 0,
\]
donde $\omega_0=\sqrt{k/m}$ y $\zeta=\frac{b}{2m\omega_0}=\frac{b}{2\sqrt{mk}}$. Si usas excitación conocida, el modelo forzado es $m\,\ddot y + b\,\dot y + k\,y = u(t)$. En la práctica conviene fijar una masa efectiva $m_{\mathrm{ef}}$ (masa añadida más ≈ un tercio de la masa del resorte) y estimar $k$ y $b$; si deseas estimar también $m$, usa datos forzados o varias corridas con masas distintas.
- En la “regresión de derivadas” se suaviza $y(t)$ y se calculan $\dot y$ y $\ddot y$ por derivación numérica (por ejemplo, con Savitzky–Golay); luego se ajusta la forma normalizada $\ddot y = -a\,\dot y - b\,y$ con $a=2\zeta\omega_0$ y $b=\omega_0^2$, resolviendo por mínimos cuadrados $\min_{a,b}\sum_i[\ddot y_i + a\,\dot y_i + b\,y_i]^2$. A partir de $a$ y $b$ se recupera $\omega_0=\sqrt{b}$ y $\zeta=a/(2\omega_0)$, y con la masa se obtiene $k=m\,\omega_0^2$ y $b_{\mathrm{visc}}=2m\zeta\omega_0$. Este método es muy simple y rápido, pero sensible al ruido introducido por las derivadas, por lo que requiere un buen suavizado.
- En el “ajuste por solución de la ODE” se toman condiciones iniciales del dato ($y(0)=y_0$ y $\dot y(0)\approx (y_1-y_0)/\Delta t$), se consideran como parámetros $\theta=(k,b)$ manteniendo $m$ fijo, y se integra la ecuación $\ddot y + \tfrac{b}{m}\dot y + \tfrac{k}{m}y=0$ sobre los tiempos medidos (por ejemplo, con Runge–Kutta). Los parámetros se estiman minimizando el costo $J(\theta)=\sum_i [y_{\mathrm{data}}(t_i)-y_{\mathrm{model}}(t_i;\theta)]^2$ —opcionalmente con pérdida robusta como Huber— y se optimiza con algoritmos tipo Levenberg–Marquardt o TRF; las incertidumbres pueden obtenerse de la jacobiana o por bootstrap. Este enfoque no depende de derivar los datos y suele ser más estable, a costa de un poco más de cómputo.
- De modo analítico, en régimen subamortiguado la solución libre es $y(t)=A\,e^{-\gamma t}\cos(\omega_d t+\phi)$, con $\gamma=\tfrac{b}{2m}$ y $\omega_d=\sqrt{\omega_0^2-\gamma^2}$. La frecuencia amortiguada $\omega_d$ se obtiene promediando periodos entre picos o ceros; el amortiguamiento $\gamma$ se estima linealizando la envolvente $\ln|y_{\text{pico}}(t)|\approx \ln A-\gamma t$; luego $\omega_0=\sqrt{\omega_d^2+\gamma^2}$ y, en consecuencia, $k=m\,\omega_0^2$ y $b=2m\gamma$. Este camino es muy transparente y rápido, aunque asume linealidad y requiere picos limpios con poco ruido.
- **Revisa los dos códigos***: [mass_spring_methods_demo.ipynb](https://github.com/sierraporta/Sensado-y-Modelado-de-Sistemas-Fisicos/blob/main/mass_spring_methods_demo.ipynb) y [mass_spring_fit.py](https://github.com/sierraporta/Sensado-y-Modelado-de-Sistemas-Fisicos/blob/main/mass_spring_fit.py).
