# Sensado-y-Modelado-de-Sistemas-Fisicos

1. El primer código [Detecting_Points.ipynb](https://github.com/sierraporta/Sensado-y-Modelado-de-Sistemas-Fisicos/blob/main/Detecting_Points.ipynb) trata de un pequeño código básico para estraer características de una imagen, especialmente rata de extraer punto de una imagen usando blur y centroides. Esto es especialmnete útil para nuestro proyecto de cálculo de gravedad usando imágenes extraidas de un video de reproduccion de una pequeña esfera cayendo libremente.
2. Módulo de ingesta de video [Video decomposing.ipynb](https://github.com/sierraporta/Sensado-y-Modelado-de-Sistemas-Fisicos/blob/main/Video_Analysis.ipynb): fotogramas, tiempos y métricas básicas. Este bloque de **utilidad** prepara cualquier video para análisis cuantitativo. Es la **primera etapa común** para los proyectos de:
- **Caída libre** (estimación de \(g\) a partir de \(y(t)=y_0+v_0 t+\tfrac12 g t^2\)),
- **Pelota rebotando** (energía y coeficiente de restitución),
- **Péndulo** (período \(T\), amortiguamiento),
- **Viscosímetro de esfera** (velocidad terminal \(v_t\) y \(\mu\)),
- y en general, **cualquier experimento** que parta de un video y necesite un *pipeline* reproducible.
 
