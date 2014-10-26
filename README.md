FGC_RayTracing
==============

Command line raytracing for a simple scene loader with YAML format

Command: <pre>raytrace -iescena1.txt -osalida.png -s640x480</pre>

Scene Example:
<pre><code>
  ---
  tipo: Camara
  posicion: 0 0 0
  direccion: 0 0 -1
  sombras: true
  aliasing: false
  ---
  tipo: Material
  nombre: Azul1
  difuso: 0.1 0.1 0.9
  especular: 1.0 1.0 1.0
  reflexion: 0.3
  refracion: 0.0
  transparencia: 0.0 # totalmente opaco
  ---
  tipo: Material
  nombre: Azul3
  difuso: 0.1 0.1 0.9
  especular: 1.0 1.0 1.0
  reflexion: 0.3
  refracion: 0.0
  transparencia: 0.0 # totalmente opaco
  ---
  tipo: Esfera
  material: Azul1
  centro: -1 1.2 -10
  radio: 1.5
  ---
  tipo: Cilindro
  material: Verde3
  centro: 0 -2
  eje: Z
  radio: 1
  min: -2
  max: 3
  ---
  tipo: Triangulo
  material: Rojo2
  punto1: 4 5 6
  punto2: -1.5 2 -3
  punto3: 1.0 -1 -1
  ---
  tipo: Luz
  posicion: -5 -5 -5
  color: 1 0 0.5
  ---
  tipo: Cubo
  material: Verde3
  carasx: -1 1
  carasy: -1 1
  carasz: -2 -3
</code></pre>
