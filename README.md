# TDA

INF421(Conception et analyse d'algorithmes) Projet : Topological Data Analysis

## SphereMin

Etant donné un nuage des points, ce programme contient deux fonctions (en 2D et en 3D) qui ont pour le but de calculer
la boule minimum qui contient ces points.

```python
draw_circle_points(SphereMin(points, []), points)
```

Il suffit de lancer directement le programme pour calculer le centre ainsi que le rayon de cette boule et la visualiser.

## Čech

Ce programme consiste à afficher et renvoyer des simplexes de Čech (de dimension inférieure à k et de valeur de
filtration plus petite que l) ainsi que leur valeur de filtration.

```python
Cech(points, k, l, dim)
```

## Alpha

Ce programme consiste à afficher et renvoyer des simplexes de Alpha (de dimension inférieure à k et de valeur de
filtration plus petite que l) ainsi que leur valeur de filtration.

```python
Alpha(points, k, l, dim)
```

## Draw

Ce programme dessine la complexe de Čech et d'Alpha pour les nuages de points 2D. Il suffit de lancer le programme et
faire déplacer la barre.

```python
show(Alpha(points, k, l, dim))
show(Cech(points, k, l, dim))
```

## Graph Optimization

Ce programme peut optimiser un graphe et écrire le graphe optimisé ainsi que le graph original dans un fichier pour que
l'on puisse tester directement sur [Ripser](https://geometrica.saclay.inria.fr/team/Marc.Glisse/tmp/ripser/) (Il faut
choisir "sparse graph").

Il faut d'abord définir le graphe. Trois exemples sont inclus. Il est aussi possible de charger un graph qui est dans un
fichier.