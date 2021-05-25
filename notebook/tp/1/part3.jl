### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 5f0b6252-2387-42e4-a46a-5b3f5b85e78f
using LinearAlgebra

# ╔═╡ c0db5a3e-f39c-11ea-12fa-3b2f78b7b2b7
md"""
Notebook au format Pluto disponible [ici](https://github.com/vlc1/Mn3.jl/blob/master/notebook/tp/1/part1.jl).

"""

# ╔═╡ 03cb7526-f392-11ea-123a-797739ba60de
md"""
# Exercice 1 -- Manipulation de matrices

1. Définir le vecteur ``U = \left [ \begin{array}{cccccc} 0 & 1 & 2 & 3 & \cdots & 49 & 50 \end{array} \right ]``. Quelle est sa taille ?
1. Définir le vecteur ``V`` contenant les cinq premiers éléments de ``U``, et le vecteur ``W`` contenant les cinq premiers et les cinq derniers éléments de ``U``.
1. Définir la matrice
```math
M =
\left [ \begin{array}{ccccccc}
1 & 2 & 3 & \cdots & 8 & 9 & 10 \\
11 & 12 & 13 & \cdots & 18 & 19 & 20 \\
21 & 22 & 23 & \cdots & 28 & 29 & 30
\end{array}
\right ].
```
4. Extraire de ``M`` les matrices
```math
N =
\left [ \begin{array}{cc}
1 & 2 \\
11 & 12 \\
21 & 22
\end{array} \right ], \quad P =
\left [ \begin{array}{ccc}
8 & 9 & 10 \\
18 & 19 & 20 \\
28 & 29 & 30
\end{array} \right ] \quad \mathrm{et} \quad Q =
\left [ \begin{array}{cc}
3 & 7 \\
23 & 27 \end{array} \right ].
```
5. Extraire de la matrice ``M`` la matrice ``R`` obtenue en prenant une colonne sur deux.
1. Définir les matrices ``A = \left [ \begin{array}{ccccc} -1 & -3 & -5 & \cdots & -99 \end{array} \right ]`` et ``B = \left [ \begin{array}{ccccc} 2 & 4 & 6 & 8 & \cdots & 100 \end{array} \right ]`` puis le vecteur ``C = \left [ \begin{array}{ccccccccc} -1 & 2 & -3 & 4 & -5 & 6 & \cdots & -99 & 100 \end{array} \right ]``.

"""

# ╔═╡ e1bea777-2e25-4aca-bd1f-251fae321923
# Q1a
U = [i for i in 0:50]

# ╔═╡ ae4090f3-a0cc-4e2f-81c5-fb17d7d6fd38
# Q1b
length(U)

# ╔═╡ c925a309-8243-4f7b-aca4-bfac2e450547
# Q2a -- option 1
# V = [U[i] for i in 1:5]
# Q2b -- option 2
V = U[1:5]

# ╔═╡ 0594bf02-2bb7-49d6-9cf3-5e008301ce91
# Q2b -- option 1
#W = vcat(U[1:5], U[47:51])
# Q2b -- option 2
W = [U[1:5]; U[47:51]]

# ╔═╡ 254eb823-e08f-4b7c-9bc0-298aa9587cc3
# Q3 -- option 1
#M = transpose(reshape(1:30, 10, 3))
# Q3 -- option 2
M = hcat(1:10:21, 2:10:22, 3:10:23, 4:10:24, 5:10:25,
	 6:10:26, 7:10:27, 8:10:28, 9:10:29, 10:10:30)

# ╔═╡ dd546a33-ddbd-4ae7-ba2c-34afc92c47a9
# Q4a
N = M[:, 1:2]

# ╔═╡ 5d3d424b-cd29-44e0-875e-4750f9fb4f3e
# Q4b -- option 1
#P = M[:, 8:10]
# Q4b -- option 2
P = M[:, end-2:end]

# ╔═╡ 40bc8fd7-f19a-4e48-9c24-dc924f776f1c
# Q4c -- option 1
#Q = M[1:2:3, 3:4:7]
# Q4c -- option 2
Q = [3 7
	23 27]

# ╔═╡ 6cfa5cf2-0dbf-44dd-81f5-bca3d7fad758
# Q5
R = M[:, 1:2:9]

# ╔═╡ 9ec04393-0ea8-4c63-a8e9-ad9e59f57fb8
# Q6a
A = [-2i + 1 for i in 1:50]

# ╔═╡ 27f4425e-3639-461b-a031-2113f706a6c1
# Q6b -- option 1
#B = [2i for i in 1:50]
# Q6b -- option 2
B = [i for i in 2:2:100]

# ╔═╡ bc9e70aa-86ee-4bf2-a23b-daad5b740259
# Q6c
C = reshape([A'; B'], 100)

# ╔═╡ 703524b0-6cca-43ea-af2b-5cd78c762b71
md"""
Concepts utilisés :

* Les "deux points" : `:`, permettent de prendre l'intégralité d'une ligne ou d'une colonne. Exemples : pour prendre la deuxième ligne de `M`
```julia
M[2, :]
```
ou la troisième colonne de `M`
```julia
M[:, 3]
```
* Le "point virgule" : `;`, permet de définir un vecteur ou effectuer une concaténation verticale. Pour définir un vecteur, soit
```julia
A = [1; π; 2 // 3]
```
ou bien de manière équivalent
```julia
A = [1
π
2 // 3]
```
* Les `Range`s permettent de définir un ensemble d'entier équidistants. Exemples : les entiers compris entre 1 et 10,
```julia
1:10
```
(pas unitaire) ou les entiers pairs entre 2 et 6
```julia
2:2:6
```
les multiples de 4 entre 8 et 20
```julia
8:4:20
```

"""

# ╔═╡ ca8cb49a-f39c-11ea-307e-97825e864b0f
md"""
# Exercice 2 -- Matrices et systèmes linéaires

1. Écrire une fonction, n'utilisant aucune boucle (`for`, `while`...) qui prend comme paramètre un entier $n$ et qui construit la matrice suivante :
$$\left [ \begin{array}{ccccccc}
1 & 1 & 0 & \cdots & 0 & 0 & 0 \\
\frac{1}{n} & 2 & \frac{n-1}{n} & \cdots & 0 & 0 & 0 \\
0 & \frac{2}{n} & 3 & \cdots & 0 & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
0 & 0 & 0 & \cdots & n - 1 & \frac{2}{n} & 0 \\
0 & 0 & 0 & \cdots & \frac{n - 1}{n} & n & \frac{1}{n} \\
0 & 0 & 0 & \cdots & 0 & 1 & n + 1
\end{array} \right ].$$
"""

# ╔═╡ 46497466-f39e-11ea-2634-cd5985460790
Markdown.MD(Markdown.Admonition("hint", "Indice", [md"Renseignez-vous sur le type `Tridiagonal` de la bibliothèque standart `LinearAlgebra` grâce au *Live docs*."]))

# ╔═╡ 836e978f-ab3c-411d-adcb-a4dce984bd2d
function foo(n)
	l = [i / n for i in 1:n]
	d = [Float64(i) for i in 1:n + 1]
	u = [(n + 1 - i) / n for i in 1:n]
	Tridiagonal(l, d, u)
end

# ╔═╡ f5c49974-3ebe-4853-abd7-c2295278f430
foo(10)

# ╔═╡ 7fb06e78-f39f-11ea-38d9-ab51b94af7c5
md"""
2. Résoudre numériquement le système
$$\left \{ \begin{aligned}
x + 2y + 3z + 4t & = 1, \\
2x + 3y + 4z + t & = -2, \\
-2x + 4y -5z + 2t & = 0, \\
8x + y - z + 3t & = 1.
\end{aligned} \right .$$

"""

# ╔═╡ 93a8d1be-e2e0-4096-bebf-a3cdd6ef314c
D = [1 2 3 4
	2 3 4 1
	-2 4 -5 2
	8 1 -1 3]

# ╔═╡ 7bbb37cc-8d23-49ca-9518-3a3ac60c2ac6
b = [1; -2; 0; 1]

# ╔═╡ 8e647f3e-083b-42e5-a57e-f3910d83bc5d
x = D \ b

# ╔═╡ f88252ce-4803-4fc8-99b1-5a50fd840af8
D * x

# ╔═╡ Cell order:
# ╟─c0db5a3e-f39c-11ea-12fa-3b2f78b7b2b7
# ╟─03cb7526-f392-11ea-123a-797739ba60de
# ╠═e1bea777-2e25-4aca-bd1f-251fae321923
# ╠═ae4090f3-a0cc-4e2f-81c5-fb17d7d6fd38
# ╠═c925a309-8243-4f7b-aca4-bfac2e450547
# ╠═0594bf02-2bb7-49d6-9cf3-5e008301ce91
# ╠═254eb823-e08f-4b7c-9bc0-298aa9587cc3
# ╠═dd546a33-ddbd-4ae7-ba2c-34afc92c47a9
# ╠═5d3d424b-cd29-44e0-875e-4750f9fb4f3e
# ╠═40bc8fd7-f19a-4e48-9c24-dc924f776f1c
# ╠═6cfa5cf2-0dbf-44dd-81f5-bca3d7fad758
# ╠═9ec04393-0ea8-4c63-a8e9-ad9e59f57fb8
# ╠═27f4425e-3639-461b-a031-2113f706a6c1
# ╠═bc9e70aa-86ee-4bf2-a23b-daad5b740259
# ╟─703524b0-6cca-43ea-af2b-5cd78c762b71
# ╟─ca8cb49a-f39c-11ea-307e-97825e864b0f
# ╟─46497466-f39e-11ea-2634-cd5985460790
# ╠═5f0b6252-2387-42e4-a46a-5b3f5b85e78f
# ╠═836e978f-ab3c-411d-adcb-a4dce984bd2d
# ╠═f5c49974-3ebe-4853-abd7-c2295278f430
# ╟─7fb06e78-f39f-11ea-38d9-ab51b94af7c5
# ╠═93a8d1be-e2e0-4096-bebf-a3cdd6ef314c
# ╠═7bbb37cc-8d23-49ca-9518-3a3ac60c2ac6
# ╠═8e647f3e-083b-42e5-a57e-f3910d83bc5d
# ╠═f88252ce-4803-4fc8-99b1-5a50fd840af8
