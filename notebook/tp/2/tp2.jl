### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 35a850c8-d24c-4593-b551-f3b9bb49e82d
begin
	using NLsolve
	solve(f!, x₀) = getproperty(nlsolve(f!, x₀), :zero)
end

# ╔═╡ ec1b0cec-af5f-492e-ba55-269c2df5b458
using LinearAlgebra

# ╔═╡ e8b460cf-3d35-4bbb-8b3a-5a2c5f8a7a21
using Plots

# ╔═╡ d61e1ea6-f924-11ea-00dc-794c93177d22
md"""
Version [Pluto](https://github.com/vlc1/Mn3.jl/blob/master/notebook/tp/2/tp2.jl) de ce notebook.

"""

# ╔═╡ 3e5d9fa5-acae-42ab-b40e-e1a197edcad2
md"""
# Programmation

Quelques points vont être abordés en début de séance :

1. Deux nouveaux *packages* : `Plots` et `NLsolve` ;
1. Fonctions et types d'arguments ;
1. La notation `.` (*broadcast) ;
1. Premier (`first`) et dernier (`last`) éléments d'un tableau.

"""

# ╔═╡ 4ae18622-f7ec-11ea-2f71-d5b166ff50fb
md"""
# Recherche de la racine d'une fonction

Nous avons vu lors de la deuxième séance que les éléments
```math
\begin{aligned}
y_1 & \simeq y \left ( t_1 \right ), \\
y_2 & \simeq y \left ( t_2 \right ), \\
& \ldots
\end{aligned}
```
de la solution numérique du problème de Cauchy
```math
\left \{ \begin{aligned}
\dot{y} \left ( t \right ) & = f \left [ t, y \left ( t \right ) \right ], \\
y \left ( 0 \right ) & = y_0
\end{aligned} \right .
```
sont définis implicitement, c'est à dire comme racines de fonctions.

L'objectif de cette première partie est de se familiariser avec le *package* `NLsolve.jl` que nous utiliserons afin de résoudre des équations non-linéaires.

Les cellules suivantes décrivent comment obtenir la racine de la function
```math
\left ( \begin{matrix}
x \\
y \end{matrix} \right ) \mapsto \left ( \begin{matrix}
\left ( x + 3 \right ) \left ( y ^ 3 - 7 \right ) + 18 \\
\sin \left [ y \exp \left ( x \right ) - 1 \right ]
\end{matrix} \right )
```
à partir de la donnée initiale
```math
\left ( x_0, y_0 \right ) = \left ( 0.1, 1.2 \right ).
```

"""

# ╔═╡ 88dd8a5e-651a-4939-9231-6696b78f024c
function example!(res, x)
    res[1] = (x[1] + 3) * (x[2] ^ 3 - 7) + 18
    res[2] = sin(x[2] * exp(x[1]) - 1)
	nothing
end

# ╔═╡ 7dccedba-aa05-4821-b1ed-ed873ad9cccb
solve(example!, [0.1; 1.2])

# ╔═╡ 073d29d8-7055-4d0c-8565-8732017b89d0
md"""
1. **Cas scalaire** -- Modifier l'exemple précédent afin de résoudre l'équation de Kepler
```math
10 - x + e \sin \left ( x \right ) = 0.
```
On définira dans un premier temps la fonction `kepler!` définie ci-dessous.

"""

# ╔═╡ bea20ce9-88ff-45e4-bb5a-494c4d2d19c2
# Q1 -- À MODIFIER
function kepler!(res, x)
	res[1] = x[1] - 1
	nothing
end

# ╔═╡ 2c4151ac-a476-4c1f-a379-f4c3a0174b3d
if norm(solve(kepler!, [0.0]) - [5.08912]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `kepler!` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `kepler!`.
	"""
end

# ╔═╡ b72fd704-46f9-4d88-86bd-7add43c47d0e
md"""
2. **Cas vectoriel** -- Modifier la fonction `system!` ci-dessous afin de résoudre le système d'équations
```math
\left \{ \begin{aligned}
x + y + z ^ 2 & = 12, \\
x ^ 2 - y + z & = 2, \\
2x - y ^ 2 + z & = 1.
\end{aligned} \right .
```

"""

# ╔═╡ d264f45c-7a66-48d5-b6c3-b297073f5ce8
# Q2 -- À MODIFIER
function system!(res, x)
	res[1] = x[1] + 1
	res[2] = x[2] + 2
	res[3] = x[3] + 3
	nothing
end

# ╔═╡ 42e83fa4-ba23-4e80-a8e6-58816319b134
if norm(solve(system!, [0.0; 0.0; 0.0]) - [1.0; 2.0; 3.0]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `system!` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `system!`.
	"""
end

# ╔═╡ 7dc345ee-f7ec-11ea-138f-a1c6b81e0260
md"""
# Modèle et solution exacte

On se concentre pour l'instant sur le modèle linéaire homogène pour lequel le second membre de l'EDO s'écrit
```math
f \colon \left ( t, y \right ) \mapsto \lambda y.
```

L'équation différentielle à résoudre s'écrit alors,
```math
\left \{ \begin{aligned}
\dot{y} \left ( t \right ) & = \lambda y \left ( t \right ), \\
y \left ( 0 \right ) & = y_0
\end{aligned} \right .
```
et la solution exacte est donnée sous la forme :
```math
y \colon t \mapsto \exp \left ( \lambda t \right ) y_0.
```

3. Implémenter la fonction ``f``, appelée ci-dessous `linear`, dans le cas ``\lambda = -1``.

"""

# ╔═╡ 58162516-f7ec-11ea-3095-85bde6d71604
# Q3 -- À MODIFIER
linear(t, y) = zero(y)

# ╔═╡ dcf06721-1996-42ce-8c77-93f6f5195c7d
if norm(linear(nothing, [π; 2 // 3]) + [π; 2 // 3]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `linear` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `linear`.
	"""
end

# ╔═╡ 56b5b3ec-7ec2-483c-8120-156665b2c49f
md"""
4. Implémenter la fonction `solution` qui correspond à la solution analytique dans le cas ``\lambda = -1`` et ``y_0 = 1``.

"""

# ╔═╡ 3ab81476-f7f8-11ea-3633-09930c9cdffe
# Q4 -- À MODIFIER
solution(t, y = ones(1)) = y

# ╔═╡ 2f1bd88d-9d9c-4fd6-aef5-5d061994759f
if norm(solution(1.0) - [0.36788]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `solution` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `solution`.
	"""
end

# ╔═╡ 8a4674da-f7ec-11ea-2faf-4b332a41d7fc
md"""
# Schéma numérique

On rappelle que lors du cours précédent, tois schémas numériques ont été présentés, à savoir :
```math
y_{n + 1} - y_n - \tau f \left ( t_n, y_n \right ) = 0 \quad \text{(Euler explicite)},
```
```math
y_{n + 1} - y_n - \tau f \left ( t_{n + 1}, y_{n + 1} \right ) = 0 \quad \text{(Euler implicite)},
```
et
```math
y_{n + 1} - y_n - \tau f \left ( \frac{t_n + t_{n + 1}}{2}, \frac{y_n + y_{n + 1}}{2} \right ) = 0 \quad \text{(point milieu)}.
```

5. En vous inspirant de l'implémentation `explicit!` du schéma explicite d'Euler présentée ci-dessous, implémenter les fonctions `implicit!` (a) et `midpoint!` (b) dont la racine est ``y_{n + 1}``. On préservera le nombre et l'ordre des paramètres, au nombre de 3, à savoir

* `y` -- la solution précédente (``y _ n``) ;
* `τ` -- le pas de temps (``\tau``) ;
* `f` -- le modèle (``f``) ;
* `t` -- l'instant précédent, (``t _ n``).

"""

# ╔═╡ 639dbbd7-e0eb-4234-a37a-254c6f751a74
function explicit!(res, x, y, τ, f, t)
	res .= x .- y .- τ .* f(t, y)
	nothing
end

# ╔═╡ ba074cc6-64cb-4f0a-9996-e6306727663c
# Q5a -- À MODIFIER
function implicit!(res, x, y, τ, f, t)
	res .= x .- y
	nothing
end

# ╔═╡ ec517a8d-bdef-41e0-bb4e-fac9df1f3c24
# Q5b -- À MODIFIER
function midpoint!(res, x, y, τ, f, t)
	res .= x .- y
	nothing
end

# ╔═╡ c87b46c0-f7ec-11ea-2918-d306ffd1c2bd
md"""
# Intégration temporelle

Il reste à présent à assembler à implémenter la bouche d'intégration temporelle. Étant donnés

* Un modèle `f` ;
* Un pas de temps `τ` ;
* Et un instant `s`

la fonction `integrate` implémentée ci-dessous retourne deux vecteurs, le premier contenant les instants
```math
t_0 \quad t_1 \quad \cdots \quad t_N = s
```
et le second la solution numérique, à savoir
```math
y_0 \quad y_1 \quad \cdots \quad y_N.
```

"""

# ╔═╡ 1887218b-4a77-4e76-8a97-54d90e69b419
function integrate(scheme!, f, τ, s)
	t, y = 0.0, ones(1)
    T, Y = [t], [y]

	while t < (1 - √eps(t)) * s

		y = solve(y) do res, x
			scheme!(res, x, y, τ, f, t)
		end
        t += τ

        push!(Y, y)
        push!(T, t)
	end

	T, Y
end

# ╔═╡ a691aa87-e01d-40f6-84ae-29b9ab13fb2b
md"""
La solution numérique peut être obtenue et visualisée comme suit.

"""

# ╔═╡ b9a3b65e-9c81-499d-acc5-85afa8d6703b
begin
	local T, Y = integrate(explicit!, linear, 0.1, 1.0)
	local fig = plot()
	scatter!(fig, T, first.(Y), label = "num")
	plot!(fig, t -> first.(solution.(t)), label = "exact")
end

# ╔═╡ 6be4f6f8-58e3-4e23-b34e-514e1045d08e
md"""
6. Utiliser les fonctions `linear` et `solution` définie précédemment dans l'implémentation de la fonction `error` ci-dessous, qui calcule l'erreur
```math
y_N - y \left ( t_N \right )
```
en fonction du schéma (`scheme!`) et du pas en temps (`τ`).

"""

# ╔═╡ 7fafbd51-99a9-4fea-bebf-6160e62a3ef4
# Q6 -- À MODIFIER
function error(scheme!, τ, s)
	T, Y = integrate(explicit!, linear, τ, s)
	norm(zero(Y))
end

# ╔═╡ d764d914-7a16-434a-8954-1f7233ee601c
md"""
7. Calculer (en utilisant la fonction `error`) et reporter les erreurs à l'instant `s = 1.0` dans le tableau ci-dessous. Commenter.

|             | `explicit!` | `implicit!` | `midpoint!` |
|:-----------:|:-----------:|:-----------:|:-----------:|
| `τ = 0.125` |             |             |             |
| `τ = 0.25`  |             |             |             |
| `τ = 0.5`   |             |             |             |
| `τ = 1.0`   |             |             |             |

8. On se place maintenant sur un horizon temporel plus long (`s = 10.0`). Augmenter la taille du pas de temps et commenter.

"""

# ╔═╡ d4637d80-f8c1-11ea-1f7f-df462373ca2d
md"""

Tout l'intérêt de l'utilisation du package `NLsolve.jl` est que notre implémentation fonctionne pour les problèmes scalaires **non-linéaires**, ainsi que pour les cas **vectoriels**.

# Au delà du cas linéaire

9. Utiliser ou imaginer un modèle scalaire non-linéaire en modifier la fonction `nonlinear` ci-dessous, et visualiser votre solution numérique pour chacun des trois schémas sur le même graphique. Vous pourrez par exemple utiliser la question 3 de l'exercice vu en TD :
```math
f \colon \left ( t, y \right ) \mapsto 2t - y ^ 2.
```

"""

# ╔═╡ 91a712b2-f8bd-11ea-3b8c-1bfbd521d29a
# Q9 -- À MODIFIER
nonlinear(t, y) = zero(y)

# ╔═╡ 580e8356-4fd2-47e1-a5a4-7063998b4ecb
begin
	local T, Y = integrate(explicit!, nonlinear, 0.1, 1.0)
	local fig = plot()
	scatter!(fig, T, first.(Y), label = "num")
	plot!(fig, t -> first.(solution.(t)), label = "exact")
end

# ╔═╡ 46222752-f91a-11ea-372e-2dde47c81add
md"""
# Au delà du cas scalaire

On se propose de résoudre l'[équation de prédation de Lotka-Volterra](https://fr.wikipedia.org/wiki/%C3%89quations_de_pr%C3%A9dation_de_Lotka-Volterra) :

> En mathématiques, les équations de prédation de Lotka-Volterra, que l'on désigne aussi sous le terme de "modèle proie-prédateur", sont un couple d'équations différentielles non-linéaires du premier ordre, et sont couramment utilisées pour décrire la dynamique de systèmes biologiques dans lesquels un prédateur et sa proie interagissent. Elles ont été proposées indépendamment par Alfred James Lotka en 1925 et Vito Volterra en 1926.

Le système d'équations s'écrit :
```math
\left \{ \begin{aligned}
\dot{x} \left ( t \right ) & = x \left ( t \right ) \left [ \alpha - \beta y \left ( t \right ) \right ], \\
\dot{y} \left ( t \right ) & = y \left ( t \right ) \left [ \delta x \left ( t \right ) - \gamma \right ]
\end{aligned} \right .
```
où

* ``t`` est le temps ;
* ``x \left ( t \right )`` est l'effectif des proies à l'instant ``t`` ;
* ``y \left ( t \right )`` est l'effectif des prédateurs à l'instant ``t``.

Les paramètres suivants enfin caractérisent les interactions entre les deux espèces :
```math
\alpha = 0.1, \quad \beta = 0.003, \quad \gamma = 0.06, \quad \delta = 0.0012.
```

10. Implémenter la fonction `lotka` correspond.
11. Résoudre et visualiser les solutions numérique correspondant à chacun des trois schéma.

"""

# ╔═╡ Cell order:
# ╟─d61e1ea6-f924-11ea-00dc-794c93177d22
# ╠═35a850c8-d24c-4593-b551-f3b9bb49e82d
# ╠═ec1b0cec-af5f-492e-ba55-269c2df5b458
# ╠═e8b460cf-3d35-4bbb-8b3a-5a2c5f8a7a21
# ╟─3e5d9fa5-acae-42ab-b40e-e1a197edcad2
# ╟─4ae18622-f7ec-11ea-2f71-d5b166ff50fb
# ╠═88dd8a5e-651a-4939-9231-6696b78f024c
# ╠═7dccedba-aa05-4821-b1ed-ed873ad9cccb
# ╟─073d29d8-7055-4d0c-8565-8732017b89d0
# ╠═bea20ce9-88ff-45e4-bb5a-494c4d2d19c2
# ╟─2c4151ac-a476-4c1f-a379-f4c3a0174b3d
# ╟─b72fd704-46f9-4d88-86bd-7add43c47d0e
# ╠═d264f45c-7a66-48d5-b6c3-b297073f5ce8
# ╟─42e83fa4-ba23-4e80-a8e6-58816319b134
# ╟─7dc345ee-f7ec-11ea-138f-a1c6b81e0260
# ╠═58162516-f7ec-11ea-3095-85bde6d71604
# ╟─dcf06721-1996-42ce-8c77-93f6f5195c7d
# ╟─56b5b3ec-7ec2-483c-8120-156665b2c49f
# ╠═3ab81476-f7f8-11ea-3633-09930c9cdffe
# ╟─2f1bd88d-9d9c-4fd6-aef5-5d061994759f
# ╟─8a4674da-f7ec-11ea-2faf-4b332a41d7fc
# ╠═639dbbd7-e0eb-4234-a37a-254c6f751a74
# ╠═ba074cc6-64cb-4f0a-9996-e6306727663c
# ╠═ec517a8d-bdef-41e0-bb4e-fac9df1f3c24
# ╟─c87b46c0-f7ec-11ea-2918-d306ffd1c2bd
# ╠═1887218b-4a77-4e76-8a97-54d90e69b419
# ╟─a691aa87-e01d-40f6-84ae-29b9ab13fb2b
# ╠═b9a3b65e-9c81-499d-acc5-85afa8d6703b
# ╟─6be4f6f8-58e3-4e23-b34e-514e1045d08e
# ╠═7fafbd51-99a9-4fea-bebf-6160e62a3ef4
# ╟─d764d914-7a16-434a-8954-1f7233ee601c
# ╟─d4637d80-f8c1-11ea-1f7f-df462373ca2d
# ╠═91a712b2-f8bd-11ea-3b8c-1bfbd521d29a
# ╠═580e8356-4fd2-47e1-a5a4-7063998b4ecb
# ╟─46222752-f91a-11ea-372e-2dde47c81add
