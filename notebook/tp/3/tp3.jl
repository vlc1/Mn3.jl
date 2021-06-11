### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ a72fac18-0fd9-49b0-8cce-80b5863c8c8a
using LinearAlgebra, NLsolve, Plots

# ╔═╡ 1a4c8af6-9ba2-416e-ae49-1d3e7715f441
md"""
# Prérequis

!!! warning "Remarque"

	Les cellules suivantes ont été copiées du TP précédent. Merci de ne pas les modifier.

	Ce TP contient trois exercices indépendants. Cependant, dans un exercice donné, les **questions doivent être traitées séquentiellement**.

"""

# ╔═╡ 82394153-ce46-4000-9c00-374105096879
explicit!(res, x, y, τ, f, t) = res .= x - y - τ * f(t, y)

# ╔═╡ bf92265c-377a-49b3-8715-25e9d0dbda28
implicit!(res, x, y, τ, f, t) = res .= x - y - τ * f(t + τ, x)

# ╔═╡ 5016a984-930e-452a-ad9c-fedd124610e3
midpoint!(res, x, y, τ, f, t) = res .= x - y - τ * f(t + τ / 2, (x + y) / 2)

# ╔═╡ bb5e9ba4-45e2-4910-a446-22aab6ee2f75
function cauchy(scheme!, f, τ, s, y₀, t₀ = zero(τ))
	t, y = t₀, y₀
    T, Y = [t], [y]

	while t < (1 - √eps(t)) * s
		y = getproperty(
			nlsolve(y) do res, x
				scheme!(res, x, y, τ, f, t)
			end,
			:zero)
        t += τ

        push!(Y, y)
        push!(T, t)
	end

	T, Y
end

# ╔═╡ da87e5ef-1f16-46e8-aa15-e1a21b85c7b9
md"""
# Méthode de Runge-Kutta explicite

La fonction `rk2` définie ci-dessous correspond à la méthode de Runge-Kutta explicite (à 2 étapes) donnée par les coefficients suivants :
```math
A_2 = \left ( \begin{matrix}
0 & 0 \\
\frac{1}{2} & 0
\end{matrix} \right ), \quad b_2 = \left ( \begin{matrix}
0 \\
1
\end{matrix} \right )
\quad \mathrm{et}
\quad c_2 = \left ( \begin{matrix}
0 \\
\frac{1}{2}
\end{matrix} \right ).
```

"""

# ╔═╡ e61c6973-177f-4885-a83c-b638a8cad271
function rk2!(res, x, y, τ, f, t)
	k₁ = f(t, y)
	k₂ = f(t + τ / 2, y + τ * k₁ / 2)
	res .= x - y - τ * k₂
end

# ╔═╡ 69538fce-c689-431f-950a-b6d5d1d7155b
md"""

1. En vous inspirant de l'exercice réalisé en TD et de la fonction `rk2` définie ci-dessus, finalisez l'implémentation `rk4` de la méthode à 4 étapes donnée par les coefficients
```math
A_4 = \left ( \begin{matrix}
0 & 0 & 0 & 0 \\
\frac{1}{2} & 0 & 0 & 0 \\
0 & \frac{1}{2} & 0 & 0 \\
0 & 0 & 1 & 0
\end{matrix} \right ), \quad b_4 = \left ( \begin{matrix}
\frac{1}{6} \\
\frac{1}{3} \\
\frac{1}{3} \\
\frac{1}{6}
\end{matrix} \right )
\quad \mathrm{et}
\quad c_4 = \left ( \begin{matrix}
0 \\
\frac{1}{2} \\
\frac{1}{2} \\
1
\end{matrix} \right ),
```

"""

# ╔═╡ 818dbcb1-df69-4ad5-ac31-3ad66e6d101a
# Q1 -- À MODIFIER
function rk4!(res, x, y, τ, f, t)
	k₁ = f(t, y)
	k₂ = f(t, y)
	k₃ = f(t, y)
	k₄ = f(t, y)
	res .= x - y
end

# ╔═╡ 396887c8-4862-4400-9426-f34ef3b308ec
if norm(rk4!(zeros(1), ones(1), ones(1), √2, (t, y) -> t .+ y, π) - [-14.2794]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `rk4!` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `rk4!`.
	"""
end

# ╔═╡ c64309e8-c39b-11eb-2841-37d321df8cd0
md"""
# L'équation de Lotka-Volterra

On se propose de résoudre l'[équation de prédation de Lotka-Volterra](https://fr.wikipedia.org/wiki/%C3%89quations_de_pr%C3%A9dation_de_Lotka-Volterra) :

> En mathématiques, les équations de prédation de Lotka-Volterra, que l'on désigne aussi sous le terme de "modèle proie-prédateur", sont un couple d'équations différentielles non-linéaires du premier ordre, et sont couramment utilisées pour décrire la dynamique de systèmes biologiques dans lesquels un prédateur et sa proie interagissent. Elles ont été proposées indépendamment par Alfred James Lotka en 1925 et Vito Volterra en 1926.

Le système d'équations s'écrit :
```math
\left \{ \begin{aligned}
\dot{y}_1 \left ( t \right ) & = y_1 \left ( t \right ) \left [ \alpha - \beta y_2 \left ( t \right ) \right ], \\
\dot{y}_2 \left ( t \right ) & = y_2 \left ( t \right ) \left [ \delta y_1 \left ( t \right ) - \gamma \right ]
\end{aligned} \right .
```
où

* ``t`` est la variable indépendante (le temps) ;
* ``y_1 \left ( t \right )`` est l'effectif des proies à l'instant ``t`` ;
* ``y_2 \left ( t \right )`` est l'effectif des prédateurs à l'instant ``t``.

Les paramètres suivants enfin caractérisent les interactions entre les deux espèces :
```math
\alpha = 0.1, \quad \beta = 0.003, \quad \gamma = 0.06, \quad \delta = 0.0012.
```

2. Implémenter la fonction `lotka` correspondant.

"""

# ╔═╡ 06256187-e634-4410-9a14-2dd9f5194ecf
# Q2 -- À MODIFIER
function lotka(t, y; α = 0.1, β = 0.003, γ = 0.06, δ = 0.0012)
	rhs = similar(y)
	rhs[1] = -y[1]
	rhs[2] = -y[2]
	rhs
end

# ╔═╡ 77343961-9f9f-4f51-a6cc-88c31cd68203
if norm(lotka(π, [√2; 1/3]) - [0.140007; -0.0194343]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `lotka` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `lotka`.
	"""
end

# ╔═╡ 08fc9e11-ec90-4521-b305-08d2bbeb83cd
md"""
3. Résoudre et visualiser les solutions numériques avec les schémas `explicit!`, `implicit!` et `midpoint!` pour les paramètres suivants : ``\tau = 1``,  ``s = 1000`` et ``y \left ( 0 \right ) = (10, 20)``. Commenter vos résultats.

"""

# ╔═╡ ad08bb68-0e9e-4a39-b84a-ed41392dd425
begin
	local T, Y = cauchy(midpoint!, lotka, 1.0, 1000.0, [10.; 20.])
	local fig = plot()
	plot!(fig, T, first.(Y), label = "Proies")
	plot!(fig, T, last.(Y), label = "Prédateurs")
end

# ╔═╡ d9972d53-ffed-4fbc-b8ee-a80501b5bb3a
md"""
# Le problème de Blasius

Comme vu en TD, l'équation de Blasius peut se réécrire sous la forme suivante :
```math
\left \{ \begin{aligned}
\dot{y}_1 & = y_2, \\
\dot{y}_2 & = y_3, \\
\dot{y}_3 & = -y_1 y_3.
\end{aligned} \right .
```

4. Finaliser l'implémentation de la fonction `blasius` ci-dessous pour qu'elle correspondre au second membre du système décrit ci-dessus.

"""

# ╔═╡ 6d558f15-1ffa-4524-a29a-6e12caf5400e
# Q4 -- À MODIFIER
function blasius(t, y)
	rhs = similar(y)
	rhs[1] = -y[1]
	rhs[2] = -y[2]
	rhs[3] = -y[3]
	rhs
end

# ╔═╡ 4bf1aea3-8fd7-44c6-92fb-fe649350e3f8
if norm(blasius(nothing, [√2; 1/3; π]) - [0.333333; 3.14159; -4.44288]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `blasius` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `blasius`.
	"""
end

# ╔═╡ 5700a08d-4cab-42fd-bc38-f75633c2051a
md"""
La symétrie identifiée en TD peut être exploitée en résolvant dans un premier temps cette équation avec la condition initiale
```math
\left \{ \begin{aligned}
y_1 \left ( 0 \right ) & = 0, \\
y_2 \left ( 0 \right ) & = 0, \\
y_3 \left ( 0 \right ) & = 1.
\end{aligned} \right .
```

On a vu en TD que le paramètre ``\alpha \equiv y_2 \left (\infty \right )`` permet de calculer la fonction de Blasius (``u``) à partir de la solution du problème de Cauchy présenté ci-dessus.

5. Implémenter la fonction `alpha` ci-dessous qui, étant donnés un schéma numérique, un pas de temps `τ` et un horizon temporel `s`, retourne ``\alpha``.

"""

# ╔═╡ 6f3aa892-c589-45ee-a4e8-56f66e471219
# Q5 -- À MODIFIER
function alpha(scheme!, τ, s)
	one(τ)
end

# ╔═╡ d713869e-d1ba-43c1-b924-e68e0e67015a
if abs(alpha(midpoint!, 0.1, 1000.) - 1.6537) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `alpha` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `alpha`.
	"""
end

# ╔═╡ c5b7bb06-b553-426f-bd23-afcc9bf0f476
md"""
6. Pour le schéma de votre choix, assurez-vous que votre choix de ``\tau`` et ``s`` garantisse au moins 4 chiffres significatifs de ``\alpha`` en remplissant **quelques cellules** du tableau suivant.

|   ``\alpha``    | ``s = 10`` | ``s = 100`` | ``s = 1000`` | ``s = 10000`` |
|:---------------:|:----------:|:-----------:|:------------:|:-------------:|
| ``\tau = 0.01`` |            |             |              |               |
| ``\tau = 0.1``  |            |             |              |               |
| ``\tau = 1.0``  |            |             |              |               |

7. Définissez la fonction de Blasius (`solution`).

"""

# ╔═╡ cdf255c7-494d-41b6-9328-0e040b83fde0
# Q7 -- À MODIFIER
function solution(x; scheme! = midpoint!, τ = 0.1, s = 1000.)
	zero(τ)
end

# ╔═╡ c293dd14-ca7c-4a08-885f-c1adb19fbe6e
if abs(solution(1.0) - 0.382) ≤ 1e-3
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

# ╔═╡ db09724a-c5e4-4486-a8a7-196c8afea296
md"""
8. Commenter l'allure du graphe suivant.

"""

# ╔═╡ 71ca93bd-debc-449f-aba5-71585ea771fe
plot(0:0.1:3.0, solution.(0:0.1:3.0), label = "Blasius")

# ╔═╡ Cell order:
# ╟─1a4c8af6-9ba2-416e-ae49-1d3e7715f441
# ╠═a72fac18-0fd9-49b0-8cce-80b5863c8c8a
# ╠═82394153-ce46-4000-9c00-374105096879
# ╠═bf92265c-377a-49b3-8715-25e9d0dbda28
# ╠═5016a984-930e-452a-ad9c-fedd124610e3
# ╠═bb5e9ba4-45e2-4910-a446-22aab6ee2f75
# ╟─da87e5ef-1f16-46e8-aa15-e1a21b85c7b9
# ╠═e61c6973-177f-4885-a83c-b638a8cad271
# ╟─69538fce-c689-431f-950a-b6d5d1d7155b
# ╠═818dbcb1-df69-4ad5-ac31-3ad66e6d101a
# ╟─396887c8-4862-4400-9426-f34ef3b308ec
# ╟─c64309e8-c39b-11eb-2841-37d321df8cd0
# ╠═06256187-e634-4410-9a14-2dd9f5194ecf
# ╟─77343961-9f9f-4f51-a6cc-88c31cd68203
# ╟─08fc9e11-ec90-4521-b305-08d2bbeb83cd
# ╠═ad08bb68-0e9e-4a39-b84a-ed41392dd425
# ╟─d9972d53-ffed-4fbc-b8ee-a80501b5bb3a
# ╠═6d558f15-1ffa-4524-a29a-6e12caf5400e
# ╟─4bf1aea3-8fd7-44c6-92fb-fe649350e3f8
# ╟─5700a08d-4cab-42fd-bc38-f75633c2051a
# ╠═6f3aa892-c589-45ee-a4e8-56f66e471219
# ╟─d713869e-d1ba-43c1-b924-e68e0e67015a
# ╟─c5b7bb06-b553-426f-bd23-afcc9bf0f476
# ╠═cdf255c7-494d-41b6-9328-0e040b83fde0
# ╟─c293dd14-ca7c-4a08-885f-c1adb19fbe6e
# ╟─db09724a-c5e4-4486-a8a7-196c8afea296
# ╠═71ca93bd-debc-449f-aba5-71585ea771fe
