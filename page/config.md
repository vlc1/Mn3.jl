<!--
Add here global page variables to use throughout your website.
-->
@def prepath = "Mn3.jl"

+++
author = "Vincent Le Chenadec"
mintoclevel = 2

# Add here files or directories that should be ignored by Franklin, otherwise
# these files might be copied and, if markdown, processed by Franklin which
# you might not want. Indicate directories by ending the name with a `/`.
# Base files such as LICENSE.md and README.md are ignored by default.
ignore = ["node_modules/"]

# RSS (the website_{title, descr, url} must be defined to get RSS)
generate_rss = true
website_title = "Méthodes Numériques pour l'Energétique"
website_descr = "Enseignement des Méthodes Numériques pour l'Energétique à l'ESIEE en première année filière énergie en apprentissage"
website_url   = "https://vlc1.github.io/Mn3.jl/"
+++

<!--
Add here global latex commands to use throughout your pages.
-->
\newcommand{\R}{\mathbb R}
\newcommand{\scal}[1]{\langle #1 \rangle}
