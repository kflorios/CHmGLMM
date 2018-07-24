# CHmGLMM
Multivariate Generalized Linear Mixed Models (mGLMMs) with CLIC Heuristic Averaging

The R package CHmGLMM implements the method in the technical report

K. Florios, I. Moustaki, D. Rizopoulos, V.G.S. Vasdekis, 
A proposal for creating a weighted composite likelihood estimator using CLIC, 
Technical Report, Athens University of Economics and Business, 2015.

Abstract of report:
Composite likelihood estimation has been proposed in the literature for handling
intractable likelihoods. In particular, pairwise likelihood estimation has been recently
proposed to estimate models with latent variables and random effects that involve high
dimensional integrals. Pairwise estimators are asymptotically consistent and normally
distributed but not the most efficient among consistent estimators.
Vasdekis et al. (2014) proposed a weighted estimator (WAVE) that is found to
be more efficient than the unweighted pairwise estimator produced by separate max-
imizations of pairwise likelihoods. Florios et al. (2015), proposed a modification to
that weighted estimator (DWAVE) that lead to simpler computations and studied its
performance through simulations and a real application.
In this paper, we propose an even simpler weighted estimator (CWAVE), based
on the concept of the Composite Likelihood Information Criterion (CLIC) (Varin &
Vidoni, 2005) which seems to combine the strengths of the AVE estimator (Fieuws &
Verbeke, 2006) for the random effects parameters and the DWAVE/WAVE estimators
for the fixed effects parameters. The new estimator CWAVE performs very well in both
fixed and random effects parameters identification, with high coverage, especially when
T is large, regardless of the size of N.

# Install
install.packages("devtools")

require(devtools)

install_github("kflorios/CHmGLMM")
