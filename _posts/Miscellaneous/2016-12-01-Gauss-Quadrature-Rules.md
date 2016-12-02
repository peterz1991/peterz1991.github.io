---
layout: post
title:  "Gauss Quadrature Rules"
date:   2016-12-01
comments: true
categories: Miscellaneous
---
* content
{:toc}

This note is served as a summarization of Gauss quadrature rules. There are lots of theoretical results either in the textbooks or papers. The many focus of this note, however, is to implement it with actual code from scratch. I will also try to provide some theoretical background when necessary and provide the sources where the readers can find the proofs. 

### 1. Big Picture

In this section, I will set up a toy problem and the process to tackle it with Gauss quadrature rules. In the later sections, I will try to introduce more applications. 

__Goal:__ Our target is to estimate 

$$\int_{-1}^{1} f(x) dx, $$

using numerical approximations (assume that the function is 'nice').

We all know how to do this with Riemann sum:


&nbsp;
$$\int_{-1}^{1} f(x) dx \approx \sum_{i=1}^{n} f(x_i^*)(x_{i+1}- x_i),$$
&nbsp;


where $-1=x_0\leq x_1 \leq \cdots \leq x_n = 1$ is a partition of the interval [-1,1], $x_i^*$'s nodes selected in between $x_i$ and $x_{i+1}$.
Gauss Quadrature rules aim to do the same thing. It tries to do the approximation in the following way:

&nbsp;
$$\int_{-1}^{1} f(x) dx = \sum_{i=1}^{n}w_i g (x_i) + R[f], \tag{Gauss Quadrature}$$
&nbsp;

where $g(x)$ will be a modified version of $f(x)$ according to the ONP, $R[f]$ is the residual of the approximation which depends on the regularity of $f(x)$. Of course, no matter what kind of approximations we do, the errors should always depend on both the approximation rules and $f(x)$ itself. The __art__ part of Gauss quadrature rules is to find the $w_i$ (called weight) and $x_i$ (called node) such that $R[f]$ is as small as possible.

To make this part short, I will list several key facts of Gauss Quadrature rules. Later on I will try to provide some proofs.

__Key Fact:__

1.  Gauss Quadrature rules are indeed a class of rules. Each of them is associated with a choice of __orthonormal polynomials (ONPs)__. For example, Gauss quadrature rule with Legendre polynomials, Gauss quadrature rule with Chebyshev polynomials of first kind, etc.

2.  The __nodes $x_i$__ of a specific Gauss Quadrature rule are selected as either __zeros__ or the local __extrema__ of the associated ONPs. 

3.  Every ONP satisfies a __three-term-recurrence-relation (3TRR)__, which will implicit indicate that the ONP can be defined through a special matrix $J$ called __Jacobi matrix__. 

4.  Given a specific ONP, the eigenvalues of the Jacobi matrix $J$ are the zeros of it. In other words, it gives the nodes $x_i$.

5.  $J$ is _symmetric_ and has eigenvalue decomposition  $J = Q\Lambda Q^T$. The weights $w_i = c q_{1i}^2$, here $q_{1i}$ is the $(1,i)$th entry of matrix Q and $c$ is a constant depends on the ONP that is used. 

Now we summarize the whole procedure of using Gauss Quadrature rules to approximate the integral.

__Workflow:__

<div> <center><img src="/Img/Gauss_Quadrature.png" alt=" " align="middle"></center> </div>
