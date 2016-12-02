---
layout: post
title:  "Notes on Nonconformist Image Processing with Graph Laplacian"
date:   2016-05-23
comments: true
categories: Image_Processing
---

* content
{:toc}


_This is my notes on SIAM Imaging Conference, May 23, 8:15 am talk: Non-conformist image processing with graph Laplacian operators, by Peyman Milanfar, Google Research_

_The video talk and the slides can be found <a href="https://www.pathlms.com/siam/courses/2426/sections/3234" target="_blank">here</a>._

### Global vs. Local filters

Usually, we represent the filtering in the following way: 

$$z = Wy,$$ 

where y is a vectorized image, W is _a data-dependent matrix_, it can represent local filters or global filters:

- Local filter: Sparse, but high-rank
- Global filter: Dense, but low-rank

>__Q__: Remember if W just represent a kernel such as Gaussian, [-1,1], it will not be data dependent, how can we get a data dependent kernel then?


# 1 Laplacian Operators

In imaging:

- (Anisotropic) Diffusion
- Curvature flow
- Adaptive sharpening
- De-blurring
- ...

A Laplacian operation is just filtering/centered average:

$$\Delta z(x) = \frac{z(x+h)-2z(x) + z(x-h)}{h^2} = \frac{2}{h^2}[\frac{z(x+h) + z(x-h)}{2} - z(x).$$

Quick explanation of some terms:

- Conformity: uniform, homogeneous
- Gradient (Dirichlet) Energy: 

$$\min_z \int \|\nabla z(x) \|^2 dx.$$


## 1.1 Discrete (Graph) Laplacian

Given a point $z_i$, the Laplacian transform applies on $z_i$ is defined as:

$$L(z_i) = \sum_j K_{ij} (z_i-z_j)$$

Notice that, this is just a __local average with certain weights__ $$K_{ij}$$. Usually, the $$K_{ij}$$ depends on the _Gaussian distance_ of $$z_i$$ and $$z_j$$.

- Gaussian distance (Gaussian kernel) applied to the pixel value

$$K(z_i,z_j) = \exp(-\|z_i-z_j\|^2/h_z^2)$$

according to the author, this is called _photometric Gaussian Kernel_.
Another is _Spatial Gaussian Kernel_ (applied to the physical location of the ith pixel):

$$K(x_i,x_j) = \exp(-\|x_i-x_j\|^2/h_x^2)$$

Now 

$$K_{ij} = K(z_i,z_j)K(x_i,x_j).$$

It's called _bilateral kernel_.

Also notice that 

$$L(z_i) = 0 \to z_i = \frac{\sum_jK_{ij}z_j}{\sum_j K_{ij}}.$$ 

(maybe unrelated, but same formula importance sampling technique in Bayesian?) Anyway, it says the unbiased estimate of $z_i$ will just be an expectation based on the weights.

I haven't checked the proofs but intuitively it's not hard to imagine that with more and more samples (denser and denser grids), the discrete Laplacian will converge to the continuous one (Ref: Lafon '04', Belkin et al '05', Hein et al '05', Singer '06').

>__Q__: Wait, we know the $$K$$ matrix now, but where is $$L$$ (the Laplacian matrix)?

_Laplacian matrix_ $L$ is defined exactly as in graph theory: 

First, we define the degree of a vertex (pixel), 

$$d_i = \sum_j K_{ij},$$ 

then the graph Laplacian is then 

$$L = D - K,$$ 

$$D = diag\{d_i\}_{i = 1}^n$$.
Here we consider K is not _normalized_, otherwise the sum of the weights $$\sum_j K_{ij}$$ will be __1__. Hence we get the definition $L = I - K$, which can be found in most of the graph theory textbooks.

In pure matrix perspective, the graph Laplacian has many interesting properties, the following is a table illustrate the many definitions of graph Laplacian and their properties:

| Graph Laplacian           | Symmetric | DC eigenvector | Spectral Range | Reference                      |
|---------------------------|-----------|----------------|----------------|--------------------------------|
| $L = D - K$               | Yes       | Yes            | [0, n]         | Un-normalized Laplacian        |
| $I - D^{-1/2} K D^{-1/2}$ | Yes       | No             | [0, 2]         | Normalized, Chung '97'         |
| $I - D^{-1}K$             | No        | Yes            | [0, 1]         | Random Walk Laplacian          |
| $I - C^{-1/2}KC^{-1/2}$   | Yes       | Yes            | [0, 1]         |'Sinkhorn' Laplacian, M. '13'  |

A naive version of normalized Laplacian is $$L = \alpha (D- K)$$ And the corresponding filter $W$ is defined as $W = I - L = I - \alpha (D - K)$.

It's kind of confusing if you think about the _regular_ graphs where each vertex has the same degree, the above definition of $W$ gives back of $K$!
You may also wonder why we need $W$ here, indeed, _suppose_  you are interested in minimizing the discrete Dirichlet energy $\|z^TLz$ (it's not hard to show that this is the discrete version of the Dirichlet energy I mentioned before), then by gradient descent you will contour the update

$$z_{k+1} = z_k - L z_k = (I - L) z_k : = W z_k.$$

Therefore, given a Laplacian matrix and an initial state of the image $z_0$, to minimize the Dirichlet energy, we just need to calculate $W$ first, and iterative update $z_k$ by $z_{k+1} = Wz_k$.

BTW, the Laplacian matrix is __positive semi-definite__, so the above target function is non-negative. There are many ways to prove this, one quick way is that $L$ by definition is _diagonal dominant_. The following theorems help me a lot and turn out to be very useful in many situations:

<div> <center><img src="/Img/Diag_dominant.png" alt=" " align="middle" height="150" width="600"></center> </div>
<div> <center><img src="/Img/Diag_dominant2.png" alt=" " height="300" width="600"></center> </div>
It can be shown that a Hermitian diagonally dominant matrix $A$ with real non-negative diagonal entries is positive semidefinite.

# 2 Application in Image Processing

## 2.1 De-blurring: adaptive sharpening

We can use the following model to do de-blurring, given blur matrix $A$:

$$(y-Az)^T(I+\beta L)(y - Az) + \eta z^TLz$$

Usually, I'd like to write the _fidelity_ as $$\|y-Az\|_{I+\beta L}^2$$. Basically, the difference between traditional norm $$\|\cdot\|_2$$ and matrix $$\|\cdot\|_\Sigma$$ is that $$\Sigma$$ may indicate correlation/covariance/weighted balance of the residual. It's more  flexible.

Most of the time we don't have the knowledge of PSF (Point spread function), then we may simply solve the above optimization problem with $$A = I$$, $$\eta = 0$$. Then
the solution will be $$z = (I + \beta L) y$$, for some $$\beta > 0$$ ( a hand tunning parameter). This is called _adaptive sharpening "nonlinear unsharp mask"_.

The following image is cropped from 
>Motion De-blurring With Graph Laplacian Regularization, Kheradmand & Milanfar, 2015.

<div> <center><figure> <img src="/Img/Laplacian_deblurring.png" alt="De-blurring with Laplacian" height="500" width="500"> </figure> </center> </div>

<figcaption>Fig1. Real motion deblurring example: (a) input blurred noisy image, (b) Output of hyper-Laplacian algorithm (c) output of "High-quality motion deblurring from a single image" and (d) output of proposed algorithm (η = 0.25, β = 2.5).</figcaption> 


The code is <a href="https://users.soe.ucsc.edu/~aminkh/KernelRestoration.html" target="_blank">here</a>.

### Let's stop talking:

Here is some Matlab code on de-blurring with Laplacian from a NIPS paper:

>Fast Image Deconvolution using Hyper-Laplacian Priors, Krishnan & Fergus, NIPS, 2009. 

Code is <a href="http://cs.nyu.edu/~dilip/research/fast-deconvolution/" target="_blank">here</a>.

I will also write some code and will upload it later in github.

## 2.2 Linear Embedding Using Laplacian

___I guess the this is the most interesting part of this post: Laplacian for manifold and hierarchical decomposition.___

- Dimension Reduction
- Nystrom Approximation
- Multi-Laplacian Scale Decomposition



### 2.2.1 Dimension Reduction (DR)

>A serious question before any fancy technique, __what is dimension reduction__? Think it again not just intuitively but also try to write some math definitions.

Recall the way we use PCA to do DR, we do the eigen-decomposition and then take the largest n eigenvalues and corresponding eigenvectors to approximate the original data matrix. Here the idea is the same, but little tricky.

We think the data points formulate the vertices of a 'large' graph.

The Laplacian DR algorithm basically goes with three steps:

 1.  Given data points $x_1, x_2,...,x_k \in R^l$, build the graph: connect 'near' points with edges.
 2.  Calculate the 'Gaussian' distance matrix W along the edges:
 
$$W_{ij} = e^{-\frac{\|x_i-x_j\|^2}{t}}$$

 3.  Solve the following optimization problem:
 
$$\min_y y^TLy,  $$

$$s.t. \hspace{0.2in} y^TDy = 1.$$

This will reduce to a generalized eigenvalue problem:

$$Ly = \lambda D y,$$

where $D_{ii} = \sum_{j} W{ji}$, $L = D - W$, $L, D, W\in R^{k \times k}$.

Suppose we solved the optimization problem and get the eigenvector $y_0, y_1, ..., y_{k-1}$, corresponding eigenvalues are $0 = \lambda_0 \leq \lambda_1 \leq \lambda_2 \leq ... \leq \lambda_{k-1}$.

We can simply take the first few (say $m$) eigenvectors starting from $y_1$, and use them to 'describe' the original data:

$$x_i \to (y_1(i), y_2(i),..., y_m(i))$$

Just as we pick up the first few PCA components! Because these contains the key information!

 We will discuss it further in the next post. 








