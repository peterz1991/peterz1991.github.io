---
layout: post
title:  "Several theorems on Matrices"
date:   2016-04-24
comments: true
categories: Miscellaneous
---
* content
{:toc}

_Things to know for matrix theory._

Other than the classical theorems or factorizations from a matrix theory class, I will list some things both interesting and good to know.

- ## [Perron–Frobenius theorem](https://en.wikipedia.org/wiki/Perron%E2%80%93Frobenius_theorem, "wiki page"):

Let all the entries of matrix $$A \in R^{n\times n}$$ are positive (i.e. $$A$$ is positive square real matrix), then the following holds:

1.  $$A$$ has a real positive eigenvalue r that is __strictly larger__ than any other eigenvalues (they can be complex), that is, $$ r > \vert \lambda \vert $$, where $$\lambda$$ is eigenvalue of A.

2.  $$r$$ is unique, that is, the multiplicity is 1.
3.  There exists an eigenvector associated with $$r$$ such that all the entries of it are **positive**.
4.  An interesting inequality: 
      
$$\min_i \sum_{j} a_{ij} \leq r \leq \max_i \sum_{j} a_{ij},$$

namely, min of column sum $$< r <$$ max of column sum.

&nbsp;

- ## [Kronecker Product]():
    
It has been a long time that I really want to learn and prove the properties of the Kronecker product, as it pops up in so many places in image processing and finite element methods. Here I give several key properties of it, suppose we have

$$A_r\in R^{m \times m}, A_c \in R^{n \times n}, X \in R^{m \times n}$$

The following holds:

1.  Matrix Multiplication/Vectorization: 
    
    $$(A_r \otimes A_c)vec(X) = vec(A_cXA_r^T).$$

2.  Transpose and Inversion:

    $$(A_r \otimes A_c)' = A_r^T \otimes A_c^T , \hspace{0.2in} (A_r \otimes A_c)^{-1} = A_r^{-1} \otimes A_c^{-1}.$$

3.  Singular value decomposition: 

    $$(U_r\Sigma_rV_r^T) \otimes (U_c\Sigma_c V_c^T) = (U_r\otimes U_c) (\Sigma_r\otimes \Sigma_c)(V_r\otimes V_c)^T.$$

&nbsp;
- ## [Gerschgorin circle theorem](https://en.wikipedia.org/wiki/Gershgorin_circle_theorem "wiki page"):

&nbsp;

- ## [Schur Complement](https://en.wikipedia.org/wiki/Schur_complement "wiki page"):

&nbsp;

- ## [Sherman–Morrison–Woodbury formula](https://en.wikipedia.org/wiki/Woodbury_matrix_identity "wiki page"):
