---
layout: post
title:  "Projection On A Set"
date:   2017-03-04
comments: true
categories: Optimization
---

* content
{:toc}

I am recently reading the code of TFOCS (Templates for First Order Conic Solvers). You can find it [here](https://github.com/cvxr/TFOCS "Github Repo") on Github. One of the key parts in TFOCS is to find the projection of a given point $x_0$ (in $R^n$ or $R^{m \times n}$) onto a convex set $\mathcal{C}$. I will try to make a neat summarization and give concrete (and useful) examples of projections on different convex sets. 

&nbsp;

### 1. General Definition