For Fem, integration is important.
I want to create a method to integrate an arbitrary polynomial over an element.
Start with the terminal state of a recursive function.
$$\int_0^1x_0^{n}dx_0 = \frac{1}{n+1}$$
Summation of integration allows me to split the method
$$\int_\Omega\sum_if_i(\vec{x})d\vec{x} = \sum_i\int_\Omega f_i(\vec{x})d\vec{x}$$
Any polynomial function can be written in a similar form
$$\sum_i^ma_i\prod_j^nx_j^{b_{ij}}$$
### Arbitrary simplex element
The domain of an element needs describing to start with.
$$\Omega = (\vec{x} | x_i\in [0, (1-\sum_j^{i-1}x_j)])$$
Integrating over an element is starting to look cumbersome
$$\prod_j^n\big[ \int_{x_j=0}^{(1-\Sigma_{i}^{j-1} x_i)} \big] \big(f(\vec{x})\big) \prod_j^n \big[dx_j\big]$$
But recursion should be possible...
$$\prod_j^n\big[ \int_{x_j=0}^{(1-\Sigma_{i}^{j-1} x_i)} \big] \big( 
\prod_j^n x_j^{b_j}
\big) \prod_j^n \big[dx_j\big] =
\prod_j^{n-1}\big[ \int_{x_j=0}^{(1-\Sigma_{i}^{j-1} x_i)} \big] \big(
\prod_j^{n-1} x_j^{b_j}
\big)
\int_{x_n=0}^{(1-\Sigma_{i}^{n-1} x_i)} x_n^{b_n}
dx_n
\prod_j^{n-1} \big[dx_j\big]$$
The internal integral can be freely evaluated.
$$\int_{x_n=0}^{(1-\Sigma_{i}^{n-1} x_i)} x_n^{b_n}
dx_n = 
\frac{1}{b_n+1}(1-\sum_i^{n-1}x_i)^{b_n+1} =
\frac{1}{b_n+1}\sum_k\big[\prod_i^{n-1}\big({b_n+1\choose c_{ki}}x_i^{c_{ki}}\big)(-1)^{(b_n+1-\Sigma_i^{n-1}c_{ki})}\big]$$
I do not know how to properly define k but the idea is clear. Find all combinations for the polynomial products.
