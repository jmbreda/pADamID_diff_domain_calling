# pADamID_diff_domain_calling
Calling domains with significantly different pADamID signal

## Model

The sign of a pADamID differential signal in successive genomic bins can be seen a coin toss.

The likelihood of sampling $k$ heads in $n$ coin toss with a coin having probability $p$ to hit a head and $(1-p)$ to hit a tail follows a binomial distribution:
$$P(k,n|p) = {n \choose k} p^k (1-p)^{n-k}$$

A coin is called *unbiased* or *fair* if $p=0.5$

Similarily, in a genomic region with a differential pADamID signal with probability $p$ to be positive, the number of positive bins $k$ in a window of size $n$ is binomial. A region that does not show strong evidence against the null hypothesis $p=0.5$ is neutral to the condition effect.

Previous domain calling methods used for LAD calling in pADamID signals defined

$$p = \frac{\sum_{i\in \text{bins}} \mathbb{I}_{x_i > 0.2}}{\sum_{i\in\text{bins}} 1}$$

with $x_i$ the differential pADamID signal in bin $i$.

A given window is significantly positive if

$$\sum_{l=k}^n P(n,l|p) < 10^{-10}$$

Being interested in the inference of $p$ rather than the data $k,n$, we would like to know $P(p\mid k,n)$ rather than $P(k,n\mid p)$

$$\rightarrow \text{Bayes' theorem} : P(A\mid B) = \frac{P(B \mid A) P(A)}{P(B)} $$

In our case:
$$P(p\mid k,n) = P(p) \frac{P(n,k\mid p)}{P(k,n)} $$
$$P(p\mid k,n) = P(p) \frac{P(n,k\mid p)}{\int_0^1 dp P(k,n\mid p)P(p) dp} $$

using a uniform prior $P(p)\sim U[0,1]$
$$P(p\mid k,n) = \frac{P(n,k\mid p)}{\int_0^1 dp P(k,n\mid p) dp} $$
$$P(p\mid k,n) \sim Beta(k+1,n-k+1)$$

A region significantly positive/negative shows a strong evidence $p\neq 0.5$

For a region to be significatly positive/negative, I ask that at least 99.9 \% of the posterior distribution is above/below $p=0.5$.



$$\text{Positive region:} \int_{p=0.5}^1 P(p\mid k,n) > 0.999 ~~~~~~ \text{Negative region:} \int_{p=0}^{0.5} P(p\mid k,n) > 0.999$$

![pADamID domains](/Fig/beta_cartoon.pdf?raw=true "Beta distribution")
